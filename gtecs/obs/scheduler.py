"""Class for monitoring the obs database and finding pointings to observe."""

import os
import threading
import time
import traceback
from copy import copy

import Pyro4

from astropy import units as u
from astropy.coordinates import AltAz, get_sun
from astropy.time import Time

from gtecs.common.logging import get_logger

from . import database as db
from . import params
from .astronomy import horizon_limit
from .scheduling import PointingQueue
from .slack import send_slack_msg


Pyro4.config.SERIALIZER = 'pickle'  # IMPORTANT - Can serialize Pointing objects
Pyro4.config.SERIALIZERS_ACCEPTED.add('pickle')


@Pyro4.expose
class Scheduler:
    """Database scheduler daemon class."""

    def __init__(self):
        # get a logger for the scheduler
        self.log = get_logger('scheduler', params.LOG_PATH,
                              log_stdout=True,
                              log_to_file=params.FILE_LOGGING,
                              log_to_stdout=params.STDOUT_LOGGING)
        self.log.info('Scheduler started')

        # scheduler variables
        self.running = False
        self.loop_time = 0
        self.check_time = 0
        self.check_period = params.SCHEDULER_CHECK_PERIOD
        self.force_check_flag = False
        self.update_lock = False
        self.update_queue = []
        self.query_lock = False

        self.readout_time = params.READOUT_TIME
        self.template_requirement = params.TEMPLATE_REQUIREMENT
        self.write_file = params.WRITE_QUEUE_FILE
        self.queue_file = os.path.join(params.QUEUE_PATH, 'queue_info')
        self.write_html = params.WRITE_QUEUE_PAGE

        # Get site data (only once, on init)
        self.site_data = db.get_site_info()
        self.telescopes = sorted([telescope_id
                                  for telescopes in [self.site_data[site_id]['telescopes'].keys()
                                                     for site_id in self.site_data]
                                  for telescope_id in telescopes])

        self.old_pointings = {telescope_id: [None] for telescope_id in self.telescopes}
        self.latest_pointings = {telescope_id: [None] for telescope_id in self.telescopes}

        self.is_night = {telescope_id: False for telescope_id in self.telescopes}
        self.pilot_query_time = {telescope_id: None for telescope_id in self.telescopes}
        self.database_update_time = 0

    def __del__(self):
        self.shutdown()

    def run(self, host, port, timeout=5):
        """Run the scheduler as a Pyro daemon."""
        self.running = True

        # Start thread
        t = threading.Thread(target=self._monitor_thread)
        t.daemon = True
        t.start()

        # Check the Pyro address is available
        try:
            pyro_daemon = Pyro4.Daemon(host, port)
        except Exception:
            raise
        else:
            pyro_daemon.close()

        # Start the daemon
        with Pyro4.Daemon(host, port) as pyro_daemon:
            self._uri = pyro_daemon.register(self, objectId='scheduler')
            Pyro4.config.COMMTIMEOUT = timeout

            # Start request loop
            self.log.info('Pyro daemon registered to {}'.format(self._uri))
            pyro_daemon.requestLoop(loopCondition=self.is_running)

        # Loop has closed
        self.log.info('Pyro daemon successfully shut down')
        time.sleep(1.)

    def is_running(self):
        """Check if the daemon is running or not.

        Used for the Pyro loop condition, it needs a function so you can't just
        give it self.running.
        """
        return self.running

    @property
    def uri(self):
        """Return the Pyro URI."""
        if hasattr(self, '_uri'):
            return self._uri
        else:
            return None

    def shutdown(self):
        """Shut down the running threads."""
        self.running = False

    # Internal thread
    def _monitor_thread(self):
        """Monitor the database and find priority Pointings."""
        self.log.info('Database monitor thread started')

        while self.running:
            self.loop_time = time.time()

            # TODO: send database reports to Slack, once a day?

            # Only check on the given period, or if forced
            if self.force_check_flag or (self.loop_time - self.check_time) > self.check_period:
                self.update_lock = True
                check_time = Time(self.loop_time, format='unix')
                msg = 'Checking database' + (' (forced update)' if self.force_check_flag else '')
                self.log.debug(msg)

                # Calculate night status for telescope sites
                self._check_night_times(check_time)

                # Check the connection to the pilots
                self._pilot_heartbeat(check_time)

                # Caretaker step to monitor and update the database
                database_updated = self._caretaker(check_time)
                if database_updated:
                    # If we just updated we need to wait to be sure that the timestamps are >1s
                    # in the past when we then get the Pointings queue.
                    # This is because in the database the time precision is only to the second,
                    # annoyingly, and we use the same timestamp for all updates.
                    check_time += 1 * u.s

                # Get the queue and find highest priority Pointings for each telescope
                self.log.debug('Updating queue')
                new_pointings = {telescope_id: [None] for telescope_id in self.telescopes}
                for site_id in self.site_data:
                    new_pointings.update(self._get_pointings(site_id, check_time))
                self.log.debug('Queue update complete')
                self.old_pointings = self.latest_pointings
                self.latest_pointings = new_pointings

                # Write status log line
                # TODO: account for horizons?
                try:
                    old_strs = []
                    new_strs = []
                    for telescope_id in self.telescopes:
                        old_pointing = self.old_pointings[telescope_id][0]
                        old_str = '{}:'.format(telescope_id)
                        if old_pointing is None:
                            old_str += 'None'
                        else:
                            pointing_str = str(old_pointing.db_id)
                            if not old_pointing.valid:
                                pointing_str = '*' + pointing_str
                            if old_pointing.too:
                                pointing_str = pointing_str + '!'
                            if old_pointing.current_telescope == telescope_id:
                                pointing_str = pointing_str + '°'
                            old_str += pointing_str
                        old_strs.append(old_str)
                        new_pointing = self.latest_pointings[telescope_id][0]
                        new_str = '{}:'.format(telescope_id)
                        if new_pointing is None:
                            new_str += 'None'
                        else:
                            pointing_str = str(new_pointing.db_id)
                            if not new_pointing.valid:
                                pointing_str = '*' + pointing_str
                            if new_pointing.too:
                                pointing_str = pointing_str + '!'
                            if new_pointing.current_telescope == telescope_id:
                                pointing_str = pointing_str + '°'
                            new_str += pointing_str
                        new_strs.append(new_str)
                    old_str = ' '.join(old_strs)
                    new_str = ' '.join(new_strs)
                    if new_str != old_str:
                        self.log.info('Telescope pointings: ' + new_str)
                except Exception:
                    self.log.error('Could not write current status')

                # Set these after the check has finished, so the check_queue() function has to wait
                self.check_time = self.loop_time
                self.force_check_flag = False
                self.update_lock = False

            time.sleep(0.1)

        self.log.info('Database monitor thread stopped')
        return

    # Internal functions
    def _check_night_times(self, check_time=None):
        """Calculate if it is night for each telescope at each site."""
        if check_time is None:
            check_time = Time.now()

        for site_id in self.site_data:
            altaz_frame = AltAz(obstime=check_time, location=self.site_data[site_id]['location'])
            sun_coords = get_sun(check_time)
            sun_altaz = sun_coords.transform_to(altaz_frame)
            for telescope_id in self.site_data[site_id]['telescopes']:
                self.is_night[telescope_id] = sun_altaz.alt.degree < params.SCHEDULER_SUNALT_LIMIT

    def _pilot_heartbeat(self, check_time):
        """Keep track of if we loose connection to any of the telescopes during the night."""
        with db.open_session() as session:
            for telescope_id in self.telescopes:
                # Check how long it's been since we had a query from this telescope.
                if self.pilot_query_time[telescope_id] is None:
                    continue
                if (check_time - self.pilot_query_time[telescope_id]).sec > 60:
                    # It's been more than a minute since we last had a query from this pilot.
                    if self.is_night[telescope_id]:
                        # Only log during the night time.
                        self.log.debug('No contact from Telescope {} since {}'.format(
                                       telescope_id, self.pilot_query_time[telescope_id].iso))

                    # If there is still a running pointing, mark it as interrupted.
                    telescope = db.get_telescope_by_id(session, telescope_id)
                    if telescope.current_pointing is not None:
                        pointing = telescope.current_pointing
                        running_time = (check_time - Time(pointing.running_time)).sec
                        obs_time = pointing.get_obstime(self.readout_time)
                        self.log.info('Pointing {} has been running for {:.1f}s/{:.1f}s'.format(
                                      pointing.db_id, running_time, obs_time))
                        self.update_queue.append((pointing.db_id, 'interrupted'))
                        # No force check needed here, we're already in one!

                # TODO: Send Slack alerts? Could be when a telescope first connects in the
                #       evening, and then if the connection is lost before sunrise.

    def _caretaker(self, check_time):
        """Monitor the database and update any Pointings."""
        database_updated = False
        with db.open_session() as session:
            # Process any Pointings that need updating
            while len(self.update_queue) > 0:
                pointing = self.update_queue.pop()
                try:
                    pointing_id = pointing[0]
                    new_status = pointing[1]
                    if new_status == 'running':
                        telescope_id = pointing[2]
                        db.mark_pointing_running(pointing_id, telescope_id)
                        self.log.info(f'Marked Pointing {pointing_id} as running'
                                      f' on Telescope {telescope_id}')
                        database_updated = True
                    elif new_status == 'completed':
                        db.mark_pointing_completed(pointing_id, schedule_next=True)
                        self.log.info(f'Marked Pointing {pointing_id} as completed')
                        database_updated = True
                    elif new_status == 'interrupted':
                        db.mark_pointing_interrupted(pointing_id, schedule_next=True)
                        self.log.info(f'Marked Pointing {pointing_id} as interrupted')
                        database_updated = True
                except Exception:
                    self.log.error('Error marking Pointing: {}'.format(pointing))
                    self.log.debug('', exc_info=True)

            # Create new Pointings for any Targets that are still unscheduled (e.g. expired)
            unscheduled_targets = db.get_targets(session, status='unscheduled', time=check_time)
            if len(unscheduled_targets) > 0:
                self.log.info('Rescheduling {} unscheduled Targets: {}'.format(
                              len(unscheduled_targets), [t.db_id for t in unscheduled_targets]))
                pointings = [target.get_next_pointing() for target in unscheduled_targets]
                pointings = [p for p in pointings if p is not None]
                db.insert_items(session, pointings)
                database_updated = True

        return database_updated

    def _get_pointings(self, site_id, check_time):
        """Calculate what to observe for telescopes at the given site."""
        try:
            pointings = {}
            telescopes = self.site_data[site_id]['telescopes']

            if (params.SCHEDULER_SKIP_DAYTIME and
                    not all(self.is_night[telescope_id] for telescope_id in sorted(telescopes))):
                # It's still daytime, no reason to check the queue
                for telescope_id in sorted(telescopes):
                    self.log.debug(f'Telescope {telescope_id}: Daytime')
                return {telescope_id: [None for _ in telescopes[telescope_id]['horizon']]
                        for telescope_id in telescopes}

            # Import the queue for this site from the database
            queue = PointingQueue.from_database(site_id, check_time)

            # Loop through each telescope
            # TODO: write queue files / HTML page for each
            for telescope_id in sorted(telescopes):
                telescope_pointings = []

                # Start with the first horizon, find what to do
                pointing, reason = queue.get_pointing(
                    telescope_id,
                    horizon=0,
                    readout_time=self.readout_time,
                    template_requirement=self.template_requirement,
                    return_reason=True,
                )
                telescope_pointings.append(copy(pointing))
                self.log.debug(f'Telescope {telescope_id}: {reason}')

                # Now check if it's above any other horizons, and if not get the next best Pointing
                # TODO: A better solution might be for get_pointing() to take
                #       multiple horizon numbers and evaluate each pointing based on both,
                #       or automatically do it for every telescope horizon when evaluating...
                for i in range(len(telescopes[telescope_id]['horizon']) - 1):
                    if pointing is None:
                        telescope_pointings.append(None)
                        continue
                    horizon = telescopes[telescope_id]['horizon'][i + 1]
                    if pointing.altaz.alt < horizon_limit(pointing.altaz.az, horizon):
                        # The Pointing is below this (presumably higher) horizon.
                        # This should be fairly rare, if it's just for wind shielding.
                        # Need to find the next best Pointing, recalculating with the new horizon.
                        pointing, reason = queue.get_pointing(
                            telescope_id,
                            horizon=i + 1,
                            readout_time=self.readout_time,
                            template_requirement=self.template_requirement,
                            return_reason=True,
                        )
                        telescope_pointings.append(copy(pointing))
                        self.log.debug(f'Telescope {telescope_id}-{i + 1}: {reason}')
                    else:
                        # This Pointing is fine
                        telescope_pointings.append(copy(pointing))

                # Add to site dict
                pointings[telescope_id] = telescope_pointings

            return pointings
        except Exception:
            self.log.error('Failed to schedule Pointings for Site {}'.format(site_id))
            self.log.debug('', exc_info=True)
            return {telescope_id: [None for _ in telescopes[telescope_id]['horizon']]
                    for telescope_id in telescopes}

    # Functions
    def update_schedule(self, telescope_id, current_pointing_id, current_status, horizon=0,
                        return_new=True, force_update=False):
        """Get what Pointing to observe for the given telescope.

        This function should be called by the pilot for each telescope. By including the
        current pointing and status the scheduler will update the database based on what the
        pilot is observing.

        Parameters
        ----------
        telescope_id : int
            The ID number of the Telescope requesting a Pointing.
        current_pointing_id : int or None
            The ID number of the current Pointing the telescope is observing, if any.
        current_status : str or None
            The status of the current Pointing to update the observing database.
            Valid statuses are 'running', 'completed' or 'interrupted',
            or None if current_pointing is None.

        horizon : int, optional
            Which horizon number to consider when scheduling, based on the horizons for the
            telescope defined in the database.
            default = 0
        return_new : bool, optional
            If False, only update the current Pointing and don't return a new one.
            This should be used when the telescope is shutting down, so we don't mark the final
            Pointing as running.
            default = True
        force_update : bool, optional
            If True force the scheduler to recalculate at the current time.
            Otherwise the pointing from the most recent check will be returned (~5s cadence).
            default = False

        """
        # TODO: Require some sort of unique API key from each telescope for security?

        # Prevent multiple queries at the same time, or queries while the schedule is updating
        while (self.query_lock is True or self.update_lock is True):
            time.sleep(0.1)
        self.query_lock = True

        try:
            msg = f'Query from Telescope {telescope_id}-{horizon}'
            msg += f' (CP={current_pointing_id}, {current_status})'
            self.log.info(msg)
            self.pilot_query_time[telescope_id] = Time.now()

            if current_pointing_id is not None:
                if current_status == 'completed':
                    # The Telescope has successfully finished this Pointing.
                    self.update_queue.append((current_pointing_id, 'completed'))
                    force_update = True
                elif current_status == 'interrupted':
                    # This Pointing was interrupted before it could finish.
                    self.update_queue.append((current_pointing_id, 'interrupted'))
                    force_update = True

            # If we don't want a new Pointing then we can return here.
            if not return_new:
                self.log.info('Returning (no Pointing requested)')
                self.query_lock = False
                return None

            # Wait for scheduler to update.
            # This is particularly important if we just marked a pointing, as during the 1s wait
            # a new check might have started which won't yet take that into account.
            # So we have to wait for the update lock to clear, then set the force check flag and
            # wait for that check to happen to mark the Pointings and update the queue.
            while self.update_lock is True:
                time.sleep(0.1)
            if force_update:
                self.force_check_flag = True
                while self.force_check_flag:
                    time.sleep(0.1)

            # Get the latest Pointing from the scheduler
            new_pointing = self.latest_pointings[telescope_id][horizon]

            # If it's none then return here
            if new_pointing is None:
                self.log.info('Returning None')
                self.query_lock = False
                return None

            # We have a new Pointing, but it might be already running.
            if current_pointing_id is not None and new_pointing.db_id == current_pointing_id:
                self.log.info(f'Returning Pointing {new_pointing.db_id}')
            else:
                self.log.info(f'Returning Pointing {new_pointing.db_id} (NEW)')

                if current_pointing_id is not None and current_status == 'running':
                    # It's interrupting a running Pointing, so mark it as interrupted.
                    self.update_queue.append((current_pointing_id, 'interrupted'))
                    self.force_check_flag = True

            # Mark the new Pointing as running, if it isn't already
            if new_pointing.current_telescope is None:
                self.update_queue.append((new_pointing.db_id, 'running', telescope_id))
                self.force_check_flag = True

            # Get the pointing info, and add anything else useful
            pointing_info = db.get_pointing_info(new_pointing.db_id)
            pointing_info['obstime'] = new_pointing.get_obstime(self.readout_time)

        finally:
            self.query_lock = False

        return pointing_info

    def check_queue(self, force_update=False):
        """Return information on the Pointings currently in the queue for all Telescopes.

        Note this returns immediately with what's been previously calculated by the scheduler,
        unless the `force_update` flag is set to True.

        Parameters
        ----------
        force_update : bool, optional
            If True force the scheduler to recalculate at the current time.
            Otherwise the pointing from the most recent check will be returned (~5s cadence).
            default = False

        """
        if force_update or self.force_check_flag:  # we might already be waiting for an update
            self.force_check_flag = True
            while self.force_check_flag:
                time.sleep(0.1)

        pointings = {telescope_id: [db.get_pointing_info(pointing.db_id)
                                    if pointing is not None else None
                                    for pointing in pointings]
                     for telescope_id, pointings in self.latest_pointings.items()}
        return pointings

    def get_pointing_info(self, pointing_id):
        """Get info from the database for a given pointing.

        This can be used by the remote clients which can't otherwise access the database.

        """
        return db.get_pointing_info(pointing_id)


def run():
    """Start the scheduler."""
    try:
        send_slack_msg('Scheduler started')
        scheduler = Scheduler()
        scheduler.run(params.PYRO_HOST, params.PYRO_PORT, params.PYRO_TIMEOUT)
    except Exception:
        print('Error detected, shutting down')
        traceback.print_exc()
    except KeyboardInterrupt:
        print('Interrupt detected, shutting down')
    finally:
        try:
            scheduler.shutdown()
        except UnboundLocalError:
            # class was never created
            pass
        time.sleep(1)  # wait to stop threads
        send_slack_msg('Scheduler shutdown')
        print('Scheduler done')
