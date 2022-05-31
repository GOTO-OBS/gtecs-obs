"""Class for monitoring the obs database and finding pointings to observe."""

import os
import threading
import time
import traceback

import Pyro4

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
                # Caretaker step to deal with any unscheduled Targets with expired Pointings
                self._caretaker()

                # Get the queue and find highest priority Pointings for each telescope
                check_time = Time(self.loop_time, format='unix')
                msg = 'Checking queue'
                if self.force_check_flag:
                    msg += ' (forced update)'
                self.log.debug(msg)
                new_pointings = {telescope_id: [None] for telescope_id in self.telescopes}
                for site_id in self.site_data:
                    new_pointings.update(self._get_pointings(site_id, check_time))
                self.log.debug('Queue check complete')
                self.old_pointings = self.latest_pointings
                self.latest_pointings = new_pointings

                # Write status log line
                # TODO: account for horizons?
                try:
                    old_strs = ['{}:{}'.format(telescope_id,
                                               self.old_pointings[telescope_id][0].db_id
                                               if self.old_pointings[telescope_id][0] is not None
                                               else 'None')
                                for telescope_id in self.old_pointings]
                    old_str = ' '.join(old_strs)
                    new_strs = ['{}:{}'.format(telescope_id,
                                               self.latest_pointings[telescope_id][0].db_id
                                               if self.latest_pointings[telescope_id][0] is not None
                                               else 'None')
                                for telescope_id in self.latest_pointings]
                    new_str = ' '.join(new_strs)
                    if new_str != old_str:
                        self.log.info(new_str)
                except Exception:
                    self.log.error('Could not write current status')

                # Set these after the check has finished, so the check_queue() function has to wait
                self.check_time = self.loop_time
                self.force_check_flag = False

            time.sleep(0.1)

        self.log.info('Database monitor thread stopped')
        return

    # Internal functions
    def _caretaker(self):
        """Schedule new Pointings for any unscheduled Targets."""
        with db.open_session() as session:
            unscheduled_targets = db.get_targets(session, status='unscheduled')
            if len(unscheduled_targets) > 0:
                self.log.info('Rescheduling {} unscheduled Targets'.format(
                              len(unscheduled_targets)))
                pointings = [target.get_next_pointing() for target in unscheduled_targets]
                pointings = [p for p in pointings if p is not None]
                db.insert_items(session, pointings)

    def _get_pointings(self, site_id, check_time):
        """Calculate what to observe for telescopes at the given site."""
        try:
            pointings = {}
            location = self.site_data[site_id]['location']
            telescopes = self.site_data[site_id]['telescopes']

            if params.SCHEDULER_SKIP_DAYTIME:
                # Check if the sun is up, if so we can just return nothing to observe
                # TODO: If we ever want to display the queue on a webpage we'll still need
                #       to calculate every time, even during the day.
                altaz_frame = AltAz(obstime=check_time, location=location)
                sun_coords = get_sun(check_time)
                altaz_coords = sun_coords.transform_to(altaz_frame)
                sunalt = altaz_coords.alt.degree
                if sunalt > params.SCHEDULER_SUNALT_LIMIT:
                    # It's still daytime, no reason to check the queue
                    for telescope_id in sorted(telescopes):
                        self.log.debug(f'Telescope {telescope_id}: Daytime (sunalt={sunalt:.1f})')
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
                telescope_pointings.append(pointing)
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
                        telescope_pointings.append(pointing)
                        self.log.debug(f'Telescope {telescope_id}-{i + 1}: {reason}')
                    else:
                        # This Pointing is fine
                        telescope_pointings.append(pointing)

                # Add to site dict
                pointings[telescope_id] = telescope_pointings

            return pointings
        except Exception:
            self.log.error('Failed to schedule Pointings for Site {}'.format(site_id))
            self.log.debug('', exc_info=True)
            return {telescope_id: [None for _ in telescopes[telescope_id]['horizon']]
                    for telescope_id in telescopes}

    # Functions
    def check_queue(self, telescope_id, horizon=0, force_update=False):
        """Return the ID of the current highest priority pointing for the given telescope.

        Note this returns immediately with what's been previously calculated by the scheduler,
        unless the `force_update` flag is set to True.

        Parameters
        ----------
        telescope_id : int, or list of int
            The ID number(s) for the Telescope(s) to query the queue for.

        horizon : int, optional
            Which horizon number to apply, based on horizons defined in the database.
            default = 0
        force_update : bool, optional
            If True force the scheduler to recalculate at the current time.
            Otherwise the pointing from the most recent check will be returned (~5s cadence).
            default = False

        """
        if force_update or self.force_check_flag:  # we might already be waiting for an update
            self.force_check_flag = True
            while self.force_check_flag:
                time.sleep(0.1)

        telescope_ids = telescope_id
        if isinstance(telescope_id, int):
            telescope_ids = [telescope_id]

        pointings = []
        for telescope_id in telescope_ids:
            msg = f'Fetching pointing for Telescope {telescope_id}:'
            pointing = self.latest_pointings[telescope_id][horizon]
            if pointing is None:
                msg += ' returns None'
                self.log.debug(msg)
                pointings.append(None)
            else:
                pointing_info = db.get_pointing_info(pointing.db_id)
                # Add any useful scheduling info
                pointing_info['obstime'] = pointing.get_obstime(self.readout_time)
                msg += f' returns Pointing {pointing_info["id"]}'
                self.log.debug(msg)
                pointings.append(pointing_info)

        if len(pointings) == 1:
            return pointings[0]
        return pointings

    def get_pointing_info(self, pointing_id):
        """Get info from the database for a given pointing.

        This can be used by the remote clients which can't otherwise access the database.

        """
        return db.get_pointing_info(pointing_id)

    def mark_pointing_running(self, pointing_id, telescope_id):
        """Mark the given pointing as running on the given telescope.

        This can be used by the remote clients which can't otherwise access the database.

        """
        db.mark_pointing_running(pointing_id, telescope_id)
        self.log.debug(f'Marked Pointing {pointing_id} as running on Telescope {telescope_id}')
        self.force_check_flag = True

    def mark_pointing_completed(self, pointing_id):
        """Mark the given pointing as completed.

        This can be used by the remote clients which can't otherwise access the database.

        """
        db.mark_pointing_completed(pointing_id, schedule_next=True)
        self.log.debug(f'Marked Pointing {pointing_id} as completed')
        self.force_check_flag = True

    def mark_pointing_interrupted(self, pointing_id):
        """Mark the given pointing as interrupted.

        This can be used by the remote clients which can't otherwise access the database.

        """
        db.mark_pointing_interrupted(pointing_id, schedule_next=True, delay=None)
        self.log.debug(f'Marked Pointing {pointing_id} as interrupted')
        self.force_check_flag = True


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
