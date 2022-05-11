"""Class for monitoring the obs database and finding pointings to observe."""

import threading
import time
import traceback

import Pyro4

from astropy.time import Time

from gtecs.common.logging import get_logger

from . import database as db
from . import params
from .astronomy import above_horizon
from .scheduling import check_queue
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
        self.write_html = params.WRITE_QUEUE_PAGE

        # Get telescope data (only once, on init)
        self.tel_data = db.get_telescope_info()
        self.latest_pointings = {tel_id: [None] for tel_id in self.tel_data}

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

            # Only check on the given period, or if forced
            if self.force_check_flag or (self.loop_time - self.check_time) > self.check_period:
                self.check_time = self.loop_time
                self.force_check_flag = False

                # Caretaker step to deal with any unscheduled Targets with expired Pointings
                self._caretaker()

                # Get the queue and find highest priority Pointings for each telescope
                check_time = Time(self.check_time, format='unix')
                new_pointings = self._get_pointings(check_time)

                # TODO: send database reports to Slack every day?

                # Write debug log line (TODO: improve, only log if something changed)
                try:
                    self.log.debug(new_pointings)
                except Exception:
                    self.log.error('Could not write current status')

                self.latest_pointings = new_pointings

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

    def _get_pointings(self, check_time):
        """Calculate what to observe for each telescope in the database."""
        tel_pointings = {tel_id: [None] for tel_id in self.tel_data}
        for telescope_id in self.tel_data:
            try:
                pointings = []

                # Start with the first horizon, find the highest priority Pointing
                pointing, reason = check_queue(telescope_id,
                                               time=check_time,
                                               horizon=self.tel_data[telescope_id]['horizon'][0],
                                               readout_time=self.readout_time,
                                               template_requirement=self.template_requirement,
                                               write_file=self.write_file,
                                               write_html=self.write_html,
                                               return_reason=True,
                                               )
                pointings.append(pointing)
                self.log.debug(f'Telescope {telescope_id}: {reason}')

                # Now check if it's above each horizon, and if not find the next best Pointing
                # TODO: AltAz is stored on the pointing, there should be a quicker way to check
                for i, horizon in enumerate(self.tel_data[telescope_id]['horizon'][1:]):
                    if pointing is not None and not above_horizon(
                            pointing.ra, pointing.dec,
                            self.tel_data[telescope_id]['location'],
                            check_time,
                            horizon):
                        # The Pointing is below this (presumably higher) horizon.
                        # This should be fairly rare, if it's just for wind shielding.
                        # There might be a more efficient way to do this,
                        # here we just recalculate using the higher horizon.
                        # A better solution might be for check_queue() to take multiple
                        # horizons and evaluate each pointing based on both, but that's
                        # probably a waste of time (even if altaz is cached).
                        pointing, reason = check_queue(telescope_id,
                            time=check_time,
                            horizon=horizon,
                            readout_time=self.readout_time,
                            template_requirement=self.template_requirement,
                            write_file=self.write_file,
                            write_html=self.write_html,
                            return_reason=True,
                            )
                        pointings.append(pointing)
                        self.log.debug(f'Telescope {telescope_id}-{i+1}: {reason}')
                    else:
                        # This Pointing is fine
                        pointings.append(pointing)

                tel_pointings[telescope_id] = pointings
            except Exception:
                self.log.error('Failed to schedule Pointing for Telescope {}'.format(telescope_id))
                self.log.debug('', exc_info=True)
                tel_pointings[telescope_id] = [None]
        return tel_pointings

    # Functions
    def check_queue(self, telescope_id, horizon=0, force_update=False):
        """Return the ID of the current highest priority pointing for the given telescope.

        Note this returns immediately with what's been previously calculated by the scheduler,
        unless the `force_update` flag is set to True.

        Parameters
        ----------
        telescope_id : int
            The ID number for the Telescope to query the queue for.

        horizon : int, optional
            Which horizon number to apply, based on horizons defined in the database.
            default = 0
        force_update : bool, optional
            If True force the scheduler to recalculate at the current time.
            Otherwise the pointing from the most recent check will be returned (~5s cadence).
            default = False

        """
        if force_update or self.force_check_flag:  # we might be waiting for a forced update
            self.force_check_flag = True
            while self.check_time < self.loop_time:
                time.sleep(0.01)
        self.log.info('CHECKING QUEUE FOR TELESCOPE {}'.format(telescope_id))
        pointing = self.latest_pointings[telescope_id][horizon]
        if pointing is None:
            return None
        else:
            pointing_info = db.get_pointing_info(pointing.db_id)
            # Add any useful scheduling info
            pointing_info['obstime'] = pointing.get_obstime(self.readout_time)
            return pointing_info

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
        self.force_check_flag = True

    def mark_pointing_completed(self, pointing_id):
        """Mark the given pointing as completed.

        This can be used by the remote clients which can't otherwise access the database.

        """
        db.mark_pointing_completed(pointing_id, schedule_next=True)
        self.force_check_flag = True

    def mark_pointing_interrupted(self, pointing_id):
        """Mark the given pointing as interrupted.

        This can be used by the remote clients which can't otherwise access the database.

        """
        db.mark_pointing_interrupted(pointing_id, schedule_next=True, delay=None)
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
