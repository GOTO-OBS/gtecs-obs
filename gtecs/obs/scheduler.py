"""Class for monitoring the observation database and finding pointings to observe."""

import json
import sys
import threading
import time
import traceback
from copy import copy

import Pyro4
from astropy import units as u
from astropy.coordinates import AltAz, get_sun
from astropy.time import Time
from flask import Flask, abort, request
from gtecs.common import logging
from gtecs.common.daemon import BaseDaemon
from gtecs.common.pyro import get_uri, wait_for_uri
from gtecs.common.system import execute_long_command, get_pid, kill_process, make_pid_file
from gtecs.common.systemd import (
    SystemctlError,
    get_service,
    get_service_pid,
    kill_service,
    restart_service,
    start_service,
    stop_service,
)
from webtest.http import StopableWSGIServer

from . import database as db
from . import params
from .astronomy import horizon_limit
from .scheduling import PointingQueue
from .slack import send_slack_msg


@Pyro4.expose
class SchedulerDaemon(BaseDaemon):
    """Database scheduler daemon class."""

    name = 'scheduler'

    def __init__(self):
        """Initialize the scheduler daemon."""
        super().__init__()

        # Scheduler properties
        self.loop_time = 0
        self.check_time = 0
        self.check_period = params.SCHEDULER_CHECK_PERIOD
        self.force_check_flag = False
        self.update_lock = False
        self.update_queue = []
        self.query_lock = False

        self.readout_time = params.READOUT_TIME
        self.template_requirement = params.TEMPLATE_REQUIREMENT
        # # These aren't implemented at the moment...
        # self.write_file = params.WRITE_QUEUE_FILE
        # self.queue_file = os.path.join(params.FILE_PATH, 'queue_info')
        # self.write_html = params.WRITE_QUEUE_PAGE
        # self.html_path = params.HTML_PATH

        # Get site data (only once, on init)
        self.site_data = db.get_site_info()
        self.telescopes = sorted(
            [
                telescope_id
                for telescopes in [
                    self.site_data[site_id]['telescopes'].keys() for site_id in self.site_data
                ]
                for telescope_id in telescopes
            ]
        )

        self.old_pointings = {telescope_id: [None] for telescope_id in self.telescopes}
        self.latest_pointings = {telescope_id: [None] for telescope_id in self.telescopes}

        self.is_night = {telescope_id: False for telescope_id in self.telescopes}
        self.pilot_query_time = {telescope_id: None for telescope_id in self.telescopes}
        self.database_update_time = 0

        # Web server params
        self.server_host = params.SERVER_HOST
        self.server_port = params.SERVER_PORT
        self.server_path = params.SERVER_PATH
        self.server_api_key = params.SERVER_API_KEY

    def run(self, host, port, timeout=5):
        """Run the scheduler as a Pyro daemon."""
        # Start monitor thread
        t = threading.Thread(target=self._monitor_thread)
        t.daemon = True
        t.start()

        # Start server thread
        t2 = threading.Thread(target=self._server_thread)
        t2.daemon = True
        t2.start()

        # Start the daemon, will block until shutdown
        self.run_daemon(host, port, timeout=timeout)

        # If we reached here the loop has exited
        # Sleep for a second to make sure the threads have time to close
        time.sleep(1.0)

    # Internal threads
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
                try:
                    self._pilot_heartbeat(check_time)
                except Exception:
                    self.log.error('Error checking pilot heartbeats')
                    self.log.debug('', exc_info=True)

                # Caretaker step to monitor and update the database
                try:
                    database_updated = self._caretaker(check_time)
                except Exception:
                    database_updated = False
                    self.log.error('Error checking database')
                    self.log.debug('', exc_info=True)
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

    def _server_thread(self):
        """Access the scheduler through a Flask server."""
        self.log.info('Web server thread started')

        # Create the Flask application
        app = Flask('scheduler')

        @app.route(self.server_path)
        def homepage():
            return 'This is the homepage for the GOTO central scheduler'

        @app.route(self.server_path + '/update_schedule/<int:telescope_id>', methods=['GET'])
        def update_schedule(telescope_id):
            try:
                # Check API key
                if 'api_key' not in request.args:
                    return '<h1>400 Bad Request</h1> No API key provided', 400
                if request.args['api_key'] != self.server_api_key:
                    return '<h1>403 Forbidden</h1> Invalid API key', 403

                # Get function args
                current_pointing_id = request.args['current_pointing_id']
                if current_pointing_id == 'None':
                    current_pointing_id = None
                else:
                    current_pointing_id = int(current_pointing_id)
                current_status = request.args['current_status']
                if current_status == 'None':
                    current_status = None
                if 'horizon' in request.args:
                    horizon = int(request.args['horizon'])
                else:
                    horizon = 0
                if 'return_new' in request.args:
                    return_new = bool(int(request.args['return_new']))
                else:
                    return_new = True
                if 'force_update' in request.args:
                    force_update = bool(int(request.args['force_update']))
                else:
                    force_update = False

                # Call the internal function
                pointing_info = self.update_schedule(
                    telescope_id,
                    current_pointing_id,
                    current_status,
                    horizon=horizon,
                    return_new=return_new,
                    force_update=force_update,
                )

                # Return the info as a JSON dict
                response = app.response_class(
                    response=json.dumps(pointing_info), mimetype='application/json'
                )
                return response
            except Exception:
                self.log.exception(f'HTTP request failed (from {request.path})')
                abort(500)

        @app.route(self.server_path + '/check_queue', methods=['GET'])
        def check_queue():
            try:
                # Check API key
                if 'api_key' not in request.args:
                    return '<h1>400 Bad Request</h1> No API key provided', 400
                if request.args['api_key'] != self.server_api_key:
                    return '<h1>403 Forbidden</h1> Invalid API key', 403

                # Get function args
                if 'force_update' in request.args:
                    force_update = bool(int(request.args['force_update']))
                else:
                    force_update = False

                # Call the internal function
                pointing_info = self.check_queue(force_update)

                # Return the info as a JSON dict
                response = app.response_class(
                    response=json.dumps(pointing_info), mimetype='application/json'
                )
                return response
            except Exception:
                self.log.exception(f'HTTP request failed (from {request.path})')
                abort(500)

        @app.route(self.server_path + '/pointing_info/<int:pointing_id>', methods=['GET'])
        def get_pointing_info(pointing_id):
            try:
                # Check API key
                if 'api_key' not in request.args:
                    return '<h1>400 Bad Request</h1> No API key provided', 400
                if request.args['api_key'] != self.server_api_key:
                    return '<h1>403 Forbidden</h1> Invalid API key', 403

                # Call the internal function
                pointing_info = self.get_pointing_info(pointing_id)

                # Return the info as a JSON dict
                response = app.response_class(
                    response=json.dumps(pointing_info), mimetype='application/json'
                )
                return response
            except Exception:
                self.log.exception(f'HTTP request failed (from {request.path})')
                abort(500)

        # Create and start the server
        server = StopableWSGIServer.create(app, host=self.server_host, port=self.server_port)
        server.wait()
        addr = f'http://{self.server_host}:{self.server_port}{self.server_path}/'
        self.log.info(f'Server started at {addr}')
        while self.running:
            time.sleep(1)

        # Stop the server nicely when shutting down
        self.log.info('Stopping server')
        server.shutdown()

        self.log.info('Web server thread stopped')
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
        with db.session_manager() as session:
            for telescope_id in self.telescopes:
                # Check how long it's been since we had a query from this telescope.
                if self.pilot_query_time[telescope_id] is None:
                    continue
                if (check_time - self.pilot_query_time[telescope_id]).sec > 60:
                    # It's been more than a minute since we last had a query from this pilot.
                    if self.is_night[telescope_id]:
                        # Only log during the night time.
                        self.log.debug(
                            'No contact from Telescope {} since {}'.format(
                                telescope_id, self.pilot_query_time[telescope_id].iso
                            )
                        )

                    # If there is still a running pointing, mark it as interrupted.
                    telescope = db.get_telescope_by_id(session, telescope_id)
                    if telescope.current_pointing is not None:
                        pointing = telescope.current_pointing
                        running_time = (check_time - Time(pointing.running_time)).sec
                        obs_time = pointing.get_obstime(self.readout_time)
                        self.log.info(
                            'Pointing {} has been running for {:.1f}s/{:.1f}s'.format(
                                pointing.db_id, running_time, obs_time
                            )
                        )
                        self.update_queue.append((pointing.db_id, 'interrupted'))
                        # No force check needed here, we're already in one!

                # TODO: Send Slack alerts? Could be when a telescope first connects in the
                #       evening, and then if the connection is lost before sunrise.

    def _caretaker(self, check_time):
        """Monitor the database and update any Pointings."""
        database_updated = False
        with db.session_manager() as session:
            # Process any Pointings that need updating
            while len(self.update_queue) > 0:
                pointing = self.update_queue.pop(0)  # First-in, first-out
                try:
                    pointing_id = pointing[0]
                    new_status = pointing[1]
                    if new_status == 'running':
                        telescope_id = pointing[2]
                        db.mark_pointing_running(pointing_id, telescope_id)
                        self.log.info(
                            f'Marked Pointing {pointing_id} as running on Telescope {telescope_id}'
                        )
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
                self.log.info(
                    'Rescheduling {} unscheduled Targets: {}'.format(
                        len(unscheduled_targets), [t.db_id for t in unscheduled_targets]
                    )
                )
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

            if params.SCHEDULER_SKIP_DAYTIME and not all(
                self.is_night[telescope_id] for telescope_id in sorted(telescopes)
            ):
                # It's still daytime, no reason to check the queue
                for telescope_id in sorted(telescopes):
                    self.log.debug(f'Telescope {telescope_id}: Daytime')
                return {
                    telescope_id: [None for _ in telescopes[telescope_id]['horizon']]
                    for telescope_id in telescopes
                }

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
            return {
                telescope_id: [None for _ in telescopes[telescope_id]['horizon']]
                for telescope_id in telescopes
            }

    # Functions
    def update_schedule(
        self,
        telescope_id,
        current_pointing_id,
        current_status,
        horizon=0,
        return_new=True,
        force_update=False,
    ):
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
            The status of the current Pointing to update the observation database.
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
        while self.query_lock is True or self.update_lock is True:
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

        pointing_info = {
            telescope_id: [
                db.get_pointing_info(pointing.db_id) if pointing is not None else None
                for pointing in pointings
            ]
            for telescope_id, pointings in self.latest_pointings.items()
        }
        return pointing_info

    def get_pointing_info(self, pointing_id):
        """Get info from the database for a given pointing.

        This can be used by the remote clients which can't otherwise access the database.

        """
        return db.get_pointing_info(pointing_id)


def query(command, args):
    """Process a query."""
    # Primary control commands
    if command == 'run':
        service = get_service('scheduler')

        # Check the systemd service is already running
        if service is not None and service['active'] and service['pid'] != os.getpid():
            print(f'The scheduler service is already running (PID={service["pid"]})')
            print('Use "scheduler stop" to stop it before running the scheduler locally')
            return

        # Check if the Pyro URI is already registered
        uri = get_uri(params.PYRO_URI)
        if uri is not None:
            pid = get_service_pid('scheduler', 'scheduler')
            print(f'Scheduler is already running ({uri}, PID={pid})')
            return

        exit_code = 0
        with make_pid_file('scheduler'):
            try:
                send_slack_msg('Scheduler started')
                daemon = SchedulerDaemon()
                daemon.run(params.PYRO_HOST, params.PYRO_PORT, params.PYRO_TIMEOUT)
            except Exception:
                print('Error detected, shutting down')
                traceback.print_exc()
                exit_code = 1
            except (KeyboardInterrupt, SystemExit):
                # SIGINT raises a KeyboardInterrupt, and make_pid_file makes SIGTERM raise
                # a SystemExit (which is how systemd stops the service),
                # so both are normal ways to stop.
                print('Interrupt detected, shutting down')
            finally:
                try:
                    daemon.shutdown()
                except UnboundLocalError:
                    # class was never created
                    pass
                time.sleep(1)  # wait to stop threads
                send_slack_msg('Scheduler shutdown')
                print('Scheduler done')

            sys.exit(exit_code)
        return

    if command == 'start':
        service = get_service('scheduler')

        # Check if the systemd service is set up
        if service is None:
            print('ERROR: Cannot connect to the systemd service, or it is not installed')
            print('Use "scheduler run" to run the scheduler in this terminal instead')
            return

        # Check if the Pyro URI is already registered
        uri = get_uri(params.PYRO_URI)
        if uri is not None:
            pid = get_service_pid('scheduler', 'scheduler')
            print(f'Scheduler is already running ({uri}, PID={pid})')
            return

        start_service('scheduler', sudo=params.SYSTEMD_SUDO)

        # Check if the service is now running
        uri = wait_for_uri(params.PYRO_URI)
        if uri is None:
            print('ERROR: The scheduler service started, but it is not responding')
            print('Check the service with "systemctl status scheduler"')
        else:
            pid = get_service_pid('scheduler', 'scheduler')
            print(f'Scheduler started ({uri}, PID={pid})')
        return

    if command in ['shutdown', 'stop']:
        service = get_service('scheduler')

        # Check if the systemd service is set up
        if service is None:
            print('ERROR: Cannot connect to the systemd service, or it is not installed')
            return

        if service['active']:
            # The service is running through systemd, so shut it down through systemd
            stop_service('scheduler', sudo=params.SYSTEMD_SUDO)
            print('Scheduler shutdown')
            return

        # It might might be running independently in a terminal,
        # so try shutting down though Pyro
        try:
            with Pyro4.Proxy(params.PYRO_URI) as proxy:
                proxy.shutdown()
                print('Scheduler shutdown')
        except Pyro4.errors.CommunicationError:
            print('Scheduler is not running')
        return

    if command == 'restart':
        service = get_service('scheduler')

        # Check if the systemd service is set up
        if service is None:
            print('ERROR: Cannot connect to the systemd service, or it is not installed')
            return

        restart_service('scheduler', sudo=params.SYSTEMD_SUDO)

        # Check if the service is now running
        uri = wait_for_uri(params.PYRO_URI)
        if uri is None:
            print('ERROR: The scheduler service restarted, but it is not responding')
            print('Check the service with "systemctl status scheduler"')
        else:
            pid = get_service_pid('scheduler', 'scheduler')
            print(f'Scheduler restarted ({uri}, PID={pid})')
        return

    if command == 'kill':
        service = get_service('scheduler')
        if service is not None and service['active']:
            # The service is running through systemd, so we need to kill it through systemd
            kill_service('scheduler', sudo=params.SYSTEMD_SUDO)
            print('Scheduler killed (note systemd may restart it)')
            return

        if get_pid('scheduler'):
            # It might might be running independently in a terminal,
            # so try killing it directly
            kill_process('scheduler')
            return

        print('Scheduler is not running')
        return

    if command in ['ping', 'status']:
        uri = get_uri(params.PYRO_URI)
        if uri is not None:
            pid = get_service_pid('scheduler', 'scheduler')
            print(f'Scheduler is running ({uri}, PID={pid})')
            with Pyro4.Proxy(params.PYRO_URI) as proxy:
                start_time = proxy.start_time
            print(f'Uptime: {time.time() - start_time:.1f}s')
        else:
            print('Scheduler is not running')

        service = get_service('scheduler')
        if service is not None:
            print('Service state: {} ({})'.format(service['active_state'], service['sub_state']))
        return

    if command in ['log', 'tail']:
        log_file = 'scheduler.log'
        log_path = logging.get_log_path() / log_file
        tail_command = 'tail {} {}'.format(log_path, ' '.join(args))

        execute_long_command(tail_command)
        return

    # Scheduler commands
    if command == 'info':
        if len(args) == 0:
            force = False
        elif len(args) == 1 and args[0] in ['f', 'F', '-f', '-F', '--force-update']:
            force = True
        else:
            raise ValueError('Invalid arguments')

        # do one query for all telescopes
        with Pyro4.Proxy(params.PYRO_URI) as proxy:
            pointings = proxy.check_queue(force_update=force)

        telescope_data = db.get_telescope_info()
        for telescope_id in telescope_data:
            # Print Pointing (or None) for each telescope
            print(telescope_data[telescope_id]['name'], end=': ')
            try:
                if telescope_id not in pointings:
                    print('ERROR')
                elif pointings[telescope_id][0] is None:
                    print('None')
                else:
                    pointing = pointings[telescope_id][0]
                    print(f'Pointing ID={pointing["id"]} ({pointing["name"]})')
                # TODO: Print higher horizon if shielding?
                # TODO: Print what it's currently observing?
                # TODO: Verbose Pointing info?
            except Exception:
                print('ERROR')
        return

    # Help or unknown command - print instructions
    if command in ['help', '?']:
        print('The scheduler monitors the database and calculates which pointings to observe.\n')
    else:
        print(f'ERROR: Unrecognized command "{command}"\n')
    print('Usage: scheduler [command] [options]\n')
    print('Scheduler commands:')
    print('   info [-f]              print the latest scheduled pointings [force-update]')
    print('Control commands:')
    print('   start                  start the scheduler systemd service')
    print('   stop/shutdown          stop the scheduler systemd service')
    print('   restart                restart the scheduler systemd service')
    print('   kill                   kill the scheduler systemd service (systemd may restart it)')
    print('   run                    run the scheduler daemon locally in this terminal')
    print('   ping/status            check if the scheduler daemon or service is running')
    print('   log [tail args]        print log file contents (alias for tail)')
    print('   help                   print these instructions')


def main():
    """Main entry point for the script."""
    if len(sys.argv) > 1:
        command, *args = sys.argv[1:]
    else:
        # Interpret no arguments as a request for help
        command, args = 'help', []

    try:
        query(command, args)
    except (SystemctlError, ConnectionError) as err:
        print(f'ERROR: {err}')
        sys.exit(1)
