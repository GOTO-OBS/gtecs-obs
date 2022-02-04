"""Robotic queue scheduler functions."""

import json
import logging
import os
import warnings

from astroplan import (AltitudeConstraint, AtNightConstraint,
                       Constraint, MoonIlluminationConstraint, MoonSeparationConstraint,
                       Observer, TimeConstraint)
from astroplan.constraints import _get_altaz

from astropy import units as u
from astropy.coordinates import EarthLocation, SkyCoord
from astropy.time import Time

from erfa import ErfaWarning

from gtecs.obs import database as db

import numpy as np

from scipy import interpolate

from . import params
from .astronomy import time_to_set
from .html import write_queue_page

warnings.simplefilter('ignore', ErfaWarning)


MOON_PHASES = {'B': 1, 'G': 0.65, 'D': 0.25}


def apply_constraints(constraints, observer, targets, times):
    """Determine if the targets are observable for given times, observer, and constraints.

    Similar to Astroplan's ``is_observable``, but returns full constraint
    results in a 3D array.

    Parameters
    ----------
    constraints : list or `~astroplan.constraints.Constraint`
        Observational constraint(s)

    observer : `~astroplan.Observer`
        The observer who has constraints ``constraints``

    targets : {list, `~astropy.coordinates.SkyCoord`, `~astroplan.FixedTarget`}
        Target or list of targets

    times : `~astropy.time.Time`
        Array of times on which to test the constraint

    Returns
    -------
    constraint_array : 3D array
        first dimension is given by the number of times
        second dimension is given by the number of targets
        third dimension is given by the number of constraints

    """
    if not hasattr(constraints, '__len__'):
        constraints = [constraints]

    applied_constraints = np.array([constraint(observer, targets, times)
                                    for constraint in constraints])
    return np.transpose(applied_constraints)


class ArtificialHorizonConstraint(Constraint):
    """Ensure altitude is above artificial horizon."""

    def __init__(self, az, alt):
        self.get_alt_limit = interpolate.interp1d(az, alt,
                                                  bounds_error=False,
                                                  fill_value='extrapolate')

    def compute_constraint(self, times, observer, targets):
        """Compute the constraint."""
        altaz = _get_altaz(times, observer, targets)['altaz']
        alt_limit = self.get_alt_limit(altaz.az) * u.deg
        return altaz.alt > alt_limit


class Pointing:
    """A class to contain infomation on each pointing."""

    def __init__(self, db_id, name, ra, dec, rank, weight, num_obs, too,
                 mintime, maxsunalt, minalt, maxmoon, minmoonsep,
                 start, stop, current, exp_sets):
        self.db_id = int(db_id)
        self.name = name
        self.ra = float(ra)
        self.dec = float(dec)
        self.rank = int(rank)
        self.weight = float(weight)
        self.num_obs = int(num_obs)
        self.too = bool(int(too))
        self.mintime = mintime
        self.maxsunalt = float(maxsunalt)
        self.minalt = float(minalt)
        self.maxmoon = maxmoon
        self.minmoonsep = minmoonsep
        self.start = start
        self.stop = stop
        self.current = bool(current)
        self.exp_sets = exp_sets

    def __eq__(self, other):
        try:
            return self.db_id == other.db_id
        except AttributeError:
            return False

    def __ne__(self, other):
        return self != other

    def __repr__(self):
        template = ('Pointing(db_id={}, name={}, ra={}, dec={}, rank={}, weight={}, ' +
                    'num_obs={}, too={}, maxsunalt={}, minalt={}, mintime={}, maxmoon={}, ' +
                    'minmoonsep={}, start={}, stop={}, ' +
                    'current={})')
        return template.format(
            self.db_id, self.name, self.ra, self.dec, self.rank, self.weight,
            self.num_obs, self.too, self.maxsunalt, self.minalt, self.mintime,
            self.maxmoon, self.minmoonsep, self.start, self.stop,
            self.current)

    @classmethod
    def from_database(cls, db_pointing):
        """Import a pointing from the database."""
        # Get DB Pointing properties
        db_id = db_pointing.db_id
        rank = db_pointing.rank
        start_time = db_pointing.start_time
        stop_time = db_pointing.stop_time
        current = db_pointing.status == 'running'

        # Get linked Target properties
        db_target = db_pointing.target
        name = db_target.name
        ra = db_target.ra
        dec = db_target.dec
        weight = db_target.weight
        num_completed = db_target.num_completed
        too = db_target.too
        min_time = db_target.min_time
        max_sunalt = db_target.max_sunalt
        min_alt = db_target.min_alt
        max_moon = db_target.max_moon
        min_moonsep = db_target.min_moonsep

        # Get linked ExposureSet properties
        exp_sets = [(es.num_exp, es.exptime) for es in db_pointing.exposure_sets]

        # Create Pointing object
        pointing = cls(db_id=db_id,
                       name=name,
                       ra=ra,
                       dec=dec,
                       rank=rank,
                       weight=weight,
                       num_obs=num_completed,
                       too=too,
                       mintime=min_time,
                       maxsunalt=max_sunalt,
                       minalt=min_alt,
                       maxmoon=max_moon,
                       minmoonsep=min_moonsep,
                       start=start_time,
                       stop=stop_time,
                       current=current,
                       exp_sets=exp_sets,
                       )
        return pointing

    def get_obstime(self, readout_time=15):
        """Get the expected time needed to observe this Pointing."""
        return sum(es[0] * (es[1] + readout_time) for es in self.exp_sets)


class PointingQueue:
    """A class to represent a queue of pointings."""

    def __init__(self, pointings=None):
        if pointings is None:
            pointings = []
        self.pointings = pointings

        # Can't do much with no Pointings
        if len(self.pointings) == 0:
            return

        # Create pointing data arrays
        self.targets = SkyCoord([float(p.ra) for p in self.pointings],
                                [float(p.dec) for p in self.pointings],
                                unit=u.deg, frame='icrs')

    def __len__(self):
        return len(self.pointings)

    def __repr__(self):
        template = ('PointingQueue(length={})')
        return template.format(len(self.pointings))

    @classmethod
    def from_database(cls, time, observer):
        """Create a Pointing Queue from current Pointings in the database.

        Parameters
        ----------
        time : `~astropy.time.Time`
            The time to fetch the queue for.
        observer : `astroplan.Observer`
            The observer of the pointings

        Returns
        -------
        pointings : list of `Pointing`
            An list containing Pointings from the database.

        """
        with db.open_session() as session:
            current, pending = db.get_filtered_queue(session,
                                                     time=time,
                                                     location=observer.location,
                                                     altitude_limit=params.HARD_ALT_LIM,
                                                     hourangle_limit=params.HARD_HA_LIM)

            pointings = [Pointing.from_database(db_pointing) for db_pointing in pending]
            if current:
                pointings.append(Pointing.from_database(current))
            queue = cls(np.array(pointings))
        return queue

    def get_current_pointing(self):
        """Return the current pointing from the queue."""
        for p in self.pointings:
            if p.current:
                return p
        return None

    def calculate_times(self, time, readout_time):
        """Calculate the finish times for all pointings."""
        # Note we use the mintime if a Pointing has it, we don't care if it's invalid after then
        obstime_arr = [p.get_obstime(readout_time)
                       if p.mintime is None else p.mintime
                       for p in self.pointings]
        self.finish_times = time + u.Quantity(obstime_arr, unit=u.s)

    def apply_constraints(self, time, observer, horizon):
        """Check if the pointings are valid, both at start and end of expected observing period."""
        # Create Constraints
        self.constraints = {}

        # SunAlt
        maxsunalts = u.Quantity([float(p.maxsunalt) for p in self.pointings], unit=u.deg)
        self.constraints['SunAlt'] = AtNightConstraint(maxsunalts)

        # MinAlt
        minalts = u.Quantity([float(p.minalt) for p in self.pointings], unit=u.deg)
        self.constraints['MinAlt'] = AltitudeConstraint(minalts, None)

        # ArtHoriz
        self.constraints['ArtHoriz'] = ArtificialHorizonConstraint(*horizon)

        # Moon
        moonphases = [MOON_PHASES[p.maxmoon] for p in self.pointings]
        self.constraints['Moon'] = MoonIlluminationConstraint(None, moonphases)

        # MoonSep
        minmoonseps = u.Quantity([float(p.minmoonsep) for p in self.pointings], unit=u.deg)
        self.constraints['MoonSep'] = MoonSeparationConstraint(minmoonseps, None)

        # Time
        starts = Time([p.start for p in self.pointings], scale='utc', format='datetime')
        # NB the stop time can be None for non-expiring pointings,
        #    but the constraint needs a Time so just say a long time from now.
        longtime = Time.now() + 10 * u.year
        stops = Time([p.stop if p.stop else longtime.datetime for p in self.pointings],
                     scale='utc', format='datetime')
        self.constraints['Time'] = TimeConstraint(starts, stops)

        # Apply constraints at the start time (now)
        start_cons = ['SunAlt', 'MinAlt', 'ArtHoriz', 'Moon', 'MoonSep']
        start_cons_valid_arr = apply_constraints([self.constraints[name] for name in start_cons],
                                                 observer, self.targets, time)

        # Apply constraints at the finish time (now + expected observing time)
        end_cons = ['SunAlt', 'MinAlt', 'ArtHoriz']
        end_cons_valid_arr = apply_constraints([self.constraints[name] for name in end_cons],
                                               observer, self.targets, self.finish_times)

        # Apply time constraint
        time_cons = ['Time']
        time_cons_valid_arr = apply_constraints([self.constraints[name] for name in time_cons],
                                                observer, self.targets, time)

        # Save constraint results on Pointings and calculate if they are valid
        self.all_constraint_names = start_cons + time_cons + [name + '_end' for name in end_cons]
        for i, pointing in enumerate(self.pointings):
            # start constraints (applies to everything)
            pointing.constraint_names = list(start_cons)
            pointing.valid_arr = start_cons_valid_arr[i]

            # time and end constraints (doesn't apply to current pointing)
            if not pointing.current:
                pointing.constraint_names += time_cons
                pointing.valid_arr = np.concatenate((pointing.valid_arr,
                                                     time_cons_valid_arr[i]))
                pointing.constraint_names += [name + '_end' for name in end_cons]
                pointing.valid_arr = np.concatenate((pointing.valid_arr,
                                                    end_cons_valid_arr[i]))

            # save all constraint names on all
            pointing.all_constraint_names = self.all_constraint_names

            # finally find out if each pointing is valid or not
            pointing.valid = np.all(pointing.valid_arr)
            pointing.valid_time = time

    def calculate_tiebreakers(self, time, observer):
        """Calculate the tiebreaker values for every pointing."""
        # Find weight values (0 to 1)
        weights = np.array([float(p.weight) for p in self.pointings])
        weight_arr = 1 - weights
        bad_weight_mask = np.logical_or(weight_arr < 0, weight_arr > 1)
        weight_arr[bad_weight_mask] = 1

        # Find airmass values (0 to 1)
        # airmass at start
        altaz_start = _get_altaz(time, observer, self.targets)['altaz']
        secz_start = altaz_start.secz

        # airmass at end
        altaz_end = _get_altaz(self.finish_times, observer, self.targets)['altaz']
        secz_end = altaz_end.secz

        # take average
        secz_arr = (secz_start + secz_end) / 2.
        airmass_arr = secz_arr.value / 10.
        bad_airmass_mask = np.logical_or(airmass_arr < 0, airmass_arr > 1)
        airmass_arr[bad_airmass_mask] = 1

        # Find time to set values (0 to 1)
        tts_arr = time_to_set(observer, self.targets, time).to(u.hour).value
        tts_arr = tts_arr / 24.
        bad_tts_mask = np.logical_or(tts_arr < 0, tts_arr > 1)
        tts_arr[bad_tts_mask] = 1

        # Construct the tiebreaker value
        total_weight = params.WEIGHTING_WEIGHT + params.AIRMASS_WEIGHT + params.TTS_WEIGHT
        w_weight = params.WEIGHTING_WEIGHT / total_weight
        a_weight = params.AIRMASS_WEIGHT / total_weight
        t_weight = params.TTS_WEIGHT / total_weight
        tiebreak_arr = (weight_arr * w_weight + airmass_arr * a_weight + tts_arr * t_weight)

        # Save values on the pointings
        for i, pointing in enumerate(self.pointings):
            pointing.altaz_start = altaz_start[i]
            pointing.altaz_end = altaz_end[i]
            pointing.weight = weight_arr[i]
            pointing.airmass = airmass_arr[i]
            pointing.tts = tts_arr[i]
            pointing.tiebreaker = tiebreak_arr[i]

    def calculate(self, time, observer, horizon, readout_time):
        """Calculate pointing validities at the given time."""
        # Calculate the expected needed time to observe each Pointing
        self.calculate_times(time, readout_time)

        # Apply constraints to all Pointings
        self.apply_constraints(time, observer, horizon)

        # Calculate tiebreaker values for all Pointings
        self.calculate_tiebreakers(time, observer)

    def get_highest_priority_pointing(self):
        """Return the pointing with the highest priority."""
        # If there are no pointings, return None
        if len(self.pointings) == 0:
            return None

        # Sort the pointings
        #   - First is by validity
        #   - Then by rank
        #   - Then by ToO flag
        #   - Then by number of times already observed
        #   - Finally use the tiebreaker
        pointings = list(self.pointings)  # make a copy
        pointings.sort(key=lambda p: (not p.valid, p.rank, not p.too, p.num_obs, p.tiebreaker))

        return pointings[0]

    def get_highest_priority_pointings(self, time, observer, horizon, readout_time, number=1):
        """Return the top X highest priority pointings."""
        # If there are no pointings, return None
        if len(self.pointings) == 0:
            return [None] * number

        # Sort the pointings
        #   - First is by validity
        #   - Then by rank
        #   - Then by ToO flag
        #   - Then by number of times already observed
        #   - Finally use the tiebreaker
        pointings = list(self.pointings)  # make a copy
        pointings.sort(key=lambda p: (not p.valid, p.rank, not p.too, p.num_obs, p.tiebreaker))

        # Return the top X pointings as requested
        if len(pointings) < number:
            pointings += [None] * (number - len(pointings))
        return pointings[0:number]

    def write_to_file(self, time, observer, filename):
        """Write any time-dependent pointing infomation to a file."""
        # The queue should already have priorities calculated
        if not hasattr(self.pointings[0], 'valid') or not hasattr(self.pointings[0], 'tiebreaker'):
            raise ValueError('Queue has not yet had priorities calculated')

        # get and sort pointings
        pointings = list(self.pointings)  # make a copy
        pointings.sort(key=lambda p: (not p.valid, p.rank, not p.too, p.num_obs, p.tiebreaker))

        # now save as json file
        with open(filename, 'w') as f:
            json.dump(str(time), f)
            f.write('\n')
            json.dump(self.all_constraint_names, f)
            f.write('\n')
            for p in pointings:
                valid_nonbool = [int(b) for b in p.valid_arr]
                json.dump([p.db_id,
                           p.name,
                           1,
                           p.rank,
                           int(p.too),
                           p.num_obs,
                           p.tiebreaker,
                           list(zip(p.constraint_names, valid_nonbool)),
                           ],
                          f)
                f.write('\n')


def what_to_do_next(current_pointing, highest_pointing, log=None):
    """Decide whether to slew to a new target, remain on the current target or park the telescope.

    Parameters
    ----------
    current_pointing : `Pointing`
        The current pointing.
        `None` if the telescope is idle.
    highest_pointing : `Pointing`
        The current highest priority pointing from the queue.
        `None` if the queue is empty.

    log: `logging.Logger`, optional
        log object to direct output to

    Returns
    -------
    new_pointing : `Pointing`
        The new pointing to send to the pilot
        Could be either:
          current_pointing (remain on current target),
          highest_pointing (slew to new target) or
          `None`           (nothing to do, park).

    """
    # Create a logger if one isn't given
    if log is None:
        logging.basicConfig(level=logging.DEBUG)
        log = logging.getLogger('scheduler')

    # Deal with either being missing (telescope is idle or queue is empty)
    if current_pointing is None and highest_pointing is None:
        reason = 'Not doing anything; Nothing to do => Do nothing'
        new_pointing = None
    elif current_pointing is None:
        if highest_pointing.valid:
            reason = 'Not doing anything; HP valid => Do HP'
            new_pointing = highest_pointing
        else:
            reason = 'Not doing anything; HP invalid => Do nothing'
            new_pointing = None
    elif highest_pointing is None:
        if not current_pointing.valid:
            reason = 'CP invalid; Nothing to do => Do nothing'
            new_pointing = None
        else:
            reason = 'CP valid; Nothing to do => Do CP'  # TODO - is that right?
            new_pointing = current_pointing

    elif current_pointing == highest_pointing:
        if not current_pointing.valid or not highest_pointing.valid:
            reason = 'CP==HP and invalid => Do nothing'
            new_pointing = None  # it's either finished or is now illegal
        else:
            reason = 'CP==HP and valid => Do HP'
            new_pointing = highest_pointing

    elif not current_pointing.valid:  # current pointing is illegal (finished)
        if highest_pointing.valid:  # new pointing is legal
            reason = 'CP invalid; HP valid => Do HP'
            new_pointing = highest_pointing
        else:
            reason = 'CP invalid; HP invalid => Do nothing'
            new_pointing = None
    else:  # telescope is observing legally
        if not highest_pointing.valid:  # no legal pointings
            reason = 'CP valid; HP invalid => Do nothing'
            new_pointing = None
        else:  # both are legal
            if highest_pointing.too:  # slew to a ToO, unless now is also a ToO
                if not current_pointing.too:
                    reason = 'CP < HP; CP is not ToO and HP is => Do HP'
                    new_pointing = highest_pointing
                else:
                    reason = 'CP < HP; CP is is ToO and HP is ToO => Do CP'
                    new_pointing = current_pointing
            else:  # stay for normal pointings
                reason = 'CP < HP; but not a ToO => Do CP'
                new_pointing = current_pointing

    # Log decision
    if current_pointing:
        current_str = '{}'.format(current_pointing.db_id)
        if not current_pointing.valid:
            current_str = '*' + current_str
        if current_pointing.too:
            current_str = current_str + '!'
    else:
        current_str = 'None'

    if highest_pointing:
        highest_str = '{}'.format(highest_pointing.db_id)
        if not highest_pointing.valid:
            highest_str = '*' + highest_str
        if highest_pointing.too:
            highest_str = highest_str + '!'
    else:
        highest_str = 'None'

    if new_pointing:
        new_str = '{}'.format(new_pointing.db_id)
        if not new_pointing.valid:
            new_str = '*' + new_str
        if new_pointing.too:
            new_str = new_str + '!'
    else:
        new_str = 'None'

    log.debug('CP={} HP={}: NP={} ({})'.format(current_str, highest_str, new_str, reason))

    return new_pointing


def check_queue(time=None, location=None, horizon=None, readout_time=10,
                write_file=True, write_html=False, log=None):
    """Check the queue and decide what to do.

    Check the current pointings in the queue, find the highest priority at
    the given time and decide whether to slew to it, stay on the current target
    or park the telescope.

    Parameters
    ----------
    time : `~astropy.time.Time`, optional
        The time to calculate the priorities at.
        Default is `astropy.time.Time.now()`.

    location : `~astropy.coordinates.EarthLocation`, optional
        The location of the observer on Earth.
        Default is `EarthLocation.of_site('lapalma')`.

    horizon : float, or tuple of (azs, alts), optional
        The horizon limits at the given site, either a flat value or varying with azimuth.
        Default is a flat horizon of 30 deg.

    readout_time : float, optional
        The time per exposure to add when calculating the expected time to observe a pointing.
        Default is 10 seconds.

    write_file : bool, optional
        Should the scheduler write out the queue to a file?
        Default is True.

    write_html : bool, optional
        Should the scheduler write the HTML queue webpage?
        Default is False.

    log: `logging.Logger`, optional
        log object to direct output to

    Returns
    -------
    new_pointing : `Pointing`
        The pointing to send to the pilot.
        Could be a new pointing, the current pointing or 'None' (park).

    """
    # Use current time if not given
    if time is None:
        time = Time.now()

    # Create an observer and load horizon file if not given
    if location is None:
        location = EarthLocation.of_site('lapalma')
    observer = Observer(location)
    if horizon is None:
        horizon = 30
    if isinstance(horizon, (int, float)):
        horizon = ([0, 90, 180, 270, 360], [horizon, horizon, horizon, horizon, horizon])

    # Import the queue from the database
    queue = PointingQueue.from_database(time, observer)
    if len(queue) == 0:
        return None

    # Calculate pointing validity at the given time
    queue.calculate(time, observer, horizon, readout_time)

    # Get the current pointing and the highest priority pointing
    current_pointing = queue.get_current_pointing()
    highest_pointing = queue.get_highest_priority_pointing()

    # Write out the queue file and web pages
    if write_file:
        queue_file = os.path.join(params.QUEUE_PATH, 'queue_info')
        queue.write_to_file(time, observer, queue_file)
    if write_html:
        # TODO this could be run from elsewhere
        write_queue_page(observer)

    # Work out what to do next
    new_pointing = what_to_do_next(current_pointing, highest_pointing, log)
    return new_pointing
