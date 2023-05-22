"""Robotic queue scheduler functions."""

import json
import warnings

from astroplan import (AltitudeConstraint, AtNightConstraint,
                       Constraint, MoonIlluminationConstraint, MoonSeparationConstraint,
                       Observer, TimeConstraint)
from astroplan.constraints import _get_altaz

from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time

from erfa import ErfaWarning

import numpy as np

from scipy import interpolate

from . import database as db
from . import params
from .astronomy import time_to_set

warnings.simplefilter('ignore', ErfaWarning)


MOON_PHASES = {'B': 1, 'G': 0.65, 'D': 0.25}  # Should be in params?


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


class TelescopeConstraint:
    """Ensure a Pointing is valid for the given telescope."""

    def __init__(self, telescope_id):
        self.telescope_id = telescope_id

    def compute_constraint(self, pointing):
        """Compute the constraint."""
        if (pointing.valid_telescopes is not None and
                self.telescope_id not in pointing.valid_telescopes):
            # Only valid for certain telescope(s), not including this one
            return False
        if (pointing.current_telescope is not None and
                pointing.current_telescope != self.telescope_id):
            # Valid for this telescope, but currently being observed by a different one
            return False
        return True


class TemplateConstraint:
    """Ensure a Pointing has a valid template."""

    def __init__(self, telescope_id, telescopes_at_site=None, requirement='ANY'):
        self.telescope_id = telescope_id
        self.telescopes_at_site = telescopes_at_site
        if self.telescopes_at_site is None:
            self.telescopes_at_site = [telescope_id]
        if requirement not in ['ANY', 'TELESCOPE', 'SITE']:
            raise ValueError('Invalid template requirement: {}'.format(requirement))
        self.requirement = requirement

    def compute_constraint(self, pointing):
        """Compute the constraint."""
        if pointing.template_telescopes is None:
            # Doesn't need a template
            return True
        if self.requirement == 'ANY' and len(pointing.template_telescopes) > 0:
            # At least one template has been observed for this tile
            return True
        if self.requirement == 'TELESCOPE' and self.telescope_id in pointing.template_telescopes:
            # A template has been observed by the given telescope
            return True
        if self.requirement == 'SITE' and any(t in pointing.template_telescopes
                                              for t in self.telescopes_at_site):
            # A template has been observed by a telescope at the correct site
            return True
        return False


class Pointing:
    """A class to contain information on each Pointing in the Queue."""

    def __init__(self, db_id, name, ra, dec, rank, weight, num_obs, too,
                 mintime, maxsunalt, minalt, maxmoon, minmoonsep,
                 start, stop, valid_telescopes, template_telescopes, exp_sets,
                 time=None, current_telescope=None):
        self.db_id = int(db_id)
        self.name = name
        self.ra = float(ra)
        self.dec = float(dec)
        self.rank = int(rank) if rank is not None else np.inf
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
        self.valid_telescopes = valid_telescopes
        self.template_telescopes = template_telescopes
        self.exp_sets = exp_sets
        self.current_telescope = current_telescope

        if time is None:
            time = Time.now()
        self.time = time

    def __eq__(self, other):
        try:
            return (self.db_id == other.db_id) and (self.time == other.time)
        except AttributeError:
            return False

    def __repr__(self):
        template = ('<Pointing: db_id={}, name="{}">')
        return template.format(self.db_id, self.name)

    @classmethod
    def from_database(cls, db_pointing, time=None, current_telescope=None):
        """Import a pointing from the database."""
        if time is None:
            time = Time.now()

        # Get DB Pointing properties
        db_id = db_pointing.db_id
        rank = db_pointing.rank
        start_time = db_pointing.start_time
        stop_time = db_pointing.stop_time

        # Get linked Target properties
        name = db_pointing.name
        ra = db_pointing.ra
        dec = db_pointing.dec
        weight = db_pointing.weight
        num_completed = db_pointing.target.num_completed

        # Get linked Strategy properties
        too = db_pointing.too
        min_time = db_pointing.min_time
        max_sunalt = db_pointing.max_sunalt
        min_alt = db_pointing.min_alt
        max_moon = db_pointing.max_moon
        min_moonsep = db_pointing.min_moonsep
        valid_telescopes = db_pointing.strategy.valid_telescopes

        # Get linked GridTile properties
        if db_pointing.requires_template:
            if db_pointing.grid_tile is not None:
                template_telescopes = db_pointing.grid_tile.get_template_telescopes()
            else:
                template_telescopes = []
        else:
            template_telescopes = None  # no restriction

        # Get linked ExposureSet properties
        # (from target not pointing, it's much quicker!)
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
                       valid_telescopes=valid_telescopes,
                       template_telescopes=template_telescopes,
                       exp_sets=exp_sets,
                       time=time,
                       current_telescope=current_telescope,
                       )
        return pointing

    def get_obstime(self, readout_time=15):
        """Get the expected time needed to observe this Pointing."""
        return sum(es[0] * (es[1] + readout_time) for es in self.exp_sets)


class PointingQueue:
    """A class to represent a queue of pointings."""

    def __init__(self, pointings, site_id, time=None):
        self.pointings = pointings
        self.site_id = site_id
        if time is None:
            time = Time.now()
        self.time = time

        # Create pointing data arrays
        self.targets = SkyCoord([float(p.ra) for p in self.pointings],
                                [float(p.dec) for p in self.pointings],
                                unit=u.deg, frame='icrs')

        # Get site details from the database
        self.site_data = db.get_site_info(site_id)
        self.observer = Observer(self.site_data['location'])

        # Mark that this queue has not been processed yet
        self.priorities_calculated = False

    def __len__(self):
        return len(self.pointings)

    def __repr__(self):
        template = ('PointingQueue(<{} Pointings>, site_id={}, time={})')
        return template.format(len(self.pointings), self.site_id, self.time)

    @classmethod
    def from_database(cls, site_id, time=None, altitude_limit=None, hourangle_limit=None):
        """Create a Pointing Queue from Pointings in the database.

        Parameters
        ----------
        site_id : int
            The site to fetch the queue for.

        time : `~astropy.time.Time`
            The time to fetch the queue for.
            default=Time.now()

        altitude_limit : float, optional
            Only pointings which ever rise above this altitude are returned.
            default = params.HARD_ALT_LIM
        hourangle_limit : float, optional
            Only pointings with hour angles closer to transit than this limit are returned.
            default = params.HARD_HA_LIM

        Returns
        -------
        pointings : list of `Pointing`
            An list containing Pointings from the database.

        """
        if time is None:
            time = Time.now()
        if altitude_limit is None:
            altitude_limit = params.HARD_ALT_LIM
        if hourangle_limit is None:
            hourangle_limit = params.HARD_HA_LIM

        with db.session_manager() as session:
            # Get the site
            site = db.get_site_by_id(session, site_id)

            # Get pending pointings
            pending_pointings = db.get_pending_pointings(session,
                                                         time=time,
                                                         location=site.location,
                                                         altitude_limit=altitude_limit,
                                                         hourangle_limit=hourangle_limit,
                                                         )
            pointings = [Pointing.from_database(db_pointing, time=time)
                         for db_pointing in pending_pointings]

            # Also get the current pointing for all telescopes at this site, if any
            for telescope in site.telescopes:
                current_pointing = db.get_current_pointing(session, telescope.db_id, time=time)
                if current_pointing is not None:
                    # Note we set the current telescope ID here, it's much faster than finding the
                    # status_at_time for every Pointing.
                    pointing = Pointing.from_database(current_pointing, time=time,
                                                      current_telescope=telescope.db_id)
                    pointings.append(pointing)

            # Create the Queue class
            queue = cls(pointings, site_id, time)
        return queue

    def calculate_times(self, readout_time):
        """Calculate the finish times for all pointings."""
        # Note we use the mintime if a Pointing has it, we don't care if it's invalid after then
        obstime_arr = [p.get_obstime(readout_time)
                       if p.mintime is None else p.mintime
                       for p in self.pointings]
        self.finish_times = self.time + u.Quantity(obstime_arr, unit=u.s)

    def apply_constraints(self, telescope_id, horizon_id, template_requirement):
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
        horizon = self.site_data['telescopes'][telescope_id]['horizon'][horizon_id]
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

        # Telescope
        self.constraints['Telescope'] = TelescopeConstraint(telescope_id)

        # Template
        telescopes_at_site = sorted(self.site_data['telescopes'].keys())
        self.constraints['Template'] = TemplateConstraint(telescope_id,
                                                          telescopes_at_site,
                                                          template_requirement)

        # Apply constraints at the start time (now)
        start_cons = ['SunAlt', 'MinAlt', 'ArtHoriz', 'Moon', 'MoonSep']
        start_cons_valid_arr = apply_constraints([self.constraints[name] for name in start_cons],
                                                 self.observer, self.targets, self.time)

        # Apply constraints at the finish time (now + expected observing time)
        end_cons = ['SunAlt', 'MinAlt', 'ArtHoriz']
        end_cons_valid_arr = apply_constraints([self.constraints[name] for name in end_cons],
                                               self.observer, self.targets, self.finish_times)
        end_cons = [name + '_end' for name in end_cons]

        # Apply time constraint
        time_cons = ['Time']
        time_cons_valid_arr = apply_constraints([self.constraints[name] for name in time_cons],
                                                self.observer, self.targets, self.time)

        # Calculate telescope constraints
        telescope_cons = ['Telescope', 'Template']
        telescope_cons_valid_arr = [[self.constraints['Telescope'].compute_constraint(p),
                                     self.constraints['Template'].compute_constraint(p)]
                                    for p in self.pointings]
        telescope_cons_valid_arr = np.array(telescope_cons_valid_arr)

        # Save constraint results on Pointings and calculate if they are valid
        self.all_constraint_names = start_cons + time_cons + end_cons + telescope_cons
        for i, pointing in enumerate(self.pointings):
            # start constraints (apply to everything)
            pointing.constraint_names = list(start_cons)
            pointing.valid_arr = start_cons_valid_arr[i]

            # time and end constraints (doesn't apply to Pointings already running)
            if pointing.current_telescope != telescope_id:
                pointing.constraint_names += time_cons
                pointing.valid_arr = np.concatenate((pointing.valid_arr,
                                                     time_cons_valid_arr[i]))
                pointing.constraint_names += end_cons
                pointing.valid_arr = np.concatenate((pointing.valid_arr,
                                                    end_cons_valid_arr[i]))

            # telescope constraints (apply to everything)
            pointing.constraint_names += telescope_cons
            pointing.valid_arr = np.concatenate((pointing.valid_arr,
                                                 telescope_cons_valid_arr[i]))

            # save all constraint names on all
            pointing.all_constraint_names = self.all_constraint_names
            pointing.constraint_dict = {name: valid for name, valid in
                                        zip(pointing.constraint_names, pointing.valid_arr)}

            # finally find out if each pointing is valid or not
            pointing.valid = np.all(pointing.valid_arr)
            pointing.valid_time = self.time

    def calculate_tiebreakers(self):
        """Calculate the tiebreaker values for every pointing."""
        # Find weight values (0 to 1)
        weights = np.array([float(p.weight) for p in self.pointings])
        weight_arr = 1 - weights
        bad_weight_mask = np.logical_or(weight_arr < 0, weight_arr > 1)
        weight_arr[bad_weight_mask] = 1

        # Find airmass values (0 to 1)
        altaz = _get_altaz(self.time, self.observer, self.targets)['altaz']
        altaz_end = _get_altaz(self.finish_times, self.observer, self.targets)['altaz']
        secz_arr = (altaz.secz + altaz_end.secz) / 2.
        airmass_arr = secz_arr.value / 10.
        bad_airmass_mask = np.logical_or(airmass_arr < 0, airmass_arr > 1)
        airmass_arr[bad_airmass_mask] = 1

        # Find time to set values (0 to 1)
        tts_arr = time_to_set(self.observer, self.targets, self.time).to(u.hour).value
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
            pointing.altaz = altaz[i]
            pointing.altaz_end = altaz_end[i]
            pointing.weight = weight_arr[i]
            pointing.airmass = airmass_arr[i]
            pointing.tts = tts_arr[i]
            pointing.tiebreaker = tiebreak_arr[i]

    def calculate_priorities(self, telescope_id, horizon, readout_time, template_requirement):
        """Calculate pointing validities and priorities for the given observer."""
        if len(self.pointings) == 0:
            return

        # Calculate the expected needed time to observe each Pointing
        self.calculate_times(readout_time)

        # Apply constraints to all Pointings
        self.apply_constraints(telescope_id, horizon, template_requirement)

        # Calculate tiebreaker values for all Pointings
        self.calculate_tiebreakers()

        # Mark the queue as processed
        self.priorities_calculated = True

    def get_sorted_pointings(self):
        """Return the pointings sorted by priority.

        Sorting metrics
            - First is by validity
            - Then by rank
            - Then by ToO flag
            - Then by number of times already observed
            - Finally use the tiebreaker
        """
        if not self.priorities_calculated:
            raise ValueError('Queue has not yet been processed')
        pointings = self.pointings.copy()
        pointings.sort(key=lambda p: (not p.valid,
                                      p.rank,
                                      not p.too,
                                      p.num_obs,
                                      p.tiebreaker))
        return pointings

    def get_highest_priority_pointing(self):
        """Return the pointing with the highest priority."""
        if len(self.pointings) == 0:
            return None
        return self.get_sorted_pointings()[0]

    def get_current_pointing(self, telescope_id):
        """Return the current pointing from the queue (there should only be one)."""
        for p in self.pointings:
            if p.current_telescope == telescope_id:
                return p
        return None

    def chose_pointing(self, highest_pointing, current_pointing, return_reason=False):
        """Decide what to do based on the current queue."""
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

        if not return_reason:
            return new_pointing

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

        log_str = 'CP={} HP={}: NP={} ({})'.format(current_str, highest_str, new_str, reason)
        return new_pointing, log_str

    def get_pointing(self, telescope_id, horizon=0, readout_time=10, template_requirement='ANY',
                     return_reason=False):
        """Decide what to do based on the current queue.

        Options are to slew to a new target, remain on the current target, or park the telescope.

        Parameters
        ----------
        telescope_id : int
            The telescope ID to calculate the queue for.

        horizon : float, default=0
            Which telescope horizon to use.
        readout_time : float, default=10
            Readout time in seconds, to use when estimating time to observe each pointing.
        template_requirement : str
            How restrictive to consider the templates required for pointings.
            One of 'ANY', 'SITE', or 'TELESCOPE'

        return_reason: bool, default=False
            if True, return a string explaining why the result was chosen

        Returns
        -------
        new_pointing : `Pointing` or None
            The new pointing to send to the pilot
            Could be either:
            current_pointing (remain on current target),
            highest_pointing (slew to new target) or
            `None`           (nothing to do, park).
        reason : str
            Returns if return_reason is True.

        """
        if telescope_id not in self.site_data['telescopes']:
            raise ValueError('Telescope ID {} is not defined for Site {}'.format(
                telescope_id, self.site_id))
        if len(self.site_data['telescopes'][telescope_id]['horizon']) < horizon + 1:
            raise IndexError('Telescope {} only has {} horizons, cannot get horizon {}'.format(
                telescope_id, len(self.site_data['telescopes'][telescope_id]['horizon']), horizon))

        # Calculate priorities and get pointings
        self.calculate_priorities(telescope_id, horizon, readout_time, template_requirement)
        highest_pointing = self.get_highest_priority_pointing()
        current_pointing = self.get_current_pointing(telescope_id)

        # Decide what to do and return
        return self.chose_pointing(highest_pointing, current_pointing, return_reason)

    def write_to_file(self, filename):
        """Write any time-dependent pointing information to a file."""
        if len(self.pointings) == 0:
            return

        # The queue should already have priorities calculated
        if not hasattr(self.pointings[0], 'valid') or not hasattr(self.pointings[0], 'tiebreaker'):
            raise ValueError('Queue has not yet had priorities calculated')

        # get and sort pointings
        pointings = list(self.pointings)  # make a copy
        pointings.sort(key=lambda p: (not p.valid, p.rank, not p.too, p.num_obs, p.tiebreaker))

        # now save as json file
        with open(filename, 'w') as f:
            json.dump(str(self.time), f)
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
