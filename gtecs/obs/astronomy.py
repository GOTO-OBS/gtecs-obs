"""Astronomy utilities."""

import warnings

from astropy import units as u
from astropy.coordinates import AltAz, GCRS, Longitude, SkyCoord, get_sun
from astropy.time import Time

import numpy as np
from numpy.polynomial.polynomial import polyval

from scipy import interpolate


def _equation_of_time(time):
    """Find the difference between apparent and mean solar time.

    Parameters
    ----------
    time : `~astropy.time.Time`
        times (array)

    Returns
    -------
    ret1 : `~astropy.units.Quantity`
        the equation of time

    """
    # Julian centuries since J2000.0
    T = (time - Time("J2000")).to(u.year).value / 100

    # obliquity of ecliptic (Meeus 1998, eq 22.2)
    poly_pars = (84381.448, 46.8150, 0.00059, 0.001813)
    eps = u.Quantity(polyval(T, poly_pars), u.arcsec)
    y = np.tan(eps / 2)**2

    # Sun's mean longitude (Meeus 1998, eq 25.2)
    poly_pars = (280.46646, 36000.76983, 0.0003032)
    L0 = u.Quantity(polyval(T, poly_pars), u.deg)

    # Sun's mean anomaly (Meeus 1998, eq 25.3)
    poly_pars = (357.52911, 35999.05029, 0.0001537)
    M = u.Quantity(polyval(T, poly_pars), u.deg)

    # eccentricity of Earth's orbit (Meeus 1998, eq 25.4)
    poly_pars = (0.016708634, -0.000042037, -0.0000001267)
    e = polyval(T, poly_pars)

    # equation of time, radians (Meeus 1998, eq 28.3)
    eot = (y * np.sin(2 * L0) - 2 * e * np.sin(M) + 4 * e * y * np.sin(M) * np.cos(2 * L0) -
           0.5 * y**2 * np.sin(4 * L0) - 5 * e**2 * np.sin(2 * M) / 4) * u.rad
    return eot.to(u.hourangle)


def _astropy_time_from_lst(time, lst, location, prev_next):
    """Convert a Local Sidereal Time to an astropy Time object.

    The local time is related to the LST through the RA of the Sun.
    This routine uses this relationship to convert a LST to an astropy
    time object.

    Returns
    -------
    ret1 : `~astropy.time.Time`
        time corresponding to LST

    """
    # now we need to figure out time to return from LST
    sun_ra = get_sun(time).ra

    # calculate Greenwich Apparent Solar Time, which we will use as ~UTC for now
    good_mask = ~np.isnan(lst)
    solar_time = lst[good_mask] - sun_ra + 12 * u.hourangle - location.lon

    # assume this is on the same day as supplied time, and fix later
    first_guess = Time(u.d * int(time.mjd) + u.hour * solar_time.wrap_at('360d').hour,
                       format='mjd')

    # Equation of time is difference between GAST and UTC
    eot = _equation_of_time(first_guess)
    first_guess = first_guess - u.hour * eot.value

    if prev_next == 'next':
        # if 'next', we want time to be greater than given time
        mask = first_guess < time
        rise_set_time = first_guess + mask * u.sday
    else:
        # if 'previous', we want time to be less than given time
        mask = first_guess > time
        rise_set_time = first_guess - mask * u.sday

    retvals = -999 * np.ones_like(lst.value)
    retvals[good_mask] = rise_set_time.jd
    return Time(retvals, format='jd')


def _rise_set_trig(time, target, location, prev_next, rise_set):
    """Crude time at next rise/set of ``target`` using spherical trig.

    This method is ~15 times faster than `_calcriseset`,
    and inherently does *not* take the atmosphere into account.

    The time returned should not be used in calculations; the purpose
    of this routine is to supply a guess to `_calcriseset`.

    Parameters
    ----------
    time : `~astropy.time.Time` or other (see below)
        Time of observation. This will be passed in as the first argument to
        the `~astropy.time.Time` initializer, so it can be anything that
        `~astropy.time.Time` will accept (including a `~astropy.time.Time`
        object)

    target : `~astropy.coordinates.SkyCoord`
        Position of target or multiple positions of that target
        at multiple times (if target moves, like the Sun)

    location : `~astropy.coordinates.EarthLocation`
        Observatory location

    prev_next : str - either 'previous' or 'next'
        Test next rise/set or previous rise/set

    rise_set : str - either 'rising' or 'setting'
        Compute prev/next rise or prev/next set

    Returns
    -------
    ret1 : `~astropy.time.Time`
        Time of rise/set

    """
    dec = target.transform_to(GCRS).dec
    cos_ha = -np.tan(dec) * np.tan(location.lat.radian)
    # find the absolute value of the hour angle
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        ha = Longitude(np.fabs(np.arccos(cos_ha)))
    # if rise, hour angle is -ve and vice versa
    if rise_set == 'rising':
        ha = -1 * ha
    # LST = HA + RA
    lst = ha + target.ra

    return _astropy_time_from_lst(time, lst, location, prev_next)


def time_to_set(observer, targets, now):
    """Time until ``target``s next set below the horizon."""
    # Create grid of times
    set_times = _rise_set_trig(now, targets, observer.location, 'next', 'setting')
    seconds_until_set = (set_times - now).to(u.s)
    return seconds_until_set


def above_horizon(ra_deg, dec_deg, location, time=None, horizon=30):
    """Check if the given coordinates are above the artificial horizon.

    Parameters
    ----------
    ra_deg : float or numpy.ndarray
        right ascension in degrees
    dec_deg : float or numpy.ndarray
        declination in degrees
    location : `astropy.coordinates.EarthLocation`
        location of the Observer
    time : `~astropy.time.Time`, optional
        time(s) to calculate at
        default is `Time.now()`
    horizon : float or tuple of (azs, alts), optional
        artificial horizon, either a flat value or varying with azimuth.
        default is a flat horizon of 30 deg

    """
    # Create a flat horizon if not given
    if isinstance(horizon, (int, float)):
        horizon = ([0, 90, 180, 270, 360], [horizon, horizon, horizon, horizon, horizon])

    # Get current altaz
    coord = SkyCoord(ra_deg * u.deg, dec_deg * u.deg)
    frame = AltAz(obstime=time, location=location)
    altaz = coord.transform_to(frame)

    # Check if the altitude is above the horizon
    get_alt_limit = interpolate.interp1d(*horizon,
                                         bounds_error=False,
                                         fill_value='extrapolate')
    alt_limit = get_alt_limit(altaz.az) * u.deg
    return altaz.alt > alt_limit
