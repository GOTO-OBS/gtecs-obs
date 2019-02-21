#!/usr/bin/env python
"""Function to create a random pointing for testing."""

from astropy.time import Time

from .models import ExposureSet, Pointing

__all__ = ['make_random_pointing']


def make_random_pointing(user_id, numexps=None, time=None):
    """Make a random pointing for testing.

    It should be observable from La Palma at the time of creation.
    However not all pointings will be valid in the queue immedietly due to
    random start and stop times.

    Parameters
    ----------
    numexps : int
        if None, a random number of exposure_sets between 1 and 5 will be added

    time : `~astropy.time.Time`
        The time to centre the pointings around
        if None, the current time is used

    Returns
    -------
    pointing : `Pointing`
        the new Pointing

    """
    import numpy as np
    from astropy.coordinates import EarthLocation
    from astropy import units as u

    lapalma = EarthLocation(lat=28 * u.deg, lon=-17 * u.deg)
    if time is None:
        time = Time.now()
    # LST in degrees
    lst = time.sidereal_time('mean', longitude=lapalma.lon).deg
    t1 = time + np.random.randint(-5, 2) * u.day
    t2 = t1 + np.random.randint(1, 10) * u.day
    p = Pointing(objectName='randObj',
                 ra=np.random.uniform(lst - 3, lst + 3),
                 dec=np.random.uniform(10, 89),
                 rank=np.random.randint(1, 100),
                 minAlt=30,
                 maxSunAlt=-15,
                 minTime=np.random.uniform(100, 3600),
                 maxMoon=np.random.choice(['D', 'G', 'B']),
                 minMoonSep=30,
                 ToO=np.random.randint(0, 2),
                 startUTC=t1,
                 stopUTC=t2,
                 status='pending',
                 user_id=user_id,
                 )

    if numexps is None:
        numexps = np.random.randint(1, 6)
    for _ in range(numexps):
        exposure_set = ExposureSet(expTime=np.random.uniform(10., 360.),
                                   numexp=1,
                                   typeFlag="SCIENCE",
                                   filt=np.random.choice(['L', 'R', 'G', 'B']),
                                   binning=1,
                                   raoff=0,
                                   decoff=0,
                                   )
        p.exposure_sets.append(exposure_set)

    return p
