#!/usr/bin/env python
"""Function to create a random pointing for testing."""

from astropy.time import Time

from .models import ExposureSet, Pointing

__all__ = ['make_random_pointing']


def make_random_pointing(user_id, num_expsets=None, time=None):
    """Make a random pointing for testing.

    It should be observable from La Palma at the time of creation.
    However not all pointings will be valid in the queue immedietly due to
    random start and stop times.

    Parameters
    ----------
    num_expsets : int
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
    p = Pointing(object_name='randObj',
                 ra=np.random.uniform(lst - 3, lst + 3),
                 dec=np.random.uniform(10, 89),
                 rank=np.random.randint(1, 100),
                 min_alt=30,
                 max_sunalt=-15,
                 min_time=np.random.uniform(100, 3600),
                 max_moon=np.random.choice(['D', 'G', 'B']),
                 min_moonsep=30,
                 too=np.random.randint(0, 2),
                 start_time=t1,
                 stop_time=t2,
                 status='pending',
                 user_id=user_id,
                 )

    if num_expsets is None:
        num_expsets = np.random.randint(1, 6)
    for _ in range(num_expsets):
        exposure_set = ExposureSet(exptime=np.random.uniform(10., 360.),
                                   num_exp=1,
                                   imgtype="SCIENCE",
                                   filt=np.random.choice(['L', 'R', 'G', 'B']),
                                   binning=1,
                                   ra_offset=0,
                                   dec_offset=0,
                                   )
        p.exposure_sets.append(exposure_set)

    return p
