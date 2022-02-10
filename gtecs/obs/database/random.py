"""Function to create a random pointing for testing."""

from astropy import units as u
from astropy.coordinates import EarthLocation
from astropy.time import Time

import numpy as np

from .models import ExposureSet, Strategy, Target


__all__ = ['make_random_target']


def make_random_target(user, num_expsets=None, time=None):
    """Make a random target for testing.

    It should be observable from La Palma at the time of creation.
    However not all pointings will be valid in the queue immediately due to
    random start and stop times.

    Parameters
    ----------
    user : `User`
        the User to associate the Target with

    num_expsets : int, optional
        the number of Exposure Sets attached to the Target
        if None, a random number of exposure_sets between 1 and 5 will be added
    time : `~astropy.time.Time`, optional
        the time to centre the pointings around
        if None, the current time is used

    Returns
    -------
    target : `Target`
        the new target

    """
    lapalma = EarthLocation(lat=28 * u.deg, lon=-17 * u.deg)
    if time is None:
        time = Time.now()
    # LST in degrees
    lst = time.sidereal_time('mean', longitude=lapalma.lon).deg
    t1 = time + np.random.randint(-2, 2) * u.day
    t2 = t1 + np.random.randint(2, 6) * u.day

    strategy = Strategy(num_todo=np.random.randint(1, 5),
                        too=np.random.randint(0, 2),
                        max_moon=np.random.choice(['D', 'G', 'B']),
                        )
    target = Target(name='random_target',
                    ra=np.random.uniform(lst - 3, lst + 3),
                    dec=np.random.uniform(10, 89),
                    start_rank=np.random.randint(1, 100),
                    start_time=t1,
                    stop_time=t2,
                    user=user,
                    strategy=strategy
                    )

    if num_expsets is None:
        num_expsets = np.random.randint(1, 6)
    for _ in range(num_expsets):
        exposure_set = ExposureSet(num_exp=1,
                                   exptime=np.random.uniform(10., 360.),
                                   filt=np.random.choice(['L', 'R', 'G', 'B']),
                                   )
        target.exposure_sets.append(exposure_set)

    return target
