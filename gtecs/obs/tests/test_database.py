"""Demonstrate the use of the database modules.

Assume we start from a clean database.
"""

import time

from astropy import units as u
from astropy.time import Time

from gtecs.obs import database as db


# add a user
with db.open_session() as session:
    # create a test user
    user = db.User(username='test_user', password='password', full_name='GOTO Test Observer')
    print(user)
    session.add(user)
    session.commit()
    print(user, end='\n\n')

    # check the password was stored correctly
    if not db.validate_user('test_user', 'password'):
        raise ValueError('Failed to validate user')

print('-------')
with db.open_session() as session:
    # create a site
    si = db.Site.from_name('LaPalma')

    # and two telescopes at that site
    t1 = db.Telescope(name='North Telescope', site=si)
    t2 = db.Telescope(name='South Telescope', site=si)

    # and create a grid
    g = db.Grid(name='testgrid',
                ra_fov=10, dec_fov=10,
                ra_overlap=0.5, dec_overlap=0.5,
                algorithm='minverlap')
    t1.grid = g
    t2.grid = g

    # and a few grid tiles
    gt1 = db.GridTile(name='T0001', ra=100, dec=20, grid=g)
    gt2 = db.GridTile(name='T0002', ra=100, dec=40, grid=g)

    # create an event
    e = db.Event(name='event', source='made-up', type='FAKE')

    # and a survey
    s = db.Survey(name='GOTO event survey', event=e)
    s.event = e

    # add them all
    db.insert_items(session, [si, t1, t2, g, gt1, gt2, e, s])
    session.commit()

    # print them to check __repr__
    print(si, end='\n\n')
    print(t1)
    print(t2, end='\n\n')
    print(g, end='\n\n')
    print(gt1)
    print(gt2, end='\n\n')
    print(e, end='\n\n')
    print(s, end='\n\n')
time.sleep(2)

print('-------')
with db.open_session() as session:
    user = db.get_user(session, username='test_user')

    # let's make an Target with ExposureSets in different filters
    L = db.ExposureSet(num_exp=3, exptime=120, filt='L')
    R = db.ExposureSet(num_exp=3, exptime=120, filt='R')
    G = db.ExposureSet(num_exp=3, exptime=120, filt='G')
    B = db.ExposureSet(num_exp=3, exptime=120, filt='B')
    strategy = db.Strategy(num_todo=5,
                           wait_time=[1 * u.hour, 120 * u.min],
                           valid_time=[None, 1 * u.hour, -1, 60],
                           )
    target = db.Target(name='M31',
                       ra=10.685,
                       dec=41.2875,
                       rank=9,
                       start_time=Time('2020-01-01 00:00'),
                       stop_time=None,
                       creation_time=Time('2020-01-01 00:00'),  # Need to fake
                       user=user,
                       strategy=strategy,
                       exposure_sets=[L, R, G, B],
                       )
    session.add(target)
    session.commit()
    print(target, end='\n\n')
    print(target.pointings, end='\n\n')
    print(target.strategy, end='\n\n')
    print(target.strategy.time_blocks, end='\n\n')
    print(L)
    print(R)
    print(G)
    print(B, end='\n\n')

print('-------')
with db.open_session() as session:
    user = db.get_user(session, username='test_user')

    # new session, add random targets
    targets = [db.make_random_target(user) for _ in range(5)]
    db.insert_items(session, targets)
    for target in targets:
        print(target)
        print(target.strategy)
        print(target.pointings)

    # make more targets, but only valid for a given telescope
    more_targets = [db.make_random_target(user) for _ in range(5)]
    for i, target in enumerate(more_targets):
        # alternate between tel 1 and tel 2
        target.strategy.tel_mask = i % 2 + 1
    db.insert_items(session, more_targets)
    for target in more_targets:
        print(target)
        print(target.strategy)
        print(target.pointings)

    # now check how many pointings we have (should be 1 for each target)
    n_pointings = len(db.get_pointings(session))
    print('{} pointings in database'.format(n_pointings))

    # check how many are pending for each telescope
    for telescope_id in [1, 2]:
        n_queue = len(db.get_pending_pointings(session, telescope_id))
        print('{} pointings in queue for telescope {}'.format(n_queue, telescope_id))
time.sleep(2)

print('-------')
with db.open_session() as session:
    def print_summary(target, now):
        """Print a summary of the database."""
        print(now, end=' ')
        print(f'Target is {target.status.upper()}', end=', ')
        print(f'completed {target.num_completed} pointings')
        print('                    ', end='')
        print('Pointings:', [p.status_at_time(now) for p in target.pointings])
        # print('                    Time blocks:')
        # for block in target.time_blocks:
        #     print('                    \t', (block.block_num, block.valid_time, block.wait_time,
        #                                      block.current))

    target = db.get_target_by_id(session, 1)
    print(target)
    print(target.pointings, end='\n\n')

    # lets pretend we're observing with a telescope to check the triggers
    t = db.get_telescope_by_id(session, 1)
    now = Time(target.start_time)
    print_summary(target, now)
    print()

    while True:
        next_pointing = target.pointings[-1]
        if next_pointing.status_at_time(now) == 'upcoming':
            # skip ahead until the pointing is pending
            print(now, 'Skipping ahead to Pointing start time')
            now = Time(next_pointing.start_time) + 30 * u.s
            print_summary(target, now)
            time.sleep(1)
            print()

        print(now, f'Marking Pointing {next_pointing.db_id} as running')
        next_pointing.mark_running(telescope=t, time=now)
        session.commit()
        print_summary(target, now)
        time.sleep(1)
        print()

        # skip ahead by some exposure time
        now = now + 120 * u.s

        print(now, f'Marking Pointing {next_pointing.db_id} as completed')
        next_pointing.mark_finished(completed=True, time=now)
        session.commit()
        print_summary(target, now)
        time.sleep(1)
        print()

        # create next pointing
        next_pointing = target.get_next_pointing(time=now)
        if next_pointing is None:
            break

        print(now, 'Creating new Pointing')
        session.add(next_pointing)
        session.commit()
        print_summary(target, now)
        time.sleep(1)
        print()

        # skip ahead by a small amount
        now = now + 30 * u.s

    # check the list of target to schedule is empty
    targets_to_schedule = db.get_targets(session, status='unscheduled')
    print('There are {} Targets to schedule'.format(len(targets_to_schedule)))

    n_pointings = len(db.get_pointings(session))
    print('{} pointings in database'.format(n_pointings))
time.sleep(2)

print('-------')
# Now let's run through a target again, but from the point of view of the pilot
now = Time('2020-01-01 00:00')

# First create a new, basic Target
with db.open_session() as session:
    user = db.get_user(session, username='test_user')

    # Let's make an Target with ExposureSets in different filters
    exps = db.ExposureSet(num_exp=3, exptime=120, filt='L')
    strategy = db.Strategy(num_todo=5,
                           wait_time=1 * u.hour,
                           valid_time=None,
                           )
    target = db.Target(name='TEST',
                       ra=0,
                       dec=50,
                       rank=9,
                       start_time=now,
                       stop_time=now + 1 * u.day,
                       creation_time=now,  # Need to fake
                       user=user,
                       strategy=strategy,
                       exposure_sets=[exps],
                       )
    session.add(target)
    session.commit()

    target_id = target.db_id
    pointing_id = target.pointings[-1].db_id

# ~~~~~~~~~~~
# The first pointing should have been created and be pending
with db.open_session() as session:
    target = db.get_target_by_id(session, target_id)
    print(now, target.status_at_time(now), [p.status_at_time(now) for p in target.pointings])
    print()

    # Check the pointing was created
    if len(target.pointings) != 1 or target.pointings[0].status_at_time(now) != 'pending':
        raise ValueError
    if target.pointings[0].start_time != Time('2020-01-01 00:00'):
        raise ValueError

# ~~~~~~~~~~~
# Let's go on 10 minutes, then mark it as running
# now += 10 * u.minute  # We can't keep adding, the floating point errors add up too quickly
now = Time('2020-01-01 00:10')
db.mark_pointing_running(pointing_id, telescope_id=1, time=now)

with db.open_session() as session:
    target = db.get_target_by_id(session, target_id)
    print(now, target.status_at_time(now), [p.status_at_time(now) for p in target.pointings])
    print()

    # Check it marked correctly
    if target.pointings[0].status_at_time(now) != 'running':
        raise ValueError
    if target.pointings[0].running_time != Time('2020-01-01 00:10'):
        raise ValueError

# ~~~~~~~~~~~
# Wait 5 minutes, mark it as completed
now = Time('2020-01-01 00:15')
db.mark_pointing_completed(pointing_id, time=now)

with db.open_session() as session:
    target = db.get_target_by_id(session, target_id)
    print(now, target.status_at_time(now), [p.status_at_time(now) for p in target.pointings])
    print()

    # Check it marked correctly
    if target.pointings[0].status_at_time(now) != 'completed':
        raise ValueError
    if target.pointings[0].finished_time != Time('2020-01-01 00:15'):
        raise ValueError
    # A new pointing should have been created, starting 60 minutes from now
    if len(target.pointings) != 2 or target.pointings[1].status_at_time(now) != 'upcoming':
        raise ValueError
    if target.pointings[1].start_time != Time('2020-01-01 01:15'):
        raise ValueError

    # Get the new pointing ID
    pointing_id = target.pointings[1].db_id

# ~~~~~~~~~~~
# Let's go on 65 minutes until after it's valid, then mark it as running
now = Time('2020-01-01 01:20')
db.mark_pointing_running(pointing_id, telescope_id=1, time=now)

with db.open_session() as session:
    target = db.get_target_by_id(session, target_id)
    print(now, target.status_at_time(now), [p.status_at_time(now) for p in target.pointings])
    print()

    # Check it marked correctly
    if target.pointings[1].status_at_time(now) != 'running':
        raise ValueError
    if target.pointings[1].running_time != Time('2020-01-01 01:20'):
        raise ValueError

# ~~~~~~~~~~~
# Wait 3 minutes, mark it as interrupted
now = Time('2020-01-01 01:23')
db.mark_pointing_interrupted(pointing_id, time=now)

with db.open_session() as session:
    target = db.get_target_by_id(session, target_id)
    print(now, target.status_at_time(now), [p.status_at_time(now) for p in target.pointings])
    print()

    # Check it marked correctly
    if target.pointings[1].status_at_time(now) != 'interrupted':
        raise ValueError
    if target.pointings[1].finished_time != Time('2020-01-01 01:23'):
        raise ValueError
    # A new pointing should have been created, starting at the previous start time
    if len(target.pointings) != 3 or target.pointings[2].status_at_time(now) != 'pending':
        raise ValueError
    if target.pointings[2].start_time != Time('2020-01-01 01:15'):
        raise ValueError

    # Get the new pointing ID
    pointing_id = target.pointings[2].db_id

# ~~~~~~~~~~~
# Let's wait for 2 minutes, then mark it as running
now = Time('2020-01-01 01:25')
db.mark_pointing_running(pointing_id, telescope_id=1, time=now)

with db.open_session() as session:
    target = db.get_target_by_id(session, target_id)
    print(now, target.status_at_time(now), [p.status_at_time(now) for p in target.pointings])
    print()

    # Check it marked correctly
    if target.pointings[2].status_at_time(now) != 'running':
        raise ValueError
    if target.pointings[2].running_time != Time('2020-01-01 01:25'):
        raise ValueError

# ~~~~~~~~~~~
# Wait another 5 minutes, mark it as completed
now = Time('2020-01-01 01:30')
db.mark_pointing_completed(pointing_id, time=now)

with db.open_session() as session:
    target = db.get_target_by_id(session, target_id)
    print(now, target.status_at_time(now), [p.status_at_time(now) for p in target.pointings])
    print()

    # Check it marked correctly
    if target.pointings[2].status_at_time(now) != 'completed':
        raise ValueError
    if target.pointings[2].finished_time != Time('2020-01-01 01:30'):
        raise ValueError
    # A new pointing should have been created, starting 60 minutes from now
    if len(target.pointings) != 4 or target.pointings[3].status_at_time(now) != 'upcoming':
        raise ValueError
    if target.pointings[3].start_time != Time('2020-01-01 02:30'):
        raise ValueError

# ~~~~~~~~~~~
# Now let's wait 5 minutes and mark the previous pointing as failed, and delay the new one by 1h
now = Time('2020-01-01 01:35')
db.mark_pointing_failed(pointing_id, time=now, delay=3600)

with db.open_session() as session:
    target = db.get_target_by_id(session, target_id)
    print(now, target.status_at_time(now), [p.status_at_time(now) for p in target.pointings])
    print()

    # Check it marked correctly
    if target.pointings[2].status_at_time(now) != 'failed':
        raise ValueError
    if target.pointings[2].validated_time != Time('2020-01-01 01:35'):
        raise ValueError
    # The previously pending pointing should have been deleted
    if target.pointings[3].status_at_time(now) != 'deleted':
        raise ValueError
    # And a new pointing should have been created, with the start time 1h after the last finished
    if len(target.pointings) != 5 or target.pointings[4].status_at_time(now) != 'upcoming':
        raise ValueError
    if target.pointings[4].start_time != Time('2020-01-01 02:30'):
        raise ValueError

# ~~~~~~~~~~~
# Finally let's wait for a day so the last Target should have expired
now = Time('2020-01-02 01:35')

with db.open_session() as session:
    target = db.get_target_by_id(session, target_id)
    print(now, target.status_at_time(now), [p.status_at_time(now) for p in target.pointings])
    print()

    # Check it marked correctly
    if target.pointings[4].status_at_time(now) != 'expired':
        raise ValueError
    # it's now too late to do any more Pointings for this Target
    if target.status_at_time(now) != 'expired':
        raise ValueError
    if target.num_completed != 1:
        raise ValueError
