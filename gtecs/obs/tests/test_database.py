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
    if not db.validate_user(session, username='test_user', password='password'):
        raise ValueError('Failed to validate user')

print('-------')
with db.open_session() as session:
    # create a site
    si = db.Site.from_name('LaPalma')

    # and two telescopes at that site
    t1 = db.Telescope(name='North Telescope')
    t2 = db.Telescope(name='South Telescope')
    t1.site = si
    t2.site = si

    # and create a grid
    g = db.Grid(name='testgrid',
                ra_fov=10, dec_fov=10,
                ra_overlap=0.5, dec_overlap=0.5,
                algorithm='NA')
    t1.grid = g
    t2.grid = g

    # and a few grid tiles
    gt1 = db.GridTile(name='T0001', ra=100, dec=20)
    gt1.grid = g
    gt2 = db.GridTile(name='T0002', ra=100, dec=40)
    gt2.grid = g

    # create an event
    e = db.Event(name='event', ivorn='ivo://goto', source='made-up', event_type='FAKE')

    # and a survey
    s = db.Survey(name='GOTO event survey')
    s.event = e
    s.grid = g

    # and some survey tiles
    st1 = db.SurveyTile(weight=0.5)
    st1.survey = s
    st1.grid_tile = gt1
    st2 = db.SurveyTile(weight=0.5)
    st2.survey = s
    st2.grid_tile = gt2

    # add them all
    db.insert_items(session, [si, t1, t2, g, gt1, gt2, e, s, st1, st2])
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
    print(st1)
    print(st2, end='\n\n')
time.sleep(2)

print('-------')
with db.open_session() as session:
    user = db.get_user(session, username='test_user')

    # let's make an Target
    target = db.Target(object_name='M31',
                       ra=10.685,
                       dec=41.2875,
                       start_rank=9,
                       min_time=3600,
                       num_todo=5,
                       valid_time=[60, 120],
                       wait_time=60,
                       user=user,
                       )
    # and add RGBL exposure set
    L = db.ExposureSet(num_exp=3, exptime=120, filt='L')
    R = db.ExposureSet(num_exp=3, exptime=120, filt='R')
    G = db.ExposureSet(num_exp=3, exptime=120, filt='G')
    B = db.ExposureSet(num_exp=3, exptime=120, filt='B')
    target.exposure_sets = [L, R, G, B]
    session.add(target)
    session.commit()
    print(target, end='\n\n')
    print(target.pointings, end='\n\n')
    print(target.time_blocks, end='\n\n')
    print(L)
    print(R)
    print(G)
    print(B, end='\n\n')

print('-------')
with db.open_session() as session:
    user = db.get_user(session, username='test_user')

    # new session, add random targets
    targets = [db.make_random_target(user) for i in range(10)]
    db.insert_items(session, targets)
    for target in targets:
        print(target, target.pointings)

    more_targets = [db.make_random_target(user) for i in range(10)]
    db.insert_items(session, more_targets)
    for target in more_targets:
        print(target, target.pointings)

    # now check how many pointings we have (should be 1 for each target)
    n_pointings = len(db.get_pointings(session))
    print('{} pointings in database'.format(n_pointings))

    # get the queue and see how many are pending
    current, pending = db.get_queue(session)
    n_queue = len(pending)
    print('{} pointings in queue'.format(n_queue))
time.sleep(2)

print('-------')
with db.open_session() as s:
    def print_summary(target, now):
        """Print a summary of the database."""
        print(now, f'Target is {target.status.upper()} ({target.num_remaining} remaining)')
        print('                    Pointings:', [p.status_at_time(now) for p in target.pointings])
        # print('                    Time blocks:')
        # for block in target.time_blocks:
        #     print('                    \t', (block.block_num, block.valid_time, block.wait_time,
        #                                      block.current))

    target = db.get_target_by_id(s, 1)
    print(target)
    print(target.pointings)
    print(target.time_blocks, end='\n\n')

    # lets pretend we're observing with a telescope to check the triggers
    t = db.get_telescope_by_id(s, 1)
    now = Time(target.start_time)
    print_summary(target, now)
    print()

    while True:
        next_pointing = target.pointings[-1]
        if next_pointing.status == 'upcoming':
            # skip ahead until the pointing is pending
            print(now, 'Skipping ahead to Pointing start time')
            now = Time(next_pointing.start_time) + 30 * u.s
            print_summary(target, now)
            time.sleep(1)
            print()

        print(now, f'Marking Pointing {next_pointing.db_id} as running')
        next_pointing.mark_running(telescope=t, time=now)
        s.commit()
        print_summary(target, now)
        time.sleep(1)
        print()

        # skip ahead by some exposure time
        now = now + 120 * u.s

        print(now, f'Marking Pointing {next_pointing.db_id} as completed')
        next_pointing.mark_finished(completed=True, time=now)
        s.commit()
        print_summary(target, now)
        time.sleep(1)
        print()

        # create next pointing
        next_pointing = target.get_next_pointing(time=now)
        if next_pointing is None:
            break

        print(now, 'Creating new Pointing')
        s.add(next_pointing)
        s.commit()
        print_summary(target, now)
        time.sleep(1)
        print()

        # skip ahead by a small amount
        now = now + 30 * u.s

    # check the list of target to schedule is empty
    targets_to_schedule = db.get_targets(s, status='unscheduled')
    print('There are {} Targets to schedule'.format(len(targets_to_schedule)))

    n_pointings = len(db.get_pointings(session))
    print('{} pointings in database'.format(n_pointings))
    current, pending = db.get_queue(session)
    n_queue = len(pending)
    print('{} pointings in queue'.format(n_queue))
