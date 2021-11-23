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
    user = db.User(username='goto_test', password='password', full_name='GOTO Test Observer')
    print(user)
    session.add(user)
    session.commit()
    print(user, end='\n\n')

    # check the password was stored correctly
    if not db.validate_user(session, username='goto_test', password='password'):
        raise ValueError('Failed to validate user')

with db.open_session() as session:
    # create an event
    e = db.Event(name='event', ivorn='ivo://goto', source='made-up', event_type='FAKE')

    # create a grid
    g = db.Grid(name='testgrid',
                ra_fov=10, dec_fov=10,
                ra_overlap=0.5, dec_overlap=0.5,
                algorithm='NA')

    # and a few grid tiles
    gt1 = db.GridTile(name='T0001', ra=100, dec=20)
    gt1.grid = g
    gt2 = db.GridTile(name='T0002', ra=100, dec=40)
    gt2.grid = g

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

    # add them
    db.insert_items(session, [e, g, gt1, gt2, s, st1, st2])
    session.commit()

    # print them to check __repr__
    print(e, end='\n\n')
    print(g, end='\n\n')
    print(gt1)
    print(gt2, end='\n\n')
    print(s, end='\n\n')
    print(st1)
    print(st2, end='\n\n')

    # let's make a Pointing
    user = db.get_user(session, username='goto_test')
    p = db.Pointing(object_name='IP Peg',
                    ra=350.785625,
                    dec=18.416472,
                    rank=9,
                    min_time=3600,
                    start_time=Time.now(),
                    stop_time=Time.now() + 3 * u.day,
                    user=user,
                    )
    e = db.ExposureSet(num_exp=20, exptime=20, filt='G')
    p.exposure_sets.append(e)
    session.add(p)
    session.commit()
    print(p, end='\n\n')

    # create an Mpointing
    mp = db.Mpointing(object_name='M31',
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
    mp.exposure_sets = [L, R, G, B]
    session.add(mp)
    session.commit()
    print(mp, end='\n\n')
    print(L)
    print(R)
    print(G)
    print(B, end='\n\n')

# new session, add random pointings
with db.open_session() as session:
    user = db.get_user(session, username='goto_test')

    pointings = [db.make_random_pointing(user) for i in range(10)]
    db.insert_items(session, pointings)

    more_pointings = [db.make_random_pointing(user) for i in range(10)]
    db.insert_items(session, more_pointings)

    # now check how many we have
    n_points = len(db.get_pointings(session))
    current, pending = db.get_queue(session)
    expired = db.get_expired_pointings(session)
    n_queue = len(pending)
    print('{} points in database'.format(n_points))
    print('{} points in queue'.format(n_queue))
    print('{} expired pointings\n'.format(len(expired)))

time.sleep(5)
print('-------')

# now let's go through a caretaker step
with db.open_session() as s:
    def print_summary(mp, now):
        """Print a summary of the database."""
        print(now, f'Mpointing is {mp.status.upper()} ({mp.num_remaining} remaining)')
        print('                    Pointings:', [p.status_at_time(now) for p in mp.pointings])
        # print('                    Time blocks:')
        # for block in mp.time_blocks:
        #     print('                    \t', (block.block_num, block.valid_time, block.wait_time,
        #                                      block.current))

    mp = db.get_mpointing_by_id(s, 1)
    print(mp)

    # lets pretend we're observing the Mpointing and check the triggers
    now = Time(mp.start_time)
    print_summary(mp, now)
    print()

    while True:
        # schedule next pointing
        next_pointing = mp.get_next_pointing()
        if next_pointing is None:
            break

        print(now, 'Creating Pointing')
        s.add(next_pointing)
        s.commit()

        time.sleep(1)
        print()

        # skip ahead until the pointing is pending
        now = Time(next_pointing.start_time) + 30 * u.s

        print(now, f'Marking Pointing {next_pointing.db_id} as running')
        next_pointing.mark_running(time=now)
        s.commit()
        print_summary(mp, now)
        time.sleep(1)
        print()

        # skip ahead by some exposure time
        now = now + 120 * u.s

        print(now, f'Marking Pointing {next_pointing.db_id} as completed')
        next_pointing.mark_finished(completed=True, time=now)
        s.commit()
        print_summary(mp, now)
        time.sleep(1)
        print()

    # check the list of mps to schedule is empty
    mps_to_schedule = db.get_mpointings(s, status='unscheduled')
    print('There are {} Mpointings to schedule'.format(len(mps_to_schedule)))

    n_points = len(db.get_pointings(s))
    current, pending = db.get_queue(s)
    expired = db.get_expired_pointings(s)
    marked = db.get_pointings(s, status='expired')
    n_queue = len(pending)

    print('{} points in database'.format(n_points))
    print('{} points in queue'.format(n_queue))
    print('{} unmarked expired pointings'.format(len(expired)))
    print('{} marked expired pointings\n'.format(len(marked)))
