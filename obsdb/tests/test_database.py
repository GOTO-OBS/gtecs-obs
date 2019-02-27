"""Demonstrate the use of the database modules.

Assume we start from a clean database.
"""

import sys
import time

from astropy import units as u
from astropy.time import Time

import obsdb as db


# add a user
with db.open_session() as session:
    # check if the database is empty
    users = session.query(db.User).all()
    if users:
        print('Error: Database is not empty')
        sys.exit(1)

    # create a test user
    user = db.User(username='goto_test', password='password', full_name='GOTO Test Observer')
    print(user)
    session.add(user)
    session.commit()
    print(user, end='\n\n')

    # check the password was stored correctly
    assert db.validate_user(session, username='goto_test', password='password')

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

    # check the survey tile weight updated correctly
    assert st1.initial_weight == st1.current_weight == 0.5

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
                    min_alt=30,
                    max_sunalt=-15,
                    min_time=3600,
                    max_moon='G',
                    min_moonsep=30,
                    too=0,
                    start_time=Time.now(),
                    stop_time=Time.now() + 3 * u.day,
                    user=user,
                    )
    e = db.ExposureSet(imgtype='SCIENCE',
                       filt='G',
                       exptime=20,
                       num_exp=20,
                       binning=2)
    p.exposure_sets.append(e)
    session.add(p)
    session.commit()
    print(p, end='\n\n')

    # create an Mpointing
    mp = db.Mpointing(object_name='M31',
                      ra=10.685,
                      dec=41.2875,
                      start_rank=9,
                      min_alt=30,
                      max_sunalt=-15,
                      min_time=3600,
                      too=0,
                      max_moon='B',
                      min_moonsep=30,
                      num_todo=5,
                      valid_time=[60, 120],
                      wait_time=60,
                      user=user,
                      )
    # and add RGBL exposure set
    L = db.ExposureSet(imgtype='SCIENCE', filt='L', exptime=120, num_exp=3, binning=2)
    R = db.ExposureSet(imgtype='SCIENCE', filt='R', exptime=120, num_exp=3, binning=2)
    G = db.ExposureSet(imgtype='SCIENCE', filt='G', exptime=120, num_exp=3, binning=2)
    B = db.ExposureSet(imgtype='SCIENCE', filt='B', exptime=120, num_exp=3, binning=2)
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
    npoints = len(db.get_pointings(session))
    current, pending = db.get_queue(session)
    expired = db.get_expired_pointings(session)
    nqueue = len(pending)
    print("{} points in database".format(npoints))
    print("{} points in queue".format(nqueue))
    print("{} expired pointings\n".format(len(expired)))

    # clean expired pointings
    print('Cleaning expired pointings')
    db.bulk_update_status(session, expired, 'expired')

    # now check how many we have
    npoints = len(db.get_pointings(session))
    current, pending = db.get_queue(session)
    expired = db.get_expired_pointings(session)
    marked = db.get_pointings(session, status='expired')
    nqueue = len(pending)
    print("{} points in database".format(npoints))
    print("{} points in queue".format(nqueue))
    print("{} unmarked expired pointings".format(len(expired)))
    print("{} marked expired pointings\n".format(len(marked)))

time.sleep(5)

# now let's go through a caretaker step
with db.open_session() as s:
    # first remove expired pointings
    expired = db.get_expired_pointings(s)
    if len(expired) > 0:
        db.bulk_update_status(s, expired, 'expired')
    print('Marked {} pointings as expired'.format(len(expired)))

    # which Mpointings need pointings submitting?
    mps_to_schedule = db.get_mpointings(s, status='unscheduled')
    print('There are {} Mpointings to schedule'.format(len(mps_to_schedule)))

    pointings_to_add = []
    for mp in mps_to_schedule:
        print('Scheduling mp #{}'.format(mp.db_id))
        next_pointing = mp.get_next_pointing()
        if next_pointing:
            pointings_to_add.append(next_pointing)
    db.insert_items(s, pointings_to_add)

    def summary(mp):
        """Print a summary of the database."""
        print('{}, num_remaining = {}'.format(mp.status.upper(), mp.num_remaining))
        print([p.status for p in mp.pointings])
        print([(b.block_num, b.valid_time, b.wait_time, b.current) for b in mp.time_blocks], '\n')

    # check that scheduling worked the way we expected
    mp = db.get_mpointing_by_id(s, 1)
    summary(mp)

    # now lets pretend we're observing the Mpointing and check the triggers
    while True:
        print('marking job as running')
        mp.pointings[-1].status = 'running'
        s.commit()
        time.sleep(1)
        summary(mp)

        print('marking as completed')
        mp.pointings[-1].status = 'completed'
        s.commit()
        time.sleep(1)
        summary(mp)

        # schedule next
        next_pointing = mp.get_next_pointing()
        if next_pointing:
            s.add(next_pointing)
        else:
            break
        s.commit()
        time.sleep(1)
        summary(mp)

    # check the list of mps to schedule is empty
    mps_to_schedule = db.get_mpointings(s, status='unscheduled')
    print('There are {} Mpointings to schedule'.format(len(mps_to_schedule)))

    npoints = len(db.get_pointings(s))
    current, pending = db.get_queue(s)
    expired = db.get_expired_pointings(s)
    marked = db.get_pointings(s, status='expired')
    nqueue = len(pending)

    print("{} points in database".format(npoints))
    print("{} points in queue".format(nqueue))
    print("{} unmarked expired pointings".format(len(expired)))
    print("{} marked expired pointings\n".format(len(marked)))
