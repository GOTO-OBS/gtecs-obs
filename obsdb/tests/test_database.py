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
    try:
        user_id = db.get_user_id(session, 'goto')
        if user_id:
            print('Error: Database is not empty')
            sys.exit()
    except ValueError:
        pass

    db.add_user(session, 'goto', 'password', "GOTO Test Observer")

    # create an event
    e = db.Event(ivo='ivo://gotoTest', name='gotoVar3', source='made-up')

    # and an event tile
    et = db.EventTile(ra=100, dec=20, probability=0.1)
    et.event = e

    # and a survey
    s = db.Survey(name='GOTO survey')

    # and a survey tile
    st = db.SurveyTile(ra=22, dec=-2, name='Tile1')
    st.survey = s

    # and an event tile linked to that survey tile
    et2 = db.EventTile(probability=0.2)
    et2.event = e
    et2.survey_tile = st

    # add them
    db.insert_items(session, [e, et, s, st, et2])
    session.commit()

    # check the second event tile updated correctly
    assert et2.ra == st.ra
    assert et2.dec == st.dec

# OK, new Session. Let's make a Pointing
with db.open_session() as session:
    user_id = db.get_user_id(session, 'goto')
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
                    user_id=user_id)
    e = db.ExposureSet(imgtype='SCIENCE',
                       filt='G',
                       exptime=20,
                       num_exp=20,
                       binning=2)
    p.exposure_sets.append(e)
    session.add(p)

# new session, add random pointings
with db.open_session() as session:
    pointings = [db.make_random_pointing(user_id) for i in range(10)]
    db.insert_items(session, pointings)

    more_pointings = [db.make_random_pointing(user_id) for i in range(10)]
    db.insert_items(session, more_pointings)

# now check how many we have
with db.open_session() as session:
    npoints = len(db.get_pointings(session))
    current, pending = db.get_queue(session)
    expired = db.get_expired_pointing_ids(session)
    nqueue = len(pending)
print("{} points in database".format(npoints))
print("{} points in queue".format(nqueue))
print("{} expired pointings\n".format(len(expired)))

# clean expired pointings
print('Cleaning expired pointings')
with db.open_session() as session:
    db.bulk_update_pointing_status(session, expired, 'expired')

# now check how many we have
with db.open_session() as session:
    npoints = len(db.get_pointings(session))
    current, pending = db.get_queue(session)
    expired = db.get_expired_pointing_ids(session)
    marked = db.get_pointings(session, status='expired')
    nqueue = len(pending)
print("{} points in database".format(npoints))
print("{} points in queue".format(nqueue))
print("{} unmarked expired pointings".format(len(expired)))
print("{} marked expired pointings\n".format(len(marked)))

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
                  user_id=user_id,
                  valid_time=[60, 120],
                  wait_time=60)
# and add RGBL exposure set
L = db.ExposureSet(imgtype='SCIENCE', filt='L', exptime=120, num_exp=3, binning=2)
R = db.ExposureSet(imgtype='SCIENCE', filt='R', exptime=120, num_exp=3, binning=2)
G = db.ExposureSet(imgtype='SCIENCE', filt='G', exptime=120, num_exp=3, binning=2)
B = db.ExposureSet(imgtype='SCIENCE', filt='B', exptime=120, num_exp=3, binning=2)
mp.exposure_sets = [L, R, G, B]
with db.open_session() as s:
    s.add(mp)

time.sleep(5)

# now let's go through a caretaker step
with db.open_session() as s:
    # first remove expired pointings
    expired = db.get_expired_pointing_ids(s)
    if len(expired) > 0:
        db.bulk_update_pointing_status(s, expired, 'expired')
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


# check that scheduling worked the way we expected
def summary(mp):
    """Print a summary of the database."""
    print('{}, num_remaining = {}'.format(mp.status.upper(), mp.num_remaining))
    print([p.status for p in mp.pointings])
    print([(b.block_num, b.valid_time, b.wait_time, b.current) for b in mp.observing_blocks], '\n')


with db.open_session() as s:
    mp = db.get_mpointing_by_id(s, 1)
    summary(mp)

# now lets pretend we're observing the Mpointing and check the triggers
s = db.load_session()
mp = db.get_mpointing_by_id(s, 1)
keepGoing = True
while keepGoing:
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
        keepGoing = False
    s.commit()
    time.sleep(1)
    summary(mp)
s.close()

with db.open_session() as s:
    # check the list of mps to schedule is empty
    mps_to_schedule = db.get_mpointings(s, status='unscheduled')
    print('There are {} Mpointings to schedule'.format(len(mps_to_schedule)))

    npoints = len(db.get_pointings(s))
    current, pending = db.get_queue(s)
    expired = db.get_expired_pointing_ids(s)
    marked = db.get_pointings(s, status='expired')
    nqueue = len(pending)
print("{} points in database".format(npoints))
print("{} points in queue".format(nqueue))
print("{} unmarked expired pointings".format(len(expired)))
print("{} marked expired pointings\n".format(len(marked)))
