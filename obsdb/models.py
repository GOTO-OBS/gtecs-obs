#!/usr/bin/env python
"""Python classes mapping on to database tables."""

import datetime

from astropy import units as u
from astropy.time import Time

from sqlalchemy import (Column, DateTime, Enum, Float, ForeignKey, Integer, String)
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship, validates


__all__ = ['User', 'Event', 'EventTile', 'Survey', 'SurveyTile', 'ExposureSet',
           'Pointing', 'TimeBlock', 'Mpointing', 'ImageLog']


Base = declarative_base()

pointing_status_list = Enum('pending', 'running', 'completed',
                            'aborted', 'interrupted', 'expired', 'deleted')

mpointing_status_list = Enum('unscheduled', 'scheduled', 'completed',
                             'aborted', 'expired', 'deleted')


class User(Base):
    """A class to represent a database User.

    Like all SQLAlchemy model classes, this object links to the
    underlying database. You can create an User, and set its attributes
    without a database session. Accessing some attributes may require
    an active database session, and some properties (like the db_id)
    will be None until the User is added to the database.

    The constructor must use keyword arguments and the arguments below
    should be supplied or set before insertion into the database.
    See Examples for details.

    Args
    ----
        username : string
            a short user name
        password : string
            password for authentication. stored as a hash in the database.
        full_name : string
            the user's full name.

    A User also has the following properties which are
    populated through database queries, but not needed for
    object creation:

    Attributes
    ----------
        db_id : int
            primary database key
        mpointings : list of `Mpointing`
            a list of any `Mpointing`s associated with this User
        pointings : list of `Pointing`
            a list of any `Pointing`s associated with this User

    Examples
    --------
        >>> from obsdb import *
        >>>
        >>> bob = User(username='bob', password='1234', full_name="Bob Marley")
        >>> session = load_session()
        >>> session.add(bob)  # add bob to database
        >>> session.commit()  # commit changes (otherwise DB not changed)
        >>> bob
        User(db_id=25, username=bob, full_name=Bob Marley)
        >>> bob.pointings
        []
        >>> pointing = make_random_pointing(25)  # make a pointing for bob
        >>> session.add(pointing)
        >>> session.commit()
        >>> bob.pointings
        [Pointing(db_id=None, status='pending', object_name=randObj, ra=352.133, dec=28.464,
        rank=84, min_alt=30, max_sunalt=-15, min_time=1575.8310236, max_moon=D, min_moonsep=30,
        too=True, start_time=2018-01-16 16:46:12, stop_time=2018-01-18 16:46:12, started_time=None,
        stopped_time=None, user_id=25, mpointing_id=None, time_block_id=None, event_id=None,
        event_tile_id=None, survey_id=None, survey_tile_id=None)]
        >>> session.close()

    """

    # Set corresponding SQL table name
    __tablename__ = 'users'

    # Primary key
    db_id = Column('id', Integer, primary_key=True)

    # Columns
    username = Column('username', String)
    password = Column('password', String)
    full_name = Column('full_name', String)

    def __repr__(self):
        return "User(db_id={}, username={}, full_name={})".format(
            self.db_id, self.username, self.full_name
        )


class Pointing(Base):
    """A class to represent an Pointing.

    Like all SQLAlchemy model classes, this object links to the
    underlying database. You can create an Pointing, and set its attributes
    without a database session. Accessing some attributes may require
    an active database session, and some properties (like the db_id)
    will be None until the Pointing is added to the database.

    The constructor must use keyword arguments and the arguments below
    should be supplied or set before insertion into the database.
    See Examples for details.

    Args
    ----
        user_id : int
            unique key identifying user to whom this Pointing belongs
        object_name : String
            object name
        ra : float, optional
            J2000 right ascension in decimal degrees
            if ra is not given and this Pointing is linked to a SurveyTile
            then the ra will be extracted from the SurveyTile
        dec : float, optional
            J2000 declination in decimal degrees
            if dec is not given and this Pointing is linked to a SurveyTile
            then the dec will be extracted from the SurveyTile
        rank : Integer
            rank to use for pointing
        min_alt : float
            minimum altitude to observer at
        min_time : float
            minimum time needed to schedule pointing
        max_sunalt : float
            altitude constraint on Sun
        max_moon : string
            Moon constraint. one of 'D', 'G', 'B'.
        min_moonsep : float
            distance constraint from the Moon, degrees
        too : int
            0 or 1 to indicate if this is a target of opportunity or not
        start_time : string, `astropy.time.Time` or datetime.datetime
            UTC time from which pointing is considered valid and can be started
        stop_time : string, `astropy.time.Time` or datetime.datetime, or None
            the latest UTC time at which pointing may be started
            can be None, if so the pointing will stay in the queue indefinitely
            (it can't be marked as expired) and will only leave when observed

        status : string, optional
            status of pointing, default 'pending'
        event_tile_id : int, optional
            unique key linking to a `EventTile`
        survey_tile_id : int, optional
            unique key linking to a `SurveyTile`
        event_id : int, optional
            unique key linking to an `Event`
        survey_id : int, optional
            unique key linking to a `Survey`
        mpointing_id : int, optional
            unique key linking to an `Mpointing`

    An Pointing also has the following properties which are
    populated through database queries, but not needed for
    object creation:

    Attributes
    ----------
        db_id : int
            primary database key
        started_time : datetime.datetime, or None
            if the pointing has started (been marked running)
            this will give the time it was updated
        stopped_time : datetime.datetime, or None
            if the pointing has finished (either completed or cancelled for
            some reason) this will give the time it was updated
        exposure_sets : list of `ExposureSet`
            the `ExposureSet` objects associated with this `Pointing`, if any
        event_tile : `EventTile`
            the `EventTile` associated with this `Pointing`, if any
        survey_tile : `SurveyTile`
            the `SurveyTile` associated with this `Pointing`, if any
        event : `Event`
            the `Event` associated with this `Pointing`, if any
        survey : `Survey`
            the `Survey` associated with this `Pointing`, if any
        mpointing : `Mpointing`
            the `Mpointing` associated with this `Pointing`, if any

    Examples
    --------
        >>> from obsdb import *
        >>> from astropy import units as u
        >>> from astropy.time import Time
        >>> session = load_session()

        Create a pointing:

        >>> p = Pointing(object_name='IP Peg', ra=350.785625, dec=18.416472, rank=9, min_alt=30,
        ... max_sunalt=-15, min_time=3600, max_moon='G', min_moonsep=30, too=0,
        ... start_time=Time.now(), stop_time=Time.now()+3*u.day, user_id=24)
        >>> p
        Pointing(db_id=None, status='None', object_name=IP Peg, ra=350.785625, dec=18.416472,
        rank=9, min_alt=30, max_sunalt=-15, min_time=3600, max_moon=G, min_moonsep=None, too=False,
        start_time=2018-01-12 17:39:54, stop_time=2018-01-15 17:39:54, started_time=None,
        stopped_time=None, user_id=24, mpointing_id=None, time_block_id=None, event_id=None,
        event_tile_id=None, survey_id=None, survey_tile_id=None)

        We can insert it into the database and the status and db_id will be set:

        >>> session.add(p)
        >>> session.commit()
        >>> p.status, p.db_id
        ('pending', 17073)

        At the moment, this pointing has no exposure sets. We can either add these to the
        `exposure_sets` attribute directly:

        >>> e1 = ExposureSet(imgtype='SCIENCE', filt='L', exptime=20, num_exp=20, binning=2)
        >>> p.exposure_sets.append(e1)

        or create `ExposureSet` instances with the `pointing_id` attribute set, and the database
        will take care of the rest:

        >>> e2 = ExposureSet(imgtype='SCIENCE', filt='G', exptime=20, num_exp=20, binning=2,
        ... pointing_id=17073)
        >>> e3 = ExposureSet(imgtype='SCIENCE', filt='R', exptime=20, num_exp=20, binning=2,
        ... pointing_id=17073)
        >>> insert_items(session, [e2, e3])
        >>> session.commit()
        >>> p.exposure_sets
        [ExposureSet(db_id=126601, ra_offset=0.0, dec_offset=0.0, imgtype=SCIENCE, filt=L,
        exptime=20.0, num_exp=20, binning=2, ut_mask=None, pointing_id=17073, mpointing_id=None),
        ExposureSet(db_id=126602, ra_offset=0.0, dec_offset=0.0, imgtype=SCIENCE, filt=G,
        exptime=20.0, num_exp=20, binning=2, ut_mask=None, pointing_id=17073, mpointing_id=None),
        ExposureSet(db_id=126603, ra_offset=0.0, dec_offset=0.0, imgtype=SCIENCE, filt=R,
        exptime=20.0, num_exp=20, binning=2, ut_mask=None, pointing_id=17073, mpointing_id=None)]

        >>> session.close()

    """

    # Set corresponding SQL table name
    __tablename__ = 'pointings'

    # Primary key
    db_id = Column('id', Integer, primary_key=True)

    # Columns
    status = Column('status', pointing_status_list, default='pending')
    object_name = Column('object', String)
    ra = Column('ra', Float)
    dec = Column('decl', Float)
    rank = Column('rank', Integer)
    min_alt = Column('min_alt', Float)
    max_sunalt = Column('max_sunalt', Float, default=-15)
    min_time = Column('min_time', Float)
    max_moon = Column('max_moon', String(1))
    min_moonsep = Column('min_moonsep', Float, default=30)
    too = Column('too', Integer)
    start_time = Column('start_time', DateTime)
    stop_time = Column('stop_time', DateTime)
    started_time = Column('started_time', DateTime, default=None)
    stopped_time = Column('stopped_time', DateTime, default=None)

    # Foreign keys
    user_id = Column('user_id', Integer, ForeignKey('users.id'), nullable=False)
    mpointing_id = Column('mpointing_id', Integer, ForeignKey('mpointings.id'), nullable=True)
    time_block_id = Column('time_block_id', Integer, ForeignKey('time_blocks.id'), nullable=True)
    event_id = Column('event_id', Integer, ForeignKey('events.id'), nullable=True)
    event_tile_id = Column('event_tile_id', Integer, ForeignKey('event_tiles.id'), nullable=True)
    survey_id = Column('survey_id', Integer, ForeignKey('surveys.id'), nullable=True)
    survey_tile_id = Column('survey_tile_id', Integer, ForeignKey('survey_tiles.id'), nullable=True)

    # Foreign relationships
    user = relationship('User', backref='pointings', uselist=False)
    mpointing = relationship('Mpointing', backref='pointings', uselist=False)
    time_block = relationship('TimeBlock', back_populates='pointings', uselist=False)
    event = relationship('Event', backref='pointings')
    event_tile = relationship('EventTile', back_populates='pointings', uselist=False)
    survey = relationship('Survey', backref='pointings')
    survey_tile = relationship('SurveyTile', back_populates='pointings', uselist=False)

    def __repr__(self):
        template = ("Pointing(db_id={}, status='{}', " +
                    "object_name={}, ra={}, dec={}, rank={}, " +
                    "min_alt={}, max_sunalt={}, min_time={}, max_moon={}, min_moonsep={}, " +
                    "too={}, start_time={}, stop_time={}, started_time={}, stopped_time={}, " +
                    "user_id={}, mpointing_id={}, time_block_id={}, " +
                    "event_id={}, event_tile_id={}, survey_id={}, survey_tile_id={})")
        return template.format(
            self.db_id, self.status, self.object_name, self.ra, self.dec, self.rank,
            self.min_alt, self.max_sunalt, self.min_time, self.max_moon, self.min_moonsep,
            bool(self.too), self.start_time, self.stop_time, self.started_time, self.stopped_time,
            self.user_id, self.mpointing_id, self.time_block_id,
            self.event_id, self.event_tile_id, self.survey_id, self.survey_tile_id
        )

    @validates('start_time', 'stop_time')
    def munge_times(self, key, field):
        """Use validators to allow various types of input for UTC.

        Also enforce write_time > start_time.
        """
        if key == 'stop_time' and field is None:
            value = None
        elif isinstance(field, datetime.datetime):
            value = field.strftime("%Y-%m-%d %H:%M:%S")
        elif isinstance(field, Time):
            field.precision = 0  # no D.P on seconds
            value = field.iso
        else:
            # just hope the string works!
            value = str(field)

        if (key == 'start_time' and self.stop_time is not None):
            if Time(value) >= Time(self.stop_time):
                raise AssertionError("stop_time must be later than start_time")
        elif key == 'stop_time' and value is not None and self.start_time is not None:
            if Time(self.start_time) >= Time(value):
                raise AssertionError("stop_time must be later than start_time")

        return value


class ExposureSet(Base):
    """A class to represent an Exposure Set: a set of repeated identical exposures.

    Like all SQLAlchemy model classes, this object links to the
    underlying database. You can create an ExposureSet, and set its attributes
    without a database session. Accessing some attributes may require
    an active database session, and some properties (like the db_id)
    will be None until the ExposureSet is added to the database.

    The constructor must use keyword arguments and the arguments below
    should be supplied or set before insertion into the database.
    See Examples for details.

    Args
    ----
        num_exp : int
            number of exposures within the set
        exptime : float
            exposure time in seconds
        filt : string
            filter to use
        binning : int
            binning to apply
        imgtype : string
            indicates the type of exposure set.
            one of SCIENCE, FOCUS, STD, FLAT, BIAS, DARK

        ut_mask : int, optional
            if set, this is a binary mask which will determine which unit
            telescopes carry out the exposure. A value of 5 (binary 0101) will
            be exposed by cameras 1 and 3.
        ra_offset : float, optional
            the size of the random offset to apply between each exposure
            if not set, no offset will be made
        dec_offset : float, optional
            the size of the random offset to apply between each exposure
            if not set, no offset will be made

        mpointing_id : int, optional
            unique key linking to an `Mpointing`
        pointing_id : int, optional
            unique key linking to an `Pointing`

    An ExposureSet also has the following properties which are
    populated through database queries, but not needed for
    object creation:

    Attributes
    ----------
        db_id : int
            primary database key
        mpointing : `Mpointing`
            the `Mpointing` associated with this `ExposureSet`, if any
        pointing : `Pointing`
            the `Pointing` associated with this `ExposureSet`, if any

    """

    # Set corresponding SQL table name
    __tablename__ = 'exposure_sets'

    # Primary key
    db_id = Column('id', Integer, primary_key=True)

    # Columns
    num_exp = Column('num_exp', Integer)
    exptime = Column('exptime', Float)
    filt = Column('filter', String(2))
    binning = Column('binning', Integer)
    imgtype = Column('imgtype', Enum('SCIENCE', 'FOCUS', 'DARK', 'BIAS', 'FLAT', 'STD'))
    ut_mask = Column('ut_mask', Integer, nullable=True)
    ra_offset = Column('ra_offset', Float, server_default='0.0')
    dec_offset = Column('dec_offset', Float, server_default='0.0')

    # Foreign keys
    pointing_id = Column('pointing_id', Integer, ForeignKey('pointings.id'), nullable=False)
    mpointing_id = Column('mpointing_id', Integer, ForeignKey('mpointings.id'), nullable=False)

    # Foreign relationships
    pointing = relationship('Pointing', backref='exposure_sets', uselist=False)
    mpointing = relationship('Mpointing', backref='exposure_sets', uselist=False)

    def __repr__(self):
        template = ("ExposureSet(db_id={}, ra_offset={}, dec_offset={}, imgtype={}, " +
                    "filt={}, exptime={}, num_exp={}, binning={}, ut_mask={}, " +
                    "pointing_id={}, mpointing_id={})")
        return template.format(
            self.db_id, self.ra_offset, self.dec_offset, self.imgtype, self.filt,
            self.exptime, self.num_exp, self.binning, self.ut_mask,
            self.pointing_id, self.mpointing_id
        )


class Mpointing(Base):
    """A class to represent an Mpointing.

    Like all SQLAlchemy model classes, this object links to the
    underlying database. You can create an Mpointing, and set its attributes
    without a database session. Accessing some attributes may require
    an active database session, and some properties (like the db_id)
    will be None until the Mpointing is added to the database.

    The constructor must use keyword arguments and the arguments below
    should be supplied or set before insertion into the database.
    See Examples for details.

    Args
    ----
        user_id : int
            unique key identifying user to whom this Mpointing belongs
        object_name : String
            object name
        ra : float, optional
            J2000 right ascension in decimal degrees
            if ra is not given and this Mpointing is linked to a SurveyTile
            then the ra will be extracted from the SurveyTile
        dec : float, optional
            J2000 declination in decimal degrees
            if dec is not given and this Mpointing is linked to a SurveyTile
            then the dec will be extracted from the SurveyTile
        start_rank : Integer
            rank to use for first pointing in series
        min_alt : float
            minimum altitude to observer at
        min_time : float
            minimum time needed to schedule pointing
        max_sunalt : float
            altitude constraint on Sun
        max_moon : string
            Moon constraint. one of 'D', 'G', 'B'.
        min_moonsep : float
            distance constraint from the Moon, degrees
        too : int
            0 or 1 to indicate if this is a too or not
        num_todo : int
            number of (sucsessful) observations required.
            less than zero means repeat infinitely.
        valid_time : float or list of float
            the amount of time the pointing(s) should be valid in the queue.
            if num_todo is greater than times given the list will be looped.
        wait_time : float or list of float
            time to wait between pointings in minutes.
            if num_todo is greater than times given the list will be looped.

        status : string, optional
            status of mpointing, default 'unscheduled'
        start_time : string, `astropy.time.Time` or datetime.datetime, optional
            UTC time from which Mpointing is considered valid and can be started
            if not given then set to now, so the Mpointing will start immediately
        stop_time : string, `astropy.time.Time` or datetime.datetime, optional
            the latest UTC time after which pointings must stop
            if not given the Mpointing will continue creating pointings until
            it is completed
        event_tile_id : int, optional
            unique key linking to `EventTile`
        survey_tile_id : int, optional
            unique key linking to `SurveyTile`
        event_id : int, optional
            unique key linking to `Event`
        survey_id : int, optional
            unique key linking to `Survey`

    An Mpointing also has the following properties which are
    populated through database queries, but not needed for
    object creation:

    Attributes
    ----------
        db_id : int
            primary database key
        rank : int
            rank for next pointing to be scheduled
        num_completed : int
            number of successfully completed pointings
        num_remaining : int
            number of pointings still to do (same as num_todo - num_completed)
        exposure_sets : list of `ExposureSet`
            the `ExposureSet` objects associated with this `Mpointing`, if any
        time_blocks : list of `TimeBlock`
            the `TimeBlock` objects associated with this `Mpointing`, if any
        survey_tile : `SurveyTile`
            the `SurveyTile` associated with this `Mpointing`, if any
        event_tile : `EventTile`
            the `EventTile` associated with this `Mpointing`, if any
        event : `Event`
            the `Event` associated with this `Mpointing`, if any
        survey : `Survey`
            the `Survey` associated with this `Mpointing`, if any

    Examples
    --------
        >>> from obsdb import *
        >>> from astropy.time import Time

        make an Mpointing - starting at midnight, 5 pointings that stay in the queue
        for 5 minutes and the next one is valid 10 minutes after the previous.

        >>> mp = Mpointing(object_name='M31', ra=22, dec=-5, start_rank=9, min_alt=30,
        ... min_time=3600, max_sunalt=-15, too=0, max_moon='B', min_moonsep=30, num_todo=5,
        ... user_id=24, valid_time=5, wait_time=10, start_time=Time('2018-01-01 00:00:00'))
        >>> mp
        Mpointing(db_id=None, status='unscheduled', num_todo=5, num_completed=0,
        num_remaining=5, infinite=False, object_name=M31, ra=22, dec=-5, rank=9, start_rank=9,
        min_alt=30, max_sunalt=-15, min_time=3600, max_moon=B, min_moonsep=30, too=False,
        start_time=2018-01-01 00:00:00, stop_time=None, user_id=24, event_id=None,
        event_tile_id=None, survey_id=None, survey_tile_id=None)

        Note that the db_id is None, because that will be filled out when we add it to the
        database.

        Looking at the TimeBlocks you can see that we only need one, and it has been generated

        >>> mp.time_blocks
        [TimeBlock(db_id=None, block_num=1, valid_time=5, wait_time=10, current=1,
        mpointing_id=None)]

        This block will be repeated 5 times, as requested.

        ~~~~~~~~~~~~~~~~~~~~~

        For a more complicated example, give a list to wait_time to have the intervals between
        pointings increase.

        >>> mp = Mpointing(object_name='M31', ra=22, dec=-5, start_rank=9, min_alt=30,
        ... min_time=3600, max_sunalt=-15, too=0, max_moon='B', min_moonsep=30, num_todo=5,
        ... user_id=24, valid_time=5, wait_time=[10,20,30], start_time=Time('2018-01-01 00:00:00'))
        >>> mp
        Mpointing(db_id=None, status='unscheduled', num_todo=5, num_completed=0,
        num_remaining=5, infinite=False, object_name=M31, ra=22, dec=-5, rank=9, start_rank=9,
        min_alt=30, max_sunalt=-15, min_time=3600, max_moon=B, min_moonsep=30, too=False,
        start_time=2018-01-01 00:00:00, stop_time=None, user_id=24, event_id=None,
        event_tile_id=None, survey_id=None, survey_tile_id=None)
        >>> mp.time_blocks
        [TimeBlock(db_id=None, block_num=1, valid_time=5, wait_time=10, current=1,
        mpointing_id=None),
        TimeBlock(db_id=None, block_num=2, valid_time=5, wait_time=20, current=None,
        mpointing_id=None),
        TimeBlock(db_id=None, block_num=3, valid_time=5, wait_time=30, current=None,
        mpointing_id=None)]

        Note this time that we only created 3 blocks, but are asking for 5 observations.
        That's not a problem, as the blocks will simply repeat.
        In this case the sequence will look like this:
        [note we set the Mpointing start time to midnight]
        Pointing 1: start_time=00:00, stop_time=00:05 (valid for 5 minutes)
        Pointing 2: start_time=00:15, stop_time=00:20 (wait for 10, then valid for 5)
        Pointing 3: start_time=00:40, stop_time=00:45 (wait for 20, then valid for 5)
        Pointing 4: start_time=01:15, stop_time=01:20 (wait for 30, then valid for 5)
        Pointing 5: start_time=01:30, stop_time=01:35 (wait for 10, then valid for 5)

        We can see what would happen by manually running through the pointings and pretending to be
        the caretaker.

        First add the Mpointing to the database:

        >>> s = load_session()
        >>> s.add(mp)
        >>> s.commit()

        Then create the first pointing and add it to the database.

        >>> p = mp.get_next_pointing()
        >>> s.add(p)
        >>> s.commit()
        >>> mp.pointings
        [Pointing(db_id=17100, status='pending', object_name=M31, ra=22.0, dec=-5.0, rank=9,
        min_alt=30.0, max_sunalt=-15.0, min_time=3600.0, max_moon=B, min_moonsep=30.0, too=False,
        start_time=2018-01-01 00:00:00, stop_time=2018-01-01 00:05:00, started_time=None,
        stopped_time=None, user_id=24, mpointing_id=2, time_block_id=1, event_id=None,
        event_tile_id=None, survey_id=None, survey_tile_id=None)]

        Note that the db_id, mpointing_id and time_block_id have been filled out,
        and this Pointing has time_block_id=1.
        Also that the start time is midnight and the stop time is 5 past, as expected.
        See what happens if we mark the pointing as completed:

        >>> mp.pointings[0].status = 'completed'
        >>> s.commit()
        >>> mp
        Mpointing(db_id=2, status='unscheduled', num_todo=5, num_completed=1, num_remaining=4,
        infinite=False, object_name=M31, ra=22.0, dec=-5.0, rank=19, start_rank=9, min_alt=30.0,
        max_sunalt=-15.0, min_time=3600.0, max_moon=B, min_moonsep=30.0, too=False,
        start_time=2018-01-01 00:00:00, stop_time=None, user_id=24, event_id=None,
        event_tile_id=None, survey_id=None, survey_tile_id=None)

        See that the num_completed atribute has gone up, the num_remaining has gone down
        and the mpointing status has changed to 'unscheduled'.

        Create the next pointing and add it:

        >>> p = mp.get_next_pointing()
        >>> s.add(p)
        >>> s.commit()
        >>> p
        Pointing(db_id=17102, status='pending', object_name=M31, ra=22.0, dec=-5.0, rank=19,
        min_alt=30.0, max_sunalt=-15.0, min_time=3600.0, max_moon=B, min_moonsep=30.0, too=False,
        start_time=2018-01-01 00:15:00, stop_time=2018-01-01 00:20:00, started_time=None,
        stopped_time=None, user_id=24, mpointing_id=2, time_block_id=2, event_id=None,
        event_tile_id=None, survey_id=None, survey_tile_id=None)
        >>> mp
        Mpointing(db_id=2, status='scheduled', num_todo=5, num_completed=1, num_remaining=4,
        infinite=False, object_name=M31, ra=22.0, dec=-5.0, rank=19, start_rank=9, min_alt=30.0,
        max_sunalt=-15.0, min_time=3600.0, max_moon=B, min_moonsep=30.0, too=False,
        start_time=2018-01-01 00:00:00, stop_time=None, user_id=24, event_id=None,
        event_tile_id=None, survey_id=None, survey_tile_id=None)

        See that the Mpointing is back to scheduled.
        Also, note that the new pointing has start time of 00:15 and stop of 00:20. That's as we
        expected, because it's linked to the second time block not the first
        (see time_block_id=2).

        If you look at the time blocks you can see the next one is now marked as current.

        >>> mp.time_blocks
        [TimeBlock(db_id=1, block_num=1, valid_time=5.0, wait_time=10.0, current=0,
        mpointing_id=2),
        TimeBlock(db_id=2, block_num=2, valid_time=5.0, wait_time=20.0, current=1,
        mpointing_id=2),
        TimeBlock(db_id=3, block_num=3, valid_time=5.0, wait_time=30.0, current=0,
        mpointing_id=2)]

        Mark this one as completed, and you'll see the Mpointing is updated.

        >>> p.status = 'completed'
        >>> s.commit()
        >>> mp.num_completed
        2
        >>> mp.status
        'unscheduled'

        Let's run through the remaining pointings:

        >>> p = mp.get_next_pointing() # Pointing 3 of 5
        >>> s.add(p)
        >>> s.commit()
        >>> p.time_block_id
        3
        >>> mp.status
        'scheduled'

        >>> p.status = 'completed'
        >>> s.commit()
        >>> mp.num_completed
        3
        >>> mp.status
        'unscheduled'

        >>> p = mp.get_next_pointing() # Pointing 4 of 5
        >>> s.add(p)
        >>> s.commit()
        >>> p.time_block_id
        1
        >>> mp.status
        'scheduled'

        >>> p.status = 'completed'
        >>> s.commit()
        >>> mp.num_completed
        4
        >>> mp.status
        'unscheduled'

        >>> p = mp.get_next_pointing() # Pointing 5 of 5
        >>> s.add(p)
        >>> s.commit()
        >>> p.time_block_id
        2
        >>> mp.status
        'scheduled'

        >>> p.status = 'completed'
        >>> s.commit()
        >>> mp.num_completed
        5
        >>> mp.status
        'completed'
        >>> mp
        Mpointing(db_id=2, status='completed', num_todo=5, num_completed=5, num_remaining=0,
        infinite=False, object_name=M31, ra=22.0, dec=-5.0, rank=59, start_rank=9, min_alt=30.0,
        max_sunalt=-15.0, min_time=3600.0, max_moon=B, min_moonsep=30.0, too=False,
        start_time=2018-01-01 00:00:00, stop_time=None, user_id=24, event_id=None,
        event_tile_id=None, survey_id=None, survey_tile_id=None)

        And the Mpointing is completed.

        ~~~~~~~~~~~~~~~~~~~~~

        To be useful, an Mpointing should have a list of `ExposureSet`s associated with it. We
        can either add these to the `exposure_sets` attribute directly:

        >>> e1 = ExposureSet(imgtype='SCIENCE', filt='L', exptime=20, num_exp=20, binning=2)
        >>> mp.exposure_sets.append(e1)

        or create `ExposureSet` instances with the `mpointing_id` attribute set, and the database
        will take care of the rest:

        >>> e2 = ExposureSet(imgtype='SCIENCE', filt='G', exptime=20, num_exp=20, binning=2,
        ... mpointing_id=1)
        >>> e3 = ExposureSet(imgtype='SCIENCE', filt='R', exptime=20, num_exp=20, binning=2,
        ... mpointing_id=1)
        >>> insert_items(session, [e2, e3])
        >>> session.commit()
        >>> mp.exposure_sets
        [ExposureSet(db_id=126598, ra_offset=0.0, dec_offset=0.0, imgtype=SCIENCE, filt=L,
        exptime=20.0, num_exp=20, binning=2, ut_mask=None, pointing_id=None, mpointing_id=1),
        ExposureSet(db_id=126599, ra_offset=0.0, dec_offset=0.0, imgtype=SCIENCE, filt=G,
        exptime=20.0, num_exp=20, binning=2, ut_mask=None, pointing_id=None, mpointing_id=1),
        ExposureSet(db_id=126600, ra_offset=0.0, dec_offset=0.0, imgtype=SCIENCE, filt=R,
        exptime=20.0, num_exp=20, binning=2, ut_mask=None, pointing_id=None, mpointing_id=1)]

        These exposure sets will be copied to the Pointings when they're created.

        >>> session.close()

    """

    # Set corresponding SQL table name
    __tablename__ = 'mpointings'

    # Primary key
    db_id = Column('id', Integer, primary_key=True)

    # Columns
    status = Column('status', mpointing_status_list, default='unscheduled')
    object_name = Column('object', String)
    ra = Column('ra', Float)
    dec = Column('decl', Float)
    rank = Column('rank', Integer)
    start_rank = Column('start_rank', Integer)
    num_todo = Column('num_todo', Integer)
    num_completed = Column('num_completed', Integer)
    infinite = Column('infinite', Integer, default=False)
    min_alt = Column('min_alt', Float)
    max_sunalt = Column('max_sunalt', Float)
    min_time = Column('min_time', Float)
    max_moon = Column('max_moon', String(1))
    min_moonsep = Column('min_moonsep', Float)
    too = Column('too', Integer)
    start_time = Column('start_time', DateTime)
    stop_time = Column('stop_time', DateTime)

    # Foreign keys
    user_id = Column('user_id', Integer, ForeignKey('users.id'), nullable=False)
    event_id = Column('event_id', Integer, ForeignKey('events.id'), nullable=True)
    event_tile_id = Column('event_tile_id', Integer, ForeignKey('event_tiles.id'), nullable=True)
    survey_id = Column('survey_id', Integer, ForeignKey('surveys.id'), nullable=True)
    survey_tile_id = Column('survey_tile_id', Integer, ForeignKey('survey_tiles.id'), nullable=True)

    # Foreign relationships
    user = relationship('User', backref='mpointings', uselist=False)
    event = relationship('Event', backref='mpointings')
    survey = relationship('Survey', backref='mpointings')
    event_tile = relationship('EventTile', back_populates='mpointing', uselist=False)
    survey_tile = relationship('SurveyTile', back_populates='mpointing', uselist=False)
    time_blocks = relationship('TimeBlock', back_populates='mpointing', viewonly=True)

    def __repr__(self):
        template = ("Mpointing(db_id={}, status='{}', num_todo={}, num_completed={}," +
                    "num_remaining={}, infinite={}, " +
                    "object_name={}, ra={}, dec={}, rank={}, start_rank={}, " +
                    "min_alt={}, max_sunalt={}, min_time={}, max_moon={}, min_moonsep={}, " +
                    "too={}, start_time={}, stop_time={}, " +
                    "user_id={}, event_id={}, event_tile_id={}, survey_id={}, survey_tile_id={})")
        return template.format(
            self.db_id, self.status, self.num_todo, self.num_completed,
            self.num_remaining, bool(self.infinite),
            self.object_name, self.ra, self.dec, self.rank, self.start_rank,
            self.min_alt, self.max_sunalt, self.min_time, self.max_moon, self.min_moonsep,
            bool(self.too), self.start_time, self.stop_time,
            self.user_id, self.event_id, self.event_tile_id, self.survey_id, self.survey_tile_id
        )

    def __init__(self, object_name=None, ra=None, dec=None,
                 start_rank=None, min_alt=None, min_time=None,
                 max_moon=None, min_moonsep=None, max_sunalt=None, too=None, start_time=None,
                 stop_time=None, num_todo=None, valid_time=None, wait_time=None,
                 status='unscheduled', **kwargs):
        self.ra = ra
        self.dec = dec
        self.object_name = object_name
        self.start_rank = start_rank
        self.max_moon = max_moon
        self.min_moonsep = min_moonsep
        self.min_alt = min_alt
        self.min_time = min_time
        self.max_sunalt = max_sunalt
        self.too = too
        self.start_time = start_time if start_time is not None else Time.now()
        self.stop_time = stop_time
        self.rank = self.start_rank
        self.status = status
        self.infinite = False
        self.num_todo = num_todo
        self.num_completed = 0

        # now add time blocks
        if valid_time is not None and wait_time is not None:

            # first convert to lists
            if not isinstance(valid_time, list):
                valid_times = [valid_time]
            else:
                valid_times = valid_time
            if not isinstance(wait_time, list):
                wait_times = [wait_time]
            else:
                wait_times = wait_time

            # check if infinite
            if num_todo < 0:
                self.infinite = True
                self.num_todo = -1

            # create TimeBlock objects
            for i in range(max(len(valid_times), len(wait_times))):
                valid = valid_times[i % len(valid_times)]
                wait = wait_times[i % len(wait_times)]

                # check if non-expiring
                if valid < 0:
                    valid = -1

                block = TimeBlock(block_num=i + 1, valid_time=valid, wait_time=wait)
                self.time_blocks.append(block)
            if len(self.time_blocks):
                self.time_blocks[0].current = 1

        if 'user' in kwargs:
            self.user = kwargs['user']
        if 'user_id' in kwargs:
            self.user_id = kwargs['user_id']
        if 'survey_id' in kwargs:
            self.survey_id = kwargs['survey_id']
        if 'survey' in kwargs:
            self.survey = kwargs['survey']
        if 'survey_tile' in kwargs:
            self.survey_tile = kwargs['survey_tile']
        if 'survey_tile_id' in kwargs:
            self.survey_tile_id = kwargs['survey_tile_id']
        if 'event' in kwargs:
            self.event = kwargs['event']
        if 'event_id' in kwargs:
            self.event_id = kwargs['event_id']
        if 'event_tile' in kwargs:
            self.event_tile = kwargs['event_tile']
        if 'event_tile_id' in kwargs:
            self.event_tile_id = kwargs['event_tile_id']

    @validates('start_time', 'stop_time')
    def munge_times(self, key, field):
        """Use validators to allow various types of input for UTC.

        Also enforce stop_time > start_time
        NB stop_time is None be default
        """
        if key == 'stop_time' and field is None:
            value = None
        elif isinstance(field, datetime.datetime):
            value = field.strftime("%Y-%m-%d %H:%M:%S")
        elif isinstance(field, Time):
            field.precision = 0  # no D.P on seconds
            value = field.iso
        else:
            # just hope the string works!
            value = str(field)

        if (key == 'start_time' and self.stop_time is not None):
            if Time(value) >= Time(self.stop_time):
                raise AssertionError("stop_time must be later than start_time")
        elif key == 'stop_time' and value is not None and self.start_time is not None:
            if Time(self.start_time) >= Time(value):
                raise AssertionError("stop_time must be later than start_time")

        return value

    @property
    def num_remaining(self):
        """Return the number of observations remaining.

        Returns
        -------
        num_remaining : int
            defined as num_todo - num_completed
            if Mpointing.infinite == True returns -1

        """
        if self.infinite:
            return -1
        else:
            return self.num_todo - self.num_completed

    def get_current_block(self):
        """Return the current time block.

        Assumes this object is still associated to an active session.

        Returns
        -------
        current : `obsdb.TimeBlock`
            The current time block (may be None).

        """
        current_block = [block for block in self.time_blocks if block.current]
        if len(current_block) == 0:
            return None
        elif len(current_block) > 1:
            raise ValueError('Multiple time blocks marked as current')
        else:
            return current_block[0]

    def get_last_block(self):
        """Return the last time block executed.

        Assumes this object is still associated to an active session.

        Returns
        -------
        last : `obsdb.TimeBlock`
            The last block done (may be None).

        """
        current_block = self.get_current_block()
        if not current_block:
            return None

        current_num = current_block.block_num
        if current_num == 1:
            last_num = len(self.time_blocks)
        else:
            last_num = current_num - 1

        last_block = [block for block in self.time_blocks
                      if block.block_num == last_num][0]
        return last_block

    def get_next_block(self):
        """Return the next time block to be executed.

        Assumes this object is still associated to an active session.

        Returns
        -------
        next : `obsdb.TimeBlock`
            The next block to do after the current one (may be None).

        """
        current_block = self.get_current_block()
        if not current_block and self.num_completed > 0:
            # no current block (for some reason, aside from just starting)
            return None

        elif not current_block or len(current_block.pointings) == 0:
            # just starting
            next_num = 1

        elif self.num_remaining == 0:
            # current block is the last one
            return None

        else:
            current_num = current_block.block_num
            latest_pointing = current_block.pointings[-1]
            if latest_pointing.status in ['completed', 'expired']:
                next_num = (current_num % len(self.time_blocks)) + 1
            else:
                next_num = current_num

        next_block = [block for block in self.time_blocks
                      if block.block_num == next_num][0]
        return next_block

    def get_next_pointing(self):
        """Retrieve the next pointing which needs to be scheduled.

        The start and stop times of the pointing are determined from
        the status of the previous pointing.

        Assumes this object is still associated with an active session.

        Returns
        -------
        pointing : `Pointing`
            the next pointing that should be sent to the database. Is None
            if there is no suitable pointing remaining, or if there is
            a pointing from this Mpointing already scheduled

        """
        # already scheduled or finished, return None
        if self.status != 'unscheduled':
            return None

        # As a safety check, see if the mpointing has any other pointings that
        # are already pending or running.
        # If so the mpointing status should be 'scheduled' and therefore be
        # caught above, but strange things can happen.
        # This should prevent the cases where Mpointings have multiple
        # pointings in the queue, hopefully!
        if len(self.pointings) > 0:
            statuses = [p.status for p in self.pointings]
            if 'pending' in statuses or 'running' in statuses:
                return None

        # all the observations have been completed
        if self.num_remaining == 0:
            return None

        current_block = self.get_current_block()
        next_block = self.get_next_block()

        # no current block or pointings, should only happen the first time
        if not current_block or len(current_block.pointings) == 0:
            start_time = self.start_time
        else:
            latest_pointing = current_block.pointings[-1]
            # decide to go onto the next block:
            #  - if the last block's pointing was completed
            #  - if the last block's pointing's valid time has expired
            if latest_pointing.status in ['completed', 'expired']:
                if latest_pointing.stop_time:
                    start_time = Time(latest_pointing.stop_time) + current_block.wait_time * u.min
                else:
                    # non-expiring pointings have no stop_time
                    start_time = Time.now() + current_block.wait_time * u.minute
            else:
                # the current block wasn't completed, and there's still time left
                # (e.g. aborted, interrupted)
                # need to re-insert the current block with a new pointing
                start_time = latest_pointing.start_time

        if next_block.valid_time < 0:
            # non-expiring pointings
            stop_time = None
        else:
            stop_time = Time(start_time) + next_block.valid_time * u.minute
            if self.stop_time and stop_time > self.stop_time:
                # force pointings to stop by an Mpointing's stop_time, if given
                stop_time = self.stop_time

            if start_time >= stop_time:
                # can happen if the Mpointing has a stop_time
                return None

        # now create a pointing
        p = Pointing(object_name=self.object_name,
                     ra=self.ra,
                     dec=self.dec,
                     rank=self.rank,
                     min_alt=self.min_alt,
                     max_sunalt=self.max_sunalt,
                     min_time=self.min_time,
                     max_moon=self.max_moon,
                     min_moonsep=self.min_moonsep,
                     too=self.too,
                     start_time=start_time,
                     stop_time=stop_time,
                     status='pending',
                     user_id=self.user_id,
                     mpointing_id=self.db_id,
                     time_block=next_block,
                     event_id=self.event_id,
                     event_tile_id=self.event_tile_id,
                     survey_id=self.survey_id,
                     survey_tile_id=self.survey_tile_id,
                     )
        # add the exposures
        p.exposure_sets = self.exposure_sets
        return p


class TimeBlock(Base):
    """A class to represent a block of time.

    Specifically, these represent a time period made up of a valid time and a wait time.
    The valid time is the time that the pointing will be valid for, and the wait time will be the
    time after the pointing is observed/invalid to wait until the next valid period.

    Like all SQLAlchemy model classes, this object links to the
    underlying database. You can create a TimeBlock, and set its attributes
    without a database session. Accessing some attributes may require
    an active database session, and some properties (like the db_id)
    will be None until the TimeBlock is added to the database.

    The constructor must use keyword arguments and the arguments below
    should be supplied or set before insertion into the database.
    See Examples for details.

    Note you shouldn't ever have to manually create TimeBlocks.
    They're used internally by Mpointings to keep track of Pointings and all
    the ones an Mpointing requires are created when the Mpointing is.

    Args
    ----
        block_num : int
            an integer indicating which block in a sequence this is
        valid_time : float
            amount of time a pointing in this block should stay valid in the queue, in minutes.
        wait_time : float
            time to wait after this block before allowing the next pointing, in minutes

        mpointing_id : int, optional
            unique key linking this TimeBlock to an Mpointing
        current : bool, optional
            True if this TimeBlock is the one that is currently linked to
            a Pointing in the queue

    An TimeBlock also has the following properties which are
    populated through database queries, but not needed for
    object creation:

    Attributes
    ----------
        pointings : list of `Pointing`
            a list of any `Pointing`s associated with this block, if any

    Examples
    --------
        >>> from obsdb import *

        make an TimeBlock
        (an Mpointing with mpointing_id = 1 is already in the database)

        >>> b = TimeBlock(block_num=1, valid_time=60, wait_time=120, mpointing_id=1)
        >>> session = load_session()
        >>> session.add(b)
        TimeBlock(db_id=7, block_num=1, valid_time=60.0, wait_time=120.0, current=False,
        mpointing_id=1)

    """

    # Set corresponding SQL table name
    __tablename__ = 'time_blocks'

    # Primary key
    db_id = Column('id', Integer, primary_key=True)

    # Columns
    block_num = Column('block_num', Integer)
    valid_time = Column('valid_time', Integer)
    wait_time = Column('wait_time', Integer)
    current = Column('current', Integer, default=False)

    # Foreign keys
    mpointing_id = Column('mpointing_id', Integer, ForeignKey('mpointings.id'), nullable=False)

    # Foreign relationships
    mpointing = relationship('Mpointing', back_populates='time_blocks', uselist=False)
    pointings = relationship('Pointing', back_populates='time_block')

    def __repr__(self):
        template = ("TimeBlock(db_id={}, block_num={}, valid_time={}, " +
                    "wait_time={}, current={}, mpointing_id={})")
        return template.format(self.db_id, self.block_num, self.valid_time,
                               self.wait_time, bool(self.current), self.mpointing_id)


class Event(Base):
    """A class to represent a transient Event.

    Like all SQLAlchemy model classes, these objects link to the
    underlying database. You can create an object, set its attributes
    without a database session. Accessing some attributes may require
    an active database session, and some properties (like the db_id)
    will be None until the Event is added to the database.

    The constructor must use keyword arguments and the arguments below
    should be supplied or set before insertion into the database.
    See Examples for details.

    Args
    ----
        ivo : string
            ivorn for the event
        name : string
            a human-readable identifier for the event
        source : string
            the event's origin, e.g. LVC
        skymap : string, optional
            the location of the source skymap file

    An Event also has the following properties which are
    populated through database queries, but not needed for
    object creation:

    Attributes
    ----------
        db_id : int
            primary database key
        event_tiles : list of `EventTile`
            a list of any `EventTiles` associated with this Event
        mpointings : list of `Mpointing`
            a list of any `Mpointing`s associated with this Event
        pointings : list of `Pointing`
            a list of any `Pointing`s associated with this Event

    Examples
    --------
        >>> e = Event(ivo='ivo://pt5mTest', name='pt5mVar3', source='pt5m')

    """

    # Set corresponding SQL table name
    __tablename__ = 'events'

    # Primary key
    db_id = Column('id', Integer, primary_key=True)

    # Columns
    name = Column('name', String)
    source = Column('source', String)
    ivo = Column('ivorn', String, unique=True)
    skymap = Column('skymap', String)

    # Foreign relationships
    event_tiles = relationship('EventTile', back_populates='event')

    def __repr__(self):
        return "Event(db_id={}, name={}, source={}, ivo={}, skymap={})".format(
            self.db_id, self.name, self.source, self.ivo, self.skymap
        )


class EventTile(Base):
    """A class to represent a Tile from an Event (e.g. a LVC skymap).

    Like all SQLAlchemy model classes, this object links to the
    underlying database. You can create an EventTile, and set its attributes
    without a database session. Accessing some attributes may require
    an active database session, and some properties (like the db_id)
    will be None until the User is added to the database.

    The constructor must use keyword arguments and the arguments below
    should be supplied or set before insertion into the database.
    See Examples for details.

    Args
    ----
        probability : float
            contained target probability within this tile
        ra : float, optional
            J2000 right ascension in decimal degrees
            if ra is not given and this EventTile is linked to a SurveyTile
            then the ra will be extracted from the SurveyTile
        dec : float, optional
            J2000 declination in decimal degrees
            if dec is not given and this EventTile is linked to a SurveyTile
            then the dec will be extracted from the SurveyTile
        event_id : int, optional
            the ID number of the Event this tile is associated with
        survey_tile_id : int, optional
            the SurveyTile this tile is associated with, if any

    An EventTile also has the following properties which are
    populated through database queries, but not needed for
    object creation:

    Attributes
    ----------
        db_id : int
            primary database key
        unobserved_probability : float
            the total probability in this tile that hasn't been observed
            should be updated whenever any overlapping tile is observed
        mpointing : `Mpointing`
            the `Mpointing` associated with this EventTile if any
        pointings : list of `Pointing`
            a list of any `Pointing`s associated with this tile, if any

    Examples
    --------
        >>> from obsdb import *

        make a LIGO event to associate our tile with

        >>> e = Event(ivo='ivo://GW150914', name='GW150914', source='LVC')
        >>> session = load_session()
        >>> session.add(e)  # add event
        >>> session.commit()  # commit changes (otherwise DB not changed)
        >>> e.db_id
        1

        construct without event_id

        >>> event_tile = EventTile(ra=122.34, dec=22.01, probability=0.01)
        >>> event_tile
        EventTile(db_id=None, ra=122.34, dec=22.01, probability=0.01, event_id=None,
        survey_tile_id=None)

        set the event_id

        >>> event_tile.event_id = 1
        >>> event_tile
        EventTile(db_id=None, ra=122.34, dec=22.01, probability=0.01, event_id=None,
        survey_tile_id=None)

        add to the database

        >>> session.add(event_tile)
        >>> session.commit()
        >>> event_tile  # note how db_id is populated now tile is in DB
        EventTile(db_id=1, ra=122.34, dec=22.01, probability=0.01, event_id=1, survey_tile_id=None)
        >>> e.event_tiles  # and event knows about all associated tiles
        [EventTile(db_id=1, ra=122.34, dec=22.01, probability=0.01, event_id=1,
        survey_tile_id=None)]
        >>> session.close()

        make a survey and a survey tile, and add them to the database

        >>> s = Survey(name='GOTO survey')
        >>> st = SurveyTile(ra=22, dec=-2, name='Tile1')
        >>> st.survey = s
        >>> session.add(s, st)
        >>> session.commit()

        make a new event tile that has no ra and dec,
        but is linked to the survey tile

        >>> et2 = EventTile(probability=0.01)
        >>> et2.event = e
        >>> et2.survey_tile = st
        >>> et2
        EventTile(db_id=None, ra=None, dec=None, probability=0.01, event_id=None,
        survey_tile_id=None)

        add to the database

        >>> session.add(event_tile)
        >>> session.commit()
        >>> et2  # note how the ra and dec have been copied from the SurveyTile
        EventTile(db_id=2, ra=22.0, dec=-2.0, probability=0.01, event_id=1, survey_tile_id=1)
        >>> st.event_tiles # and the SurveyTile knows about the linked EventTiles
        [EventTile(db_id=2, ra=22.0, dec=-2.0, probability=0.01, event_id=1, survey_tile_id=1)]

    """

    # Set corresponding SQL table name
    __tablename__ = 'event_tiles'

    # Primary key
    db_id = Column('id', Integer, primary_key=True)

    # Columns
    ra = Column('ra', Float)
    dec = Column('decl', Float)
    probability = Column('probability', Float)
    unobserved_probability = Column('unobserved_probability', Float)

    # Foreign keys
    event_id = Column('event_id', Integer, ForeignKey('events.id'), nullable=False)
    survey_tile_id = Column('survey_tile_id', Integer, ForeignKey('survey_tiles.id'), nullable=True)

    # Foreign relationships
    pointings = relationship('Pointing', back_populates='event_tile')
    mpointing = relationship('Mpointing', back_populates='event_tile')
    event = relationship('Event', back_populates='event_tiles', uselist=False)
    survey_tile = relationship('SurveyTile', back_populates='event_tiles', uselist=False)

    def __repr__(self):
        template = ("EventTile(db_id={}, ra={}, dec={}, " +
                    "probability={}, event_id={}, survey_tile_id={})")
        return template.format(
            self.db_id, self.ra, self.dec, self.probability, self.event_id,
            self.survey_tile_id)


class Survey(Base):
    """A class to represent an observation survey.

    Like all SQLAlchemy model classes, these objects link to the
    underlying database. You can create an object, set its attributes
    without a database session. Accessing some attributes may require
    an active database session, and some properties (like the db_id)
    will be None until the Survey is added to the database.

    The constructor must use keyword arguments and the arguments below
    should be supplied or set before insertion into the database.
    See Examples for details.

    Args
    ----
        name : string
            a human-readable identifier for the survey

    A Survey also has the following properties which are
    populated through database queries, but not needed for
    object creation:

    Attributes
    ----------
        db_id : int
            primary database key
        survey_tiles : list of `SurveyTile`
            a list of any `SurveyTiles` associated with this Survey
        mpointings : list of `Mpointing`
            a list of any `Mpointing`s associated with this Survey
        pointings : list of `Pointing`
            a list of any `Pointing`s associated with this Survey

    Examples
    --------
        >>> s = Survey(name='GOTO4-allsky')

    """

    # Set corresponding SQL table name
    __tablename__ = 'surveys'

    # Primary key
    db_id = Column('id', Integer, primary_key=True)

    # Columns
    name = Column('name', String)

    # Foreign relationships
    survey_tiles = relationship('SurveyTile', back_populates='survey')

    def __repr__(self):
        return "Survey(db_id={}, name={})".format(
            self.db_id, self.name
        )


class SurveyTile(Base):
    """A class to represent a Survey Tile.

    Like all SQLAlchemy model classes, this object links to the
    underlying database. You can create an SurveyTile, and set its attributes
    without a database session. Accessing some attributes may require
    an active database session, and some properties (like the db_id)
    will be None until the User is added to the database.

    The constructor must use keyword arguments and the arguments below
    should be supplied or set before insertion into the database.
    See Examples for details.

    Args
    ----
        ra : float
            J2000 right ascension in decimal degrees
        dec : float
            J2000 declination in decimal degrees
        survey_id : int, optional
            the ID number of the Survey this tile is associated with
        name : str, optional
            a human-readable identifier for the tile

    A SurveyTile also has the following properties which are
    populated through database queries, but not needed for
    object creation:

    Attributes
    ----------
        db_id : int
            primary database key
        event_tiles : list of `EventTile`
            a list of any `EventTiles` associated with this SurveyTiles, if any
        mpointing : `Mpointing`
            the `Mpointing` associated with this SurveyTile if any
        pointings : list of `Pointing`
            a list of any `Pointing`s associated with this tile, if any

    Examples
    --------
        >>> from obsdb import *

        make a survey to associate our tile with

        >>> s = Survey(name='GOTO4-allsky')
        >>> session = load_session()
        >>> session.add(s)  # add survey
        >>> session.commit()  # commit changes (otherwise DB not changed)
        >>> s.db_id
        1

        construct without survey_id

        >>> tile = SurveyTile(ra=100, dec=20)
        >>> tile
        SurveyTile(db_id=None, ra=100.00, dec=20.00)

        set the survey_id

        >>> event_tile.event_id = 1
        >>> event_tile
        EventTile(db_id=None, ra=122.34, dec=22.01, probability=0.01, event_id=None)

        add to the database, demonstrating the open_session context manager,
        which will handle committing changes for you:

        >>> with open_session as session():
        >>>     insert_items([tile])
        >>>     tile  # note how db_id is populated now tile is in DB
        SurveyTile(db_id=1, ra=100.00, dec=20.00, survey_id=1)
        >>> with open_session as session():
        >>>     s.survey_tiles  # and survey knows about all associated tiles
        [SurveyTile(db_id=1, ra=100.00, dec=20.00, survey_id=1)]

    """

    # Set corresponding SQL table name
    __tablename__ = 'survey_tiles'

    # Primary key
    db_id = Column('id', Integer, primary_key=True)

    # Columns
    name = Column('name', String)
    ra = Column('ra', Float)
    dec = Column('decl', Float)

    # Foreign keys
    survey_id = Column('survey_id', Integer, ForeignKey('surveys.id'), nullable=False)

    # Foreign relationships
    survey = relationship('Survey', back_populates='survey_tiles', uselist=False)
    event_tiles = relationship('EventTile', back_populates='survey_tile')
    mpointing = relationship('Mpointing', back_populates='survey_tile')
    pointings = relationship('Pointing', back_populates='survey_tile')

    def __repr__(self):
        template = ("SurveyTile(db_id={}, ra={}, dec={}, " +
                    ", name={}, survey_id={})")
        return template.format(
            self.db_id, self.ra, self.dec, self.name, self.survey_id)


class ImageLog(Base):
    """A class to store a record of a FITS image file created by the camera daemon.

    The ImageLog is a simple way to link pointings in the database to physical
    FITS files. An ImageLog should be created by the camera daemon each time
    an exposure is taken, even if taken manually rather than originating from
    the database.

    Like all SQLAlchemy model classes, this object links to the
    underlying database. You can create an ImageLog, and set its attributes
    without a database session. Accessing some attributes may require
    an active database session, and some properties (like the db_id)
    will be None until the log is added to the database.

    The constructor must use keyword arguments and the arguments below
    should be supplied or set before insertion into the database.
    See Examples for details.

    Args
    ----
        filename : string
            the name of the image file
        run_number : int
            the run ID number for this exposure
        ut : int
            the unit telescope this frame was captured on
        ut_mask : int
            a binary mask for which unit telescopes carried out this exposure.
            A value of 5 (binary 0101) will have been exposed by cameras 1 and 3.
        start_time : string, `astropy.time.Time` or datetime.datetime
            the time that the exposure began
        write_time : string, `astropy.time.Time` or datetime.datetime
            the time that the image file was written

        set_position : int, optional
            position of this exposure in a set, if it's in one
            if not, it will default to 1
        set_total : int, optional
            total number of exposures in this set, if any
            if not given, it will default to 1
        exposure_set_id : int, optional
            the unique key of the exposure set this frame was part of
        pointing_id : int, optional
            the unique key of the pointing which generated this frame
        mpointing_id : int, optional, optional
            unique key of the mpointing which generated the pointing which
            generated this frame

    An ImageLog also has the following properties which are
    populated through database queries, but not needed for
    object creation:

    Attributes
    ----------
        db_id : int
            primary database key
        exposure_set : `ExposureSet`
            the `ExposureSet` object associated with this `ImageLog`, if any
        pointing : `Pointing`
            the `Pointing` associated with this `ImageLog`, if any
        mpointing : `Mpointing`
            the `Mpointing` associated with this `ImageLog`, if any

    """

    # Set corresponding SQL table name
    __tablename__ = 'image_logs'

    # Primary key
    db_id = Column('id', Integer, primary_key=True)

    # Columns
    filename = Column('filename', String)
    run_number = Column('run_number', Integer)
    ut = Column('ut', Integer)
    ut_mask = Column('ut_mask', Integer)
    start_time = Column('start_time', DateTime)
    write_time = Column('write_time', DateTime)
    set_position = Column('set_position', Integer, default=1)
    set_total = Column('set_total', Integer, default=1)

    # Foreign keys
    exposure_set_id = Column('exposure_set_id', Integer, ForeignKey('exposure_sets.id'), nullable=True)
    pointing_id = Column('pointing_id', Integer, ForeignKey('pointings.id'), nullable=True)
    mpointing_id = Column('mpointing_id', Integer, ForeignKey('mpointings.id'), nullable=True)

    # Foreign relationships
    exposure_set = relationship('ExposureSet', backref='image_logs', uselist=False)
    pointing = relationship('Pointing', backref='image_logs', uselist=False)
    mpointing = relationship('Mpointing', backref='image_logs', uselist=False)

    def __repr__(self):
        template = ("ImageLog(db_id={}, filename={}, run_number={}, " +
                    "ut={}, ut_mask={}, start_time={}, write_time={}, " +
                    "exposure_set_id={}, pointing_id={}, mpointing_id={})")
        return template.format(
            self.db_id, self.filename, self.run_number, self.ut, self.ut_mask,
            self.start_time, self.write_time, self.exposure_set_id, self.pointing_id,
            self.mpointing_id
        )

    @validates('start_time', 'write_time')
    def munge_times(self, key, field):
        """Use validators to allow various types of input for UTC.

        Also enforce write_time > start_time.
        """
        if isinstance(field, datetime.datetime):
            value = field.strftime("%Y-%m-%d %H:%M:%S")
        elif isinstance(field, Time):
            field.precision = 0  # no D.P on seconds
            value = field.iso
        else:
            # just hope the string works!
            value = str(field)

        if (key == 'start_time' and self.write_time is not None):
            if Time(value) >= Time(self.write_time):
                raise AssertionError("write_time must be later than start_time")
        elif key == 'write_time' and self.write_time is not None:
            if Time(self.start_time) >= Time(value):
                raise AssertionError("write_time must be later than start_time")

        return value
