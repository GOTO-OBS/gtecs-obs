#!/usr/bin/env python
"""Python classes mapping on to database tables."""

import datetime
import hashlib

from astropy import units as u
from astropy.time import Time

from sqlalchemy import (Boolean, Column, DateTime, Enum, Float, ForeignKey, Integer, String)
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship, validates


__all__ = ['User', 'Pointing', 'ExposureSet', 'Mpointing', 'TimeBlock',
           'Grid', 'GridTile', 'Survey', 'SurveyTile', 'Event', 'ImageLog']


Base = declarative_base()

pointing_status_list = Enum('pending', 'running', 'completed',
                            'aborted', 'interrupted', 'expired', 'deleted')

mpointing_status_list = Enum('unscheduled', 'scheduled', 'completed',
                             'aborted', 'expired', 'deleted')


class User(Base):
    """A class to represent a database User.

    Every Pointing in the database must be associated with a User.

    Like all SQLAlchemy model classes, this object links to the
    underlying database. You can create an instance, and set its attributes
    without a database session. Accessing some attributes may require
    an active database session, and some properties (like the db_id)
    will be None until the instance is added to the database.

    Parameters
    ----------
    username : string
        a short user name
    password : string
        password for authentication
        the password is stored as a hash in the database, so the unshased string isn't stored
    full_name : string
        the user's full name.

    Attributes
    ----------
    db_id : int
        primary database key
        only populated when the instance is added to the database
    password_hash : str
        the hashed password as stored in the database

    When created the instance can be linked to the following other tables,
    otherwise they are populated when it is added to the database:

    Relationships
    -------------
    mpointings : list of `Mpointing`, optional
        the Mpointings associated with this User, if any
    pointings : list of `Pointing`, optional
        the Pointings associated with this User, if any

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
    [Pointing(db_id=None, status='pending', object_name=random_object, ra=352.133, dec=28.464,
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
    username = Column(String)
    password_hash = Column('password', String)  # Don't allow access to the unhashed password
    full_name = Column(String)

    # Foreign relationships
    pointings = relationship('Pointing', back_populates='user', uselist=True)
    mpointings = relationship('Mpointing', back_populates='user', uselist=True)

    def __init__(self, username=None, password=None, full_name=None):
        self.username = username
        self.password_hash = hashlib.sha512(password.encode()).hexdigest()
        self.full_name = full_name

    def __repr__(self):
        strings = ['db_id={}'.format(self.db_id),
                   'username={}'.format(self.username),
                   'full_name={}'.format(self.full_name),
                   ]
        return 'User({})'.format(', '.join(strings))


class Pointing(Base):
    """A class to represent an Pointing.

    Pointings are the primary table of the Observation Database.
    When decideing what to observe all the Pointings are processed
    based on their status and constraints.

    Like all SQLAlchemy model classes, this object links to the
    underlying database. You can create an instance, and set its attributes
    without a database session. Accessing some attributes may require
    an active database session, and some properties (like the db_id)
    will be None until the instance is added to the database.

    Parameters
    ----------
    object_name : String
        object name
    ra : float, optional
        J2000 right ascension in decimal degrees
        if ra is not given and this Pointing is linked to a GridTile
        then the ra will be extracted from the GridTile
    dec : float, optional
        J2000 declination in decimal degrees
        if dec is not given and this Pointing is linked to a GridTile
        then the dec will be extracted from the GridTile
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
    too : bool
        indicates if this is a Target of Opportunity (ToO)
    start_time : string, `astropy.time.Time` or datetime.datetime
        UTC time from which pointing is considered valid and can be started
    stop_time : string, `astropy.time.Time` or datetime.datetime, or None
        the latest UTC time at which pointing may be started
        can be None, if so the pointing will stay in the queue indefinitely
        (it can't be marked as expired) and will only leave when observed

    status : string, optional
        status of pointing, default 'pending'

    Attributes
    ----------
    db_id : int
        primary database key
        only populated when the instance is added to the database
    started_time : datetime.datetime, or None
        if the pointing has started (been marked running)
        this will give the time it was updated
    stopped_time : datetime.datetime, or None
        if the pointing has finished (either completed or cancelled for
        some reason) this will give the time it was updated

    When created the instance can be linked to the following other tables,
    otherwise they are populated when it is added to the database:

    Relationships
    -------------
    user : `User`
        the User associated with this Pointing
        required before addition to the database
        can also be added with the user_id parameter

    exposure_sets : list of `ExposureSet`, optional
        the Exposure Sets associated with this Pointing, if any
    mpointing : `Mpointing`, optional
        the Mpointing associated with this Pointing, if any
        can also be added with the mpointing_id parameter
    grid_tile : `GridTile`, optional
        the GridTile associated with this Pointing, if any
        can also be added with the grid_tile_id parameter
    survey_tile : `SurveyTile`, optional
        the SurveyTile associated with this Pointing, if any
        can also be added with the survey_tile_id parameter
    event : `Event`, optional
        the Event associated with this Pointing, if any
        can also be added with the event_id parameter

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
    status = Column(pointing_status_list, default='pending')
    object_name = Column('object', String)  # object is a built in class in Python
    ra = Column(Float)
    dec = Column('decl', Float)  # dec is reserved in SQL so can't be a column name
    rank = Column(Integer)
    min_alt = Column(Float)
    max_sunalt = Column(Float, default=-15)
    min_time = Column(Float)
    max_moon = Column(String(1))
    min_moonsep = Column(Float, default=30)
    too = Column(Boolean, default=False)
    start_time = Column(DateTime)
    stop_time = Column(DateTime)
    started_time = Column(DateTime, default=None)
    stopped_time = Column(DateTime, default=None)

    # Foreign keys
    user_id = Column(Integer, ForeignKey('users.id'), nullable=False)
    mpointing_id = Column(Integer, ForeignKey('mpointings.id'), nullable=True)
    time_block_id = Column(Integer, ForeignKey('time_blocks.id'), nullable=True)
    grid_tile_id = Column(Integer, ForeignKey('grid_tiles.id'), nullable=True)
    survey_tile_id = Column(Integer, ForeignKey('survey_tiles.id'), nullable=True)
    event_id = Column(Integer, ForeignKey('events.id'), nullable=True)

    # Foreign relationships
    user = relationship('User', back_populates='pointings', uselist=False)
    exposure_sets = relationship('ExposureSet', back_populates='pointing', uselist=True)
    mpointing = relationship('Mpointing', back_populates='pointings', uselist=False)
    time_block = relationship('TimeBlock', back_populates='pointings', uselist=False)
    grid_tile = relationship('GridTile', back_populates='pointings', uselist=False)
    survey_tile = relationship('SurveyTile', back_populates='pointings', uselist=False)
    event = relationship('Event', back_populates='pointings', uselist=False)

    def __repr__(self):
        strings = ['db_id={}'.format(self.db_id),
                   'status={}'.format(self.status),
                   'object_name={}'.format(self.object_name),
                   'ra={}'.format(self.ra),
                   'dec={}'.format(self.dec),
                   'rank={}'.format(self.rank),
                   'min_alt={}'.format(self.min_alt),
                   'max_sunalt={}'.format(self.max_sunalt),
                   'min_time={}'.format(self.min_time),
                   'max_moon={}'.format(self.max_moon),
                   'min_moonsep={}'.format(self.min_moonsep),
                   'too={}'.format(self.too),
                   'start_time={}'.format(self.start_time),
                   'stop_time={}'.format(self.stop_time),
                   'started_time={}'.format(self.started_time),
                   'stopped_time={}'.format(self.stopped_time),
                   'user_id={}'.format(self.user_id),
                   'mpointing_id={}'.format(self.mpointing_id),
                   'time_block_id={}'.format(self.time_block_id),
                   'grid_tile_id={}'.format(self.grid_tile_id),
                   'survey_tile_id={}'.format(self.survey_tile_id),
                   'event_id={}'.format(self.event_id),
                   ]
        return 'Pointing({})'.format(', '.join(strings))

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
    """A class to represent an Exposure Set.

    An Exposure Set is a set of repeated identical exposures, with the same
    exposure time, filter, binning factor etc. Each Pointing should have one or more
    Exposure Sets.

    Like all SQLAlchemy model classes, this object links to the
    underlying database. You can create an instance, and set its attributes
    without a database session. Accessing some attributes may require
    an active database session, and some properties (like the db_id)
    will be None until the instance is added to the database.

    Parameters
    ----------
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

    Attributes
    ----------
    db_id : int
        primary database key
        only populated when the instance is added to the database

    When created the instance can be linked to the following other tables,
    otherwise they are populated when it is added to the database:

    Relationships
    -------------
    pointing : `Pointing`, optional
        the Pointing associated with this Exposure Set, if any
        can also be added with the pointing_id parameter
    mpointing : `Mpointing`, optional
        the Mpointing associated with this Exposure Set, if any
        can also be added with the mpointing_id parameter

    """

    # Set corresponding SQL table name
    __tablename__ = 'exposure_sets'

    # Primary key
    db_id = Column('id', Integer, primary_key=True)

    # Columns
    num_exp = Column(Integer)
    exptime = Column(Float)
    filt = Column('filter', String(2))  # filter is a built in function in Python
    binning = Column(Integer)
    imgtype = Column(Enum('SCIENCE', 'FOCUS', 'DARK', 'BIAS', 'FLAT', 'STD'))
    ut_mask = Column(Integer, nullable=True)
    ra_offset = Column(Float, server_default='0.0')
    dec_offset = Column(Float, server_default='0.0')

    # Foreign keys
    pointing_id = Column(Integer, ForeignKey('pointings.id'), nullable=False)
    mpointing_id = Column(Integer, ForeignKey('mpointings.id'), nullable=False)

    # Foreign relationships
    pointing = relationship('Pointing', back_populates='exposure_sets', uselist=False)
    mpointing = relationship('Mpointing', back_populates='exposure_sets', uselist=False)

    def __repr__(self):
        strings = ['db_id={}'.format(self.db_id),
                   'num_exp={}'.format(self.num_exp),
                   'exptime={}'.format(self.exptime),
                   'filt={}'.format(self.filt),
                   'binning={}'.format(self.binning),
                   'imgtype={}'.format(self.imgtype),
                   'ut_mask={}'.format(self.ut_mask),
                   'ra_offset={}'.format(self.ra_offset),
                   'dec_offset={}'.format(self.dec_offset),
                   ]
        return 'ExposureSet({})'.format(', '.join(strings))


class Mpointing(Base):
    """A class to represent an Mpointing.

    Mpointings are generation engines for Pointings. Most of the parameters
    are passed directly on to the generated Pointings. Mpointings allow the
    database to regenerate Pointings at given time intervals.

    Like all SQLAlchemy model classes, this object links to the
    underlying database. You can create an instance, and set its attributes
    without a database session. Accessing some attributes may require
    an active database session, and some properties (like the db_id)
    will be None until the instance is added to the database.

    Parameters
    ----------
    object_name : String
        object name
    ra : float, optional
        J2000 right ascension in decimal degrees
        if ra is not given and this Mpointing is linked to a GridTile
        then the ra will be extracted from the GridTile
    dec : float, optional
        J2000 declination in decimal degrees
        if dec is not given and this Mpointing is linked to a GridTile
        then the dec will be extracted from the GridTile
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
    too : bool
        indicates if this is a Target of Opportunity (ToO)
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

    Attributes
    ----------
    db_id : int
        primary database key
        only populated when the instance is added to the database
    current_rank : int
        rank for next pointing to be scheduled
    initial_rank : int
        initial rank set (it will increase as pointings are observed)
    num_completed : int
        number of successfully completed pointings
    num_remaining : int
        number of pointings still to do (same as num_todo - num_completed)
    infinite : bool
        if the Mpointing will continue infnitely (set if num_todo is < 0)

    When created the instance can be linked to the following other tables,
    otherwise they are populated when it is added to the database:

    Relationships
    -------------
    user : `User`
        the User associated with this Mpointing
        required before addition to the database
        can also be added with the user_id parameter

    pointings : list of `Pointing`, optional
        the Pointings associated with this Mpointing, if any
    exposure_sets : list of `ExposureSet`, optional
        the Exposure Sets associated with this Mpointing, if any
    time_blocks : list of `TimeBlock`, optional
        the Time Blocks associated with this Mpointing, if any
    grid_tile : `GridTile`, optional
        the GridTile associated with this Mpointing, if any
        can also be added with the grid_tile_id parameter
    survey_tile : `SurveyTile`, optional
        the SurveyTile associated with this Mpointing, if any
        can also be added with the survey_tile_id parameter
    event : `Event`, optional
        the Event associated with this Mpointing, if any
        can also be added with the event_id parameter

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
    status = Column(mpointing_status_list, default='unscheduled')
    object_name = Column('object', String)  # object is a built in class in Python
    ra = Column(Float)
    dec = Column('decl', Float)  # dec is reserved in SQL so can't be a column name
    current_rank = Column(Integer)
    initial_rank = Column(Integer)
    num_todo = Column(Integer)
    num_completed = Column(Integer)
    infinite = Column(Boolean, default=False)
    min_alt = Column(Float)
    max_sunalt = Column(Float)
    min_time = Column(Float)
    max_moon = Column(String(1))
    min_moonsep = Column(Float)
    too = Column(Boolean, default=False)
    start_time = Column(DateTime)
    stop_time = Column(DateTime)

    # Foreign keys
    user_id = Column(Integer, ForeignKey('users.id'), nullable=False)
    grid_tile_id = Column(Integer, ForeignKey('grid_tiles.id'), nullable=True)
    survey_tile_id = Column(Integer, ForeignKey('survey_tiles.id'), nullable=True)
    event_id = Column(Integer, ForeignKey('events.id'), nullable=True)

    # Foreign relationships
    user = relationship('User', back_populates='mpointings', uselist=False)
    pointings = relationship('Pointing', back_populates='mpointing', uselist=True)
    exposure_sets = relationship('ExposureSet', back_populates='mpointing', uselist=True)
    time_blocks = relationship('TimeBlock', back_populates='mpointing', viewonly=True)
    grid_tile = relationship('GridTile', back_populates='mpointings', uselist=False)
    survey_tile = relationship('SurveyTile', back_populates='mpointings', uselist=False)
    event = relationship('Event', back_populates='mpointings', uselist=False)

    def __init__(self, object_name=None, ra=None, dec=None,
                 start_rank=None, min_alt=None, min_time=None,
                 max_moon=None, min_moonsep=None, max_sunalt=None, too=None, start_time=None,
                 stop_time=None, num_todo=None, valid_time=None, wait_time=None,
                 status='unscheduled', **kwargs):
        self.ra = ra
        self.dec = dec
        self.object_name = object_name
        self.start_rank = start_rank
        self.current_rank = start_rank
        self.initial_rank = start_rank
        self.max_moon = max_moon
        self.min_moonsep = min_moonsep
        self.min_alt = min_alt
        self.min_time = min_time
        self.max_sunalt = max_sunalt
        self.too = too
        self.start_time = start_time if start_time is not None else Time.now()
        self.stop_time = stop_time
        self.status = status
        self.infinite = False  # changed below
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
                self.time_blocks[0].current = True

        if 'user' in kwargs:
            self.user = kwargs['user']
        if 'user_id' in kwargs:
            self.user_id = kwargs['user_id']
        if 'grid_tile_id' in kwargs:
            self.grid_tile_id = kwargs['grid_tile_id']
        if 'grid_tile' in kwargs:
            self.grid_tile = kwargs['grid_tile']
        if 'survey_tile' in kwargs:
            self.survey_tile = kwargs['survey_tile']
        if 'survey_tile_id' in kwargs:
            self.survey_tile_id = kwargs['survey_tile_id']
        if 'event' in kwargs:
            self.event = kwargs['event']
        if 'event_id' in kwargs:
            self.event_id = kwargs['event_id']

    def __repr__(self):
        strings = ['db_id={}'.format(self.db_id),
                   'status={}'.format(self.status),
                   'num_todo={}'.format(self.num_todo),
                   'num_completed={}'.format(self.num_completed),
                   'num_remaining={}'.format(self.num_remaining),
                   'infinite={}'.format(self.infinite),
                   'object_name={}'.format(self.object_name),
                   'ra={}'.format(self.ra),
                   'dec={}'.format(self.dec),
                   'current_rank={}'.format(self.current_rank),
                   'initial_rank={}'.format(self.initial_rank),
                   'min_alt={}'.format(self.min_alt),
                   'max_sunalt={}'.format(self.max_sunalt),
                   'min_time={}'.format(self.min_time),
                   'max_moon={}'.format(self.max_moon),
                   'min_moonsep={}'.format(self.min_moonsep),
                   'too={}'.format(self.too),
                   'start_time={}'.format(self.start_time),
                   'stop_time={}'.format(self.stop_time),
                   'user_id={}'.format(self.user_id),
                   'grid_tile_id={}'.format(self.grid_tile_id),
                   'survey_tile_id={}'.format(self.survey_tile_id),
                   'event_id={}'.format(self.event_id),
                   ]
        return 'Mpointing({})'.format(', '.join(strings))

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
        current_block : `obsdb.TimeBlock`
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
        next_block : `obsdb.TimeBlock`
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
                     rank=self.current_rank,
                     min_alt=self.min_alt,
                     max_sunalt=self.max_sunalt,
                     min_time=self.min_time,
                     max_moon=self.max_moon,
                     min_moonsep=self.min_moonsep,
                     too=self.too,
                     start_time=start_time,
                     stop_time=stop_time,
                     status='pending',
                     mpointing=self,
                     user=self.user,
                     time_block=next_block,
                     grid_tile=self.grid_tile,
                     survey_tile=self.survey_tile,
                     event=self.event,
                     )
        # add the exposures
        p.exposure_sets = self.exposure_sets
        return p


class TimeBlock(Base):
    """A class to represent a block of time.

    Specifically, these represent a time period made up of a valid time and a wait time.
    The valid time is the time that the pointing will be valid for, and the wait time will be the
    time after the pointing is observed/invalid to wait until the next valid period.

    NB: you shouldn't ever have to manually create TimeBlocks.
    They're used internally by Mpointings to keep track of Pointings and all
    the ones an Mpointing requires are created when the Mpointing is.

    Like all SQLAlchemy model classes, this object links to the
    underlying database. You can create an instance, and set its attributes
    without a database session. Accessing some attributes may require
    an active database session, and some properties (like the db_id)
    will be None until the instance is added to the database.

    Parameters
    ----------
    block_num : int
        an integer indicating which block in a sequence this is
    valid_time : float
        amount of time a pointing in this block should stay valid in the queue, in minutes.
    wait_time : float
        time to wait after this block before allowing the next pointing, in minutes

    current : bool, optional
        True if this Time Block is the one that is currently linked to
        a Pointing in the queue

    Attributes
    ----------
    db_id : int
        primary database key
        only populated when the instance is added to the database

    When created the instance can be linked to the following other tables,
    otherwise they are populated when it is added to the database:

    Relationships
    -------------
    mpointing : `Mpointing`
        the Mpointing associated with this Time Block
        required before addition to the database
        can also be added with the mpointing_id parameter

    pointings : list of `Pointing`, optional
        the Pointings associated with this Time Block, if any

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
    block_num = Column(Integer)
    valid_time = Column(Integer)
    wait_time = Column(Integer)
    current = Column(Boolean, default=False)

    # Foreign keys
    mpointing_id = Column(Integer, ForeignKey('mpointings.id'), nullable=False)

    # Foreign relationships
    mpointing = relationship('Mpointing', back_populates='time_blocks', uselist=False)
    pointings = relationship('Pointing', back_populates='time_block')

    def __repr__(self):
        strings = ['db_id={}'.format(self.db_id),
                   'block_num={}'.format(self.block_num),
                   'valid_time={}'.format(self.valid_time),
                   'wait_time={}'.format(self.wait_time),
                   'current={}'.format(self.current),
                   'mpointing_id={}'.format(self.mpointing_id),
                   ]
        return 'TimeBlock({})'.format(', '.join(strings))


class Grid(Base):
    """A class to represent a Grid on the stellar sphere.

    A Grid is a way of generating Grid Tiles.
    See `gototile.grid.SkyGrid` for more details.

    Like all SQLAlchemy model classes, this object links to the
    underlying database. You can create an instance, and set its attributes
    without a database session. Accessing some attributes may require
    an active database session, and some properties (like the db_id)
    will be None until the instance is added to the database.

    Parameters
    ----------
    name : string
        a human-readable identifier for the grid
    ra_fov : float
        field of view in the RA direction (decimal degrees)
    dec_fov : float
        field of view in the Dec direction (decimal degrees)
    ra_overlap : float
        overlap in the RA direction (between 0 and 1)
    dec_overlap : float
        overlap in the Dec direction (between 0 and 1)
    algorithm : string
        the gridding algorithm used to generate the grid

    Attributes
    ----------
    db_id : int
        primary database key
        only populated when the instance is added to the database

    When created the instance can be linked to the following other tables,
    otherwise they are populated when it is added to the database:

    Relationships
    -------------
    grid_tiles : list of `GridTile`, optional
        the Grid Tiles associated with this Grid, if any
    surveys : list of `Survey`, optional
        the Surveys associated with this Grid, if any

    Examples
    --------
    TODO

    """

    # Set corresponding SQL table name
    __tablename__ = 'grids'

    # Primary key
    db_id = Column('id', Integer, primary_key=True)

    # Columns
    name = Column(String)
    ra_fov = Column(Float)
    dec_fov = Column(Float)
    ra_overlap = Column(Float)
    dec_overlap = Column(Float)
    algorithm = Column(String)

    # Foreign relationships
    grid_tiles = relationship('GridTile', back_populates='grid', uselist=True)
    surveys = relationship('Survey', back_populates='grid', uselist=True)

    def __repr__(self):
        strings = ['db_id={}'.format(self.db_id),
                   'name={}'.format(self.name),
                   'ra_fov={}'.format(self.ra_fov),
                   'dec_fov={}'.format(self.dec_fov),
                   'ra_overlap={}'.format(self.ra_overlap),
                   'dec_overlap={}'.format(self.dec_overlap),
                   'algorithm={}'.format(self.algorithm),
                   ]
        return 'Grid({})'.format(', '.join(strings))


class GridTile(Base):
    """A class to represent a Grid Tile.

    Grid Tiles are associated with Grids.
    See `gototile.grid.SkyGrid` for more details.

    Like all SQLAlchemy model classes, this object links to the
    underlying database. You can create an instance, and set its attributes
    without a database session. Accessing some attributes may require
    an active database session, and some properties (like the db_id)
    will be None until the instance is added to the database.

    Parameters
    ----------
    name : str
        a human-readable identifier for the tile
    ra : float
        J2000 right ascension in decimal degrees
    dec : float
        J2000 declination in decimal degrees

    Attributes
    ----------
    db_id : int
        primary database key
        only populated when the instance is added to the database

    When created the instance can be linked to the following other tables,
    otherwise they are populated when it is added to the database:

    Relationships
    -------------
    grid : `Grid`
        the Grid associated with this GridTile
        required before addition to the database
        can also be added with the grid_id parameter

    pointings : list of `Pointing`, optional
        the Pointings associated with this GridTile, if any
    mpointings : list of `Mpointing`, optional
        the Mpointing associated with this GridTile, if any
    survey_tiles : list of `SurveyTile`, optional
        the Survey Tiles associated with this GridTile, if any

    Examples
    --------
    TODO

    """

    # Set corresponding SQL table name
    __tablename__ = 'grid_tiles'

    # Primary key
    db_id = Column('id', Integer, primary_key=True)

    # Columns
    name = Column(String)
    ra = Column(Float)
    dec = Column('decl', Float)  # dec is reserved in SQL so can't be a column name

    # Foreign keys
    grid_id = Column(Integer, ForeignKey('grids.id'), nullable=False)

    # Foreign relationships
    grid = relationship('Grid', back_populates='grid_tiles', uselist=False)
    survey_tiles = relationship('SurveyTile', back_populates='grid_tile', uselist=True)
    mpointings = relationship('Mpointing', back_populates='grid_tile', uselist=True)
    pointings = relationship('Pointing', back_populates='grid_tile', uselist=True)

    def __repr__(self):
        strings = ['db_id={}'.format(self.db_id),
                   'name={}'.format(self.name),
                   'ra={}'.format(self.ra),
                   'dec={}'.format(self.dec),
                   'grid_id={}'.format(self.grid_id),
                   ]
        return 'GridTile({})'.format(', '.join(strings))


class Survey(Base):
    """A class to represent a Survey.

    A Survey is a grouping of tiles from a specific Grid.
    The purpose of Surveys is to add weighting to the base Grid Tiles.
    This is done using Survey Tiles.

    Like all SQLAlchemy model classes, this object links to the
    underlying database. You can create an instance, and set its attributes
    without a database session. Accessing some attributes may require
    an active database session, and some properties (like the db_id)
    will be None until the instance is added to the database.

    Parameters
    ----------
    name : string
        a human-readable identifier for the survey

    Attributes
    ----------
    db_id : int
        primary database key
        only populated when the instance is added to the database

    When created the instance can be linked to the following other tables,
    otherwise they are populated when it is added to the database:

    Relationships
    -------------
    grid : `Grid`
        the Grid associated with this Survey
        required before addition to the database
        can also be added with the grid_id parameter

    survey_tiles : list of `SurveyTile`, optional
        the Survey Tiles associated with this Survey, if any
    event : `Event`, optional
        the Event associated with this Survey, if any
        can also be added with the event_id parameter

    Examples
    --------
    TODO

    """

    # Set corresponding SQL table name
    __tablename__ = 'surveys'

    # Primary key
    db_id = Column('id', Integer, primary_key=True)

    # Columns
    name = Column(String)

    # Foreign keys
    grid_id = Column(Integer, ForeignKey('grids.id'), nullable=False)
    event_id = Column(Integer, ForeignKey('events.id'), nullable=True)

    # Foreign relationships
    survey_tiles = relationship('SurveyTile', back_populates='survey')
    grid = relationship('Grid', back_populates='surveys', uselist=False)
    event = relationship('Event', back_populates='surveys', uselist=False)

    def __repr__(self):
        strings = ['db_id={}'.format(self.db_id),
                   'name={}'.format(self.name),
                   'grid_id={}'.format(self.grid_id),
                   'event_id={}'.format(self.event_id),
                   ]
        return 'Survey({})'.format(', '.join(strings))


class SurveyTile(Base):
    """A class to represent a Survey Tile.

    Survey Tiles map onto Grid Tiles, but contain an additional weighting.

    Like all SQLAlchemy model classes, this object links to the
    underlying database. You can create an instance, and set its attributes
    without a database session. Accessing some attributes may require
    an active database session, and some properties (like the db_id)
    will be None until the instance is added to the database.

    Parameters
    ----------
    weight : float
        initial weighting for this tile (between 0 and 1)

    Attributes
    ----------
    db_id : int
        primary database key
        only populated when the instance is added to the database
    current_weight : float
        the current weighting of this tile
    initial_weight : float
        the initial weighting of this tile (it be modified as neighbours are observed)

    When created the instance can be linked to the following other tables,
    otherwise they are populated when it is added to the database:

    Relationships
    -------------
    survey : `Survey`
        the Survey associated with this SurveyTile
        required before addition to the database
        can also be added with the survey_id parameter
    grid_tile : `GridTile`
        the GridTile associated with this SurveyTile
        required before addition to the database
        can also be added with the grid_tile_id parameter

    pointings : list of `Pointing`, optional
        the Pointings associated with this SurveyTile, if any
    mpointings : list of `Mpointing`, optional
        the Mpointing associated with this SurveyTile, if any

    Examples
    --------
    TODO

    """

    # Set corresponding SQL table name
    __tablename__ = 'survey_tiles'

    # Primary key
    db_id = Column('id', Integer, primary_key=True)

    # Columns
    current_weight = Column(Float)
    initial_weight = Column(Float)

    # Foreign keys
    survey_id = Column(Integer, ForeignKey('surveys.id'), nullable=False)
    grid_tile_id = Column(Integer, ForeignKey('grid_tiles.id'), nullable=False)

    # Foreign relationships
    survey = relationship('Survey', back_populates='survey_tiles', uselist=False)
    grid_tile = relationship('GridTile', back_populates='survey_tiles', uselist=False)
    pointings = relationship('Pointing', back_populates='survey_tile', uselist=True)
    mpointings = relationship('Mpointing', back_populates='survey_tile', uselist=True)

    def __init__(self, weight=None, survey_id=None, grid_tile_id=None):
        self.current_weight = weight
        self.initial_weight = weight
        self.survey_id = survey_id
        self.grid_tile_id = grid_tile_id

    def __repr__(self):
        strings = ['db_id={}'.format(self.db_id),
                   'current_weight={}'.format(self.current_weight),
                   'initial_weight={}'.format(self.initial_weight),
                   'survey_id={}'.format(self.survey_id),
                   'grid_tile_id={}'.format(self.grid_tile_id),
                   ]
        return 'SurveyTile({})'.format(', '.join(strings))


class Event(Base):
    """A class to represent a transient Event.

    Events in the database are used to store infomation and link Surveys and Pointings.

    Like all SQLAlchemy model classes, this object links to the
    underlying database. You can create an instance, and set its attributes
    without a database session. Accessing some attributes may require
    an active database session, and some properties (like the db_id)
    will be None until the instance is added to the database.

    Parameters
    ----------
    name : string
        a human-readable identifier for the event
    ivorn : string
        unique IVORN (International Virtual Observatory Resource Name) for the event
    source : string
        the event's origin, e.g. LVC, Fermi, GAIA
    event_type : string
        the type of event, e.g. GW, GRB

    time : string, `astropy.time.Time` or datetime.datetime, optional
        time the event occured
    skymap : string, optional
        the location of the source skymap file

    Attributes
    ----------
    db_id : int
        primary database key
        only populated when the instance is added to the database

    When created the instance can be linked to the following other tables,
    otherwise they are populated when it is added to the database:

    Relationships
    -------------
    surveys : list of `Survey`, optional
        the Surveys associated with this Event, if any
    pointings : list of `Pointing`, optional
        the Pointings associated with this Event, if any
    mpointings : list of `Mpointing`, optional
        the Mpointings associated with this Event, if any

    Examples
    --------
    TODO

    """

    # Set corresponding SQL table name
    __tablename__ = 'events'

    # Primary key
    db_id = Column('id', Integer, primary_key=True)

    # Columns
    name = Column(String)
    ivorn = Column(String, unique=True)
    source = Column(String)
    event_type = Column('type', String)  # type is a built in function in Python
    time = Column(DateTime)
    skymap = Column(String)

    # Foreign relationships
    surveys = relationship('Survey', back_populates='event', uselist=True)
    pointings = relationship('Pointing', back_populates='event', uselist=True)
    mpointings = relationship('Mpointing', back_populates='event', uselist=True)

    def __repr__(self):
        strings = ['db_id={}'.format(self.db_id),
                   'name={}'.format(self.name),
                   'ivorn={}'.format(self.ivorn),
                   'source={}'.format(self.source),
                   'event_type={}'.format(self.event_type),
                   'time={}'.format(self.time),
                   'skymap={}'.format(self.skymap),
                   ]
        return 'Event({})'.format(', '.join(strings))

    @validates('time')
    def munge_times(self, key, field):
        """Use validators to allow various types of input for UTC."""
        if isinstance(field, datetime.datetime):
            value = field.strftime("%Y-%m-%d %H:%M:%S")
        elif isinstance(field, Time):
            field.precision = 0  # no D.P on seconds
            value = field.iso
        else:
            # just hope the string works!
            value = str(field)
        return value


class ImageLog(Base):
    """A class to store a record of an Image.

    The ImageLog is a simple way to link Pointings in the database to physical
    FITS files. An ImageLog should be created each time an exposure is taken,
    even if taken manually rather than originating from the database.

    Like all SQLAlchemy model classes, this object links to the
    underlying database. You can create an instance, and set its attributes
    without a database session. Accessing some attributes may require
    an active database session, and some properties (like the db_id)
    will be None until the instance is added to the database.

    Parameters
    ----------
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

    Attributes
    ----------
    db_id : int
        primary database key
        only populated when the instance is added to the database

    When created the instance can be linked to the following other tables,
    otherwise they are populated when it is added to the database:

    Relationships
    -------------
    exposure_set : `ExposureSet`, optional
        the Exposure Set associated with this ImageLog, if any
        can also be added with the exposure_set_id parameter
    pointing : `Pointing`, optional
        the Pointing associated with this ImageLog, if any
        can also be added with the pointing_id parameter
    mpointing : `Mpointing`, optional
        the Mpointing associated with this ImageLog, if any
        can also be added with the mpointing_id parameter

    """

    # Set corresponding SQL table name
    __tablename__ = 'image_logs'

    # Primary key
    db_id = Column('id', Integer, primary_key=True)

    # Columns
    filename = Column(String)
    run_number = Column(Integer)
    ut = Column(Integer)
    ut_mask = Column(Integer)
    start_time = Column(DateTime)
    write_time = Column(DateTime)
    set_position = Column(Integer, default=1)
    set_total = Column(Integer, default=1)

    # Foreign keys
    exposure_set_id = Column(Integer, ForeignKey('exposure_sets.id'), nullable=True)
    pointing_id = Column(Integer, ForeignKey('pointings.id'), nullable=True)
    mpointing_id = Column(Integer, ForeignKey('mpointings.id'), nullable=True)

    # Foreign relationships
    exposure_set = relationship('ExposureSet', backref='image_logs', uselist=False)
    pointing = relationship('Pointing', backref='image_logs', uselist=False)
    mpointing = relationship('Mpointing', backref='image_logs', uselist=False)

    def __repr__(self):
        strings = ['db_id={}'.format(self.db_id),
                   'filename={}'.format(self.filename),
                   'run_number={}'.format(self.run_number),
                   'ut={}'.format(self.ut),
                   'ut_mask={}'.format(self.ut_mask),
                   'start_time={}'.format(self.start_time),
                   'write_time={}'.format(self.write_time),
                   'exposure_set_id={}'.format(self.exposure_set_id),
                   'pointing_id={}'.format(self.pointing_id),
                   'mpointing_id={}'.format(self.mpointing_id),
                   ]
        return 'ImageLog({})'.format(', '.join(strings))

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
