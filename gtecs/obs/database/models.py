"""Python classes mapping on to database tables."""

import datetime
import hashlib

from astropy import units as u
from astropy.coordinates import EarthLocation
from astropy.time import Time

from gototile.grid import SkyGrid

from sqlalchemy import Boolean, Column, DateTime, Enum, Float, ForeignKey, Integer, String, Text
from sqlalchemy import exists, func, select, text
from sqlalchemy.dialects.mysql import TIMESTAMP
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.ext.hybrid import hybrid_method, hybrid_property
from sqlalchemy.orm import column_property, relationship, validates
from sqlalchemy.sql import and_, case, or_


__all__ = ['User', 'Pointing', 'ExposureSet', 'Mpointing', 'TimeBlock',
           'Site', 'Telescope',
           'Grid', 'GridTile', 'Survey', 'SurveyTile', 'Event',
           'ImageLog',
           'TRIGGERS']


class Utf8Base:
    """A simple base to force the table to use UTF8 encoding.

    See https://stackoverflow.com/questions/63488890/
    """

    __table_args__ = {'mysql_default_charset': 'utf8'}


Base = declarative_base(cls=Utf8Base)


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
        the password is stored as a hash in the database, this unhashed string isn't stored
    full_name : string
        the user's full name.

    When created the instance can be linked to the following other tables as parameters,
    otherwise they are populated when it is added to the database:

    Primary relationships
    -------------
    mpointings : list of `Mpointing`, optional
        the Mpointings created by this User, if any

    Attributes
    ----------
    db_id : int
        primary database key
        only populated when the instance is added to the database
    password_hash : str
        the hashed password as stored in the database

    """

    # Set corresponding SQL table name
    __tablename__ = 'users'

    # Primary key
    db_id = Column('id', Integer, primary_key=True)

    # Columns
    username = Column(String(255), nullable=False, unique=True)
    password_hash = Column('password', String(255), nullable=False)
    full_name = Column(Text, nullable=False)

    # Foreign relationships
    mpointings = relationship('Mpointing', back_populates='user')

    def __init__(self, **kwargs):
        # Use hashed password instead of plain
        kwargs['password_hash'] = hashlib.sha512(kwargs.pop('password').encode()).hexdigest()

        # Init base class
        super().__init__(**kwargs)

    def __repr__(self):
        strings = ['db_id={}'.format(self.db_id),
                   'username={}'.format(self.username),
                   'full_name={}'.format(self.full_name),
                   ]
        return 'User({})'.format(', '.join(strings))


class ExposureSet(Base):
    """A class to represent an Exposure Set.

    An Exposure Set is a set of repeated identical exposures, with the same
    exposure time, filter, binning factor etc. Each Mpointing should have one or more
    Exposure Sets defined, to be observed when its Pointings are selected by the scheduler.

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

    binning : int, optional
        binning factor to apply
        default = 1 (no binning)
    imgtype : string, optional
        indicates the type of exposure set.
        one of SCIENCE, FOCUS, STD, FLAT, BIAS, DARK
        default = 'SCIENCE'
    ut_mask : int, optional
        if set, this is a binary mask which will determine which unit
        telescopes carry out the exposure. A value of 5 (binary 0101) will
        be exposed by cameras 1 and 3.
        default = None
    ra_offset : float, optional
        the size of the random offset to apply between each exposure
        if not set, no offset will be made
        default = 0
    dec_offset : float, optional
        the size of the random offset to apply between each exposure
        if not set, no offset will be made
        default = 0

    When created the instance can be linked to the following other tables as parameters,
    otherwise they are populated when it is added to the database:

    Primary relationships
    ---------------------
    mpointing : `Mpointing`, optional
        the Mpointing this Exposure Set is linked to, if any
        can also be added with the mpointing_id parameter

    Attributes
    ----------
    db_id : int
        primary database key
        only populated when the instance is added to the database

    The following secondary relationships are not settable directly,
    but are populated through the primary relationships and are available as attributes:

    Secondary relationships
    -------------
    pointings : list of `Pointing`
        the Pointings associated with the Mpointing this Exposure Set is linked to, if any

    """

    # Set corresponding SQL table name
    __tablename__ = 'exposure_sets'

    # Primary key
    db_id = Column('id', Integer, primary_key=True)

    # Columns
    num_exp = Column(Integer, nullable=False)
    exptime = Column(Float, nullable=False)
    filt = Column('filter',   # filter is a built in function in Python
                  String(1), nullable=False)
    binning = Column(Integer, nullable=False, default=1)
    imgtype = Column(Enum('SCIENCE', 'FOCUS', 'DARK', 'BIAS', 'FLAT', 'STD'), nullable=False,
                     default='SCIENCE')
    ut_mask = Column(Integer, default=None)
    ra_offset = Column(Float, nullable=False, default=0)
    dec_offset = Column(Float, nullable=False, default=0)

    # Foreign keys
    mpointing_id = Column(Integer, ForeignKey('mpointings.id'), nullable=True)

    # Foreign relationships
    mpointing = relationship('Mpointing', back_populates='exposure_sets')
    image_logs = relationship('ImageLog', lazy='joined', back_populates='exposure_set')

    # Secondary relationships
    pointings = relationship('Pointing',
                             secondary='mpointings',
                             primaryjoin='ExposureSet.mpointing_id == Mpointing.db_id',
                             secondaryjoin='Pointing.mpointing_id == Mpointing.db_id',
                             back_populates='exposure_sets',
                             viewonly=True,
                             uselist=True,
                             )

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
                   'mpointing_id={}'.format(self.mpointing_id),
                   ]
        return 'ExposureSet({})'.format(', '.join(strings))


class Pointing(Base):
    """A class to represent an Pointing.

    Pointings are the primary table of the Observation Database, and are usually generated
    automatically by Mpointings. When deciding what to observe the scheduler processes a queue
    made of Pointings based on their status and constraints.

    NB: you shouldn't ever have to manually create Pointings, instead you should generate them
    from an Mpointing using its get_next_pointing() method. The first Pointing is also created when
    the Mpointing is.

    Like all SQLAlchemy model classes, this object links to the
    underlying database. You can create an instance, and set its attributes
    without a database session. Accessing some attributes may require
    an active database session, and some properties (like the db_id)
    will be None until the instance is added to the database.

    Parameters
    ----------
    rank : Integer
        rank to use when scheduling (should be the current_rank of the linked mpointing)

    start_time : string, `astropy.time.Time` or datetime.datetime, optional
        UTC time from which pointing is considered valid and can be started
        default = Time.now()
    stop_time : string, `astropy.time.Time` or datetime.datetime, or None, optional
        the latest UTC time at which pointing may be started
        can be None, if so the pointing will stay in the queue indefinitely
        (it can't be marked as expired) and will only leave when observed
        default = None
    creation_time : string, `astropy.time.Time` or datetime.datetime, optional
        the time the pointing was created
        this should usually be set automatically, unless doing simulations
        default = Time.now()
    running_time : string, `astropy.time.Time` or datetime.datetime, or None, optional
        the time the pointing was marked as running (i.e. started being observed)
        this should usually be set automatically when the pointing status is set
        default = None
    finished_time : string, `astropy.time.Time` or datetime.datetime, or None, optional
        the time the pointing was marked as finished (i.e. completed or interupted)
        this should usually be set automatically when the pointing status is set
        default = None
    completed : bool or None, optional
        if the pointing was successfully completed (True) or interupted (False), or just
        hasn't started yet (None)
        default = None

    When created the instance can be linked to the following other tables as parameters,
    otherwise they are populated when it is added to the database:

    Primary relationships
    ---------------------
    mpointing : `Mpointing`
        the Mpointing associated with this Pointing
        can also be added with the mpointing_id parameter

    time_block : `TimeBlock`, optional
        the TimeBlock this Pointing as generated from, if any
        can also be added with the time_block_id parameter
    telescope : `Telescope`, optional
        the Telescope this Pointing was observed by, if any
        can also be added with the telescope_id parameter

    Attributes
    ----------
    db_id : int
        primary database key
        only populated when the instance is added to the database
    status : string
        status of the pointing

    The following secondary relationships are not settable directly,
    but are populated through the primary relationships and are available as attributes:

    Secondary relationships
    -----------------------
    exposure_sets : list of `ExposureSet`
        the ExposureSets linked to the Mpointing this Pointing was generated from
    grid_tile : `GridTile`
        the GridTile the Mpointing this Pointing was generated from covers, if any
    survey_tile : `SurveyTile`
        the SurveyTile the Mpointing this Pointing was generated from covers, if any
    event : `Event`
        the Event the Mpointing this Pointing was generated from was part of, if any

    Methods
    -------
    status_at_time(time) : str
        status of the pointing at the given time
    mark_deleted(time=None)
        mark the pointing as deleted
        if no time is given, default to Time.now()
    mark_running(self, telescope=None, time=None)
        mark the pointing as running on the given telescope
        if no time is given, default to Time.now()
    mark_finished(self, completed=True, time=None)
        mark the pointing as finished, either completed (True) or interupted (False)
        if no time is given, default to Time.now()

    """

    # Set corresponding SQL table name
    __tablename__ = 'pointings'

    # Primary key
    db_id = Column('id', Integer, primary_key=True)

    # Columns
    # # Basic properties
    rank = Column(Integer, nullable=False)
    # # Constraints
    start_time = Column(DateTime, nullable=False, index=True, server_default=func.now())
    stop_time = Column(DateTime, nullable=True, index=True, default=None)
    # # Status
    creation_time = Column(DateTime, nullable=False, server_default=func.now())
    running_time = Column(DateTime, nullable=True, default=None)
    finished_time = Column(DateTime, nullable=True, default=None)
    completed = Column(Boolean, nullable=False, default=False)
    # # Update timestamp (TODO: Do we still need this? Or should more tables have it?)
    ts = Column(TIMESTAMP(fsp=3), nullable=False,
                server_default=text('CURRENT_TIMESTAMP(3) ON UPDATE CURRENT_TIMESTAMP(3)'))

    # Foreign keys
    mpointing_id = Column(Integer, ForeignKey('mpointings.id'), nullable=False)
    time_block_id = Column(Integer, ForeignKey('time_blocks.id'), nullable=False)
    telescope_id = Column(Integer, ForeignKey('telescopes.id'), nullable=True)

    # Foreign relationships
    mpointing = relationship('Mpointing', lazy='joined', back_populates='pointings')
    time_block = relationship('TimeBlock', lazy='joined', back_populates='pointings')
    telescope = relationship('Telescope', lazy='joined', back_populates='pointings')
    image_logs = relationship('ImageLog', lazy='joined', back_populates='pointing')

    # Secondary relationships
    exposure_sets = relationship('ExposureSet',
                                 secondary='mpointings',
                                 primaryjoin='Pointing.mpointing_id == Mpointing.db_id',
                                 secondaryjoin='ExposureSet.mpointing_id == Mpointing.db_id',
                                 back_populates='pointings',
                                 viewonly=True,
                                 uselist=True,
                                 )
    grid_tile = relationship('GridTile',
                             lazy='joined',
                             secondary='mpointings',
                             primaryjoin='Pointing.mpointing_id == Mpointing.db_id',
                             secondaryjoin='GridTile.db_id == Mpointing.grid_tile_id',
                             back_populates='pointings',
                             viewonly=True,
                             uselist=False,
                             )
    survey_tile = relationship('SurveyTile',
                               lazy='joined',
                               secondary='mpointings',
                               primaryjoin='Pointing.mpointing_id == Mpointing.db_id',
                               secondaryjoin='SurveyTile.db_id == Mpointing.survey_tile_id',
                               back_populates='pointings',
                               viewonly=True,
                               uselist=False,
                               )
    event = relationship('Event',
                         lazy='joined',
                         secondary='mpointings',
                         primaryjoin='Pointing.mpointing_id == Mpointing.db_id',
                         secondaryjoin='Event.db_id == Mpointing.event_id',
                         back_populates='pointings',
                         viewonly=True,
                         uselist=False,
                         )

    def __init__(self, **kwargs):
        # Set default times
        if 'creation_time' not in kwargs or kwargs['creation_time'] is None:
            kwargs['creation_time'] = Time.now()
        if 'start_time' not in kwargs or kwargs['start_time'] is None:
            kwargs['start_time'] = kwargs['creation_time']

        # Init base class
        super().__init__(**kwargs)

    def __repr__(self):
        strings = ['db_id={}'.format(self.db_id),
                   'status={}'.format(self.status),
                   'rank={}'.format(self.rank),
                   'start_time={}'.format(self.start_time),
                   'stop_time={}'.format(self.stop_time),
                   'creation_time={}'.format(self.creation_time),
                   'running_time={}'.format(self.running_time),
                   'finished_time={}'.format(self.finished_time),
                   'completed={}'.format(self.completed),
                   'mpointing_id={}'.format(self.mpointing_id),
                   'time_block_id={}'.format(self.time_block_id),
                   'telescope_id={}'.format(self.telescope_id),
                   ]
        return 'Pointing({})'.format(', '.join(strings))

    @validates('start_time', 'stop_time', 'creation_time', 'running_time', 'finished_time')
    def validate_times(self, key, field):
        """Use validators to allow various types of input for times.

        Also enforce stop_time > start_time and finished_time > running_time.
        """
        if key not in ['start_time', 'creation_time'] and field is None:
            # start_time and creation_time are not nullable, the others are
            return None

        if isinstance(field, datetime.datetime):
            value = field.strftime('%Y-%m-%d %H:%M:%S')
        elif isinstance(field, Time):
            field.precision = 0  # no D.P on seconds
            value = field.iso
        else:
            # just hope the string works!
            value = str(field)

        # force start_time > stop_time
        if (key == 'start_time' and self.stop_time is not None and
                Time(value) >= Time(self.stop_time)):
            raise ValueError(f'start_time must be before stop_time ({self.stop_time})')
        if (key == 'stop_time' and self.start_time is not None and
                Time(self.start_time) >= Time(value)):
            raise ValueError(f'stop_time must be after start_time ({self.start_time})')

        # force running_time > finished_time
        if (key == 'running_time' and self.finished_time is not None and
                Time(value) >= Time(self.finished_time)):
            raise ValueError(f'running_time must be before finished_time ({self.finished_time})')
        if (key == 'finished_time' and self.running_time is not None and
                Time(self.running_time) >= Time(value)):
            raise ValueError(f'finished_time must be after running_time ({self.running_time})')

        return value

    @hybrid_property
    def status(self):
        """Return a string giving the current status."""
        time = Time.now()

        if self.running_time is None and self.finished_time is not None:
            # The Pointing was marked as finished before it started: deleted
            return 'deleted'
        elif self.running_time is not None and self.finished_time is None:
            # The Pointing has started but hasn't finished yet: running
            return 'running'
        elif self.finished_time is not None and self.completed:
            # The Pointing has finished and is flagged as completed: completed
            return 'completed'
        elif self.finished_time is not None and not self.completed:
            # The Pointing has finished and is not flagged as completed: interrupted
            return 'interrupted'
        elif time < Time(self.start_time):  # start_time can't be None
            # The Pointing hasn't yet reached its start time: upcoming
            return 'upcoming'
        elif self.stop_time is not None and time > Time(self.stop_time):
            # The Pointing has passed its start time: expired
            return 'expired'
        else:
            # The Pointing isn't running or finished: pending
            return 'pending'

    @status.setter
    def status(self, value):
        if value == 'deleted':
            self.mark_deleted()
        elif value == 'running':
            self.mark_running()  # Will raise an error since we don't give a Telescope
        elif value == 'completed':
            self.mark_finished(completed=True)
        elif value == 'interrupted':
            self.mark_finished(completed=False)
        elif value in ['upcoming', 'expired', 'pending']:
            raise ValueError(f'Can not set status to "{value}", change start/stop times instead')
        else:
            raise ValueError(f'Invalid status: "{value}"')

    @status.expression
    def status(self):
        time = datetime.datetime.utcnow()
        return case([(and_(self.running_time.is_(None), self.finished_time.isnot(None)),
                      'deleted'),
                     (and_(self.running_time.isnot(None), self.finished_time.is_(None)),
                      'running'),
                     (and_(self.finished_time.isnot(None), self.completed.is_(True)),
                      'completed'),
                     (and_(self.finished_time.isnot(None), self.completed.is_(False)),
                      'interrupted'),
                     (time < self.start_time,
                      'upcoming'),
                     (and_(self.stop_time.isnot(None), time > self.stop_time),
                      'expired'),
                     ],
                    else_='pending')

    @hybrid_method
    def status_at_time(self, time):
        """Return a string giving the status at the given time (for testing and simulations)."""
        if isinstance(time, (str, datetime.datetime)):
            time = Time(time)

        if time < Time(self.creation_time):
            # This shouldn't usually happen, unless we are doing simulations and set it manually,
            # then it's worth knowing if this entry wouldn't have been created yet..
            return None

        if (self.running_time is None and self.finished_time is not None and
                time >= Time(self.finished_time)):
            return 'deleted'
        elif (self.running_time is not None and time >= Time(self.running_time) and
              (self.finished_time is None or time < Time(self.finished_time))):
            return 'running'
        elif (self.finished_time is not None and time >= Time(self.finished_time) and
              self.completed):
            return 'completed'
        elif (self.finished_time is not None and time >= Time(self.finished_time) and
              not self.completed):
            return 'interrupted'
        elif time < Time(self.start_time):
            return 'upcoming'
        elif self.stop_time is not None and time >= Time(self.stop_time):
            return 'expired'
        else:
            return 'pending'

    @status_at_time.expression
    def status_at_time(self, time):
        if isinstance(time, str):
            time = Time(time)
        if isinstance(time, Time):
            time = time.datetime

        c = case([(time < self.creation_time,
                   None),
                  (and_(self.running_time.is_(None), self.finished_time.isnot(None),
                        time >= self.finished_time),
                   'deleted'),
                  (and_(self.running_time.isnot(None), time >= self.running_time,
                        or_(self.finished_time.is_(None), time < self.finished_time)),
                   'running'),
                  (and_(self.finished_time.isnot(None), time >= self.finished_time,
                        self.completed.is_(True)),
                   'completed'),
                  (and_(self.finished_time.isnot(None), time >= self.finished_time,
                        self.completed.is_(False)),
                   'interrupted'),
                  (time < self.start_time,
                   'upcoming'),
                  (and_(self.stop_time.isnot(None), time >= self.stop_time),
                   'expired'),
                  ],
                 else_='pending')
        return c

    def mark_deleted(self, time=None):
        """Mark this Pointing as deleted."""
        if time is None:
            time = Time.now()
        if isinstance(time, (str, datetime.datetime)):
            time = Time(time)

        if self.status_at_time(time) == 'deleted':
            raise ValueError(f'Pointing is already deleted (at {self.finished_time})')
        if self.status_at_time(time) == 'expired':
            raise ValueError(f'Pointing is already expired (at {self.stop_time})')
        if self.status_at_time(time) == 'running':
            raise ValueError(f'Pointing is already running (at {self.running_time})')
        if self.status_at_time(time) in ['completed', 'interrupted']:
            raise ValueError(f'Pointing is already finished (at {self.finished_time})')

        self.finished_time = time

    def mark_running(self, telescope=None, time=None):
        """Mark this Pointing as running on the given telescope."""
        if time is None:
            time = Time.now()
        if isinstance(time, (str, datetime.datetime)):
            time = Time(time)
        if telescope is None:
            raise ValueError('Pointings must be linked to a Telescope when marked running')

        if self.status_at_time(time) == 'deleted':
            raise ValueError(f'Pointing is already deleted (at {self.finished_time})')
        if self.status_at_time(time) == 'upcoming':
            raise ValueError(f'Pointing has not yet reached its start time (at {self.start_time})')
        if self.status_at_time(time) == 'expired':
            raise ValueError(f'Pointing is already expired (at {self.stop_time})')
        if self.status_at_time(time) == 'running':
            raise ValueError(f'Pointing is already running (at {self.running_time})')
        if self.status_at_time(time) in ['completed', 'interrupted']:
            raise ValueError(f'Pointing is already finished (at {self.finished_time})')

        if self.telescope is not None or self.telescope_id is not None:
            raise ValueError(f'Pointing is already linked to a Telescope: {self.telescope}')

        if (self.mpointing.valid_telescopes is not None and
                telescope.db_id not in self.mpointing.valid_telescopes):
            raise ValueError(f'Telescope ID ({telescope.db_id}) not in '
                             f'list of valid telescopes ({self.mpointing.valid_telescopes})')
        if telescope.current_pointing is not None:
            pointing_id = telescope.current_pointing.db_id
            raise ValueError(f'Telescope is already observing Pointing ID {pointing_id}')

        self.running_time = time
        self.telescope = telescope

    def mark_finished(self, completed=True, time=None):
        """Mark this Pointing as stopped, either completed (True) or interrupted (False)."""
        if time is None:
            time = Time.now()
        if isinstance(time, (str, datetime.datetime)):
            time = Time(time)

        if self.status_at_time(time) == 'deleted':
            raise ValueError(f'Pointing is already deleted (at {self.finished_time})')
        if self.status_at_time(time) == 'expired':
            raise ValueError(f'Pointing is already expired (at {self.stop_time})')
        if self.status_at_time(time) in ['completed', 'interrupted']:
            raise ValueError(f'Pointing is already finished (at {self.finished_time})')
        if self.status_at_time(time) != 'running':
            raise ValueError('Pointing is not running (use mark_deleted to stop before running)')

        self.finished_time = time

        self.completed = completed


class Mpointing(Base):
    """A class to represent an Mpointing.

    Mpointings are generation engines for Pointings, and allow the
    database to regenerate Pointings at given time intervals.

    Like all SQLAlchemy model classes, this object links to the
    underlying database. You can create an instance, and set its attributes
    without a database session. Accessing some attributes may require
    an active database session, and some properties (like the db_id)
    will be None until the instance is added to the database.

    Parameters
    ----------
    object_name : str
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
    num_todo : int
        number of (successful) observations required.
        less than zero means repeat infinitely.

    wait_time : float or list of float, optional
        time to wait between pointings in minutes.
        if num_todo is greater than times given the list will be looped.
        default = 0 (no delay)
    valid_time : float or list of float, optional
        the amount of time the pointing(s) should be valid in the queue.
        if num_todo is greater than times given the list will be looped.
        less than zero means valid indefinitely
        default = -1 (indefinitely valid)

    min_alt : float, optional
        minimum altitude to observe at, degrees
        default = 30
    max_sunalt : float, optional
        altitude constraint on Sun, degrees
        default = -15
    min_time : float, optional
        minimum time needed to schedule pointing
        default = None
    max_moon : string, optional
        Moon constraint
        one of 'D', 'G', 'B'
        default = 'B'
    min_moonsep : float, optional
        distance constraint from the Moon, degrees
        default = 30
    too : bool, optional
        indicates if this is a Target of Opportunity (ToO)
        default = False
    tel_mask : int, optional
        if set, this is a binary mask which will determine which telescopes
        can carry out the observation.
        A value of 4 (binary (0100) means only observable by Telescope 3, a value of
        5 (binary 0101) will allow either Telescope 1 or 3 to observe the Pointing.
        default = None (no restriction)
    start_time : string, `astropy.time.Time` or datetime.datetime, optional
        UTC time from which Mpointing is considered valid and can be started
        if not given then set to now, so the Mpointing will start immediately
        default = Time.now()
    stop_time : string, `astropy.time.Time` or datetime.datetime, optional
        the latest UTC time after which pointings must stop
        if not given the Mpointing will continue creating pointings until
        it is completed
        default = None
    creation_time : string, `astropy.time.Time` or datetime.datetime, optional
        the time the Mpointing was created
        this should usually be set automatically, unless doing simulations
        default = Time.now()
    deleted_time : string, `astropy.time.Time` or datetime.datetime, or None, optional
        the time the pointing was marked as deleted
        this should usually be set automatically when the pointing status is set
        default = None

    When created the instance can be linked to the following other tables as parameters,
    otherwise they are populated when it is added to the database:

    Primary relationships
    ---------------------
    user : `User`
        the User associated with this Mpointing
        required before addition to the database
        can also be added with the user_id parameter

    exposure_sets : list of `ExposureSet`, optional
        the Exposure Sets associated with this Mpointing, if any
    pointings : list of `Pointing`, optional
        the Pointings generated by this Mpointing, if any
    time_blocks : list of `TimeBlock`, optional
        the TimeBlocks generated by this Mpointing, if any
    grid_tile : `GridTile`, optional
        the GridTile this Mpointing covers, if any
        can also be added with the grid_tile_id parameter
    survey_tile : `SurveyTile`, optional
        the SurveyTile this Mpointing covers, if any
        can also be added with the survey_tile_id parameter
    event : `Event`, optional
        the Event this Mpointing is part of, if any
        can also be added with the event_id parameter

    Attributes
    ----------
    db_id : int
        primary database key
        only populated when the instance is added to the database
    status : string
        status of the mpointing
    scheduled : bool
        returns True if the Mpointing is scheduled (has a pending Pointing)
    current_rank : int
        rank for next pointing to be scheduled (it will increase as pointings are observed)
    valid_telescopes : list of int
        a list of valid Telescope IDs, if any, based on the tel_mask

    num_completed : int
        number of successfully completed pointings
    num_remaining : int
        number of pointings still to do (same as num_todo - num_completed)
    infinite : bool
        if the Mpointing will continue infnitely (set if num_todo is < 0)
    last_observed : `astropy.time.Time`
        the most recent time of completion of any Pointing linked to this Mpointing, if any

    The following secondary relationships are not settable directly,
    but are populated through the primary relationships and are available as attributes:

    Secondary relationships
    -----------------------
    grid : `Grid`
        the Grid that the GridTile this Mpointing covers, if any, is part of
    survey : `Survey`
        the Survey that the SurveyTile this Mpointing covers, if any, is part of

    Methods
    -------
    status_at_time(time) : str
        status of the pointing at the given time
    scheduled_at_time(time) : bool
        returns True if the Mpointing was scheduled (has a pending Pointing) at the given time
    num_completed_at_time(time) : int
        number of successfully completed pointings at the given time
    num_remaining_at_time(time) : int
        number of pointings still to do at the given time
    current_rank_at_time(time) : int
        rank for next pointing to be scheduled at the given time
    mark_deleted(time=None)
        mark the pointing as deleted
        if no time is given, default to Time.now()
    get_current_block() : TimeBlock, or None
        get the current TimeBlock
    get_last_block() : TimeBlock, or None
        get the previous TimeBlock
    get_next_block() : TimeBlock, or None
        get the upcoming TimeBlock
    get_next_pointing(time=None) : Pointing, or None
        generate the next Pointing for this Mpointing
        returns None if a Pointing is already scheduled, or if the Mpointing is completed
        the time is set as the Pointing creation time, and should usually be set automatically
        unless doing simulations

    """

    # Set corresponding SQL table name
    __tablename__ = 'mpointings'

    # Primary key
    db_id = Column('id', Integer, primary_key=True)

    # Columns
    # # Basic properties
    object_name = Column('object', Text, nullable=False)  # object is a built in class in Python
    ra = Column(Float, nullable=False)
    dec = Column(Float, nullable=False)
    start_rank = Column(Integer, nullable=False)
    num_todo = Column(Integer, nullable=False)
    # # Constraints
    min_alt = Column(Float, nullable=False, default=30)
    max_sunalt = Column(Float, nullable=False, default=-15)
    min_time = Column(Float, nullable=False)
    max_moon = Column(String(1), nullable=False, default='B')
    min_moonsep = Column(Float, nullable=False, default=30)
    too = Column(Boolean, nullable=False, default=False)
    tel_mask = Column(Integer, default=None)
    start_time = Column(DateTime, nullable=False, index=True, server_default=func.now())
    stop_time = Column(DateTime, nullable=True, index=True, default=None)
    # # Status
    creation_time = Column(DateTime, nullable=False, server_default=func.now())
    deleted_time = Column(DateTime, nullable=True, default=None)

    # Foreign keys
    user_id = Column(Integer, ForeignKey('users.id'), nullable=False)
    grid_tile_id = Column(Integer, ForeignKey('grid_tiles.id'), nullable=True)
    survey_tile_id = Column(Integer, ForeignKey('survey_tiles.id'), nullable=True)
    event_id = Column(Integer, ForeignKey('events.id'), nullable=True)

    # Foreign relationships
    # (remember to add to __init__)
    user = relationship('User', lazy='joined', back_populates='mpointings')
    pointings = relationship('Pointing', lazy='joined', back_populates='mpointing')
    exposure_sets = relationship('ExposureSet', lazy='joined', back_populates='mpointing')
    time_blocks = relationship('TimeBlock', lazy='joined', back_populates='mpointing')
    grid_tile = relationship('GridTile', lazy='joined', back_populates='mpointings')
    survey_tile = relationship('SurveyTile', lazy='joined', back_populates='mpointings')
    event = relationship('Event', lazy='joined', back_populates='mpointings')

    # Secondary relationships
    grid = relationship('Grid',
                        lazy='joined',
                        secondary='grid_tiles',
                        primaryjoin='Mpointing.grid_tile_id == GridTile.db_id',
                        secondaryjoin='Grid.db_id == GridTile.grid_id',
                        back_populates='mpointings',
                        viewonly=True,
                        uselist=False,
                        )
    survey = relationship('Survey',
                          lazy='joined',
                          secondary='survey_tiles',
                          primaryjoin='Mpointing.survey_tile_id == SurveyTile.db_id',
                          secondaryjoin='Survey.db_id == SurveyTile.survey_id',
                          back_populates='mpointings',
                          viewonly=True,
                          uselist=False,
                          )

    # Column properties
    last_observed = column_property(select([Pointing.finished_time])
                                    .where(and_(Pointing.mpointing_id == db_id,
                                                Pointing.status == 'completed',
                                                )
                                           )
                                    .order_by(Pointing.finished_time.desc())
                                    .limit(1)
                                    .correlate_except(Pointing)
                                    .scalar_subquery()
                                    )

    def __init__(self, **kwargs):
        # Get extra arguments (need to remove from kwargs before super init)
        valid_times = kwargs.pop('valid_time') if 'valid_time' in kwargs else -1
        if not isinstance(valid_times, list):
            valid_times = [valid_times]
        wait_times = kwargs.pop('wait_time') if 'wait_time' in kwargs else 0
        if not isinstance(wait_times, list):
            wait_times = [wait_times]

        # Set default times
        if 'creation_time' not in kwargs or kwargs['creation_time'] is None:
            kwargs['creation_time'] = Time.now()
        if 'start_time' not in kwargs or kwargs['start_time'] is None:
            kwargs['start_time'] = kwargs['creation_time']

        # Init base class
        super().__init__(**kwargs)

        # Create TimeBlocks
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

        # Now create the first Pointing
        # It's unlikely that it returns None, but in theory you could create an Mpointing with
        # num_todo=0. Also note get_next_pointing() allows creating Pointings if the Mpointing
        # status is 'upcoming', so we can create it now even if it's not valid yet.
        pointing = self.get_next_pointing(time=self.creation_time)
        if pointing is not None:
            self.pointings.append(pointing)

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
                   'start_rank={}'.format(self.start_rank),
                   'min_alt={}'.format(self.min_alt),
                   'max_sunalt={}'.format(self.max_sunalt),
                   'min_time={}'.format(self.min_time),
                   'max_moon={}'.format(self.max_moon),
                   'min_moonsep={}'.format(self.min_moonsep),
                   'too={}'.format(self.too),
                   'tel_mask={}'.format(self.tel_mask),
                   'start_time={}'.format(self.start_time),
                   'stop_time={}'.format(self.stop_time),
                   'creation_time={}'.format(self.creation_time),
                   'deleted_time={}'.format(self.deleted_time),
                   'user_id={}'.format(self.user_id),
                   'grid_tile_id={}'.format(self.grid_tile_id),
                   'survey_tile_id={}'.format(self.survey_tile_id),
                   'event_id={}'.format(self.event_id),
                   ]
        return 'Mpointing({})'.format(', '.join(strings))

    @validates('start_time', 'stop_time', 'creation_time', 'deleted_time')
    def validate_times(self, key, field):
        """Use validators to allow various types of input for times.

        Also enforce stop_time > start_time.
        """
        if key not in ['start_time', 'creation_time'] and field is None:
            # start_time and creation_time are not nullable, the others are
            return None

        if isinstance(field, datetime.datetime):
            value = field.strftime('%Y-%m-%d %H:%M:%S')
        elif isinstance(field, Time):
            field.precision = 0  # no D.P on seconds
            value = field.iso
        else:
            # just hope the string works!
            value = str(field)

        # force start_time > stop_time
        if (key == 'start_time' and self.stop_time is not None and
                Time(value) >= Time(self.stop_time)):
            raise ValueError(f'start_time must be before stop_time ({self.stop_time})')
        if (key == 'stop_time' and self.start_time is not None and
                Time(self.start_time) >= Time(value)):
            raise ValueError(f'stop_time must be after start_time ({self.start_time})')

        return value

    @hybrid_property
    def scheduled(self):
        """Return True if currently linked to a pending or running Pointing."""
        return any(p.status in ['upcoming', 'pending', 'running'] for p in self.pointings)

    @scheduled.expression
    def scheduled(self):
        return exists().where(and_(Pointing.mpointing_id == self.db_id,
                                   or_(Pointing.status == 'upcoming',
                                       Pointing.status == 'pending',
                                       Pointing.status == 'running',
                                       )
                                   ))

    @hybrid_method
    def scheduled_at_time(self, time):
        """Return True if linked to a pending or running Pointing at the given time."""
        return any(p.status_at_time(time) in ['upcoming', 'pending', 'running']
                   for p in self.pointings)

    @scheduled_at_time.expression
    def scheduled_at_time(self, time):
        return exists().where(and_(Pointing.mpointing_id == self.db_id,
                                   or_(Pointing.status_at_time(time) == 'upcoming',
                                       Pointing.status_at_time(time) == 'pending',
                                       Pointing.status_at_time(time) == 'running',
                                       )
                                   ))

    @hybrid_property
    def infinite(self):
        """Return if the Mpointing is infinite."""
        return self.num_todo < 0

    @hybrid_property
    def num_completed(self):
        """Return the number of Pointings successfully completed."""
        return sum(p.status == 'completed' for p in self.pointings)

    @num_completed.expression
    def num_completed(self):
        return select([func.count(Pointing.db_id)]).\
            where(and_(Pointing.mpointing_id == self.db_id,
                       Pointing.status == 'completed',
                       )).scalar_subquery()

    @hybrid_method
    def num_completed_at_time(self, time):
        """Return the number of Pointings successfully completed at the given time."""
        return sum(p.status_at_time(time) == 'completed' for p in self.pointings)

    @num_completed_at_time.expression
    def num_completed_at_time(self, time):
        return select([func.count(Pointing.db_id)]).\
            where(and_(Pointing.mpointing_id == self.db_id,
                       Pointing.status_at_time(time) == 'completed',
                       ))

    @hybrid_property
    def num_remaining(self):
        """Return the number of observations remaining (num_todo - num_completed).

        Returns -1 for infinite Mpointings.
        """
        if self.infinite:
            return -1
        else:
            return self.num_todo - self.num_completed

    @num_remaining.expression
    def num_remaining(self):
        return case([(self.infinite, -1)],
                    else_=self.num_todo - self.num_completed)

    @hybrid_method
    def num_remaining_at_time(self, time):
        """Return the number of observations remaining at the given time."""
        if self.infinite:
            return -1
        else:
            return self.num_todo - self.num_completed_at_time(time)

    @num_remaining_at_time.expression
    def num_remaining_at_time(self, time):
        return case([(self.infinite, -1)],
                    else_=self.num_todo - self.num_completed_at_time(time))

    @hybrid_property
    def status(self):
        """Return a string giving the current status."""
        time = Time.now()
        if self.deleted_time is not None:
            # The Mpointing has been deleted before it completed: deleted
            return 'deleted'
        elif self.num_remaining == 0:
            # The Mpointing has enough completed Pointings: completed
            return 'completed'
        elif self.start_time is not None and time < Time(self.start_time):
            # The Mpointing hasn't yet reached its start time: upcoming
            return 'upcoming'
        elif self.stop_time is not None and time >= Time(self.stop_time):
            # The Mpointing has passed its start time: expired
            return 'expired'
        elif not self.scheduled:
            # The Mpointing doesn't have a pending Pointing: unscheduled
            return 'unscheduled'
        else:
            # The Mpointing has a pending Pointing: scheduled
            return 'scheduled'

    @status.setter
    def status(self, value):
        if value == 'deleted':
            self.mark_deleted()
        elif value == 'completed':
            raise ValueError('Can not set status to "completed", mark Pointings instead')
        elif value in ['upcoming', 'expired']:
            raise ValueError(f'Can not set status to "{value}", change start/stop times instead')
        elif value == 'unscheduled':
            raise ValueError('Can not set status to "unscheduled", mark Pointings instead')
        elif value == 'scheduled':
            raise ValueError('Can not set status to "scheduled", use `get_next_pointing()`')
        else:
            raise ValueError(f'Invalid status: "{value}"')

    @status.expression
    def status(self):
        time = datetime.datetime.utcnow()
        return case([(self.deleted_time.isnot(None),
                      'deleted'),
                     (self.num_remaining == 0,
                      'completed'),
                     (and_(self.start_time.isnot(None), time < self.start_time),
                      'upcoming'),
                     (and_(self.stop_time.isnot(None), time >= self.stop_time),
                      'expired'),
                     (self.scheduled.is_(False),
                      'unscheduled'),
                     ],
                    else_='scheduled')

    @hybrid_method
    def status_at_time(self, time):
        """Return a string giving the status at the given time (for testing and simulations)."""
        # This is very useful for running through simulations.
        # Then how do you deal with Pointings which didn't exist at that time but were created
        # afterwards? Do we need a Pointing.creation_time? We could have it set on creation,
        # but for simulations we need to fake that time.
        if isinstance(time, (str, datetime.datetime)):
            time = Time(time)

        if time < Time(self.creation_time):
            # This shouldn't usually happen, unless we are doing simulations and set it manually,
            # then it's worth knowing if this entry wouldn't have been created yet..
            return None

        if (self.deleted_time is not None and time >= Time(self.deleted_time)):
            return 'deleted'
        elif self.num_remaining_at_time(time) == 0:
            return 'completed'
        elif self.start_time is not None and time < Time(self.start_time):
            return 'upcoming'
        elif self.stop_time is not None and time >= Time(self.stop_time):
            return 'expired'
        elif not self.scheduled_at_time(time):
            return 'unscheduled'
        else:
            return 'scheduled'

    @status_at_time.expression
    def status_at_time(self, time):
        if isinstance(time, str):
            time = Time(time)
        if isinstance(time, Time):
            time = time.datetime

        c = case([(time < self.creation_time,
                   None),
                  (and_(self.deleted_time.isnot(None), time >= self.deleted_time),
                   'deleted'),
                  (self.num_remaining_at_time(time) == 0,
                   'completed'),
                  (and_(self.start_time.isnot(None), time < self.start_time),
                   'upcoming'),
                  (and_(self.stop_time.isnot(None), time >= self.stop_time),
                   'expired'),
                  (self.scheduled_at_time(time).is_(False),
                   'unscheduled'),
                  ],
                 else_='scheduled')
        return c

    def mark_deleted(self, time=None):
        """Mark this Mpointing (and any pending Pointings) as deleted."""
        if time is None:
            time = Time.now()
        if isinstance(time, (str, datetime.datetime)):
            time = Time(time)

        if self.status_at_time(time) == 'deleted':
            raise ValueError(f'Mpointing is already deleted (at={self.deleted_time})')
        if self.status_at_time(time) == 'expired':
            raise ValueError(f'Mpointing is already expired (at {self.stop_time})')
        if self.status_at_time(time) == 'completed':
            raise ValueError('Mpointing is already completed')

        # Pointings can use finished_time, but we need an explicit deleted_time
        self.deleted_time = time

        for pointing in self.pointings:
            if pointing.status in ['upcoming', 'pending']:
                pointing.mark_deleted(time=time)

    @property
    def valid_telescopes(self):
        """Get a list of valid Telescope IDs, if any, based on the tel_mask."""
        if self.tel_mask is None:
            return None
        tel_mask_str = str(bin(self.tel_mask))[2:]
        valid_dict = {i + 1: bool(int(d)) for i, d in enumerate(tel_mask_str[::-1])}
        return [tel for tel in valid_dict if valid_dict[tel] is True]

    @valid_telescopes.setter
    def valid_telescopes(self, telescope_ids):
        if telescope_ids is None:
            self.tel_mask = None
        else:
            if not isinstance(telescope_ids, (list, tuple)):
                telescope_ids = [telescope_ids]
            telescope_ids = {int(t) for t in telescope_ids}
            tel_mask = sum(2 ** (t - 1) for t in telescope_ids)
            self.tel_mask = tel_mask

    @hybrid_property
    def current_rank(self):
        """Calculate current rank as the starting rank + 10 * the number of completed Pointings."""
        return self.start_rank + 10 * self.num_completed

    @hybrid_method
    def current_rank_at_time(self, time):
        """Calculate current rank as the starting rank + 10 * the number of completed Pointings."""
        return self.start_rank + 10 * self.num_completed_at_time(time)

    def get_current_block(self):
        """Return the current time block.

        Assumes this object is still associated to an active session.

        Returns
        -------
        current_block : `TimeBlock`
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
        last : `TimeBlock`
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
        next_block : `TimeBlock`
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

    def get_next_pointing(self, time=None):
        """Retrieve the next pointing which needs to be scheduled.

        The start and stop times of the pointing are determined from
        the status of the previous pointing.

        Assumes this object is still associated with an active session.

        Any time given is set as the creation time, it should only be used for simulations
        (just like status_at_time() etc).

        Returns
        -------
        pointing : `Pointing` or None
            the next pointing that should be sent to the database. Is None
            if there is no suitable pointing remaining, or if there is
            a pointing from this Mpointing already scheduled

        """
        # already scheduled or finished, return None
        if self.status not in ['upcoming', 'unscheduled']:
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
                if current_block.valid_time > 0:
                    # Want to wait after the valid period of the latest pointing, which started at
                    # it's start time.
                    # We can't just use the stop_time, because it might have been set by the
                    # Mpointing's stop_time instead (see below)
                    start_time = (Time(latest_pointing.start_time) +
                                  current_block.valid_time * u.min +
                                  current_block.wait_time * u.min)
                else:
                    # non-expiring pointing, no valid_time
                    # (unless the mpointing had one, but we don't want to wait until after that)
                    # Start the wait time after the previous one stopped
                    start_time = (Time(latest_pointing.stopped_time) +
                                  current_block.wait_time * u.minute)
            else:
                # the current block wasn't completed, and there's still time left
                # (e.g. aborted, interrupted)
                # need to re-insert the current block with a new pointing
                start_time = latest_pointing.start_time

        if next_block.valid_time < 0 and not self.stop_time:
            # If valid time is infinite (and there's no Mpointing stop_time)
            # then the pointings should never expire.
            stop_time = None
        else:
            if next_block.valid_time < 0:
                # The Mpointing stop time takes priority over infinite valid time
                stop_time = self.stop_time
            else:
                stop_time = Time(start_time) + next_block.valid_time * u.minute
                if self.stop_time and stop_time > self.stop_time:
                    # force pointings to stop by an Mpointing's stop_time, if given
                    stop_time = self.stop_time

            if start_time >= stop_time:
                # can happen if the Mpointing has a stop_time
                return None

        # mark the new block as current
        current_block.current = False
        next_block.current = True

        # now create a pointing
        p = Pointing(rank=self.current_rank,
                     start_time=start_time,
                     stop_time=stop_time,
                     mpointing=self,
                     time_block=next_block,
                     )
        # add the exposures
        p.exposure_sets = self.exposure_sets
        # set the creation time (for simulations)
        if time is not None:
            p.creation_time = time

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

    When created the instance can be linked to the following other tables as parameters,
    otherwise they are populated when it is added to the database:

    Primary relationships
    ---------------------
    mpointing : `Mpointing`
        the Mpointing associated with this Time Block
        required before addition to the database
        can also be added with the mpointing_id parameter

    pointings : list of `Pointing`, optional
        the Pointings associated with this Time Block, if any

    Attributes
    ----------
    db_id : int
        primary database key
        only populated when the instance is added to the database

    """

    # Set corresponding SQL table name
    __tablename__ = 'time_blocks'

    # Primary key
    db_id = Column('id', Integer, primary_key=True)

    # Columns
    block_num = Column(Integer, nullable=False, index=True)
    valid_time = Column(Float, nullable=False)
    wait_time = Column(Float, nullable=False)
    current = Column(Boolean, nullable=False, index=True, default=False)

    # Foreign keys
    mpointing_id = Column(Integer, ForeignKey('mpointings.id'), nullable=False)

    # Foreign relationships
    mpointing = relationship('Mpointing', back_populates='time_blocks')
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


class Site(Base):
    """A class to represent an observing Site.

    Sites are used to group Telescopes, which are then linked to Pointings to observe.

    Like all SQLAlchemy model classes, this object links to the
    underlying database. You can create an instance, and set its attributes
    without a database session. Accessing some attributes may require
    an active database session, and some properties (like the db_id)
    will be None until the instance is added to the database.

    Parameters
    ----------
    name : str
        name of this site (must be unique)
    latitude : float
        latitude of this site, in degrees
    longitude : float
        longitude of this site, in degrees
    height : float
        height of this site, in metres above sea level

    When created the instance can be linked to the following other tables as parameters,
    otherwise they are populated when it is added to the database:

    Primary relationships
    ---------------------
    telescopes : list of `Telescope`, optional
        the Telescopes located at this Site, if any

    Attributes
    ----------
    db_id : int
        primary database key
        only populated when the instance is added to the database
    location : `astropy.cooridnates.EarthLocation`
        Astropy EarthLocation class for this Site
    tel_mask : int
        the binary telescope mask for the telescopes at this site

    Class methods
    -------------
    from_location(`Astropy.coordinates.EarthLocation`) : Site
        create a Site from an Astropy EarthLocation
    from_name(string) : Site
        class method to create a Site from any name recognised by `EarthLocation.of_site()`

    """

    # Set corresponding SQL table name
    __tablename__ = 'sites'

    # Primary key
    db_id = Column('id', Integer, primary_key=True)

    # Columns
    name = Column(String(255), nullable=False, unique=True, index=True)
    latitude = Column(Float, nullable=False)
    longitude = Column(Float, nullable=False)
    height = Column(Float, nullable=False)

    # Foreign relationships
    telescopes = relationship('Telescope', back_populates='site')

    def __repr__(self):
        strings = ['db_id={}'.format(self.db_id),
                   'name="{}"'.format(self.name),
                   'latitude={}'.format(self.latitude),
                   'longitude={}'.format(self.longitude),
                   'height={}'.format(self.height),
                   ]
        return 'Site({})'.format(', '.join(strings))

    @classmethod
    def from_location(cls, location, name=None):
        """Create a Site from an Astropy `EarthLocation` class and a name string."""
        if not isinstance(location, EarthLocation):
            raise ValueError('"location" must be an `astropy.coordinates.EarthLocation`')
        if hasattr(location.info, 'name'):
            name = location.info.name
        if name is None:
            raise ValueError('Missing name for site')

        return cls(name=name,
                   latitude=location.lat.value,
                   longitude=location.lon.value,
                   height=location.height.value,
                   )

    @classmethod
    def from_name(cls, name):
        """Create a Site from a name recognised by Astropy's `EarthLocation.of_site()`."""
        location = EarthLocation.of_site(name)
        return cls.from_location(location)

    @property
    def location(self):
        """Return an Astropy EarthLocation for this Site."""
        return EarthLocation(self.latitude * u.deg, self.longitude * u.deg, self.height * u.m)

    @property
    def tel_mask(self):
        """Get the binary telescope mask for the telescopes at this site."""
        return sum(t.tel_mask for t in self.telescopes)  # Isn't binary neat!


class Telescope(Base):
    """A class to represent a Telescope to observe Pointings.

    Telescopes should be linked to a Site, and can be linked to an observing Grid.

    Like all SQLAlchemy model classes, this object links to the
    underlying database. You can create an instance, and set its attributes
    without a database session. Accessing some attributes may require
    an active database session, and some properties (like the db_id)
    will be None until the instance is added to the database.

    Parameters
    ----------
    name : str
        name of this telescope (must be unique)

    When created the instance can be linked to the following other tables as parameters,
    otherwise they are populated when it is added to the database:

    Primary relationships
    ---------------------
    site : `Site`
        the Site this Telescope is located at
        required before addition to the database
        can also be added with the site_id parameter

    grid : `Grid`, optional
        the Grid this Telescope uses to observe, if any
        can also be added with the grid_id parameter
    pointings : list of `Pointing`, optional
        the Pointings observed by this Telescope, if any

    Attributes
    ----------
    db_id : int
        primary database key
        only populated when the instance is added to the database
    status : string
        the status of this telescope (either 'observing' or 'idle')
    current_pointing : `Pointing`
        the Pointing currently being observed by this Telescope, if any
    tel_mask : int
        the binary telescope mask for this Telescope, based on the db_id
        only populated when the instance is added to the database
    """

    # Set corresponding SQL table name
    __tablename__ = 'telescopes'

    # Primary key
    db_id = Column('id', Integer, primary_key=True)

    # Columns
    name = Column(String(255), nullable=False, unique=True, index=True)

    # Foreign keys
    site_id = Column(Integer, ForeignKey('sites.id'), nullable=False)
    grid_id = Column(Integer, ForeignKey('grids.id'), nullable=True)

    # Foreign relationships
    site = relationship('Site', back_populates='telescopes')
    grid = relationship('Grid', back_populates='telescopes')
    pointings = relationship('Pointing', back_populates='telescope')
    current_pointing = relationship('Pointing', viewonly=True, uselist=False,
                                    primaryjoin='and_(Telescope.db_id==Pointing.telescope_id,'
                                                'Pointing.status=="running")')

    def __repr__(self):
        strings = ['db_id={}'.format(self.db_id),
                   'status={}'.format(self.status),
                   'name="{}"'.format(self.name),
                   'site_id={}'.format(self.site_id),
                   'grid_id={}'.format(self.grid_id),
                   ]
        return 'Telescope({})'.format(', '.join(strings))

    @hybrid_property
    def status(self):
        """Return a string giving the current status."""
        if self.current_pointing is not None:
            return 'observing'
        else:
            return 'idle'

    @status.expression
    def status(self):
        return case([(self.current_pointing.isnot(None), 'observing')],
                    else_='idle')

    @property
    def tel_mask(self):
        """Get the binary telescope mask for this Telescope."""
        if self.db_id is None:
            raise ValueError('Telescope needs to be added to database to assign ID')
        return 2 ** (self.db_id - 1)


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

    When created the instance can be linked to the following other tables as parameters,
    otherwise they are populated when it is added to the database:

    Primary relationships
    ---------------------
    grid_tiles : list of `GridTile`, optional
        the Grid Tiles that are part of this Grid, if any
    surveys : list of `Survey`, optional
        the Surveys observed using this Grid, if any
    telescopes : list of `Telescope`, optional
        the Telescopes which observe using this Grid, if any

    Attributes
    ----------
    db_id : int
        primary database key
        only populated when the instance is added to the database

    The following secondary relationships are not settable directly,
    but are populated through the primary relationships and are available as attributes:

    Secondary relationships
    -----------------------
    mpointings : list of `Mpointing`
        the Mpointings which cover any of this Grid's GridTiles, if any

    """

    # Set corresponding SQL table name
    __tablename__ = 'grids'

    # Primary key
    db_id = Column('id', Integer, primary_key=True)

    # Columns
    name = Column(String(255), nullable=False, unique=True, index=True)
    ra_fov = Column(Float, nullable=False)
    dec_fov = Column(Float, nullable=False)
    ra_overlap = Column(Float, nullable=False)
    dec_overlap = Column(Float, nullable=False)
    algorithm = Column(String(255), nullable=False)

    # Foreign relationships
    grid_tiles = relationship('GridTile', back_populates='grid')
    surveys = relationship('Survey', back_populates='grid')
    telescopes = relationship('Telescope', back_populates='grid')

    # Secondary relationships
    mpointings = relationship('Mpointing',
                              secondary='grid_tiles',
                              primaryjoin='Grid.db_id == GridTile.grid_id',
                              secondaryjoin='Mpointing.grid_tile_id == GridTile.db_id',
                              back_populates='grid',
                              viewonly=True,
                              uselist=True,
                              )

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

    def get_skygrid(self):
        """Create a GOTO-tile SkyGrid from the current database Grid."""
        fov = {'ra': self.ra_fov * u.deg, 'dec': self.dec_fov * u.deg}
        overlap = {'ra': self.ra_overlap, 'dec': self.dec_overlap}
        skygrid = SkyGrid(fov, overlap, kind=self.algorithm)
        return skygrid


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

    When created the instance can be linked to the following other tables as parameters,
    otherwise they are populated when it is added to the database:

    Primary relationships
    ---------------------
    grid : `Grid`
        the Grid this GridTile is part of
        required before addition to the database
        can also be added with the grid_id parameter

    mpointings : list of `Mpointing`, optional
        the Mpointings that cover this GridTile, if any
    survey_tiles : list of `SurveyTile`, optional
        the SurveyTiles that cover this GridTile, if any

    Attributes
    ----------
    db_id : int
        primary database key
        only populated when the instance is added to the database

    The following secondary relationships are not settable directly,
    but are populated through the primary relationships and are available as attributes:

    Secondary relationships
    -----------------------
    pointings : list of `Pointing`
        the Pointings that cover this GridTile, if any

    """

    # Set corresponding SQL table name
    __tablename__ = 'grid_tiles'

    # Primary key
    db_id = Column('id', Integer, primary_key=True)

    # Columns
    name = Column(String(255), nullable=False, index=True)
    ra = Column(Float, nullable=False)
    dec = Column(Float, nullable=False)

    # Foreign keys
    grid_id = Column(Integer, ForeignKey('grids.id'), nullable=False)

    # Foreign relationships
    grid = relationship('Grid', back_populates='grid_tiles')
    survey_tiles = relationship('SurveyTile', back_populates='grid_tile')
    mpointings = relationship('Mpointing', back_populates='grid_tile')

    # Secondary relationships
    pointings = relationship('Pointing',
                             secondary='mpointings',
                             primaryjoin='GridTile.db_id == Mpointing.grid_tile_id',
                             secondaryjoin='Pointing.mpointing_id == Mpointing.db_id',
                             back_populates='grid_tile',
                             viewonly=True,
                             uselist=True,
                             )

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

    When created the instance can be linked to the following other tables as parameters,
    otherwise they are populated when it is added to the database:

    Primary relationships
    ---------------------
    grid : `Grid`
        the Grid associated with this Survey
        required before addition to the database
        can also be added with the grid_id parameter

    survey_tiles : list of `SurveyTile`, optional
        the Survey Tiles that are part of this Survey, if any
    event : `Event`, optional
        the Event this Survey is part of, if any
        can also be added with the event_id parameter

    Attributes
    ----------
    db_id : int
        primary database key
        only populated when the instance is added to the database

    The following secondary relationships are not settable directly,
    but are populated through the primary relationships and are available as attributes:

    Secondary relationships
    -----------------------
    mpointings : list of `Mpointing`
        the Mpointings which cover any of this Survey's SurveyTiles, if any

    """

    # Set corresponding SQL table name
    __tablename__ = 'surveys'

    # Primary key
    db_id = Column('id', Integer, primary_key=True)

    # Columns
    name = Column(String(255), nullable=False, index=True)

    # Foreign keys
    grid_id = Column(Integer, ForeignKey('grids.id'), nullable=False)
    event_id = Column(Integer, ForeignKey('events.id'), nullable=True)

    # Foreign relationships
    survey_tiles = relationship('SurveyTile', back_populates='survey')
    grid = relationship('Grid', back_populates='surveys')
    event = relationship('Event', back_populates='surveys')

    # Secondary relationships
    mpointings = relationship('Mpointing',
                              secondary='survey_tiles',
                              primaryjoin='Survey.db_id == SurveyTile.survey_id',
                              secondaryjoin='Mpointing.survey_tile_id == SurveyTile.db_id',
                              back_populates='survey',
                              viewonly=True,
                              uselist=True,
                              )

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
        weighting for this tile (between 0 and 1)

    When created the instance can be linked to the following other tables as parameters,
    otherwise they are populated when it is added to the database:

    Primary relationships
    ---------------------
    Survey : `Survey`
        the Survey this SurveyTile is part of
        required before addition to the database
        can also be added with the survey_id parameter

    mpointings : list of `Mpointing`, optional
        the Mpointings that cover this SurveyTile, if any
    grid_tiles : list of `GridTile`, optional
        the GridTiles that cover this SurveyTile, if any

    Attributes
    ----------
    db_id : int
        primary database key
        only populated when the instance is added to the database

    The following secondary relationships are not settable directly,
    but are populated through the primary relationships and are available as attributes:

    Secondary relationships
    -----------------------
    pointings : list of `Pointing`
        the Pointings that cover this SurveyTile, if any

    """

    # Set corresponding SQL table name
    __tablename__ = 'survey_tiles'

    # Primary key
    db_id = Column('id', Integer, primary_key=True)

    # Columns
    weight = Column(Float, nullable=False)

    # Foreign keys
    survey_id = Column(Integer, ForeignKey('surveys.id'), nullable=False)
    grid_tile_id = Column(Integer, ForeignKey('grid_tiles.id'), nullable=False)

    # Foreign relationships
    survey = relationship('Survey', back_populates='survey_tiles')
    grid_tile = relationship('GridTile', back_populates='survey_tiles')
    mpointings = relationship('Mpointing', back_populates='survey_tile')

    # Secondary relationships
    pointings = relationship('Pointing',
                             secondary='mpointings',
                             primaryjoin='SurveyTile.db_id == Mpointing.survey_tile_id',
                             secondaryjoin='Pointing.mpointing_id == Mpointing.db_id',
                             back_populates='survey_tile',
                             viewonly=True,
                             uselist=True,
                             )

    def __repr__(self):
        strings = ['db_id={}'.format(self.db_id),
                   'weight={}'.format(self.weight),
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
        time the event occurred
    skymap : string, optional
        the location of the source skymap file

    When created the instance can be linked to the following other tables as parameters,
    otherwise they are populated when it is added to the database:

    Primary relationships
    ---------------------
    surveys : list of `Survey`, optional
        the Surveys created for this Event, if any
    mpointings : list of `Mpointing`, optional
        the Mpointings which are part of this Event, if any

    Attributes
    ----------
    db_id : int
        primary database key
        only populated when the instance is added to the database

    The following secondary relationships are not settable directly,
    but are populated through the primary relationships and are available as attributes:

    Secondary relationships
    -----------------------
    pointings : list of `Pointing`
        the Pointings which are part of this Event, if any

    """

    # Set corresponding SQL table name
    __tablename__ = 'events'

    # Primary key
    db_id = Column('id', Integer, primary_key=True)

    # Columns
    name = Column(String(255), nullable=False, index=True)
    ivorn = Column(String(255), nullable=False, unique=True)
    source = Column(String(255), nullable=False)
    # type is a built in function in Python
    event_type = Column('type', String(255), nullable=False, index=True)
    time = Column(DateTime)
    skymap = Column(String(255))

    # Foreign relationships
    surveys = relationship('Survey', back_populates='event')
    mpointings = relationship('Mpointing', back_populates='event')

    # Secondary relationships
    pointings = relationship('Pointing',
                             secondary='mpointings',
                             primaryjoin='Event.db_id == Mpointing.event_id',
                             secondaryjoin='Pointing.mpointing_id == Mpointing.db_id',
                             back_populates='event',
                             viewonly=True,
                             uselist=True,
                             )

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
    def validate_times(self, key, field):
        """Use validators to allow various types of input for times."""
        if isinstance(field, datetime.datetime):
            value = field.strftime('%Y-%m-%d %H:%M:%S')
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

    When created the instance can be linked to the following other tables as parameters,
    otherwise they are populated when it is added to the database:

    Primary relationships
    ---------------------
    exposure_set : `ExposureSet`, optional
        the Exposure Set associated with this ImageLog, if any
        can also be added with the exposure_set_id parameter
    pointing : `Pointing`, optional
        the Pointing associated with this ImageLog, if any
        can also be added with the pointing_id parameter

    Attributes
    ----------
    db_id : int
        primary database key
        only populated when the instance is added to the database

    """

    # Set corresponding SQL table name
    __tablename__ = 'image_logs'

    # Primary key
    db_id = Column('id', Integer, primary_key=True)

    # Columns
    filename = Column(String(30), nullable=False, unique=True)
    run_number = Column(Integer, nullable=False, index=True)
    ut = Column(Integer, nullable=False)
    ut_mask = Column(Integer, nullable=False)
    start_time = Column(DateTime, nullable=False)
    write_time = Column(DateTime, nullable=False, index=True)
    set_position = Column(Integer, nullable=False, default=1)
    set_total = Column(Integer, nullable=False, default=1)

    # Foreign keys
    exposure_set_id = Column(Integer, ForeignKey('exposure_sets.id'), nullable=True)
    pointing_id = Column(Integer, ForeignKey('pointings.id'), nullable=True)

    # Foreign relationships
    exposure_set = relationship('ExposureSet', back_populates='image_logs')
    pointing = relationship('Pointing', back_populates='image_logs')

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
    def validate_times(self, key, field):
        """Use validators to allow various types of input for times.

        Also enforce write_time > start_time.
        """
        if isinstance(field, datetime.datetime):
            value = field.strftime('%Y-%m-%d %H:%M:%S')
        elif isinstance(field, Time):
            field.precision = 0  # no D.P on seconds
            value = field.iso
        else:
            # just hope the string works!
            value = str(field)

        # force write_time > start_time
        if (key == 'start_time' and self.write_time is not None and
                Time(value) >= Time(self.write_time)):
            raise ValueError(f'start_time must be before write_time ({self.write_time})')
        if (key == 'write_time' and self.start_time is not None and
                Time(self.start_time) >= Time(value)):
            raise ValueError(f'write_time must be after start_time ({self.start_time})')

        return value


TRIGGERS = [
    # Mpointing trigger before insert
    # If no coordinates are given but it's linked to a grid tile then use its coordinates.
    """CREATE TRIGGER `mpointings_BEFORE_INSERT` BEFORE INSERT ON `mpointings` FOR EACH ROW
    BEGIN
    IF ((NEW.`grid_tile_id` is not NULL) and (NEW.`ra` is NULL) and (NEW.`dec` is NULL)) THEN
        SET NEW.`ra` = (SELECT `ra` FROM `grid_tiles`
                        WHERE NEW.`grid_tile_id` = `grid_tiles`.`id`);
        SET NEW.`dec` = (SELECT `dec` FROM `grid_tiles`
                            WHERE NEW.`grid_tile_id` = `grid_tiles`.`id`);
    END IF;
    END
    """,
]
