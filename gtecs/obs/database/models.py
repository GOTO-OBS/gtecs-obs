"""Python classes mapping on to database tables."""

import datetime
import hashlib

from astropy import units as u
from astropy.coordinates import EarthLocation
from astropy.time import Time

from gototile.grid import SkyGrid

from sqlalchemy import Boolean, Column, DateTime, Float, ForeignKey, Integer, String, Text
from sqlalchemy import exists, func, select, text
from sqlalchemy.dialects.mysql import TIMESTAMP
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.ext.hybrid import hybrid_method, hybrid_property
from sqlalchemy.orm import column_property, relationship, validates
from sqlalchemy.sql import and_, case, or_

from .. import params


__all__ = ['User', 'ExposureSet', 'Pointing', 'Strategy', 'Target', 'TimeBlock',
           'Site', 'Telescope', 'Grid', 'GridTile',
           'Survey', 'Event',
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
    targets : list of `Target`, optional
        the Targets created by this User, if any

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
    targets = relationship('Target', back_populates='user')

    def __init__(self, username, password, full_name, **kwargs):
        kwargs['username'] = username
        kwargs['full_name'] = full_name

        # Use hashed password instead of plain
        kwargs['password_hash'] = hashlib.sha512(password.encode()).hexdigest()

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
    exposure time, filter, binning factor etc. Each Target should have one or more
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
    ut_mask : int, optional
        if set, this is a binary mask which will determine which unit
        telescopes carry out the exposure. A value of 5 (binary 0101) will
        be exposed by cameras 1 and 3.
        default = None

    When created the instance can be linked to the following other tables as parameters,
    otherwise they are populated when it is added to the database:

    Primary relationships
    ---------------------
    target : `Target`, optional
        the Target this Exposure Set is linked to, if any
        can also be added with the target_id parameter

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
        the Pointings associated with the Target this Exposure Set is linked to, if any

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
    ut_mask = Column(Integer, default=None)

    # Foreign keys
    target_id = Column(Integer, ForeignKey('targets.id'), nullable=True)

    # Foreign relationships
    target = relationship('Target', back_populates='exposure_sets')

    # Secondary relationships
    pointings = relationship('Pointing',
                             secondary='targets',
                             primaryjoin='ExposureSet.target_id == Target.db_id',
                             secondaryjoin='Pointing.target_id == Target.db_id',
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
                   'ut_mask={}'.format(self.ut_mask),
                   'target_id={}'.format(self.target_id),
                   ]
        return 'ExposureSet({})'.format(', '.join(strings))


class Pointing(Base):
    """A class to represent a Pointing.

    Pointings are the primary table of the Observation Database, and are usually generated
    automatically by Targets. When deciding what to observe the scheduler processes a queue
    made of Pointings based on their status and constraints.

    NB: you shouldn't ever have to manually create Pointings, instead you should generate them
    from a Target using its get_next_pointing() method. The first Pointing is also created when
    the Target is.

    Like all SQLAlchemy model classes, this object links to the
    underlying database. You can create an instance, and set its attributes
    without a database session. Accessing some attributes may require
    an active database session, and some properties (like the db_id)
    will be None until the instance is added to the database.

    Parameters
    ----------
    rank : int
        rank to use when scheduling (should be the current_rank of the linked target)

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
    target : `Target`
        the Target associated with this Pointing
        can also be added with the target_id parameter

    time_block : `TimeBlock`, optional
        the TimeBlock this Pointing as generated from, if any
        can also be added with the time_block_id parameter
    strategy : `Strategy`, optional
        the Strategy defined for this Pointing, if any
        can also be added with the strategy_id parameter
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
        the ExposureSets linked to the Target this Pointing was generated from
    grid_tile : `GridTile`
        the GridTile the Target this Pointing was generated from covers, if any
    survey : `Survey`
        the Survey the Target this Pointing was generated from was part of, if any
    event : `Event`
        the Event the Target this Pointing was generated from was part of, if any

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
    target_id = Column(Integer, ForeignKey('targets.id'), nullable=False)
    time_block_id = Column(Integer, ForeignKey('time_blocks.id'), nullable=False)
    strategy_id = Column(Integer, ForeignKey('strategies.id'), nullable=False)
    telescope_id = Column(Integer, ForeignKey('telescopes.id'), nullable=True)

    # Foreign relationships
    target = relationship('Target', lazy='joined', back_populates='pointings')
    time_block = relationship('TimeBlock', lazy='joined', back_populates='pointings')
    strategy = relationship('Strategy', lazy='joined', back_populates='pointings')
    telescope = relationship('Telescope', lazy='joined', back_populates='pointings')

    # Secondary relationships
    exposure_sets = relationship('ExposureSet',
                                 secondary='targets',
                                 primaryjoin='Pointing.target_id == Target.db_id',
                                 secondaryjoin='ExposureSet.target_id == Target.db_id',
                                 back_populates='pointings',
                                 viewonly=True,
                                 uselist=True,
                                 )
    grid_tile = relationship('GridTile',
                             lazy='joined',
                             secondary='targets',
                             primaryjoin='Pointing.target_id == Target.db_id',
                             secondaryjoin='GridTile.db_id == Target.grid_tile_id',
                             back_populates='pointings',
                             viewonly=True,
                             uselist=False,
                             )
    survey = relationship('Survey',
                          lazy='joined',
                          secondary='targets',
                          primaryjoin='Pointing.target_id == Target.db_id',
                          secondaryjoin='Survey.db_id == Target.event_id',
                          back_populates='pointings',
                          viewonly=True,
                          uselist=False,
                          )
    event = relationship('Event',
                         lazy='joined',
                         secondary='targets',
                         primaryjoin='Pointing.target_id == Target.db_id',
                         secondaryjoin='Event.db_id == Target.event_id',
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
                   'target_id={}'.format(self.target_id),
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
            # The Pointing has passed its stop time: expired
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

        if (self.strategy.valid_telescopes is not None and
                telescope.db_id not in self.strategy.valid_telescopes):
            raise ValueError(f'Telescope ID ({telescope.db_id}) not in '
                             f'list of valid telescopes ({self.strategy.valid_telescopes})')
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

    def get_obstime(self, readout_time=15):
        """Get the expected time needed to observe this Pointing."""
        return sum(es.num_exp * (es.exptime + readout_time) for es in self.exposure_sets)


class Target(Base):
    """A class to represent an astronomical Target.

    Targets are generation engines for Pointings, and allow the
    database to regenerate Pointings at given time intervals.
    Targets must be linked to one or more Strategies in order to
    generate new Pointings.

    Like all SQLAlchemy model classes, this object links to the
    underlying database. You can create an instance, and set its attributes
    without a database session. Accessing some attributes may require
    an active database session, and some properties (like the db_id)
    will be None until the instance is added to the database.

    Parameters
    ----------
    name : str
        object name
    ra : float, optional
        J2000 right ascension in decimal degrees
        if ra is not given and this Target is linked to a GridTile
        then the ra will be extracted from the GridTile
    dec : float, optional
        J2000 declination in decimal degrees
        if dec is not given and this Target is linked to a GridTile
        then the dec will be extracted from the GridTile
    start_rank : int
        rank to use for first Pointing generated by this Target
        lower values are prioritised first by the scheduler

    rank_decay : bool, optional
        if True, the Target's `current_rank` will increase by 10 for each completed Pointing
        default = True
    weight : int, optional
        weighting relative to other Targets in the same survey
        default = 1

    start_time : string, `astropy.time.Time` or datetime.datetime, optional
        UTC time from which Target is considered valid and can be started.
        if not given then set to now, so the Target will start immediately.
        default = Time.now()
    stop_time : string, `astropy.time.Time` or datetime.datetime, optional
        the latest UTC time after which Pointings must stop.
        if not given the Target will continue creating Pointings until
        it is completed.
        default = None

    creation_time : string, `astropy.time.Time` or datetime.datetime, optional
        the time the Target was created
        this should usually be set automatically, unless doing simulations
        default = Time.now()
    deleted_time : string, `astropy.time.Time` or datetime.datetime, or None, optional
        the time the Target was marked as deleted
        this should usually be set automatically with `Target.mark_deleted()`
        default = None

    When created the instance can be linked to the following other tables as parameters,
    otherwise they are populated when it is added to the database:

    Primary relationships
    ---------------------
    user : `User`
        the User associated with this Target
        required before addition to the database
        can also be added with the user_id parameter

    exposure_sets : list of `ExposureSet`, optional
        the Exposure Sets associated with this Target, if any
    strategies : list of `Strategy`, optional
        the Strategies associated with this Target, if any
    pointings : list of `Pointing`, optional
        the Pointings generated by this Target, if any
    grid_tile : `GridTile`, optional
        the GridTile this Target covers, if any
        can also be added with the grid_tile_id parameter
    survey : `Survey`, optional
        the Survey this Target is part of, if any
        can also be added with the survey_id parameter
    event : `Event`, optional
        the Event this Target is part of, if any
        can also be added with the event_id parameter

    Attributes
    ----------
    db_id : int
        primary database key
        only populated when the instance is added to the database
    status : string
        status of the target
        see also the `status_at_time()` method
    scheduled : bool
        returns True if the Target is scheduled (has a pending Pointing)
        see also the `scheduled_at_time()` method
    current_rank : int
        rank for next Pointing to be scheduled (it will increase as Pointings are observed)
        see also the `current_rank_at_time()` method
    num_completed : int
        number of successfully completed Pointings generated from this Target
        (across all Strategies)
        see also the `num_completed_at_time()` method
    last_observed : `astropy.time.Time`
        the most recent time of completion of any Pointing generated from this Target, if any
    strategy : `Strategy`
        the current Strategy for this Target, if any
        (see get_current_strategy())

    The following secondary relationships are not settable directly,
    but are populated through the primary relationships and are available as attributes:

    Secondary relationships
    -----------------------
    grid : `Grid`
        the Grid that the GridTile this Target covers, if any, is part of

    Methods
    -------
    mark_deleted(time=None)
        mark the pointing as deleted
        if no time is given, default to Time.now()
    get_current_strategy(time=None) : Strategy
        get the currently valid Strategy, or at the given time
    get_next_pointing(time=None) : Pointing, or None
        generate the next Pointing for this Target
        returns None if a Pointing is already scheduled, or if the Target is completed
        the time is set as the Pointing creation time, and should usually be set automatically
        unless doing simulations

    """

    # Set corresponding SQL table name
    __tablename__ = 'targets'

    # Primary key
    db_id = Column('id', Integer, primary_key=True)

    # Columns
    # # Basic properties
    name = Column(Text, nullable=False)
    ra = Column(Float, nullable=False)
    dec = Column(Float, nullable=False)
    start_rank = Column(Integer, nullable=False)
    rank_decay = Column(Boolean, nullable=False, default=True)
    weight = Column(Integer, nullable=False, default=1)
    # # Constraints
    start_time = Column(DateTime, nullable=False, index=True, server_default=func.now())
    stop_time = Column(DateTime, nullable=True, index=True, default=None)
    # # Status
    creation_time = Column(DateTime, nullable=False, server_default=func.now())
    deleted_time = Column(DateTime, nullable=True, default=None)

    # Foreign keys
    user_id = Column(Integer, ForeignKey('users.id'), nullable=False)
    grid_tile_id = Column(Integer, ForeignKey('grid_tiles.id'), nullable=True)
    survey_id = Column(Integer, ForeignKey('surveys.id'), nullable=True)
    event_id = Column(Integer, ForeignKey('events.id'), nullable=True)

    # Foreign relationships
    # (remember to add to __init__)
    user = relationship('User', lazy='joined', back_populates='targets')
    pointings = relationship('Pointing', lazy='joined', back_populates='target')
    exposure_sets = relationship('ExposureSet', lazy='joined', back_populates='target')
    strategies = relationship('Strategy', lazy='joined', back_populates='target')
    grid_tile = relationship('GridTile', lazy='joined', back_populates='targets')
    survey = relationship('Survey', lazy='joined', back_populates='targets')
    event = relationship('Event', lazy='joined', back_populates='targets')

    # Secondary relationships
    grid = relationship('Grid',
                        lazy='joined',
                        secondary='grid_tiles',
                        primaryjoin='Target.grid_tile_id == GridTile.db_id',
                        secondaryjoin='Grid.db_id == GridTile.grid_id',
                        back_populates='targets',
                        viewonly=True,
                        uselist=False,
                        )

    # Column properties
    last_observed = column_property(select([Pointing.finished_time])
                                    .where(and_(Pointing.target_id == db_id,
                                                Pointing.status == 'completed',
                                                )
                                           )
                                    .order_by(Pointing.finished_time.desc())
                                    .limit(1)
                                    .correlate_except(Pointing)
                                    .scalar_subquery()
                                    )

    def __init__(self, **kwargs):
        # Set default times
        if 'creation_time' not in kwargs or kwargs['creation_time'] is None:
            kwargs['creation_time'] = Time.now()
        if 'start_time' not in kwargs or kwargs['start_time'] is None:
            kwargs['start_time'] = kwargs['creation_time']

        # Get coordinates from tile if not given
        if ('grid_tile' in kwargs and kwargs['grid_tile'] is not None and
                ('ra' not in kwargs or kwargs['ra'] is None) and
                ('dec' not in kwargs or kwargs['dec'] is None)):
            kwargs['ra'] = kwargs['grid_tile'].ra
            kwargs['dec'] = kwargs['grid_tile'].dec

        # Allow singular strategy
        if 'strategy' in kwargs:
            strategy = kwargs.pop('strategy')
            if 'strategies' not in kwargs:
                kwargs['strategies'] = [strategy]

        # Init base class
        super().__init__(**kwargs)

        # Create the first Pointing
        # This assumes the Target is already linked to a Strategy, which has TimeBlocks, or it
        # will return None.
        # Otherwise it's unlikely that it returns None, but in theory you could create a Strategy
        # with num_todo=0 or a stop_time in the past or some other contrivance.
        # Also note get_next_pointing() allows creating Pointings if the Target
        # status is 'upcoming', so we can create it now even if it's not valid yet.
        pointing = self.get_next_pointing(time=self.creation_time)
        if pointing is not None:
            self.pointings.append(pointing)

    def __repr__(self):
        strings = ['db_id={}'.format(self.db_id),
                   'status={}'.format(self.status),
                   'name={}'.format(self.name),
                   'ra={}'.format(self.ra),
                   'dec={}'.format(self.dec),
                   'current_rank={}'.format(self.current_rank),
                   'start_rank={}'.format(self.start_rank),
                   'weight={}'.format(self.weight),
                   'num_completed={}'.format(self.num_completed),
                   'start_time={}'.format(self.start_time),
                   'stop_time={}'.format(self.stop_time),
                   'creation_time={}'.format(self.creation_time),
                   'deleted_time={}'.format(self.deleted_time),
                   'user_id={}'.format(self.user_id),
                   'grid_tile_id={}'.format(self.grid_tile_id),
                   'survey_id={}'.format(self.survey_id),
                   'event_id={}'.format(self.event_id),
                   ]
        return 'Target({})'.format(', '.join(strings))

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
        return exists().where(and_(Pointing.target_id == self.db_id,
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
        return exists().where(and_(Pointing.target_id == self.db_id,
                                   or_(Pointing.status_at_time(time) == 'upcoming',
                                       Pointing.status_at_time(time) == 'pending',
                                       Pointing.status_at_time(time) == 'running',
                                       )
                                   ))

    @hybrid_property
    def num_completed(self):
        """Return the number of Pointings successfully completed."""
        return sum(p.status == 'completed' for p in self.pointings)

    @num_completed.expression
    def num_completed(self):
        return select([func.count(Pointing.db_id)]).\
            where(and_(Pointing.target_id == self.db_id,
                       Pointing.status == 'completed',
                       )).scalar_subquery()

    @hybrid_method
    def num_completed_at_time(self, time):
        """Return the number of Pointings successfully completed at the given time."""
        return sum(p.status_at_time(time) == 'completed' for p in self.pointings)

    @num_completed_at_time.expression
    def num_completed_at_time(self, time):
        return select([func.count(Pointing.db_id)]).\
            where(and_(Pointing.target_id == self.db_id,
                       Pointing.status_at_time(time) == 'completed',
                       ))

    @hybrid_property
    def completed(self):
        """Are all Strategies linked to this Target completed (or expired)."""
        # I guess we're completed by default if we have no Strategies, but it's not ideal...
        return all(strategy.finished for strategy in self.strategies)

    @completed.expression
    def completed(self):
        # Have to invert and return if any are incomplete
        return ~exists().where(and_(Strategy.target_id == self.db_id,
                                    Strategy.finished.is_(False)
                                    ))

    @hybrid_method
    def completed_at_time(self, time):
        """Were all Strategies linked to this Target completed (or expired) at the given time."""
        return all(strategy.finished_at_time(time) for strategy in self.strategies)

    @completed_at_time.expression
    def completed_at_time(self, time):
        return ~exists().where(and_(Strategy.target_id == self.db_id,
                                    Strategy.finished_at_time(time).is_(False)
                                    ))

    @hybrid_property
    def status(self):
        """Return a string giving the current status."""
        time = Time.now()
        if self.deleted_time is not None:
            # The Target has been deleted before it completed: deleted
            return 'deleted'
        elif self.completed:
            # The Target has enough completed Pointings: completed
            return 'completed'
        elif self.start_time is not None and time < Time(self.start_time):
            # The Target hasn't yet reached its start time: upcoming
            return 'upcoming'
        elif self.stop_time is not None and time >= Time(self.stop_time):
            # The Target has passed its stop time: expired
            return 'expired'
        elif not self.scheduled:
            # The Target doesn't have a pending Pointing: unscheduled
            return 'unscheduled'
        else:
            # The Target has a pending Pointing: scheduled
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
                     (self.completed,
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
        if isinstance(time, (str, datetime.datetime)):
            time = Time(time)

        if time < Time(self.creation_time):
            # This shouldn't usually happen, unless we are doing simulations and set it manually,
            # then it's worth knowing if this entry wouldn't have been created yet..
            return None

        if (self.deleted_time is not None and time >= Time(self.deleted_time)):
            return 'deleted'
        elif self.completed_at_time(time):
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
                  (self.completed_at_time(time),
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
        """Mark this Target (and any pending Pointings) as deleted."""
        if time is None:
            time = Time.now()
        if isinstance(time, (str, datetime.datetime)):
            time = Time(time)

        if self.status_at_time(time) == 'deleted':
            raise ValueError(f'Target is already deleted (at={self.deleted_time})')
        if self.status_at_time(time) == 'expired':
            raise ValueError(f'Target is already expired (at {self.stop_time})')
        if self.status_at_time(time) == 'completed':
            raise ValueError('Target is already completed')

        # Pointings can use finished_time, but we need an explicit deleted_time
        self.deleted_time = time

        for pointing in self.pointings:
            if pointing.status in ['upcoming', 'pending']:
                pointing.mark_deleted(time=time)

    @hybrid_property
    def current_rank(self):
        """Calculate current rank for new Pointings.

        If rank_decay is True, this is the start_rank + 10 * the number of completed Pointings.
        If False then it's just the start_rank.

        """
        if self.rank_decay:
            return self.start_rank + 10 * self.num_completed
        else:
            return self.start_rank

    @current_rank.expression
    def current_rank(self):
        return case([(self.rank_decay.is_(True),
                      self.start_rank + 10 * self.num_completed),
                     ],
                    else_=self.start_rank)

    @hybrid_method
    def current_rank_at_time(self, time):
        """Calculate current rank for new Pointings at the given time.

        If rank_decay is True, this is the start_rank + 10 * the number of completed Pointings.
        If False then it's just the start_rank.

        """
        return self.start_rank + 10 * self.num_completed_at_time(time)

    @current_rank_at_time.expression
    def current_rank_at_time(self, time):
        if isinstance(time, str):
            time = Time(time)
        if isinstance(time, Time):
            time = time.datetime

        return case([(self.rank_decay.is_(True),
                      self.start_rank + 10 * self.num_completed_at_time(time)),
                     ],
                    else_=self.start_rank)

    def get_current_strategy(self, time=None):
        """Get the currently valid Strategy, or at the given time."""
        if time is None:
            time = Time.now()

        if len(self.strategies) == 0:
            return None

        if len(self.strategies) == 1:
            strategy = self.strategies[0]
            # Check it's not finished (complete or past its stop_time)
            if not strategy.finished_at_time(time):
                return strategy
        else:
            # Loop through them until we find the first that isn't finished
            for strategy in self.strategies:
                if not strategy.finished_at_time(time):
                    return strategy

        # If we haven't found one then they're all either complete or too late (or both)
        return None

    @property
    def strategy(self):
        """The currently valid Strategy."""
        return self.get_current_strategy()

    @strategy.setter
    def strategy(self, strategy):
        """Set the given Strategy as the only strategy for this Target.

        This is a sneaky way to shortcut `Target.strategies = [strategy]
        """
        self.strategies = [strategy]

    def get_next_pointing(self, time=None):
        """Retrieve the next Pointing which needs to be scheduled.

        The start and stop times of the Pointing are determined from
        the status of the previous Pointing and the TimeBlocks attached
        to the current Strategy.

        Assumes this object is still associated with an active session.

        Any time given is used to find the valid strategy and set as the creation time,
        it should only be used for simulations (just like status_at_time() etc).

        Returns
        -------
        pointing : `Pointing` or None
            the next Pointing that should be sent to the database. Is None
            if there is no suitable Pointing remaining, or if there is
            a Pointing from this Target already scheduled.

        """
        # Already scheduled or finished, return None
        if self.status not in ['upcoming', 'unscheduled']:
            return None

        # As a safety check, see if the Target has any other Pointings that
        # are already pending or running.
        # If so the Target status should be 'scheduled' and therefore be
        # caught above, but strange things can happen.
        # This should prevent the cases where Targets have multiple
        # Pointings in the queue, hopefully!
        if len(self.pointings) > 0:
            statuses = [p.status for p in self.pointings]
            if 'pending' in statuses or 'running' in statuses:
                return None

        # Get the current Strategy (or at the given time)
        strategy = self.get_current_strategy(time)

        if strategy is None:
            # All observations were completed (or we have no linked Strategies)
            return None

        # Get the current and next TimeBlocks from the Strategy
        current_block = strategy.get_current_block()
        next_block = strategy.get_next_block()

        # Find the Pointing start_time
        if current_block is None or len(current_block.pointings) == 0:
            # No current block or Pointings, should only happen the first time
            # So use the start_time from the Target
            start_time = self.start_time
        else:
            latest_pointing = current_block.pointings[-1]
            # Decide to go onto the next block:
            #  - if the last block's pointing was completed, OR
            #  - if the last block's pointing's valid time has expired
            if latest_pointing.status in ['completed', 'expired']:
                if current_block.valid_time > 0:
                    # We want to start after the wait_time from the last Pointing's block.
                    # We can't just use the stop_time, because it might have been set by the
                    # Target's stop_time instead (see below).
                    # So we take the previous Pointing's start time + the valid_time, then
                    # add the wait_time.
                    start_time = (Time(latest_pointing.start_time) +
                                  current_block.valid_time * u.min +
                                  current_block.wait_time * u.min)
                else:
                    # The valid_time means we have non-expiring Pointings.
                    # So we just add the wait_time to the stopped_time of the previous Pointing.
                    start_time = (Time(latest_pointing.stopped_time) +
                                  current_block.wait_time * u.minute)
            else:
                # The current block wasn't completed, and there's still time left
                # (i.e. the Pointing was interrupted)
                # Need to re-insert the current block with a new pointing
                start_time = latest_pointing.start_time

        # Find the Pointing stop_time
        if next_block.valid_time < 0 and not self.stop_time:
            # The valid_time is infinite and there's no Target stop_time,
            # so the Pointings should never expire and they don't have a stop_time.
            stop_time = None
        else:
            if next_block.valid_time > 0:
                # The stop_time is calculated by adding the valid_time to the start_time.
                stop_time = Time(start_time) + next_block.valid_time * u.minute
                if self.stop_time and stop_time > self.stop_time:
                    # The Pointing stop_time would overrun past when this Target is valid for,
                    # so we force it to stop at the earlier time.
                    stop_time = self.stop_time
            else:
                # The Target stop time takes priority over infinite valid time.
                stop_time = self.stop_time

            if start_time >= stop_time:
                # This can happen if the Target has a stop_time, the Pointing would be impossible.
                return None

        # Mark the new block as current
        current_block.current = False
        next_block.current = True

        # Now create the pointing
        new_pointing = Pointing(rank=self.current_rank,
                                start_time=start_time,
                                stop_time=stop_time,
                                target=self,
                                time_block=next_block,
                                strategy=strategy,
                                )

        # set the creation time (for simulations)
        if time is not None:
            new_pointing.creation_time = time

        return new_pointing


class Strategy(Base):
    """A class to represent observing Strategies.

    Strategies contain observing parameters to define when Pointings are created and constraints
    which are applied when the validity of a Pointing is calculated by the scheduler.
    Each Target should have one or more Strategy classes defined, and the currently valid one is
    passed on to its Pointings when they are generated.

    Like all SQLAlchemy model classes, this object links to the
    underlying database. You can create an instance, and set its attributes
    without a database session. Accessing some attributes may require
    an active database session, and some properties (like the db_id)
    will be None until the instance is added to the database.

    Parameters
    ----------
    num_todo : int
        number of (successful) observations required for this Strategy to be completed.
        less than zero means repeat infinitely (or until stop_time is reached)

    stop_time : string, `astropy.time.Time` or datetime.datetime, optional
        UTC time after which this Strategy is considered invalid.
        if a Target has multiple Strategies then the next one will be selected once this time is
        passed, even if this Strategy is incomplete (i.e. hasn't reached num_todo Pointings)
        if not given then the Strategy is always valid (until it's completed)
        default = None

    wait_time : float or list of float, optional
        minimum time to wait between completed Pointings, in minutes.
        if the given `num_todo` is greater than times given the list will be looped.
        default = 0 (no delay)
    valid_time : float or list of float, optional
        the amount of time Pointing should reamin valid in the queue, in minutes.
        less than zero means valid indefinitely
        if the given `num_todo` is greater than times given the list will be looped.
        default = -1 (indefinitely valid)
    min_time : float, optional
        minimum time desired when observing a Pointing, in seconds.
        if a Pointing is interrupted after observing for this time then it is marked as complete
        instead of interrupted.
        default = None
    too : bool, optional
        indicates if Pointings should be considered Targets of Opportunity (ToO) or not.
        Pointings which are ToOs can interupt lower-(or equal-) ranked Pointings in the scheduler
        queue.
        default = False

    min_alt : float, optional
        minimum altitude to observe at, degrees
        default = 30
    max_sunalt : float, optional
        altitude constraint on Sun, degrees
        default = -15
    max_moon : string, optional
        Moon constraint
        one of 'D', 'G', 'B'
        default = 'B'
    min_moonsep : float, optional
        distance constraint from the Moon, degrees
        default = 30
    tel_mask : int, optional
        if set, this is a binary mask which will determine which telescopes
        can carry out the observation.
        A value of 4 (binary (0100) means only observable by Telescope 3, a value of
        5 (binary 0101) will allow either Telescope 1 or 3 to observe the Pointing.
        default = None (no restriction)

    When created the instance can be linked to the following other tables as parameters,
    otherwise they are populated when it is added to the database:

    Primary relationships
    ---------------------
    target : `Target`, optional
        the Target this Strategy is linked to, if any
        can also be added with the target_id parameter
    pointings : list of `Pointing`, optional
        the Pointings associated with this Strategy, if any
    time_blocks : list of `TimeBlock`, optional
        the TimeBlocks generated by this Strategy, if any

    Attributes
    ----------
    db_id : int
        primary database key
        only populated when the instance is added to the database
    num_completed : int
        number of successfully completed Pointings linked to the Strategy
        see also the `num_completed_at_time()` method
    num_remaining : int
        number of Pointings still to do (same as num_todo - num_completed)
        see also the `num_remaining_at_time()` method
    infinite : bool
        if Pointings will be regenerated infnitely (if num_todo is < 0)
    valid_telescopes : list of int
        a list of valid Telescope IDs, if any, based on the tel_mask

    Methods
    -------
    get_current_block() : TimeBlock, or None
        get the current TimeBlock
    get_last_block() : TimeBlock, or None
        get the previous TimeBlock
    get_next_block() : TimeBlock, or None
        get the upcoming TimeBlock

    """

    # Set corresponding SQL table name
    __tablename__ = 'strategies'

    # Primary key
    db_id = Column('id', Integer, primary_key=True)

    # Columns
    # # Basic properties
    num_todo = Column(Integer, nullable=False)
    stop_time = Column(DateTime, nullable=True, default=None)
    # # Pointing properties
    # wait_time - not stored in database (use TimeBlocks)
    # valid_time - not stored in database (use TimeBlocks)
    min_time = Column(Float, nullable=True, default=None)
    too = Column(Boolean, nullable=False, default=False)
    # # Observing constraints
    min_alt = Column(Float, nullable=False, default=params.DEFAULT_MIN_ALT)
    max_sunalt = Column(Float, nullable=False, default=params.DEFAULT_MAX_SUNALT)
    max_moon = Column(String(1), nullable=False, default=params.DEFAULT_MAX_MOON)
    min_moonsep = Column(Float, nullable=False, default=params.DEFAULT_MIN_MOONSEP)
    tel_mask = Column(Integer, nullable=True, default=None)

    # Foreign keys
    target_id = Column(Integer, ForeignKey('targets.id'), nullable=True)

    # Foreign relationships
    target = relationship('Target', back_populates='strategies')
    pointings = relationship('Pointing', back_populates='strategy')
    time_blocks = relationship('TimeBlock', lazy='joined', back_populates='strategy')

    def __init__(self, **kwargs):
        # Get extra arguments (need to remove from kwargs before super init)
        valid_times = kwargs.pop('valid_time') if 'valid_time' in kwargs else -1
        if not isinstance(valid_times, list):
            valid_times = [valid_times]
        wait_times = kwargs.pop('wait_time') if 'wait_time' in kwargs else 0
        if not isinstance(wait_times, list):
            wait_times = [wait_times]

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

    def __repr__(self):
        strings = ['db_id={}'.format(self.db_id),
                   'num_todo={}'.format(self.num_todo),
                   'num_completed={}'.format(self.num_completed),
                   'num_remaining={}'.format(self.num_remaining),
                   'stop_time={}'.format(self.stop_time),
                   'min_time={}'.format(self.min_time),
                   'too={}'.format(self.too),
                   'min_alt={}'.format(self.min_alt),
                   'max_sunalt={}'.format(self.max_sunalt),
                   'max_moon={}'.format(self.max_moon),
                   'min_moonsep={}'.format(self.min_moonsep),
                   'tel_mask={}'.format(self.tel_mask),
                   'target_id={}'.format(self.target_id),
                   ]
        return 'Strategy({})'.format(', '.join(strings))

    @validates('stop_time')
    def validate_times(self, key, field):
        """Use validators to allow various types of input for times."""
        if field is None:
            return None

        if isinstance(field, datetime.datetime):
            value = field.strftime('%Y-%m-%d %H:%M:%S')
        elif isinstance(field, Time):
            field.precision = 0  # no D.P on seconds
            value = field.iso
        else:
            # just hope the string works!
            value = str(field)

        return value

    @hybrid_property
    def infinite(self):
        """Return True if the number of Pointings to generate is infinite."""
        return self.num_todo < 0

    @hybrid_property
    def num_completed(self):
        """Return the number of Pointings successfully completed."""
        return sum(p.status == 'completed' for p in self.pointings)

    @num_completed.expression
    def num_completed(self):
        return select([func.count(Pointing.db_id)]).\
            where(and_(Pointing.strategy_id == self.db_id,
                       Pointing.status == 'completed',
                       )).scalar_subquery()

    @hybrid_method
    def num_completed_at_time(self, time):
        """Return the number of Pointings successfully completed at the given time."""
        return sum(p.status_at_time(time) == 'completed' for p in self.pointings)

    @num_completed_at_time.expression
    def num_completed_at_time(self, time):
        return select([func.count(Pointing.db_id)]).\
            where(and_(Pointing.strategy_id == self.db_id,
                       Pointing.status_at_time(time) == 'completed',
                       ))

    @hybrid_property
    def num_remaining(self):
        """Return the number of observations remaining (num_todo - num_completed).

        Returns -1 if num_todo is infinite.
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
    def finished(self):
        """Return True if this Strategy is either completed or expired."""
        time = Time.now()
        if self.num_remaining == 0:
            # The Strategy has enough completed Pointings: completed
            return True
        elif self.stop_time is not None and time >= Time(self.stop_time):
            # The Strategy has passed its stop time: expired
            return True
        else:
            # The Strategy still needs Pointings (and is still valid): incomplete
            return False

    @finished.expression
    def finished(self):
        time = datetime.datetime.utcnow()
        return case([(self.num_remaining == 0,
                      True),
                     (and_(self.stop_time.isnot(None), time >= self.stop_time),
                      True),
                     ],
                    else_=False)

    @hybrid_method
    def finished_at_time(self, time):
        """Return True if this Strategy is was completed or expired at the given time."""
        if self.num_remaining_at_time(time) == 0:
            return True
        elif self.stop_time is not None and time >= Time(self.stop_time):
            return True
        else:
            return False

    @finished_at_time.expression
    def finished_at_time(self, time):
        if isinstance(time, str):
            time = Time(time)
        if isinstance(time, Time):
            time = time.datetime

        c = case([(self.num_remaining_at_time(time) == 0,
                   True),
                  (and_(self.stop_time.isnot(None), time >= self.stop_time),
                   True),
                  ],
                 else_=False)
        return c

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


class TimeBlock(Base):
    """A class to represent a block of time.

    Specifically, these represent a time period made up of a valid time and a wait time.
    The valid time is the time that the pointing will be valid for, and the wait time will be the
    time after the pointing is observed/invalid to wait until the next valid period.

    NB: you shouldn't ever have to manually create TimeBlocks.
    They're used internally by Strategies to keep track of Pointings and all
    the ones a Strategy requires are created when the Strategy is.

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
    strategy : `Strategy`
        the Strategy associated with this Time Block
        required before addition to the database
        can also be added with the target_id parameter

    target : `Target`, optional
        the Target associated with this Time Block, if any
        can also be added with the target_id parameter
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
    strategy_id = Column(Integer, ForeignKey('strategies.id'), nullable=False)

    # Foreign relationships
    strategy = relationship('Strategy', back_populates='time_blocks')
    pointings = relationship('Pointing', back_populates='time_block')

    def __repr__(self):
        strings = ['db_id={}'.format(self.db_id),
                   'block_num={}'.format(self.block_num),
                   'valid_time={}'.format(self.valid_time),
                   'wait_time={}'.format(self.wait_time),
                   'current={}'.format(self.current),
                   'strategy_id={}'.format(self.strategy_id),
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
        """Create a Site from an Astropy `EarthLocation` class and optional name string."""
        if not isinstance(location, EarthLocation):
            raise ValueError('"location" must be an `astropy.coordinates.EarthLocation`')
        if name is None:
            if hasattr(location.info, 'name'):
                name = location.info.name
            else:
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
    targets : list of `Target`
        the Targets which cover any of this Grid's GridTiles, if any

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
    telescopes = relationship('Telescope', back_populates='grid')

    # Secondary relationships
    targets = relationship('Target',
                           secondary='grid_tiles',
                           primaryjoin='Grid.db_id == GridTile.grid_id',
                           secondaryjoin='Target.grid_tile_id == GridTile.db_id',
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

    @classmethod
    def from_skygrid(cls, skygrid, name=None):
        """Create a Grid from a GOTO-tile `SkyGrid` class and optional name string."""
        if not isinstance(skygrid, SkyGrid):
            raise ValueError('"location" must be a `gototile.grid.SkyGrid`')
        if name is None:
            if hasattr(skygrid, 'name'):
                name = skygrid.name
            else:
                raise ValueError('Missing name for grid')

        return cls(name=name,
                   ra_fov=skygrid.fov['ra'].value,
                   dec_fov=skygrid.fov['dec'].value,
                   ra_overlap=skygrid.overlap['ra'],
                   dec_overlap=skygrid.overlap['dec'],
                   algorithm=skygrid.algorithm,
                   )

    @property
    def skygrid(self):
        """Return a GOTO-tile SkyGrid for this Grid."""
        fov = {'ra': self.ra_fov * u.deg, 'dec': self.dec_fov * u.deg}
        overlap = {'ra': self.ra_overlap, 'dec': self.dec_overlap}
        return SkyGrid(fov, overlap, kind=self.algorithm)


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

    targets : list of `Target`, optional
        the Targets that cover this GridTile, if any

    Attributes
    ----------
    db_id : int
        primary database key
        only populated when the instance is added to the database
    last_observed : `astropy.time.Time`
        the most recent time of completion of any Target covering this GridTile, if any

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
    targets = relationship('Target', back_populates='grid_tile')

    # Secondary relationships
    pointings = relationship('Pointing',
                             secondary='targets',
                             primaryjoin='GridTile.db_id == Target.grid_tile_id',
                             secondaryjoin='Pointing.target_id == Target.db_id',
                             back_populates='grid_tile',
                             viewonly=True,
                             uselist=True,
                             )

    # Column properties
    last_observed = column_property(select([Target.last_observed])
                                    .where(Target.grid_tile_id == db_id,
                                           )
                                    .order_by(Target.last_observed.desc())
                                    .limit(1)
                                    .correlate_except(Pointing)
                                    .scalar_subquery()
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

    A Survey is a way to group Targets, usually from a paticular Event.

    Like all SQLAlchemy model classes, this object links to the
    underlying database. You can create an instance, and set its attributes
    without a database session. Accessing some attributes may require
    an active database session, and some properties (like the db_id)
    will be None until the instance is added to the database.

    Parameters
    ----------
    name : string
        a human-readable identifier for the survey

    skymap : string, optional
        the location of the source skymap for this survey

    When created the instance can be linked to the following other tables as parameters,
    otherwise they are populated when it is added to the database:

    Primary relationships
    ---------------------
    targets : list of `Target`, optional
        the Targets that are part of this Survey, if any
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
    pointings : list of `Pointing`
        the Pointings which are part of this Survey, if any

    """

    # Set corresponding SQL table name
    __tablename__ = 'surveys'

    # Primary key
    db_id = Column('id', Integer, primary_key=True)

    # Columns
    name = Column(String(255), nullable=False, index=True)
    skymap = Column(String(255), nullable=True, default=None)

    # Foreign keys
    event_id = Column(Integer, ForeignKey('events.id'), nullable=True)

    # Foreign relationships
    targets = relationship('Target', back_populates='survey')
    event = relationship('Event', back_populates='surveys')

    # Secondary relationships
    pointings = relationship('Pointing',
                             secondary='targets',
                             primaryjoin='Survey.db_id == Target.event_id',
                             secondaryjoin='Pointing.target_id == Target.db_id',
                             back_populates='survey',
                             viewonly=True,
                             uselist=True,
                             )

    def __repr__(self):
        strings = ['db_id={}'.format(self.db_id),
                   'name={}'.format(self.name),
                   'skymap={}'.format(self.skymap),
                   'event_id={}'.format(self.event_id),
                   ]
        return 'Survey({})'.format(', '.join(strings))


class Event(Base):
    """A class to represent a transient Event.

    Events can be linked to specific Surveys made of Targets, usually weighted
    based on a skymap. A specific astrophysical event (e.g. GW170817) might
    produce multiple Surveys as the skymap is updated, which will all be linked
    to the one Event.

    Like all SQLAlchemy model classes, this object links to the
    underlying database. You can create an instance, and set its attributes
    without a database session. Accessing some attributes may require
    an active database session, and some properties (like the db_id)
    will be None until the instance is added to the database.

    Parameters
    ----------
    name : string
        a unique, human-readable identifier for the event
    source : string
        the event's origin, e.g. LVC, Fermi, GAIA

    type : string, optional
        the type of event, e.g. GW, GRB
    time : string, `astropy.time.Time` or datetime.datetime, optional
        time the event occurred

    When created the instance can be linked to the following other tables as parameters,
    otherwise they are populated when it is added to the database:

    Primary relationships
    ---------------------
    surveys : list of `Survey`, optional
        the Surveys created for this Event, if any
    targets : list of `Target`, optional
        the Targets created for this Event, if any

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
        the Pointings created for this Event, if any

    """

    # Set corresponding SQL table name
    __tablename__ = 'events'

    # Primary key
    db_id = Column('id', Integer, primary_key=True)

    # Columns
    name = Column(String(255), nullable=False, unique=True, index=True)
    source = Column(String(255), nullable=False)
    type = Column(String(255), nullable=True, index=True, default=None)  # noqa: A003
    time = Column(DateTime, nullable=True, default=None)

    # Foreign relationships
    surveys = relationship('Survey', back_populates='event')
    targets = relationship('Target', back_populates='event')

    # Secondary relationships
    pointings = relationship('Pointing',
                             secondary='targets',
                             primaryjoin='Event.db_id == Target.event_id',
                             secondaryjoin='Pointing.target_id == Target.db_id',
                             back_populates='event',
                             viewonly=True,
                             uselist=True,
                             )

    def __repr__(self):
        strings = ['db_id={}'.format(self.db_id),
                   'name={}'.format(self.name),
                   'source={}'.format(self.source),
                   'type={}'.format(self.type),
                   'time={}'.format(self.time),
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


TRIGGERS = []
