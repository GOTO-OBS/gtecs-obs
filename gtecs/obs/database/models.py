"""Python classes mapping on to database tables."""

import datetime
import hashlib

from astropy import units as u
from astropy.coordinates import EarthLocation
from astropy.time import Time

from gototile.grid import SkyGrid

import numpy as np

from sqlalchemy import Boolean, Column, DateTime, Float, ForeignKey, Integer, String, Text
from sqlalchemy import event, exists, func, select, DDL
from sqlalchemy.ext.associationproxy import association_proxy
from sqlalchemy.ext.hybrid import hybrid_method, hybrid_property
from sqlalchemy.orm import column_property, declarative_base, relationship, validates
from sqlalchemy.sql import and_, case, or_

from .. import params


__all__ = ['User', 'ExposureSet', 'Pointing', 'Strategy', 'Target', 'TimeBlock',
           'Site', 'Telescope', 'Grid', 'GridTile', 'Survey',
           ]

SCHEMA = 'obs'

Base = declarative_base()


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

    # Set corresponding SQL table name and schema
    __tablename__ = 'users'
    __table_args__ = {'schema': SCHEMA}

    # Primary key
    db_id = Column('id', Integer, primary_key=True)

    # Columns
    username = Column(String(255), nullable=False, unique=True)
    password_hash = Column('password', String(255), nullable=False)
    full_name = Column(Text, nullable=False)

    # Update timestamp
    ts = Column(DateTime, nullable=False, server_default=func.now())

    # Foreign relationships
    targets = relationship(
        'Target',
        order_by='Target.db_id',
        back_populates='user',
    )

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
    dithering : int, optional
        dithering pattern to use (0=None, 1=default, 2+=custom (TODO: not implemented))
        default = 1 (default dithering)
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

    # Set corresponding SQL table name and schema
    __tablename__ = 'exposure_sets'
    __table_args__ = {'schema': SCHEMA}

    # Primary key
    db_id = Column('id', Integer, primary_key=True)

    # Columns
    num_exp = Column(Integer, nullable=False)
    exptime = Column(Float, nullable=False)
    filt = Column('filter',   # filter is a built in function in Python
                  String(1), nullable=False)
    binning = Column(Integer, nullable=False, default=1)
    dithering = Column(Integer, nullable=False, default=1)
    ut_mask = Column(Integer, default=None)

    # Foreign keys
    target_id = Column(Integer, ForeignKey(f'{SCHEMA}.targets.id'), nullable=True, index=True)

    # Update timestamp
    ts = Column(DateTime, nullable=False, server_default=func.now())

    # Foreign relationships
    target = relationship(
        'Target',
        order_by='Target.db_id',
        back_populates='exposure_sets',
    )

    # Secondary relationships
    pointings = relationship(
        'Pointing',
        order_by='Pointing.db_id',
        secondary=f'{SCHEMA}.targets',
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
                   'dithering={}'.format(self.dithering),
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
    rank : int or None
        rank to use when scheduling (lower values are prioritised first by the scheduler)
        rank=None means infinite rank, will always be at the bottom of the queue (queue-fillers)

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
    completed : bool or None, optional
        if the pointing was successfully completed (True) or interrupted (False), or just
        hasn't started yet (None)
        default = None
    finished_time : string, `astropy.time.Time` or datetime.datetime, or None, optional
        the time the pointing was marked as finished (i.e. completed or interrupted)
        this should usually be set automatically when the pointing status is set
        default = None
    valid : bool or None, optional
        if the results of the completed Pointing was confirmed as good (True) or failed (False),
        or just hasn't been reviewed yet (None)
        default = None
    validated_time : string, `astropy.time.Time` or datetime.datetime, or None, optional
        the time the resulting images from this Pointing were reviewed
        this should usually be set automatically when the validity is set
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

    The following proxy attributes are attributes of related classes, and are included
    as a useful shortcut (e.g. `pointing.ra` is the same as `pointing.target.ra`)

    Proxy attributes
    ----------------
    name : str
        object name
        (inherited from linked Target)
    ra : float, optional
        J2000 right ascension in decimal degrees
        (inherited from linked Target)
    dec : float, optional
        J2000 declination in decimal degrees
        (inherited from linked Target)
    weight : float, optional
        weighting relative to other Targets in the same survey
        (inherited from linked Target)
    is_template : bool, optional
        if True, this Pointing will count as a template for a linked GridTile
        (inherited from linked Target)

    min_time : float, optional
        minimum time desired when observing this Pointing, in seconds.
        (inherited from linked Strategy)
    too : bool, optional
        indicates if this Pointing should be considered a Targets of Opportunity (ToO) or not.
        (inherited from linked Strategy)
    requires_template : bool, optional
        if True, this Pointing won't be valid in the scheduler queue unless the linked GridTile
        has at least 1 completed Pointing with `is_template=True`
        (inherited from linked Strategy)
    min_alt : float, optional
        minimum altitude to observe at, degrees
        (inherited from linked Strategy)
    max_sunalt : float, optional
        altitude constraint on Sun, degrees
        (inherited from linked Strategy)
    max_moon : string, optional
        Moon constraint, one of 'D', 'G', 'B'
        (inherited from linked Strategy)
    min_moonsep : float, optional
        distance constraint from the Moon, degrees
        (inherited from linked Strategy)
    tel_mask : int, optional
        a binary mask which will determine which telescopes can carry out the observation
        (inherited from linked Strategy)

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
        mark the pointing as finished, either completed (True) or interrupted (False)
        if no time is given, default to Time.now()
    mark_validated(self, good=True, time=None)
        mark the pointing as validated, either good (True) or failed (False)
        if no time is given, default to Time.now()

    """

    # Set corresponding SQL table name and schema
    __tablename__ = 'pointings'
    __table_args__ = {'schema': SCHEMA}

    # Primary key
    db_id = Column('id', Integer, primary_key=True)

    # Columns
    # # Basic properties
    rank = Column(Integer, nullable=True)
    # # Constraints
    start_time = Column(DateTime, nullable=False, index=True, server_default=func.now())
    stop_time = Column(DateTime, nullable=True, index=True, default=None)
    # # Status
    creation_time = Column(DateTime, nullable=False, server_default=func.now())
    running_time = Column(DateTime, nullable=True, default=None)
    completed = Column(Boolean, nullable=False, default=False)
    finished_time = Column(DateTime, nullable=True, default=None)
    valid = Column(Boolean, nullable=True, default=None)
    validated_time = Column(DateTime, nullable=True, default=None)

    # Foreign keys
    target_id = Column(Integer, ForeignKey(f'{SCHEMA}.targets.id'), nullable=False, index=True)
    time_block_id = Column(Integer, ForeignKey(f'{SCHEMA}.time_blocks.id'), nullable=False, index=True)
    strategy_id = Column(Integer, ForeignKey(f'{SCHEMA}.strategies.id'), nullable=False, index=True)
    telescope_id = Column(Integer, ForeignKey(f'{SCHEMA}.telescopes.id'), nullable=True, index=True)

    # Update timestamp
    ts = Column(DateTime, nullable=False, server_default=func.now())

    # Foreign relationships
    target = relationship(
        'Target',
        order_by='Target.db_id',
        lazy='joined',  # SAVE TIME IN SCHEDULER
        back_populates='pointings',
    )
    time_block = relationship(
        'TimeBlock',
        order_by='TimeBlock.db_id',
        back_populates='pointings',
    )
    strategy = relationship(
        'Strategy',
        order_by='Strategy.db_id',
        lazy='joined',  # SAVE TIME IN SCHEDULER
        back_populates='pointings',
    )
    telescope = relationship(
        'Telescope',
        order_by='Telescope.db_id',
        back_populates='pointings',
    )

    # Secondary relationships
    exposure_sets = relationship(
        'ExposureSet',
        order_by='ExposureSet.db_id',
        lazy='joined',  # SAVE TIME IN SCHEDULER
        secondary=f'{SCHEMA}.targets',
        primaryjoin='Pointing.target_id == Target.db_id',
        secondaryjoin='ExposureSet.target_id == Target.db_id',
        back_populates='pointings',
        viewonly=True,
        uselist=True,
    )
    grid_tile = relationship(
        'GridTile',
        order_by='GridTile.db_id',
        secondary=f'{SCHEMA}.targets',
        primaryjoin='Pointing.target_id == Target.db_id',
        secondaryjoin='GridTile.db_id == Target.grid_tile_id',
        back_populates='pointings',
        viewonly=True,
        uselist=False,
    )
    grid_tile_id = association_proxy('grid_tile', 'db_id')
    survey = relationship(
        'Survey',
        order_by='Survey.db_id',
        secondary=f'{SCHEMA}.targets',
        primaryjoin='Pointing.target_id == Target.db_id',
        secondaryjoin='Survey.db_id == Target.survey_id',
        back_populates='pointings',
        viewonly=True,
        uselist=False,
    )
    survey_id = association_proxy('survey', 'db_id')

    # Proxy attributes
    name = association_proxy('target', 'name')
    ra = association_proxy('target', 'ra')
    dec = association_proxy('target', 'dec')
    weight = association_proxy('target', 'weight')
    is_template = association_proxy('target', 'is_template')
    min_time = association_proxy('strategy', 'min_time')
    too = association_proxy('strategy', 'too')
    requires_template = association_proxy('strategy', 'requires_template')
    min_alt = association_proxy('strategy', 'min_alt')
    max_sunalt = association_proxy('strategy', 'max_sunalt')
    max_moon = association_proxy('strategy', 'max_moon')
    min_moonsep = association_proxy('strategy', 'min_moonsep')
    tel_mask = association_proxy('strategy', 'tel_mask')

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
                   'completed={}'.format(self.completed),
                   'finished_time={}'.format(self.finished_time),
                   'valid={}'.format(self.valid),
                   'validated_time={}'.format(self.validated_time),
                   'target_id={}'.format(self.target_id),
                   'time_block_id={}'.format(self.time_block_id),
                   'telescope_id={}'.format(self.telescope_id),
                   ]
        return 'Pointing({})'.format(', '.join(strings))

    @validates('start_time', 'stop_time', 'creation_time', 'running_time', 'finished_time',
               'validated_time')
    def validate_times(self, key, field):
        """Use validators to allow various types of input for times.

        Also enforce logical orders:
            - stop_time > start_time
            - finished_time > running_time
            - validated_time > finished_time
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

        # force start_time < stop_time
        if (key == 'start_time' and self.stop_time is not None and
                Time(value) >= Time(self.stop_time)):
            raise ValueError(f'start_time must be before stop_time ({self.stop_time})')
        if (key == 'stop_time' and self.start_time is not None and
                Time(self.start_time) >= Time(value)):
            raise ValueError(f'stop_time must be after start_time ({self.start_time})')

        # force running_time < finished_time
        if (key == 'running_time' and self.finished_time is not None and
                Time(value) >= Time(self.finished_time)):
            raise ValueError(f'running_time must be before finished_time ({self.finished_time})')
        if (key == 'finished_time' and self.running_time is not None and
                Time(self.running_time) >= Time(value)):
            raise ValueError(f'finished_time must be after running_time ({self.running_time})')

        # force finished_time < validated_time
        if (key == 'finished_time' and self.validated_time is not None and
                Time(value) >= Time(self.validated_time)):
            msg = f'finished_time must be before validated_time ({self.validated_time})'
            raise ValueError(msg)
        if (key == 'validated_time' and self.finished_time is not None and
                Time(self.finished_time) >= Time(value)):
            raise ValueError(f'validated_time must be after finished_time ({self.finished_time})')

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
        elif self.validated_time is not None and not self.valid:
            # The Pointing has been completed and judged as good enough: failed
            return 'failed'
        # # We don't need a successful status, there's no real difference with completed
        # # since we consider everything as successful by default
        # elif self.validated_time is not None and self.valid:
        #     # The Pointing has been completed and judged as good enough: confirmed
        #     return 'confirmed'
        elif self.finished_time is not None and self.completed:
            # The Pointing has finished and is flagged as successful: completed
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
        elif value == 'confirmed':
            self.mark_validated(good=True)
        elif value == 'failed':
            self.mark_validated(good=False)
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
        return case(
            (and_(self.running_time.is_(None), self.finished_time.isnot(None)),
             'deleted'),
            (and_(self.running_time.isnot(None), self.finished_time.is_(None)),
             'running'),
            (and_(self.validated_time.isnot(None), self.valid.is_(False)),
             'failed'),
            (and_(self.finished_time.isnot(None), self.completed.is_(True)),
             'completed'),
            (and_(self.finished_time.isnot(None), self.completed.is_(False)),
             'interrupted'),
            (time < self.start_time,
             'upcoming'),
            (and_(self.stop_time.isnot(None), time > self.stop_time),
             'expired'),
            else_='pending')

    @hybrid_method
    def status_at_time(self, time):
        """Return a string giving the status at the given time (for testing and simulations)."""
        if time is None:
            return self.status
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
        elif (self.validated_time is not None and time >= Time(self.validated_time) and
              not self.valid):
            return 'failed'
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
        if time is None:
            return self.status
        if isinstance(time, str):
            time = Time(time)
        if isinstance(time, Time):
            time = time.datetime

        c = case(
            (time < self.creation_time,
             None),
            (and_(self.running_time.is_(None),
                  self.finished_time.isnot(None),
                  time >= self.finished_time),
             'deleted'),
            (and_(self.running_time.isnot(None),
                  time >= self.running_time,
                  or_(self.finished_time.is_(None), time < self.finished_time)),
             'running'),
            (and_(self.validated_time.isnot(None),
                  time >= self.validated_time,
                  self.valid.is_(False)),
             'failed'),
            (and_(self.finished_time.isnot(None),
                  time >= self.finished_time,
                  self.completed.is_(True)),
             'completed'),
            (and_(self.finished_time.isnot(None),
                  time >= self.finished_time,
                  self.completed.is_(False)),
             'interrupted'),
            (time < self.start_time,
             'upcoming'),
            (and_(self.stop_time.isnot(None), time >= self.stop_time),
             'expired'),
            else_='pending')
        return c

    @property
    def valid_telescopes(self):
        """Get a list of valid Telescope IDs, if any, based on the tel_mask."""
        return self.strategy.valid_telescopes

    @valid_telescopes.setter
    def valid_telescopes(self, telescope_ids):
        self.strategy.valid_telescopes = telescope_ids

    @hybrid_method
    def tel_is_valid(self, telescope_id):
        """Return True if the given Telescope ID is valid, based on the tel_mask."""
        return self.strategy.tel_is_valid(telescope_id)

    @tel_is_valid.expression
    def tel_is_valid(self, telescope_id):
        return self.strategy.has(Strategy.tel_is_valid(telescope_id).is_(True))

    def mark_deleted(self, time=None):
        """Mark this Pointing as deleted."""
        if time is None:
            time = Time.now()
        if isinstance(time, (str, datetime.datetime)):
            time = Time(time)

        current_status = self.status_at_time(time)
        if current_status == 'deleted':
            msg = f'Pointing {self.db_id} is already deleted (at {self.finished_time})'
            raise ValueError(msg)
        elif current_status == 'expired':
            msg = f'Pointing {self.db_id} is already expired (at {self.stop_time})'
            raise ValueError(msg)
        elif current_status == 'running':
            msg = f'Pointing {self.db_id} is already running (at {self.running_time})'
            raise ValueError(msg)
        elif current_status in ['completed', 'interrupted', 'failed']:
            msg = f'Pointing {self.db_id} is already {current_status} (at {self.finished_time})'
            raise ValueError(msg)

        self.finished_time = time

    def mark_running(self, telescope=None, time=None):
        """Mark this Pointing as running on the given telescope."""
        if time is None:
            time = Time.now()
        if isinstance(time, (str, datetime.datetime)):
            time = Time(time)
        if telescope is None:
            raise ValueError('Pointings must be linked to a Telescope when marked running')

        current_status = self.status_at_time(time)
        if current_status == 'deleted':
            msg = f'Pointing {self.db_id} is already deleted (at {self.finished_time})'
            raise ValueError(msg)
        elif current_status == 'upcoming':
            msg = f'Pointing {self.db_id} has not yet reached its start time (at {self.start_time})'
            raise ValueError(msg)
        elif current_status == 'expired':
            msg = f'Pointing {self.db_id} is already expired (at {self.stop_time})'
            raise ValueError(msg)
        elif current_status == 'running':
            msg = f'Pointing {self.db_id} is already running (at {self.running_time})'
            raise ValueError(msg)
        elif current_status in ['completed', 'interrupted', 'failed']:
            msg = f'Pointing {self.db_id} is already {current_status} (at {self.finished_time})'
            raise ValueError(msg)

        if self.telescope is not None:
            msg = f'Pointing {self.db_id} is already linked to Telescope {self.telescope.db_id}'
            raise ValueError(msg)
        elif self.telescope_id is not None:
            msg = f'Pointing {self.db_id} is already linked to Telescope {self.telescope_id}'
            raise ValueError(msg)

        if (self.valid_telescopes is not None and
                telescope.db_id not in self.valid_telescopes):
            msg = f'Telescope {telescope.db_id} not in list of valid telescopes'
            msg += f' ({self.valid_telescopes})'
            raise ValueError(msg)
        if telescope.current_pointing is not None:
            pointing_id = telescope.current_pointing.db_id
            msg = f'Telescope {telescope.db_id} is already observing Pointing {pointing_id}'
            raise ValueError(msg)

        self.running_time = time
        self.telescope = telescope

    def mark_finished(self, completed=True, time=None):
        """Mark this Pointing as stopped, either completed (True) or interrupted (False)."""
        if time is None:
            time = Time.now()
        if isinstance(time, (str, datetime.datetime)):
            time = Time(time)

        current_status = self.status_at_time(time)
        if current_status == 'deleted':
            msg = f'Pointing {self.db_id} is already deleted (at {self.finished_time})'
            raise ValueError(msg)
        elif current_status == 'expired':
            msg = f'Pointing {self.db_id} is already expired (at {self.stop_time})'
            raise ValueError(msg)
        elif current_status in ['completed', 'interrupted', 'failed']:
            msg = f'Pointing {self.db_id} is already {current_status} (at {self.finished_time})'
            raise ValueError(msg)
        elif current_status != 'running':
            msg = f'Pointing {self.db_id} is not running (use mark_deleted to stop before running)'
            raise ValueError(msg)

        self.finished_time = time
        self.completed = completed

    def mark_validated(self, good=True, time=None):
        """Mark this Pointing as validated, either good (True) or failed (False)."""
        if time is None:
            time = Time.now()
        if isinstance(time, (str, datetime.datetime)):
            time = Time(time)

        current_status = self.status_at_time(time)
        if current_status == 'deleted':
            msg = f'Pointing {self.db_id} is already deleted (at {self.finished_time})'
            raise ValueError(msg)
        elif current_status == 'expired':
            msg = f'Pointing {self.db_id} is already expired (at {self.stop_time})'
            raise ValueError(msg)
        elif current_status == 'running':
            msg = f'Pointing {self.db_id} is still running (at {self.running_time})'
            raise ValueError(msg)
        elif current_status == 'interrupted':
            msg = f'Pointing {self.db_id} was interrupted (at {self.finished_time})'
            raise ValueError(msg)
        elif current_status != 'completed':
            msg = f'Pointing {self.db_id} is not completed (use mark_finished before validating)'
            raise ValueError(msg)
        elif self.validated_time is not None and self.validated_time > time:
            msg = f'Pointing {self.db_id} has already been validated (at {self.validated_time})'
            raise ValueError(msg)

        self.validated_time = time
        self.valid = good

    def get_obstime(self, readout_time=15):
        """Get the expected time needed to observe this Pointing."""
        return sum(es.num_exp * (es.exptime + readout_time) for es in self.exposure_sets)


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

    wait_time : float, `astropy.units.Quantity` or list of same, optional
        minimum time to wait between completed Pointings.
        if a Quantity it should have units of time, if a float assume in seconds.
        if the given `num_todo` is greater than times given the list will be looped.
        default = 0 (no delay)
    valid_time : float, `astropy.units.Quantity`, list of same, or None, optional
        the amount of time Pointing should remain valid in the queue.
        if a Quantity it should have units of time, if a float assume in seconds.
        less than zero or None means valid indefinitely
        if the given `num_todo` is greater than times given the list will be looped.
        default = None (indefinitely valid)
    rank_change : int, list of int, or None, optional
        the amount to change the rank of Pointings as they get completed.
        the starting rank is defined by the Target.
        if the given `num_todo` is greater than changes given the list will be looped.
        default = 10 (add 10 to the start rank each time a Pointing is completed)
    min_time : float or `astropy.units.Quantity`, optional
        minimum time desired when observing a Pointing.
        if a Quantity it should have units of time, if a float assume in seconds.
        if a Pointing is interrupted after observing for this time then it is marked as complete
        instead of interrupted.
        default = None
    too : bool, optional
        indicates if Pointings should be considered Targets of Opportunity (ToO) or not.
        Pointings which are ToOs can interrupt lower-(or equal-) ranked Pointings in the scheduler
        queue.
        default = False
    requires_template : bool, optional
        if True, Pointings won't be valid in the scheduler queue unless their GridTile has at
        least 1 completed Pointing with `is_template=True`
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
        if Pointings will be regenerated indefinitely (if num_todo is < 0)
    valid_telescopes : list of int
        a list of valid Telescope IDs, if any, based on the tel_mask

    Methods
    -------
    get_current_block() : TimeBlock
        get the current TimeBlock
    get_last_block() : TimeBlock
        get the previous TimeBlock
    get_next_block() : TimeBlock
        get the upcoming TimeBlock

    """

    # Set corresponding SQL table name and schema
    __tablename__ = 'strategies'
    __table_args__ = {'schema': SCHEMA}

    # Primary key
    db_id = Column('id', Integer, primary_key=True)

    # Columns
    # # Basic properties
    num_todo = Column(Integer, nullable=False)
    stop_time = Column(DateTime, nullable=True, default=None)
    # # TimeBlock properties
    # wait_time - not stored in database (use TimeBlocks)
    # valid_time - not stored in database (use TimeBlocks)
    # rank_change - not stored in database (use TimeBlocks)
    # # Scheduling
    min_time = Column(Float, nullable=True, default=None)
    too = Column(Boolean, nullable=False, default=False)
    requires_template = Column(Boolean, nullable=False, default=False)
    # # Observing constraints
    min_alt = Column(Float, nullable=False, default=params.DEFAULT_MIN_ALT)
    max_sunalt = Column(Float, nullable=False, default=params.DEFAULT_MAX_SUNALT)
    max_moon = Column(String(1), nullable=False, default=params.DEFAULT_MAX_MOON)
    min_moonsep = Column(Float, nullable=False, default=params.DEFAULT_MIN_MOONSEP)
    tel_mask = Column(Integer, nullable=True, default=None)

    # Foreign keys
    target_id = Column(Integer, ForeignKey(f'{SCHEMA}.targets.id'), nullable=True, index=True)

    # Update timestamp
    ts = Column(DateTime, nullable=False, server_default=func.now())

    # Foreign relationships
    target = relationship(
        'Target',
        order_by='Target.db_id',
        back_populates='strategies',
    )
    pointings = relationship(
        'Pointing',
        order_by='Pointing.db_id',
        back_populates='strategy',
    )
    time_blocks = relationship(
        'TimeBlock',
        order_by='TimeBlock.db_id',
        back_populates='strategy',
    )

    # Column properties
    num_completed = column_property(select(func.count(Pointing.db_id))
                                    .where(and_(Pointing.strategy_id == db_id,
                                                Pointing.status == 'completed',
                                                )
                                           )
                                    .correlate_except(Pointing)
                                    .scalar_subquery()
                                    )

    def __init__(self, **kwargs):
        # Get extra arguments (need to remove from kwargs before super init)
        wait_times = kwargs.pop('wait_time') if 'wait_time' in kwargs else 0
        valid_times = kwargs.pop('valid_time') if 'valid_time' in kwargs else None
        rank_changes = kwargs.pop('rank_change') if 'rank_change' in kwargs else None
        try:
            len(wait_times)
        except TypeError:
            wait_times = [wait_times]
        try:
            len(valid_times)
        except TypeError:
            valid_times = [valid_times]
        try:
            len(rank_changes)
        except TypeError:
            rank_changes = [rank_changes]

        # Init base class
        super().__init__(**kwargs)

        # Create TimeBlocks (if none were given)
        if len(self.time_blocks) == 0:
            for i in range(max(len(wait_times), len(valid_times), len(rank_changes))):
                # Loop through lists to get values
                wait_time = wait_times[i % len(wait_times)]
                valid_time = valid_times[i % len(valid_times)]
                rank_change = rank_changes[i % len(rank_changes)]
                # Allow negative values as None
                if valid_time is not None and valid_time < 0:
                    valid_time = None
                # Create the block and add to list
                block = TimeBlock(
                    block_num=i + 1,
                    wait_time=wait_time,
                    valid_time=valid_time,
                    rank_change=rank_change,
                    )
                self.time_blocks.append(block)

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
            # stop_time is nullable
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

    @validates('min_time')
    def validate_durations(self, key, field):
        """Use validators to allow various types of input for time periods."""
        if field is None:
            # min_time is nullable
            return None

        if isinstance(field, u.Quantity):
            # Should have units of time
            if field.unit not in u.s.find_equivalent_units():
                raise ValueError('{} does not have valid units of time: {}'.format(key, field.unit))
            # Convert to seconds when storing in DB
            value = float(field.to(u.s).value)
        else:
            # Assume it is already in seconds
            value = float(field)

        return value

    @hybrid_property
    def infinite(self):
        """Return True if the number of Pointings to generate is infinite."""
        return self.num_todo < 0

    # Replaced with a column property:

    # @hybrid_property
    # def num_completed(self):
    #     """Return the number of Pointings successfully completed."""
    #     return sum(p.status == 'completed' for p in self.pointings)

    # @num_completed.expression
    # def num_completed(self):
    #     return select(func.count(Pointing.db_id)).\
    #         where(and_(Pointing.strategy_id == self.db_id,
    #                    Pointing.status == 'completed',
    #                    )).scalar_subquery()

    @hybrid_method
    def num_completed_at_time(self, time):
        """Return the number of Pointings successfully completed at the given time."""
        if time is None:
            return self.num_completed
        return sum(p.status_at_time(time) == 'completed' for p in self.pointings)

    @num_completed_at_time.expression
    def num_completed_at_time(self, time):
        if time is None:
            return self.num_completed
        return select(func.count(Pointing.db_id)).\
            where(and_(Pointing.strategy_id == self.db_id,
                       Pointing.status_at_time(time) == 'completed',
                       )).scalar_subquery()

    @hybrid_property
    def num_remaining(self):
        """Return the number of observations remaining (num_todo - num_completed).

        Returns -1 if num_todo is infinite.
        """
        if self.infinite:
            return -1
        elif self.num_completed is None:
            # column property is None: not added to DB yet
            return self.num_todo
        else:
            return self.num_todo - self.num_completed

    @num_remaining.expression
    def num_remaining(self):
        return case((self.infinite, -1),
                    else_=self.num_todo - self.num_completed)

    @hybrid_method
    def num_remaining_at_time(self, time):
        """Return the number of observations remaining at the given time."""
        if time is None:
            return self.num_remaining
        if self.infinite:
            return -1
        elif self.num_completed_at_time(time) is None:
            # column property is None: not added to DB yet
            return self.num_todo
        else:
            return self.num_todo - self.num_completed_at_time(time)

    @num_remaining_at_time.expression
    def num_remaining_at_time(self, time):
        if time is None:
            return self.num_remaining
        return case((self.infinite, -1),
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
        return case((self.num_remaining == 0, True),
                    (and_(self.stop_time.isnot(None), time >= self.stop_time), True),
                    else_=False)

    @hybrid_method
    def finished_at_time(self, time):
        """Return True if this Strategy is was completed or expired at the given time."""
        if time is None:
            return self.finished
        if self.num_remaining_at_time(time) == 0:
            return True
        elif self.stop_time is not None and time >= Time(self.stop_time):
            return True
        else:
            return False

    @finished_at_time.expression
    def finished_at_time(self, time):
        if time is None:
            return self.finished
        if isinstance(time, str):
            time = Time(time)
        if isinstance(time, Time):
            time = time.datetime

        c = case((self.num_remaining_at_time(time) == 0, True),
                 (and_(self.stop_time.isnot(None), time >= self.stop_time), True),
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

    @hybrid_method
    def tel_is_valid(self, telescope_id):
        """Return True if the given Telescope ID is valid, based on the tel_mask."""
        if self.tel_mask is None:
            return True
        else:
            return bool(self.tel_mask & 1 << int(telescope_id) - 1)

    @tel_is_valid.expression
    def tel_is_valid(self, telescope_id):
        return case((self.tel_mask.is_(None), True),
                    else_=self.tel_mask.op('&', precedence=2, is_comparison=True)(1)
                                       .op('<<', precedence=1)(telescope_id - 1)
                    )

    def get_block_by_number(self, block_num):
        """Return the TimeBlock with the given number."""
        if len(self.time_blocks) == 0:
            # This shouldn't happen, since TimeBlocks are created on init
            raise ValueError('Strategy has no linked TimeBlocks')

        time_block = [block for block in self.time_blocks if block.block_num == block_num]
        if len(time_block) > 1:
            raise ValueError('Multiple TimeBlocks with the same number: {}'.format(block_num))
        else:
            return time_block[0]

    def get_previous_block(self, block):
        """Return the last TimeBlock executed before the given one."""
        if block.block_num == 1:
            # This is the first block, so loop back to the end
            block_num = len(self.time_blocks)
        else:
            block_num = block.block_num - 1
        return self.get_block_by_number(block_num)

    def get_next_block(self, block):
        """Return the next TimeBlock to be executed after the given one."""
        if block.block_num == len(self.time_blocks):
            # This is the last block, so loop back to the start
            block_num = 1
        else:
            block_num = block.block_num + 1
        return self.get_block_by_number(block_num)


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
    wait_time : float or `astropy.units.Quantity`
        time to wait after this block before allowing the next pointing.
        if a Quantity it should have units of time, if a float assume in seconds.
        default = 0 (no delay)
    valid_time : float, `astropy.units.Quantity`, or None
        amount of time a pointing in this block should stay valid in the queue.
        if a Quantity it should have units of time, if a float assume in seconds.
        less than zero or None means valid indefinitely.
        default = None (indefinitely valid)
    rank_change : int, optional, or None
        the amount to change the rank of the next pointing after this block is completed.
        0 or None means no change, positive integers will increase the rank (lower priority),
        negative integers will decrease the rank (higher priority).
        default = 10 (increasing rank i.e. lower priority as each pointing is observed)

    When created the instance can be linked to the following other tables as parameters,
    otherwise they are populated when it is added to the database:

    Primary relationships
    ---------------------
    strategy : `Strategy`
        the Strategy associated with this Time Block
        required before addition to the database
        can also be added with the target_id parameter

    pointings : list of `Pointing`, optional
        the Pointings associated with this Time Block, if any

    Attributes
    ----------
    db_id : int
        primary database key
        only populated when the instance is added to the database

    """

    # Set corresponding SQL table name and schema
    __tablename__ = 'time_blocks'
    __table_args__ = {'schema': SCHEMA}

    # Primary key
    db_id = Column('id', Integer, primary_key=True)

    # Columns
    block_num = Column(Integer, nullable=False, index=True)
    wait_time = Column(Float, nullable=False, default=0)
    valid_time = Column(Float, nullable=True, default=None)
    rank_change = Column(Integer, nullable=True, default=10)

    # Foreign keys
    strategy_id = Column(Integer, ForeignKey(f'{SCHEMA}.strategies.id'), nullable=False, index=True)

    # Update timestamp
    ts = Column(DateTime, nullable=False, server_default=func.now())

    # Foreign relationships
    strategy = relationship(
        'Strategy',
        order_by='Strategy.db_id',
        back_populates='time_blocks',
    )
    pointings = relationship(
        'Pointing',
        order_by='Pointing.db_id',
        back_populates='time_block',
    )

    def __repr__(self):
        strings = ['db_id={}'.format(self.db_id),
                   'block_num={}'.format(self.block_num),
                   'wait_time={}'.format(self.wait_time),
                   'valid_time={}'.format(self.valid_time),
                   'rank_change={}'.format(self.rank_change),
                   'strategy_id={}'.format(self.strategy_id),
                   ]
        return 'TimeBlock({})'.format(', '.join(strings))

    @validates('wait_time', 'valid_time')
    def validate_durations(self, key, field):
        """Use validators to allow various types of input for time periods."""
        if key == 'valid_time' and field is None:
            # valid_time is nullable, wait_time isn't
            return None

        if isinstance(field, u.Quantity):
            # Should have units of time
            if field.unit not in u.s.find_equivalent_units():
                raise ValueError('{} does not have valid units of time: {}'.format(key, field.unit))
            # Convert to seconds when storing in DB
            value = float(field.to(u.s).value)
        else:
            # Assume it is already in seconds
            value = float(field)

        return value


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

    rank : int or None
        rank to use for Pointings generated by this Target
        (modified by any rank_change values in the associated Strategy)
        lower values are prioritised first by the scheduler
        rank=None means infinite rank, will always be at the bottom of the queue (queue-fillers)
    weight : float, optional
        weighting relative to other Targets in the same survey
        default = 1
    is_template : bool, optional
        if True, completed Pointings for this Target will count as templates for their linked
        GridTiles (this parameter is meaningless unless this Target is linked to a GridTile)
        default = False

    start_time : string, `astropy.time.Time` or datetime.datetime, optional
        UTC time from which Target is considered valid and can be started.
        if not given then set to now, so the Target will start immediately.
        default = Time.now()
    stop_time : string, `astropy.time.Time` or datetime.datetime, or None, optional
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

    # Set corresponding SQL table name and schema
    __tablename__ = 'targets'
    __table_args__ = {'schema': SCHEMA}

    # Primary key
    db_id = Column('id', Integer, primary_key=True)

    # Columns
    # # Basic properties
    name = Column(Text, nullable=False)
    ra = Column(Float, nullable=False)
    dec = Column(Float, nullable=False)
    # # Scheduling
    rank = Column(Integer, nullable=True)
    weight = Column(Float, nullable=False, default=1)
    is_template = Column(Boolean, nullable=False, default=False)
    # # Constraints
    start_time = Column(DateTime, nullable=False, index=True, server_default=func.now())
    stop_time = Column(DateTime, nullable=True, index=True, default=None)
    # # Status
    creation_time = Column(DateTime, nullable=False, server_default=func.now())
    deleted_time = Column(DateTime, nullable=True, default=None)

    # Foreign keys
    user_id = Column(Integer, ForeignKey(f'{SCHEMA}.users.id'), nullable=False, index=True)
    grid_tile_id = Column(Integer, ForeignKey(f'{SCHEMA}.grid_tiles.id'), nullable=True, index=True)
    survey_id = Column(Integer, ForeignKey(f'{SCHEMA}.surveys.id'), nullable=True, index=True)

    # Update timestamp
    ts = Column(DateTime, nullable=False, server_default=func.now())

    # Foreign relationships
    # (remember to add to __init__)
    user = relationship(
        'User',
        order_by='User.db_id',
        back_populates='targets',
    )
    pointings = relationship(
        'Pointing',
        order_by='Pointing.db_id',
        back_populates='target',
    )
    exposure_sets = relationship(
        'ExposureSet',
        order_by='ExposureSet.db_id',
        back_populates='target',
    )
    strategies = relationship(
        'Strategy',
        order_by='Strategy.db_id',
        back_populates='target',
    )
    grid_tile = relationship(
        'GridTile',
        order_by='GridTile.db_id',
        back_populates='targets',
    )
    survey = relationship(
        'Survey',
        order_by='Survey.db_id',
        back_populates='targets',
    )

    # Secondary relationships
    grid = relationship(
        'Grid',
        order_by='Grid.db_id',
        secondary=f'{SCHEMA}.grid_tiles',
        primaryjoin='Target.grid_tile_id == GridTile.db_id',
        secondaryjoin='Grid.db_id == GridTile.grid_id',
        back_populates='targets',
        viewonly=True,
        uselist=False,
    )
    grid_id = association_proxy('grid', 'db_id')

    # Column properties
    scheduled = column_property(exists()
                                .where(and_(Pointing.target_id == db_id,
                                            or_(Pointing.status == 'upcoming',
                                                Pointing.status == 'pending',
                                                Pointing.status == 'running',
                                                )
                                            )
                                       )
                                .correlate_except(Pointing)
                                )
    num_completed = column_property(select(func.count(Pointing.db_id))
                                    .where(and_(Pointing.target_id == db_id,
                                                Pointing.status == 'completed',
                                                )
                                           )
                                    .correlate_except(Pointing)
                                    .scalar_subquery()
                                    )
    last_observed = column_property(select(Pointing.finished_time)
                                    .where(and_(Pointing.target_id == db_id,
                                                Pointing.status == 'completed',
                                                )
                                           )
                                    .order_by(Pointing.finished_time.desc())
                                    .limit(1)
                                    .correlate_except(Pointing)
                                    .scalar_subquery()
                                    )
    # Have to invert and return if any are incomplete, no way to do all()
    completed = column_property(~exists()
                                .where(and_(Strategy.target_id == db_id,
                                            Strategy.finished.is_(False),
                                            )
                                       )
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
                   'rank={}'.format(self.rank),
                   'weight={}'.format(self.weight),
                   'num_completed={}'.format(self.num_completed),
                   'start_time={}'.format(self.start_time),
                   'stop_time={}'.format(self.stop_time),
                   'creation_time={}'.format(self.creation_time),
                   'deleted_time={}'.format(self.deleted_time),
                   'user_id={}'.format(self.user_id),
                   'grid_tile_id={}'.format(self.grid_tile_id),
                   'survey_id={}'.format(self.survey_id),
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

        # force start_time < stop_time
        if (key == 'start_time' and self.stop_time is not None and
                Time(value) >= Time(self.stop_time)):
            raise ValueError(f'start_time must be before stop_time ({self.stop_time})')
        if (key == 'stop_time' and self.start_time is not None and
                Time(self.start_time) >= Time(value)):
            raise ValueError(f'stop_time must be after start_time ({self.start_time})')

        return value

    # Replaced with a column property:

    # @hybrid_property
    # def scheduled(self):
    #     """Return True if currently linked to a pending or running Pointing."""
    #     return any(p.status in ['upcoming', 'pending', 'running'] for p in self.pointings)

    # @scheduled.expression
    # def scheduled(self):
    #     return exists().where(and_(Pointing.target_id == self.db_id,
    #                                or_(Pointing.status == 'upcoming',
    #                                    Pointing.status == 'pending',
    #                                    Pointing.status == 'running',
    #                                    )
    #                                ))

    @hybrid_method
    def scheduled_at_time(self, time):
        """Return True if linked to a pending or running Pointing at the given time."""
        if time is None:
            return self.scheduled
        return any(p.status_at_time(time) in ['upcoming', 'pending', 'running']
                   for p in self.pointings)

    @scheduled_at_time.expression
    def scheduled_at_time(self, time):
        if time is None:
            return self.scheduled
        return exists().where(and_(Pointing.target_id == self.db_id,
                                   or_(Pointing.status_at_time(time) == 'upcoming',
                                       Pointing.status_at_time(time) == 'pending',
                                       Pointing.status_at_time(time) == 'running',
                                       )
                                   ))

    # Replaced with a column property:

    # @hybrid_property
    # def num_completed(self):
    #     """Return the number of Pointings successfully completed."""
    #     return sum(p.status == 'completed' for p in self.pointings)

    # @num_completed.expression
    # def num_completed(self):
    #     return select(func.count(Pointing.db_id)).\
    #         where(and_(Pointing.target_id == self.db_id,
    #                    Pointing.status == 'completed',
    #                    )).scalar_subquery()

    @hybrid_method
    def num_completed_at_time(self, time):
        """Return the number of Pointings successfully completed at the given time."""
        if time is None:
            return self.num_completed
        return sum(p.status_at_time(time) == 'completed' for p in self.pointings)

    @num_completed_at_time.expression
    def num_completed_at_time(self, time):
        if time is None:
            return self.num_completed
        return select(func.count(Pointing.db_id)).\
            where(and_(Pointing.target_id == self.db_id,
                       Pointing.status_at_time(time) == 'completed',
                       )).scalar_subquery()

    # Replaced with a column property:

    # @hybrid_property
    # def completed(self):
    #     """Are all Strategies linked to this Target completed (or expired)."""
    #     # I guess we're completed by default if we have no Strategies, but it's not ideal...
    #     return all(strategy.finished for strategy in self.strategies)

    # @completed.expression
    # def completed(self):
    #     # Have to invert and return if any are incomplete, no way to do all()
    #     return ~exists().where(and_(Strategy.target_id == self.db_id,
    #                                 Strategy.finished.is_(False)
    #                                 ))

    @hybrid_method
    def completed_at_time(self, time):
        """Were all Strategies linked to this Target completed (or expired) at the given time."""
        if time is None:
            return self.completed
        return all(strategy.finished_at_time(time) for strategy in self.strategies)

    @completed_at_time.expression
    def completed_at_time(self, time):
        if time is None:
            return self.completed
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
        return case(
            (self.deleted_time.isnot(None),
             'deleted'),
            (self.completed,
             'completed'),
            (and_(self.start_time.isnot(None), time < self.start_time),
             'upcoming'),
            (and_(self.stop_time.isnot(None), time >= self.stop_time),
             'expired'),
            (self.scheduled.is_(False),
             'unscheduled'),
            else_='scheduled')

    @hybrid_method
    def status_at_time(self, time):
        """Return a string giving the status at the given time (for testing and simulations)."""
        if time is None:
            return self.status
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
        if time is None:
            return self.status
        if isinstance(time, str):
            time = Time(time)
        if isinstance(time, Time):
            time = time.datetime

        c = case(
            (time < self.creation_time,
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
            else_='scheduled')
        return c

    def mark_deleted(self, time=None):
        """Mark this Target (and any pending Pointings) as deleted."""
        if time is None:
            time = Time.now()
        if isinstance(time, (str, datetime.datetime)):
            time = Time(time)

        current_status = self.status_at_time(time)
        if current_status == 'deleted':
            msg = f'Target {self.db_id} is already deleted (at={self.deleted_time})'
            raise ValueError(msg)
        elif current_status == 'expired':
            msg = f'Target {self.db_id} is already expired (at {self.stop_time})'
            raise ValueError(msg)
        elif current_status == 'completed':
            msg = f'Target {self.db_id} is already completed'
            raise ValueError(msg)

        # Pointings can use finished_time, but we need an explicit deleted_time
        self.deleted_time = time

        for pointing in self.pointings:
            if pointing.status_at_time(time) in ['upcoming', 'pending']:
                pointing.mark_deleted(time=time)

    def undelete():
        raise NotImplementedError

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
        """Get the currently valid Strategy."""
        return self.get_current_strategy()

    @strategy.setter
    def strategy(self, strategy):
        """Set the given Strategy as the only strategy for this Target.

        This is a sneaky way to shortcut `Target.strategies = [strategy]
        """
        self.strategies = [strategy]

    def get_next_pointing(self, time=None, reschedule_delay=None):
        """Retrieve the next Pointing which needs to be scheduled.

        The start and stop times of the Pointing are determined from
        the status of the previous Pointing and the TimeBlocks attached
        to the current Strategy.

        If the previous Pointing was interrupted or failed you can add an optional
        delay (in seconds) to stop repeatedly observing poor patches of the sky.

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
        if self.status_at_time(time) not in ['upcoming', 'unscheduled']:
            return None

        # As a safety check, see if the Target has any other Pointings that
        # are already pending or running.
        # If so the Target status should be 'scheduled' and therefore be
        # caught above, but strange things can happen.
        # This should prevent the cases where Targets have multiple
        # Pointings in the queue, hopefully!
        if len(self.pointings) > 0:
            statuses = [p.status_at_time(time) for p in self.pointings]
            if 'pending' in statuses or 'running' in statuses:
                return None

        # Get the current Strategy (or at the given time)
        strategy = self.get_current_strategy(time)

        if strategy is None:
            # All observations were completed (or we have no linked Strategies)
            return None

        # Get last "finished" Pointing for this Target
        # (could be completed, interrupted, failed or expired)
        # The checks earlier should mean there aren't any pending or running Pointings,
        # but we also want to ignore "deleted" Pointings when getting the times.
        # NB This Pointing might have been observed under a different Strategy,
        # so we can't just do current_block.pointings[-1] as we used to.
        finished_statuses = ['completed', 'interrupted', 'failed', 'expired']
        finished_pointings = [p for p in self.pointings
                              if p.status_at_time(time) in finished_statuses]
        if len(finished_pointings) == 0:
            latest_pointing = None
        else:
            latest_pointing = finished_pointings[-1]

        # Find the Pointing start_time
        if latest_pointing is None:
            # We haven't observed anything yet, so use the start_time from the Target
            start_time = self.start_time
            # The "next" block should be the first one in the Strategy list
            next_block = strategy.time_blocks[0]
        elif latest_pointing.status_at_time(time) in ['interrupted', 'failed']:
            # The Pointing was interrupted before it could complete or expire,
            # or it completed but was validated as bad and needs to be re-observed.
            if reschedule_delay is None:
                # Re-create the Pointing starting from the same time and try again.
                start_time = latest_pointing.start_time
            else:
                # We want to delay the Pointing becoming valid.
                # So we'll take this delay from the time the Pointing was finished
                # (interrupted sets finish_time, failed implies it was completed).
                start_time = Time(latest_pointing.finished_time) + reschedule_delay * u.s
            # Use the block as the previous Pointing
            next_block = latest_pointing.time_block
        else:
            # We go onto the next block if the previous pointing was completed,
            #  or if the previous pointing's valid time has expired.
            if latest_pointing.time_block.valid_time is not None:
                # We want to start after "wait_time" from when the last Pointing was due to finish.
                # (i.e. the stop_time, the end of its "block").
                # But that might have been foreshortened by the Target's stop_time, so we
                # take the previous Pointing's start_time, add the valid time (if it has one) to
                # get its nominal stop_time, then add the wait time.
                start_time = (Time(latest_pointing.start_time) +
                              latest_pointing.time_block.valid_time * u.s +
                              latest_pointing.time_block.wait_time * u.s)
            else:
                # There was no set "block" for this Pointing, since it had no valid_time.
                # Note in this case it's impossible for the Pointing to be expired,
                # so it must have been completed and we want to start "wait_time" after
                # the time it was completed.
                start_time = (Time(latest_pointing.finished_time) +
                              latest_pointing.time_block.wait_time * u.s)

            # Get the next TimeBlock after the one this Pointing used
            # Since the "current" strategy might be different from the one this Pointing was
            # observed under (e.g. that Strategy might have expired or been completed) we need
            # to check if the Strategy has changed.
            if latest_pointing.strategy == strategy:
                # The Pointing was observed under the current Strategy, so we can just use
                # the next block in the list.
                next_block = strategy.get_next_block(latest_pointing.time_block)
            else:
                # The Strategy has changed, so the new Pointing should take the first block from
                # the new Strategy.
                next_block = strategy.time_blocks[0]

        # Find the Pointing stop_time
        if next_block.valid_time is not None:
            # The stop_time is calculated by adding the valid_time to the start_time.
            stop_time = Time(start_time) + next_block.valid_time * u.s
            if self.stop_time and stop_time > self.stop_time:
                # The Pointing stop_time would overrun past when this Target is valid for,
                # so we force it to stop at the earlier time.
                stop_time = self.stop_time
        else:
            if self.stop_time:
                # The Target stop time takes priority over infinite valid time.
                stop_time = self.stop_time
            else:
                # The Pointing should never expire, so it doesn't have a stop_time.
                stop_time = None

        # Sanity check that the stop_time is after the start time.
        # It can break if the Target has a stop_time, then the Pointing would be impossible.
        if stop_time is not None and start_time >= stop_time:
            # In this case there is no valid Pointing, so we have to return None.
            return None

        # Find the Pointing rank
        if latest_pointing is None:
            # This is the first Pointing for this Target, so use the base rank from the Target.
            rank = self.rank
        else:
            # For Pointings after the first take the rank from the previous Pointing.
            # If it was completed add any rank_change from that Pointing's TimeBlock
            # (the default is increasing by 10 after each Pointing).
            # If it was interrupted, failed or expired then keep the same rank and re-observe.
            # Only applies if the rank and rank change are both not None.
            if (latest_pointing.status_at_time(time) == 'completed' and
                    latest_pointing.rank is not None and
                    latest_pointing.time_block.rank_change is not None):
                rank = latest_pointing.rank + latest_pointing.time_block.rank_change
            else:
                rank = latest_pointing.rank

        # Now create the pointing
        new_pointing = Pointing(rank=rank,
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
    location : `astropy.coordinates.EarthLocation`
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

    # Set corresponding SQL table name and schema
    __tablename__ = 'sites'
    __table_args__ = {'schema': SCHEMA}

    # Primary key
    db_id = Column('id', Integer, primary_key=True)

    # Columns
    name = Column(String(255), nullable=False, unique=True, index=True)
    latitude = Column(Float, nullable=False)
    longitude = Column(Float, nullable=False)
    height = Column(Float, nullable=False)

    # Update timestamp
    ts = Column(DateTime, nullable=False, server_default=func.now())

    # Foreign relationships
    telescopes = relationship(
        'Telescope',
        order_by='Telescope.db_id',
        back_populates='site',
    )

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
        return EarthLocation.from_geodetic(lat=self.latitude * u.deg,
                                           lon=self.longitude * u.deg,
                                           height=self.height * u.m)

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

    horizon : str, optional
        local path to file containing the telescope artificial horizon
        multiple paths can be included, separated by a ";"
        default = None

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

    Methods
    -------
    get_horizon() : 2-tuple of lists, or list of same
        get this Telescope's artificial horizon(s)

    """

    # Set corresponding SQL table name and schema
    __tablename__ = 'telescopes'
    __table_args__ = {'schema': SCHEMA}
    # Primary key
    db_id = Column('id', Integer, primary_key=True)

    # Columns
    name = Column(String(255), nullable=False, unique=True, index=True)
    horizon = Column(Text, nullable=True)

    # Foreign keys
    site_id = Column(Integer, ForeignKey(f'{SCHEMA}.sites.id'), nullable=False, index=True)
    grid_id = Column(Integer, ForeignKey(f'{SCHEMA}.grids.id'), nullable=True, index=True)

    # Update timestamp
    ts = Column(DateTime, nullable=False, server_default=func.now())

    # Foreign relationships
    site = relationship(
        'Site',
        uselist=False,
        back_populates='telescopes',
    )
    grid = relationship(
        'Grid',
        uselist=False,
        back_populates='telescopes',
    )
    pointings = relationship(
        'Pointing',
        order_by='Pointing.db_id',
        back_populates='telescope',
    )
    current_pointing = relationship(
        'Pointing',
        viewonly=True,
        uselist=False,
        primaryjoin='and_(Telescope.db_id==Pointing.telescope_id, Pointing.status=="running")',
    )

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
        return case((self.current_pointing.isnot(None), 'observing'),
                    else_='idle')

    @property
    def tel_mask(self):
        """Get the binary telescope mask for this Telescope."""
        if self.db_id is None:
            raise ValueError('Telescope needs to be added to database to assign ID')
        return 2 ** (self.db_id - 1)

    def get_horizon(self):
        """Get the az/alt values of this Telescope's artificial horizon(s)."""
        if self.horizon is None:
            return None
        paths = self.horizon.split(';')
        horizons = []
        for path in paths:
            azs, alts = np.loadtxt(path, usecols=(0, 1)).T
            horizons.append((azs, alts))
        if len(horizons) == 1:
            return horizons[0]
        return horizons


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

    # Set corresponding SQL table name and schema
    __tablename__ = 'grids'
    __table_args__ = {'schema': SCHEMA}

    # Primary key
    db_id = Column('id', Integer, primary_key=True)

    # Columns
    name = Column(String(255), nullable=False, unique=True, index=True)
    ra_fov = Column(Float, nullable=False)
    dec_fov = Column(Float, nullable=False)
    ra_overlap = Column(Float, nullable=False)
    dec_overlap = Column(Float, nullable=False)
    algorithm = Column(String(255), nullable=False)

    # Update timestamp
    ts = Column(DateTime, nullable=False, server_default=func.now())

    # Foreign relationships
    grid_tiles = relationship(
        'GridTile',
        order_by='GridTile.db_id',
        back_populates='grid',
    )
    telescopes = relationship(
        'Telescope',
        order_by='Telescope.db_id',
        back_populates='grid',
    )

    # Secondary relationships
    targets = relationship(
        'Target',
        order_by='Target.db_id',
        secondary=f'{SCHEMA}.grid_tiles',
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

    Methods
    -------
    get_templates(telescope_id=None) : list of Pointings
        returns any Pointings linked to this GridTile that have been successfully observed
        (i.e. status=completed) and has is_template=True.
        if a Telescope ID is given it will only return Pointings which have been observed
        by the given Telescope.
        see also the `get_templates_at_time()` method

    has_template(telescope_id=None) : bool
        returns True if at least one template Pointing has been observed (for the given Telescope).
        equivalent to `GridTile.len(get_templates(telescope_id)) >= 1`.
        see also the `has_template_at_time()` method

    get_template_telescopes() : list of int
        returns the Telescope IDs, if any, which have a completed template Pointing for this tile.
        see also the `get_template_telescopes_at_time()` method

    """

    # Set corresponding SQL table name and schema
    __tablename__ = 'grid_tiles'
    __table_args__ = {'schema': SCHEMA}

    # Primary key
    db_id = Column('id', Integer, primary_key=True)

    # Columns
    name = Column(String(255), nullable=False, index=True)
    ra = Column(Float, nullable=False)
    dec = Column(Float, nullable=False)

    # Foreign keys
    grid_id = Column(Integer, ForeignKey(f'{SCHEMA}.grids.id'), nullable=False, index=True)

    # Update timestamp
    ts = Column(DateTime, nullable=False, server_default=func.now())

    # Foreign relationships
    grid = relationship(
        'Grid',
        order_by='Grid.db_id',
        back_populates='grid_tiles',
    )
    targets = relationship(
        'Target',
        order_by='Target.db_id',
        back_populates='grid_tile',
    )

    # Secondary relationships
    pointings = relationship(
        'Pointing',
        order_by='Pointing.db_id',
        secondary=f'{SCHEMA}.targets',
        primaryjoin='GridTile.db_id == Target.grid_tile_id',
        secondaryjoin='Pointing.target_id == Target.db_id',
        back_populates='grid_tile',
        viewonly=True,
        uselist=True,
    )

    # Column properties
    last_observed = column_property(select(Target.last_observed)
                                    .where(Target.grid_tile_id == db_id,
                                           )
                                    .order_by(Target.last_observed.desc())
                                    .limit(1)
                                    .correlate_except(Target)
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

    def get_templates(self, telescope_id=None):
        """Return any completed template Pointings for this tile."""
        if telescope_id is None:
            return [p for p in self.pointings
                    if (p.is_template is True and
                        p.status == 'completed')
                    ]
        else:
            return [p for p in self.pointings
                    if (p.is_template is True and
                        p.status == 'completed' and
                        p.telescope_id == telescope_id)
                    ]

    def get_templates_at_time(self, time, telescope_id=None):
        """Return any completed template Pointings for this tile at the given time."""
        if time is None:
            return self.get_templates(telescope_id)
        if telescope_id is None:
            return [p for p in self.pointings
                    if (p.is_template is True and
                        p.status_at_time(time) == 'completed')
                    ]
        else:
            return [p for p in self.pointings
                    if (p.is_template is True and
                        p.status_at_time(time) == 'completed' and
                        p.telescope_id == telescope_id)
                    ]

    def has_template(self, telescope_id=None):
        """Return True if this tile has any completed template Pointings."""
        return len(self.get_templates(telescope_id)) >= 1

    def has_template_at_time(self, time, telescope_id=None):
        """Return True if this tile has any completed template Pointings at the given time."""
        if time is None:
            return self.has_template
        return len(self.get_templates_at_time(time, telescope_id)) >= 1

    def get_template_telescopes(self):
        """Return a list of Telescope IDs with completed template Pointings for this tile."""
        pointings = self.get_templates(telescope_id=None)
        return sorted({p.telescope.db_id for p in pointings})

    def get_template_telescopes_at_time(self, time):
        """Return a list of Telescope IDs with completed template Pointings at the given time."""
        if time is None:
            return self.get_template_telescopes()
        pointings = self.get_templates_at_time(time, telescope_id=None)
        return sorted({p.telescope.db_id for p in pointings})


class Survey(Base):
    """A class to represent a Survey.

    A Survey is a way to group multiple Targets together,
    usually covering a particular area of the sky.

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
    targets : list of `Target`, optional
        the Targets that are part of this Survey, if any

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

    # Set corresponding SQL table name and schema
    __tablename__ = 'surveys'
    __table_args__ = {'schema': SCHEMA}

    # Primary key
    db_id = Column('id', Integer, primary_key=True)

    # Columns
    name = Column(String(255), nullable=False, index=True)

    # Update timestamp
    ts = Column(DateTime, nullable=False, server_default=func.now())

    # Foreign relationships
    targets = relationship(
        'Target',
        order_by='Target.db_id',
        back_populates='survey',
    )

    # Secondary relationships
    pointings = relationship(
        'Pointing',
        order_by='Pointing.db_id',
        secondary=f'{SCHEMA}.targets',
        primaryjoin='Survey.db_id == Target.survey_id',
        secondaryjoin='Pointing.target_id == Target.db_id',
        back_populates='survey',
        viewonly=True,
        uselist=True,
    )

    def __repr__(self):
        strings = ['db_id={}'.format(self.db_id),
                   'name={}'.format(self.name),
                   ]
        return 'Survey({})'.format(', '.join(strings))


# Registries for other entities
# Note: These are created in the database via alembic-utils
functions = {}
triggers = {}

# Define ts update function and triggers for tracking changes
ts_function = 'update_ts()'
ts_function_sql = """RETURNS TRIGGER
LANGUAGE plpgsql AS
$function$
BEGIN
    NEW.ts := now();
    RETURN NEW;
END
$function$;
"""
functions[ts_function] = {
    'schema': SCHEMA,
    'signature': ts_function,
    'definition': ts_function_sql,
}
tables_with_ts = [
    table for table in Base.metadata.tables.values()
    if table.schema == SCHEMA and 'ts' in table.columns
]
for table in tables_with_ts:
    ts_trigger = f'trig_update_ts_{table.name}'
    ts_trigger_sql = f"""BEFORE UPDATE ON {SCHEMA}.{table.name}
    FOR EACH ROW EXECUTE FUNCTION {SCHEMA}.update_ts();
    """
    triggers[ts_trigger] = {
        'schema': SCHEMA,
        'signature': ts_trigger,
        'on_entity': f"{SCHEMA}.{table.name}",
        'definition': ts_trigger_sql,
    }
