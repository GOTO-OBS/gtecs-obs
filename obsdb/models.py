#!/usr/bin/env python
"""Python classes mapping on to database tables."""

import datetime

from astropy import units as u
from astropy.time import Time

from sqlalchemy import (Column, DateTime, Enum, Float, ForeignKey, Integer, String)
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship, validates


__all__ = ['User', 'Event', 'EventTile', 'Survey', 'SurveyTile', 'ExposureSet',
           'Pointing', 'ObservingBlock', 'Mpointing', 'ImageLog']


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
        userName : string
            a short user name
        password : string
            password for authentication. stored as a hash in the database.
        fullName : string
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
        >>> bob = User(userName='bob', password='1234', fullName="Bob Marley")
        >>> session = load_session()
        >>> session.add(bob)  # add bob to database
        >>> session.commit()  # commit changes (otherwise DB not changed)
        >>> bob
        User(db_id=25, username=bob, fullName=Bob Marley)
        >>> bob.pointings
        []
        >>> pointing = make_random_pointing(25)  # make a pointing for bob
        >>> session.add(pointing)
        >>> session.commit()
        >>> bob.pointings
        [Pointing(db_id=None, status='pending', objectName=randObj, ra=352.133, dec=28.464,
        rank=84, minAlt=30, maxSunAlt=-15, minTime=1575.8310236068696, maxMoon=D, minMoonSep=30,
        ToO=True, startUTC=2018-01-16 16:46:12, stopUTC=2018-01-18 16:46:12, startedUTC=None,
        stoppedUTC=None, user_id=25, mpointing_id=None, observing_block_id=None, event_id=None,
        event_tile_id=None, survey_id=None, survey_tile_id=None)]
        >>> session.close()

    """

    # Set corresponding SQL table name
    __tablename__ = 'users'

    # Primary key
    db_id = Column('id', Integer, primary_key=True)

    # Columns
    userName = Column('username', String)
    password = Column('password', String)
    fullName = Column('full_name', String)

    def __repr__(self):
        return "User(db_id={}, userName={}, fullName={})".format(
            self.db_id, self.userName, self.fullName
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
        objectName : String
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
        minAlt : float
            minimum altitude to observer at
        minTime : float
            minimum time needed to schedule pointing
        maxSunAlt : float
            altitude constraint on Sun
        maxMoon : string
            Moon constraint. one of 'D', 'G', 'B'.
        minMoonSep : float
            distance constraint from the Moon, degrees
        ToO : int
            0 or 1 to indicate if this is a ToO or not
        startUTC : string, `astropy.time.Time` or datetime.datetime
            UTC time from which pointing is considered valid and can be started
        stopUTC : string, `astropy.time.Time` or datetime.datetime, or None
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
        startedUTC : datetime.datetime, or None
            if the pointing has started (been marked running)
            this will give the time it was updated
        stoppedUTC : datetime.datetime, or None
            if the pointing has finished (either completed or cancelled for
            some reason) this will give the time it was updated
        exposure_sets : list of `ExposureSet`
            the `ExposureSet` objects associated with this `Pointing`, if any
        eventTile : `EventTile`
            the `EventTile` associated with this `Pointing`, if any
        surveyTile : `SurveyTile`
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

        >>> p = Pointing(objectName='IP Peg', ra=350.785625, dec=18.416472, rank=9, minAlt=30,
        ... maxSunAlt=-15, minTime=3600, maxMoon='G', minMoonSep=30, ToO=0, startUTC=Time.now(),
        ... stopUTC=Time.now()+3*u.day, user_id=24)
        >>> p
        Pointing(db_id=None, status='None', objectName=IP Peg, ra=350.785625, dec=18.416472,
        rank=9, minAlt=30, maxSunAlt=-15, minTime=3600, maxMoon=G, minMoonSep=None, ToO=False,
        startUTC=2018-01-12 17:39:54, stopUTC=2018-01-15 17:39:54, startedUTC=None, stoppedUTC=None,
        user_id=24, mpointing_id=None, observing_block_id=None, event_id=None, event_tile_id=None,
        survey_id=None, survey_tile_id=None)

        We can insert it into the database and the status and db_id will be set:

        >>> session.add(p)
        >>> session.commit()
        >>> p.status, p.db_id
        ('pending', 17073)

        At the moment, this pointing has no exposure sets. We can either add these to the
        `exposure_sets` attribute directly:

        >>> e1 = ExposureSet(typeFlag='SCIENCE', filt='L', expTime=20, numexp=20, binning=2)
        >>> p.exposure_sets.append(e1)

        or create `ExposureSet` instances with the `pointing_id` attribute set, and the database
        will take care of the rest:

        >>> e2 = ExposureSet(typeFlag='SCIENCE', filt='G', expTime=20, numexp=20, binning=2,
        ... pointing_id=17073)
        >>> e3 = ExposureSet(typeFlag='SCIENCE', filt='R', expTime=20, numexp=20, binning=2,
        ... pointing_id=17073)
        >>> insert_items(session, [e2, e3])
        >>> session.commit()
        >>> p.exposure_sets
        [ExposureSet(db_id=126601, raoff=0.0, decoff=0.0, typeFlag=SCIENCE, filt=L, expTime=20.0,
        numexp=20, binning=2, utMask=None, pointing_id=17073, mpointing_id=None),
        ExposureSet(db_id=126602, raoff=0.0, decoff=0.0, typeFlag=SCIENCE, filt=G, expTime=20.0,
        numexp=20, binning=2, utMask=None, pointing_id=17073, mpointing_id=None),
        ExposureSet(db_id=126603, raoff=0.0, decoff=0.0, typeFlag=SCIENCE, filt=R, expTime=20.0,
        numexp=20, binning=2, utMask=None, pointing_id=17073, mpointing_id=None)]

        >>> session.close()

    """

    # Set corresponding SQL table name
    __tablename__ = 'pointings'

    # Primary key
    db_id = Column('id', Integer, primary_key=True)

    # Columns
    status = Column('status', pointing_status_list, default='pending')
    objectName = Column('object', String)
    ra = Column('ra', Float)
    dec = Column('decl', Float)
    rank = Column('rank', Integer)
    minAlt = Column('min_alt', Float)
    maxSunAlt = Column('max_sunalt', Float, default=-15)
    minTime = Column('min_time', Float)
    maxMoon = Column('max_moon', String(1))
    minMoonSep = Column('min_moonsep', Float, default=30)
    ToO = Column('too', Integer)
    startUTC = Column('start_time', DateTime)
    stopUTC = Column('stop_time', DateTime)
    startedUTC = Column('started_time', DateTime, default=None)
    stoppedUTC = Column('stopped_time', DateTime, default=None)

    # Foreign keys
    user_id = Column('user_id', Integer, ForeignKey('users.id'), nullable=False)
    mpointing_id = Column('mpointing_id', Integer, ForeignKey('mpointings.id'), nullable=True)
    observing_block_id = Column('observing_block_id', Integer, ForeignKey('observing_blocks.id'), nullable=True)
    event_id = Column('event_id', Integer, ForeignKey('events.id'), nullable=True)
    event_tile_id = Column('event_tile_id', Integer, ForeignKey('event_tiles.id'), nullable=True)
    survey_id = Column('survey_id', Integer, ForeignKey('surveys.id'), nullable=True)
    survey_tile_id = Column('survey_tile_id', Integer, ForeignKey('survey_tiles.id'), nullable=True)

    # Foreign relationships
    user = relationship('User', backref='pointings', uselist=False)
    mpointing = relationship('Mpointing', backref='pointings', uselist=False)
    observing_block = relationship('ObservingBlock', back_populates='pointings', uselist=False)
    event = relationship('Event', backref='pointings')
    eventTile = relationship('EventTile', back_populates='pointings', uselist=False)
    survey = relationship('Survey', backref='pointings')
    surveyTile = relationship('SurveyTile', back_populates='pointings', uselist=False)

    def __repr__(self):
        template = ("Pointing(db_id={}, status='{}', " +
                    "objectName={}, ra={}, dec={}, rank={}, " +
                    "minAlt={}, maxSunAlt={}, minTime={}, maxMoon={}, minMoonSep={}, " +
                    "ToO={}, startUTC={}, stopUTC={}, startedUTC={}, stoppedUTC={}, " +
                    "user_id={}, mpointing_id={}, observing_block_id={}, " +
                    "event_id={}, event_tile_id={}, survey_id={}, survey_tile_id={})")
        return template.format(
            self.db_id, self.status, self.objectName, self.ra, self.dec, self.rank,
            self.minAlt, self.maxSunAlt, self.minTime, self.maxMoon, self.minMoonSep,
            bool(self.ToO), self.startUTC, self.stopUTC, self.startedUTC, self.stoppedUTC,
            self.user_id, self.mpointing_id, self.observing_block_id,
            self.event_id, self.event_tile_id, self.survey_id, self.survey_tile_id
        )

    @validates('startUTC', 'stopUTC')
    def munge_times(self, key, field):
        """Use validators to allow various types of input for UTC.

        Also enforce writeUTC > startUTC.
        """
        if key == 'stopUTC' and field is None:
            value = None
        elif isinstance(field, datetime.datetime):
            value = field.strftime("%Y-%m-%d %H:%M:%S")
        elif isinstance(field, Time):
            field.precision = 0  # no D.P on seconds
            value = field.iso
        else:
            # just hope the string works!
            value = str(field)

        if (key == 'startUTC' and self.stopUTC is not None):
            if Time(value) >= Time(self.stopUTC):
                raise AssertionError("stopUTC must be later than startUTC")
        elif key == 'stopUTC' and value is not None and self.startUTC is not None:
            if Time(self.startUTC) >= Time(value):
                raise AssertionError("stopUTC must be later than startUTC")

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
        typeFlag : string
            indicates the type of exposure set.
            one of SCIENCE, FOCUS, STD, FLAT, BIAS, DARK
        filt : string
            filter to use
        expTime : float
            exposure time in seconds
        numexp : int
            number of exposures within the set
        binning : int
            binning to apply

        raoff : float, optional
            the size of the random offset to apply between each exposure
            if not set, no offset will be made
        decoff : float, optional
            the size of the random offset to apply between each exposure
            if not set, no offset will be made
        utMask : int, optional
            if set, this is a binary mask which will determine which unit
            telescopes carry out the exposure. A value of 5 (binary 0101) will
            be exposed by cameras 1 and 3.
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
    numexp = Column('num_exp', Integer)
    expTime = Column('exptime', Float)
    filt = Column('filter', String(2))
    binning = Column('binning', Integer)
    typeFlag = Column('imgtype', Enum('SCIENCE', 'FOCUS', 'DARK', 'BIAS', 'FLAT', 'STD'))
    utMask = Column('ut_mask', Integer, nullable=True)
    raoff = Column('ra_offset', Float, server_default='0.0')
    decoff = Column('dec_offset', Float, server_default='0.0')

    # Foreign keys
    pointing_id = Column('pointing_id', Integer, ForeignKey('pointings.id'), nullable=False)
    mpointing_id = Column('mpointing_id', Integer, ForeignKey('mpointings.id'), nullable=False)

    # Foreign relationships
    pointing = relationship('Pointing', backref='exposure_sets', uselist=False)
    mpointing = relationship('Mpointing', backref='exposure_sets', uselist=False)

    def __repr__(self):
        template = ("ExposureSet(db_id={}, raoff={}, decoff={}, typeFlag={}, " +
                    "filt={}, expTime={}, numexp={}, binning={}, utMask={}, " +
                    "pointing_id={}, mpointing_id={})")
        return template.format(
            self.db_id, self.raoff, self.decoff, self.typeFlag, self.filt,
            self.expTime, self.numexp, self.binning, self.utMask,
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
        objectName : String
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
        minAlt : float
            minimum altitude to observer at
        minTime : float
            minimum time needed to schedule pointing
        maxSunAlt : float
            altitude constraint on Sun
        maxMoon : string
            Moon constraint. one of 'D', 'G', 'B'.
        minMoonSep : float
            distance constraint from the Moon, degrees
        ToO : int
            0 or 1 to indicate if this is a ToO or not
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
        startUTC : string, `astropy.time.Time` or datetime.datetime, optional
            UTC time from which Mpointing is considered valid and can be started
            if not given then set to now, so the Mpointing will start immediately
        stopUTC : string, `astropy.time.Time` or datetime.datetime, optional
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
        observing_blocks : list of `ObservingBlock`
            the `ObservingBlock` objects associated with this `Mpointing`, if any
        surveyTile : `SurveyTile`
            the `SurveyTile` associated with this `Mpointing`, if any
        eventTile : `EventTile`
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

        >>> mp = Mpointing(objectName='M31', ra=22, dec=-5, start_rank=9, minAlt=30, minTime=3600,
        ... maxSunAlt=-15, ToO=0, maxMoon='B', minMoonSep=30, num_todo=5, user_id=24, valid_time=5,
        ... wait_time=10, startUTC=Time('2018-01-01 00:00:00'))
        >>> mp
        Mpointing(db_id=None, status='unscheduled', num_todo=5, num_completed=0,
        num_remaining=5, infinite=False, objectName=M31, ra=22, dec=-5, rank=9, start_rank=9,
        minAlt=30, maxSunAlt=-15, minTime=3600, maxMoon=B, minMoonSep=30, ToO=False,
        startUTC=2018-01-01 00:00:00, stopUTC=None, user_id=24, event_id=None, event_tile_id=None,
        survey_id=None, survey_tile_id=None)

        Note that the db_id is None, because that will be filled out when we add it to the
        database.

        Looking at the ObservingBlocks you can see that we only need one, and it has been generated

        >>> mp.observing_blocks
        [ObservingBlock(db_id=None, blockNum=1, valid_time=5, wait_time=10, current=1,
        mpointing_id=None)]

        This block will be repeated 5 times, as requested.

        ~~~~~~~~~~~~~~~~~~~~~

        For a more complicated example, give a list to wait_time to have the intervals between
        pointings increase.

        >>> mp = Mpointing(objectName='M31', ra=22, dec=-5, start_rank=9, minAlt=30, minTime=3600,
        ... maxSunAlt=-15, ToO=0, maxMoon='B', minMoonSep=30, num_todo=5, user_id=24, valid_time=5,
        ... wait_time=[10,20,30], startUTC=Time('2018-01-01 00:00:00'))
        >>> mp
        Mpointing(db_id=None, status='unscheduled', num_todo=5, num_completed=0,
        num_remaining=5, infinite=False, objectName=M31, ra=22, dec=-5, rank=9, start_rank=9,
        minAlt=30, maxSunAlt=-15, minTime=3600, maxMoon=B, minMoonSep=30, ToO=False,
        startUTC=2018-01-01 00:00:00, stopUTC=None, user_id=24, event_id=None, event_tile_id=None,
        survey_id=None, survey_tile_id=None)
        >>> mp.observing_blocks
        [ObservingBlock(db_id=None, blockNum=1, valid_time=5, wait_time=10, current=1,
        mpointing_id=None),
        ObservingBlock(db_id=None, blockNum=2, valid_time=5, wait_time=20, current=None,
        mpointing_id=None),
        ObservingBlock(db_id=None, blockNum=3, valid_time=5, wait_time=30, current=None,
        mpointing_id=None)]

        Note this time that we only created 3 blocks, but are asking for 5 observations.
        That's not a problem, as the blocks will simply repeat.
        In this case the sequence will look like this:
        [note we set the Mpointing start time to midnight]
        Pointing 1: startUTC=00:00, stopUTC=00:05 (valid for 5 minutes)
        Pointing 2: startUTC=00:15, stopUTC=00:20 (wait for 10, then valid for 5)
        Pointing 3: startUTC=00:40, stopUTC=00:45 (wait for 20, then valid for 5)
        Pointing 4: startUTC=01:15, stopUTC=01:20 (wait for 30, then valid for 5)
        Pointing 5: startUTC=01:30, stopUTC=01:35 (wait for 10, then valid for 5)

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
        [Pointing(db_id=17100, status='pending', objectName=M31, ra=22.0, dec=-5.0, rank=9,
        minAlt=30.0, maxSunAlt=-15.0, minTime=3600.0, maxMoon=B, minMoonSep=30.0, ToO=False,
        startUTC=2018-01-01 00:00:00, stopUTC=2018-01-01 00:05:00, startedUTC=None, stoppedUTC=None,
        user_id=24, mpointing_id=2, observing_block_id=1, event_id=None, event_tile_id=None,
        survey_id=None, survey_tile_id=None)]

        Note that the db_id, mpointing_id and observing_block_id have been filled out,
        and this Pointing has observing_block_id=1.
        Also that the start time is midnight and the stop time is 5 past, as expected.
        See what happens if we mark the pointing as completed:

        >>> mp.pointings[0].status = 'completed'
        >>> s.commit()
        >>> mp
        Mpointing(db_id=2, status='unscheduled', num_todo=5, num_completed=1, num_remaining=4,
        infinite=False, objectName=M31, ra=22.0, dec=-5.0, rank=19, start_rank=9, minAlt=30.0,
        maxSunAlt=-15.0, minTime=3600.0, maxMoon=B, minMoonSep=30.0, ToO=False,
        startUTC=2018-01-01 00:00:00, stopUTC=None, user_id=24, event_id=None, event_tile_id=None,
        survey_id=None, survey_tile_id=None)

        See that the num_completed atribute has gone up, the num_remaining has gone down
        and the mpointing status has changed to 'unscheduled'.

        Create the next pointing and add it:

        >>> p = mp.get_next_pointing()
        >>> s.add(p)
        >>> s.commit()
        >>> p
        Pointing(db_id=17102, status='pending', objectName=M31, ra=22.0, dec=-5.0, rank=19,
        minAlt=30.0, maxSunAlt=-15.0, minTime=3600.0, maxMoon=B, minMoonSep=30.0, ToO=False,
        startUTC=2018-01-01 00:15:00, stopUTC=2018-01-01 00:20:00, startedUTC=None, stoppedUTC=None,
        user_id=24, mpointing_id=2, observing_block_id=2, event_id=None, event_tile_id=None,
        survey_id=None, survey_tile_id=None)
        >>> mp
        Mpointing(db_id=2, status='scheduled', num_todo=5, num_completed=1, num_remaining=4,
        infinite=False, objectName=M31, ra=22.0, dec=-5.0, rank=19, start_rank=9, minAlt=30.0,
        maxSunAlt=-15.0, minTime=3600.0, maxMoon=B, minMoonSep=30.0, ToO=False,
        startUTC=2018-01-01 00:00:00, stopUTC=None, user_id=24, event_id=None, event_tile_id=None,
        survey_id=None, survey_tile_id=None)

        See that the Mpointing is back to scheduled.
        Also, note that the new pointing has start time of 00:15 and stop of 00:20. That's as we
        expected, because it's linked to the second observing block not the first
        (see observing_block_id=2).

        If you look at the observing blocks you can see the next one is now marked as current.

        >>> mp.observing_blocks
        [ObservingBlock(db_id=1, blockNum=1, valid_time=5.0, wait_time=10.0, current=0,
        mpointing_id=2),
        ObservingBlock(db_id=2, blockNum=2, valid_time=5.0, wait_time=20.0, current=1,
        mpointing_id=2),
        ObservingBlock(db_id=3, blockNum=3, valid_time=5.0, wait_time=30.0, current=0,
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
        >>> p.observing_block_id
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
        >>> p.observing_block_id
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
        >>> p.observing_block_id
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
        infinite=False, objectName=M31, ra=22.0, dec=-5.0, rank=59, start_rank=9, minAlt=30.0,
        maxSunAlt=-15.0, minTime=3600.0, maxMoon=B, minMoonSep=30.0, ToO=False,
        startUTC=2018-01-01 00:00:00, stopUTC=None, user_id=24, event_id=None, event_tile_id=None,
        survey_id=None, survey_tile_id=None)

        And the Mpointing is completed.

        ~~~~~~~~~~~~~~~~~~~~~

        To be useful, an Mpointing should have a list of `ExposureSet`s associated with it. We
        can either add these to the `exposure_sets` attribute directly:

        >>> e1 = ExposureSet(typeFlag='SCIENCE', filt='L', expTime=20, numexp=20, binning=2)
        >>> mp.exposure_sets.append(e1)

        or create `ExposureSet` instances with the `mpointing_id` attribute set, and the database
        will take care of the rest:

        >>> e2 = ExposureSet(typeFlag='SCIENCE', filt='G', expTime=20, numexp=20, binning=2,
        ... mpointing_id=1)
        >>> e3 = ExposureSet(typeFlag='SCIENCE', filt='R', expTime=20, numexp=20, binning=2,
        ... mpointing_id=1)
        >>> insert_items(session, [e2, e3])
        >>> session.commit()
        >>> mp.exposure_sets
        [ExposureSet(db_id=126598, raoff=0.0, decoff=0.0, typeFlag=SCIENCE, filt=L, expTime=20.0,
        numexp=20, binning=2, utMask=None, pointing_id=None, mpointing_id=1),
        ExposureSet(db_id=126599, raoff=0.0, decoff=0.0, typeFlag=SCIENCE, filt=G, expTime=20.0,
        numexp=20, binning=2, utMask=None, pointing_id=None, mpointing_id=1),
        ExposureSet(db_id=126600, raoff=0.0, decoff=0.0, typeFlag=SCIENCE, filt=R, expTime=20.0,
        numexp=20, binning=2, utMask=None, pointing_id=None, mpointing_id=1)]

        These exposure sets will be copied to the Pointings when they're created.

        >>> session.close()

    """

    # Set corresponding SQL table name
    __tablename__ = 'mpointings'

    # Primary key
    db_id = Column('id', Integer, primary_key=True)

    # Columns
    status = Column('status', mpointing_status_list, default='unscheduled')
    objectName = Column('object', String)
    ra = Column('ra', Float)
    dec = Column('decl', Float)
    rank = Column('rank', Integer)
    start_rank = Column('start_rank', Integer)
    num_todo = Column('num_todo', Integer)
    num_completed = Column('num_completed', Integer)
    infinite = Column('infinite', Integer, default=False)
    minAlt = Column('min_alt', Float)
    maxSunAlt = Column('max_sunalt', Float)
    minTime = Column('min_time', Float)
    maxMoon = Column('max_moon', String(1))
    minMoonSep = Column('min_moonsep', Float)
    ToO = Column('too', Integer)
    startUTC = Column('start_time', DateTime)
    stopUTC = Column('stop_time', DateTime)

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
    eventTile = relationship('EventTile', back_populates='mpointing', uselist=False)
    surveyTile = relationship('SurveyTile', back_populates='mpointing', uselist=False)
    observing_blocks = relationship('ObservingBlock', back_populates='mpointing', viewonly=True)

    def __repr__(self):
        template = ("Mpointing(db_id={}, status='{}', num_todo={}, num_completed={}," +
                    "num_remaining={}, infinite={}, " +
                    "objectName={}, ra={}, dec={}, rank={}, start_rank={}, " +
                    "minAlt={}, maxSunAlt={}, minTime={}, maxMoon={}, minMoonSep={}, " +
                    "ToO={}, startUTC={}, stopUTC={}, " +
                    "user_id={}, event_id={}, event_tile_id={}, survey_id={}, survey_tile_id={})")
        return template.format(
            self.db_id, self.status, self.num_todo, self.num_completed,
            self.num_remaining, bool(self.infinite),
            self.objectName, self.ra, self.dec, self.rank, self.start_rank,
            self.minAlt, self.maxSunAlt, self.minTime, self.maxMoon, self.minMoonSep,
            bool(self.ToO), self.startUTC, self.stopUTC,
            self.user_id, self.event_id, self.event_tile_id, self.survey_id, self.survey_tile_id
        )

    def __init__(self, objectName=None, ra=None, dec=None,
                 start_rank=None, minAlt=None, minTime=None,
                 maxMoon=None, minMoonSep=None, maxSunAlt=None, ToO=None, startUTC=None,
                 stopUTC=None, num_todo=None, valid_time=None, wait_time=None,
                 status='unscheduled', **kwargs):
        self.ra = ra
        self.dec = dec
        self.objectName = objectName
        self.start_rank = start_rank
        self.maxMoon = maxMoon
        self.minMoonSep = minMoonSep
        self.minAlt = minAlt
        self.minTime = minTime
        self.maxSunAlt = maxSunAlt
        self.ToO = ToO
        self.startUTC = startUTC if startUTC is not None else Time.now()
        self.stopUTC = stopUTC
        self.rank = self.start_rank
        self.status = status
        self.infinite = False
        self.num_todo = num_todo
        self.num_completed = 0

        # now add observing blocks
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

            # create ObservingBlock objects
            for i in range(max(len(valid_times), len(wait_times))):
                valid = valid_times[i % len(valid_times)]
                wait = wait_times[i % len(wait_times)]

                # check if non-expiring
                if valid < 0:
                    valid = -1

                block = ObservingBlock(blockNum=i + 1, valid_time=valid, wait_time=wait)
                self.observing_blocks.append(block)
            if len(self.observing_blocks):
                self.observing_blocks[0].current = 1

        if 'user' in kwargs:
            self.user = kwargs['user']
        if 'user_id' in kwargs:
            self.user_id = kwargs['user_id']
        if 'survey_id' in kwargs:
            self.survey_id = kwargs['survey_id']
        if 'survey' in kwargs:
            self.survey = kwargs['survey']
        if 'surveyTile' in kwargs:
            self.surveyTile = kwargs['surveyTile']
        if 'survey_tile_id' in kwargs:
            self.survey_tile_id = kwargs['survey_tile_id']
        if 'event' in kwargs:
            self.event = kwargs['event']
        if 'event_id' in kwargs:
            self.event_id = kwargs['event_id']
        if 'eventTile' in kwargs:
            self.eventTile = kwargs['eventTile']
        if 'event_tile_id' in kwargs:
            self.event_tile_id = kwargs['event_tile_id']

    @validates('startUTC', 'stopUTC')
    def munge_times(self, key, field):
        """Use validators to allow various types of input for UTC.

        Also enforce stopUTC > startUTC
        NB stopUTC is None be default
        """
        if key == 'stopUTC' and field is None:
            value = None
        elif isinstance(field, datetime.datetime):
            value = field.strftime("%Y-%m-%d %H:%M:%S")
        elif isinstance(field, Time):
            field.precision = 0  # no D.P on seconds
            value = field.iso
        else:
            # just hope the string works!
            value = str(field)

        if (key == 'startUTC' and self.stopUTC is not None):
            if Time(value) >= Time(self.stopUTC):
                raise AssertionError("stopUTC must be later than startUTC")
        elif key == 'stopUTC' and value is not None and self.startUTC is not None:
            if Time(self.startUTC) >= Time(value):
                raise AssertionError("stopUTC must be later than startUTC")

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
        """Return the current observing block.

        Assumes this object is still associated to an active session.

        Returns
        -------
        current : `obsdb.ObservingBlock`
            The current observing block (may be None).

        """
        current_block = [block for block in self.observing_blocks if block.current]
        if len(current_block) == 0:
            return None
        elif len(current_block) > 1:
            raise ValueError('Multiple observing blocks marked as current')
        else:
            return current_block[0]

    def get_last_block(self):
        """Return the last observing block executed.

        Assumes this object is still associated to an active session.

        Returns
        -------
        last : `obsdb.ObservingBlock`
            The last block done (may be None).

        """
        current_block = self.get_current_block()
        if not current_block:
            return None

        current_num = current_block.blockNum
        if current_num == 1:
            last_num = len(self.observing_blocks)
        else:
            last_num = current_num - 1

        last_block = [block for block in self.observing_blocks
                      if block.blockNum == last_num][0]
        return last_block

    def get_next_block(self):
        """Return the next observing block to be executed.

        Assumes this object is still associated to an active session.

        Returns
        -------
        next : `obsdb.ObservingBlock`
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
            current_num = current_block.blockNum
            latest_pointing = current_block.pointings[-1]
            if latest_pointing.status in ['completed', 'expired']:
                next_num = (current_num % len(self.observing_blocks)) + 1
            else:
                next_num = current_num

        next_block = [block for block in self.observing_blocks
                      if block.blockNum == next_num][0]
        return next_block

    def get_next_pointing(self):
        """Retrieve the next pointing which needs to be scheduled.

        The start and stop UTC of the pointing are determined from
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
            startUTC = self.startUTC
        else:
            latest_pointing = current_block.pointings[-1]
            # decide to go onto the next block:
            #  - if the last block's pointing was completed
            #  - if the last block's pointing's valid time has expired
            if latest_pointing.status in ['completed', 'expired']:
                if latest_pointing.stopUTC:
                    startUTC = Time(latest_pointing.stopUTC) + current_block.wait_time * u.minute
                else:
                    # non-expiring pointings have no stopUTC
                    startUTC = Time.now() + current_block.wait_time * u.minute
            else:
                # the current block wasn't completed, and there's still time left
                # (e.g. aborted, interrupted)
                # need to re-insert the current block with a new pointing
                startUTC = latest_pointing.startUTC

        if next_block.valid_time < 0:
            # non-expiring pointings
            stopUTC = None
        else:
            stopUTC = Time(startUTC) + next_block.valid_time * u.minute
            if self.stopUTC and stopUTC > self.stopUTC:
                # force pointings to stop by an Mpointing's stopUTC, if given
                stopUTC = self.stopUTC

            if startUTC >= stopUTC:
                # can happen if the Mpointing has a stopUTC
                return None

        # now create a pointing
        p = Pointing(objectName=self.objectName,
                     ra=self.ra,
                     dec=self.dec,
                     rank=self.rank,
                     minAlt=self.minAlt,
                     maxSunAlt=self.maxSunAlt,
                     minTime=self.minTime,
                     maxMoon=self.maxMoon,
                     minMoonSep=self.minMoonSep,
                     ToO=self.ToO,
                     startUTC=startUTC,
                     stopUTC=stopUTC,
                     status='pending',
                     user_id=self.user_id,
                     mpointing_id=self.db_id,
                     observing_block=next_block,
                     event_id=self.event_id,
                     event_tile_id=self.event_tile_id,
                     survey_id=self.survey_id,
                     survey_tile_id=self.survey_tile_id,
                     )
        # add the exposures
        p.exposure_sets = self.exposure_sets
        return p


class ObservingBlock(Base):
    """A class to represent a block of observing time.

    Like all SQLAlchemy model classes, this object links to the
    underlying database. You can create a ObservingBlock, and set its attributes
    without a database session. Accessing some attributes may require
    an active database session, and some properties (like the db_id)
    will be None until the ObservingBlock is added to the database.

    The constructor must use keyword arguments and the arguments below
    should be supplied or set before insertion into the database.
    See Examples for details.

    Note you shouldn't ever have to manually create ObservingBlocks.
    They're used internally by Mpointings to keep track of Pointings and all
    the ones an Mpointing requires are created when the Mpointing is.

    Args
    ----
        blockNum : int
            an integer indicating which block in a sequence this is
        valid_time : float
            amount of time a pointing in this block should stay valid in the queue, in minutes.
        wait_time : float
            time to wait after this block before allowing the next pointing, in minutes

        mpointing_id : int, optional
            unique key linking this ObservingBlock to an Mpointing
        current : bool, optional
            True if this ObservingBlock is the one that is currently linked to
            a Pointing in the queue

    An ObservingBlock also has the following properties which are
    populated through database queries, but not needed for
    object creation:

    Attributes
    ----------
        pointings : list of `Pointing`
            a list of any `Pointing`s associated with this block, if any

    Examples
    --------
        >>> from obsdb import *

        make an ObservingBlock
        (an Mpointing with mpointing_id = 1 is already in the database)

        >>> b = ObservingBlock(blockNum=1, valid_time=60, wait_time=120, mpointing_id=1)
        >>> session = load_session()
        >>> session.add(b)
        ObservingBlock(db_id=7, blockNum=1, valid_time=60.0, wait_time=120.0, current=False,
        mpointing_id=1)

    """

    # Set corresponding SQL table name
    __tablename__ = 'observing_blocks'

    # Primary key
    db_id = Column('id', Integer, primary_key=True)

    # Columns
    blockNum = Column('block_num', Integer)
    valid_time = Column('valid_time', Integer)
    wait_time = Column('wait_time', Integer)
    current = Column('current', Integer, default=False)

    # Foreign keys
    mpointing_id = Column('mpointing_id', Integer, ForeignKey('mpointings.id'), nullable=False)

    # Foreign relationships
    mpointing = relationship('Mpointing', back_populates='observing_blocks', uselist=False)
    pointings = relationship('Pointing', back_populates='observing_block')

    def __repr__(self):
        template = ("ObservingBlock(db_id={}, blockNum={}, valid_time={}, " +
                    "wait_time={}, current={}, mpointing_id={})")
        return template.format(self.db_id, self.blockNum, self.valid_time,
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
        eventTiles : list of `EventTile`
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
    eventTiles = relationship('EventTile', back_populates='event')

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
        >>> e.eventTiles  # and event knows about all associated tiles
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
        >>> et2.surveyTile = st
        >>> et2
        EventTile(db_id=None, ra=None, dec=None, probability=0.01, event_id=None,
        survey_tile_id=None)

        add to the database

        >>> session.add(event_tile)
        >>> session.commit()
        >>> et2  # note how the ra and dec have been copied from the surveyTile
        EventTile(db_id=2, ra=22.0, dec=-2.0, probability=0.01, event_id=1, survey_tile_id=1)
        >>> st.eventTiles # and the surveyTile knows about the linked eventTiles
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
    pointings = relationship('Pointing', back_populates='eventTile')
    mpointing = relationship('Mpointing', back_populates='eventTile')
    event = relationship('Event', back_populates='eventTiles', uselist=False)
    surveyTile = relationship('SurveyTile', back_populates='eventTiles', uselist=False)

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
        surveyTiles : list of `SurveyTile`
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
    surveyTiles = relationship('SurveyTile', back_populates='survey')

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
        eventTiles : list of `EventTile`
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
        >>>     s.surveyTiles  # and survey knows about all associated tiles
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
    survey = relationship('Survey', back_populates='surveyTiles', uselist=False)
    eventTiles = relationship('EventTile', back_populates='surveyTile')
    mpointing = relationship('Mpointing', back_populates='surveyTile')
    pointings = relationship('Pointing', back_populates='surveyTile')

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
        runNumber : int
            the run ID number for this exposure
        ut : int
            the unit telescope this frame was captured on
        utMask : int
            a binary mask for which unit telescopes carried out this exposure.
            A value of 5 (binary 0101) will have been exposed by cameras 1 and 3.
        startUTC : string, `astropy.time.Time` or datetime.datetime
            the time that the exposure began
        writeUTC : string, `astropy.time.Time` or datetime.datetime
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
    runNumber = Column('run_number', Integer)
    ut = Column('ut', Integer)
    utMask = Column('ut_mask', Integer)
    startUTC = Column('start_time', DateTime)
    writeUTC = Column('write_time', DateTime)
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
        template = ("ImageLog(db_id={}, filename={}, runNumber={}, " +
                    "ut={}, utMask={}, startUTC={}, writeUTC={}, " +
                    "exposure_set_id={}, pointing_id={}, mpointing_id={})")
        return template.format(
            self.db_id, self.filename, self.runNumber, self.ut, self.utMask,
            self.startUTC, self.writeUTC, self.exposure_set_id, self.pointing_id,
            self.mpointing_id
        )

    @validates('startUTC', 'writeUTC')
    def munge_times(self, key, field):
        """Use validators to allow various types of input for UTC.

        Also enforce writeUTC > startUTC.
        """
        if isinstance(field, datetime.datetime):
            value = field.strftime("%Y-%m-%d %H:%M:%S")
        elif isinstance(field, Time):
            field.precision = 0  # no D.P on seconds
            value = field.iso
        else:
            # just hope the string works!
            value = str(field)

        if (key == 'startUTC' and self.writeUTC is not None):
            if Time(value) >= Time(self.writeUTC):
                raise AssertionError("writeUTC must be later than startUTC")
        elif key == 'writeUTC' and self.writeUTC is not None:
            if Time(self.startUTC) >= Time(value):
                raise AssertionError("writeUTC must be later than startUTC")

        return value
