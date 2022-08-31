"""Database management functions."""

import os
from contextlib import contextmanager

import numpy as np

import pymysql

from sqlalchemy import create_engine
from sqlalchemy.exc import ProgrammingError
from sqlalchemy.orm import sessionmaker

from .models import Base, TRIGGERS
from .. import params

__all__ = ['create_database', 'fill_database',
           'get_engine', 'open_session', 'load_session']


# Encode Numpy floats
# https://stackoverflow.com/questions/46205532/
pymysql.converters.encoders[np.float64] = pymysql.converters.escape_float
pymysql.converters.conversions = pymysql.converters.encoders.copy()
pymysql.converters.conversions.update(pymysql.converters.decoders)


def get_engine(url=params.DATABASE_URL, db_name=params.DATABASE_NAME,
               echo=params.DATABASE_ECHO, *args, **kwargs):
    """Create the database engine."""
    if db_name:
        url = os.path.join(url, db_name)
    engine = create_engine(f'mysql+pymysql://{url}?charset=utf8',
                           echo=echo,
                           pool_pre_ping=params.DATABASE_PRE_PING,
                           *args, *kwargs)
    return engine


def create_database(overwrite=False, verbose=False):
    """Create the blank database.

    The name of the database is 'goto_obs' by default, it can be changed in
    `gtecs.obs.params.DATABASE_NAME`.

    Parameters
    ----------
    overwrite : bool, default=False
        If True and the database already exists then drop it before creating the new one.
        If False and the database already exists then an error is raised.

    verbose : bool, default=False
        If True, echo SQL output.

    """
    db_name = params.DATABASE_NAME
    engine = get_engine(db_name=None, echo=verbose)
    with engine.connect() as conn:
        create_command = f'CREATE DATABASE `{db_name}`'
        # Set default encoding to UTF8 (see https://dba.stackexchange.com/questions/76788)
        create_command += 'CHARACTER SET utf8mb4 COLLATE utf8mb4_unicode_ci'
        if not overwrite:
            try:
                conn.execute(create_command)  # will raise if it exists
            except ProgrammingError as err:
                raise ValueError(f'WARNING: Database "{db_name}" already exists!') from err
        else:
            conn.execute(f'DROP DATABASE IF EXISTS `{db_name}`')
            conn.execute(create_command)


def fill_database(verbose=False):
    """Fill a blank database with the ObsDB metadata."""
    # Create the schema from the base
    engine = get_engine(echo=verbose)
    Base.metadata.create_all(engine)

    # Create triggers
    for trigger in TRIGGERS:
        with engine.connect() as conn:
            conn.execute(trigger)


@contextmanager
def open_session(echo=params.DATABASE_ECHO):
    """Create a DB session context manager.

    Automatically takes care of commiting changes
    to DB when scope closes and rolls back on exceptions.

    Needless to say it also closes the session when it goes out of scope.

    Parameters
    ----------
    echo : bool
        If True, echo SQL output.
        Default is `gtecs.obs.params.DATABASE_ECHO`

    Returns
    -------
    session : `sqlalchemy.orm.session.Session`
        a database session context manager

    Examples
    --------
    To create a new user:

    >>> with open_session() as session:
    >>>     bob = User(username='bob', password='1234', full_name='Bob Marley')
    >>>     session.add(bob)
    >>> bob
    User(db_id=1, username=bob, full_name=Bob Marley)

    Note it was committed automatically when leaving the context.

    """
    engine = get_engine(echo=echo)
    new_session = sessionmaker(bind=engine)
    session = new_session()
    try:
        yield session
        session.commit()
    except Exception:
        session.rollback()
        raise
    finally:
        session.close()


def load_session(echo=params.DATABASE_ECHO):
    """Create a DB session.

    By making a DB session this way, you must commit and rollback changes
    yourself.

    Parameters
    ----------
    echo : bool
        If True, echo SQL output.
        Default is `gtecs.obs.params.DATABASE_ECHO`

    Returns
    -------
    session : `sqlalchemy.orm.session.Session`
        a database session

    Examples
    --------
    To create a new user:

    >>> session = load_session()
    >>> bob = User(username='bob', password='1234', full_name='Bob Marley')
    >>> try:
    >>>    session.add(new_user)
    >>> except:
    >>>    session.rollback()
    >>> finally:
    >>>    session.close()
    >>> bob
    User(db_id=1, username=bob, full_name=Bob Marley)

    """
    engine = get_engine(echo=echo)
    new_session = sessionmaker(bind=engine)
    session = new_session()
    return session
