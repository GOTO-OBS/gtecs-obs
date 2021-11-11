"""Session management functions."""

import os
from contextlib import contextmanager

import numpy as np

import pymysql

from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

from .. import params

__all__ = ['get_engine', 'open_session', 'load_session']


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
