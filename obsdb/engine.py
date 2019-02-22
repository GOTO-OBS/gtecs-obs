#!/usr/bin/env python
"""Session management functions."""

from contextlib import contextmanager

from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

from . import params

__all__ = ['open_session', 'load_session']


ENGINE = create_engine('mysql+pymysql://{}'.format(params.DATABASE_LOCATION),
                       echo=params.DATABASE_ECHO,
                       pool_pre_ping=params.DATABASE_PRE_PING)


@contextmanager
def open_session():
    """Create a DB session context manager.

    Automatically takes care of commiting changes
    to DB when scope closes and rolls back on exceptions.

    Needless to say it also closes the session when it goes out of scope.

    Returns
    -------
    session : `sqlalchemy.orm.session.Session`
        a database session context manager

    Examples
    --------
    To create a new user, and add to the database:

    >>> with open_session() as session:
    >>>     new_user = User('sl', '1234', 'Stuey')
    >>>     session.add(new_user)

    """
    new_session = sessionmaker(bind=ENGINE)
    session = new_session()
    try:
        yield session
        session.commit()
    except Exception:
        session.rollback()
        raise
    finally:
        session.close()


def load_session():
    """Create a DB session.

    By making a DB session this way, you must commit and rollback changes
    yourself.

    Returns
    -------
    session : `sqlalchemy.orm.session.Session`
        a database session

    Examples
    --------
    To create a new user, and add to the database:

    >>> session = load_session()
    >>> new_user = User('sl', '1234', 'Stuey')
    >>> try:
    >>>    session.add(new_user)
    >>> except:
    >>>    session.rollback()
    >>> finally:
    >>>    session.close()

    """
    new_session = sessionmaker(bind=ENGINE)
    session = new_session()
    return session
