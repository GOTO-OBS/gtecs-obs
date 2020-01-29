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
    To create a new user:

    >>> with open_session() as session:
    >>>     bob = User(username='bob', password='1234', full_name='Bob Marley')
    >>>     session.add(bob)
    >>> bob
    User(db_id=1, username=bob, full_name=Bob Marley)

    Note it was committed automatically when leaving the context.

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
    new_session = sessionmaker(bind=ENGINE)
    session = new_session()
    return session
