from contextlib import contextmanager

from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

__all__ = ['open_session', 'load_session']

@contextmanager
def open_session(host='localhost', echo=False):
    """
    Create a DB session context manager.

    Automatically takes care of commiting changes
    to DB when scope closes and rolls back on exceptions.

    Needless to say it also closes the session when it goes out of scope.

    Parameters
    ----------
    host : str, optional
        the host location of the database
        default is 'localhost'
    echo : bool, optional
        if true, will log all statements to the engines logger,
        which defaults to stdout

    Returns
    -------
    session : `sqlalchemy.orm.session.Session`
        a database session context manager

    Examples
    --------
    To create a new user, and add to the database:

    >>> with open_session() as session:
    >>>     newUser = User('sl', '1234', 'Stuey')
    >>>     session.add(newUser)
    """

    location_string = 'goto:gotoobs@{}/goto_obs'.format(host)
    engine_string = 'mysql+pymysql://{}'.format(location_string)
    engine = create_engine(engine_string, echo=echo, pool_pre_ping=True)

    Session = sessionmaker(bind=engine)
    session = Session()
    try:
        yield session
        session.commit()
    except:
        session.rollback()
        raise
    finally:
        session.close()


def load_session(host='localhost', echo=False):
    """
    Create a DB session.

    By making a DB session this way, you must commit and rollback changes
    yourself.

    Parameters
    ----------
    host : str, optional
        the host location of the database
        default is 'localhost'
    echo : bool, optional
        if true, will log all statements to the engines logger,
        which defaults to stdout

    Returns
    -------
    session : `sqlalchemy.orm.session.Session`
        a database session

    Examples
    --------
    To create a new user, and add to the database:

    >>> session = load_session()
    >>> newUser = User('sl', '1234', 'Stuey')
    >>> try:
    >>>    session.add(newUser)
    >>> except:
    >>>    session.rollback()
    >>> finally:
    >>>    session.close()
    """

    location_string = 'goto:gotoobs@{}/goto_obs'.format(host)
    engine_string = 'mysql+pymysql://{}'.format(location_string)
    engine = create_engine(engine_string, echo=echo, pool_pre_ping=True)

    Session = sessionmaker(bind=engine)
    session = Session()
    return session
