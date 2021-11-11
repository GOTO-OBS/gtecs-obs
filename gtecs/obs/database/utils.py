"""Utility functions for using the database."""

import hashlib

from astropy import units as u
from astropy.coordinates import Longitude
from astropy.time import Time

from sqlalchemy import or_
from sqlalchemy.exc import ProgrammingError

from .engine import get_engine, open_session
from .models import Base, Event, ExposureSet, Grid, GridTile, Mpointing, Pointing, TRIGGERS, User
from .. import params


__all__ = ['create_database', 'fill_database',
           'get_user', 'validate_user',
           'get_filtered_queue', 'get_queue',
           'get_pointings', 'get_pointing_by_id',
           'get_mpointings', 'get_mpointing_by_id',
           'get_exposure_set_by_id',
           'get_current_grid', 'get_grid_tile_by_name',
           'get_events', 'delete_event_pointings',
           'get_expired_pointings', 'get_expired_mpointings',
           'insert_items',
           'update_pointing_status', 'bulk_update_status',
           'mark_completed', 'mark_aborted', 'mark_interrupted', 'mark_running',
           ]


def create_database(overwrite=False, verbose=False):
    """Create the blank goto_obs database.

    Parameters
    ----------
    overwrite : bool, default=False
        If True and a database named 'goto_obs' exists then drop it before creating the new one.
        If False and the database already exists then an error is raised.

    verbose : bool, default=False
        If True, echo SQL output.

    """
    engine = get_engine(db_name=None, echo=verbose)
    with engine.connect() as conn:
        if not overwrite:
            try:
                conn.execute('CREATE DATABASE `goto_obs`')  # will raise if it exists
            except ProgrammingError as err:
                raise ValueError('WARNING: goto_obs database already exists!') from err
        else:
            conn.execute('DROP DATABASE IF EXISTS `goto_obs`')
            conn.execute('CREATE DATABASE `goto_obs`')


def fill_database(verbose=False):
    """Fill a blank database with the ObsDB metadata."""
    # Create the schema from the base
    engine = get_engine(echo=verbose)
    Base.metadata.create_all(engine)

    # Create triggers
    for trigger in TRIGGERS:
        with engine.connect() as conn:
            conn.execute(trigger)


def get_user(session, username):
    """Return the `User` for a given username.

    Parameters
    ----------
    session : `sqlalchemy.Session.session`
        a session object - see `load_session` or `open_session`
    username : string
        short name of user

    Returns
    -------
    user : `User`
        the User class for the given username

    Raises
    ------
    ValueError : if no matching User is found in the database

    """
    user = session.query(User).filter(User.username == username).one_or_none()
    if not user:
        raise ValueError('No matching User found')
    return user


def validate_user(session, username, password):
    """Check user exists and password is correct.

    Parameters
    ----------
    session : `sqlalchemy.Session.session`
        a session object - see `load_session` or `open_session`
    username : string
        short name of user
    password : string
        plain text password (stored in DB using a sha512 hash for security)

    Returns
    -------
    passed : bool
        True if user exists and password is correct

    Raises
    ------
    ValueError : if username is not found in DB

    """
    user = get_user(session, username)
    return user.password_hash == hashlib.sha512(password.encode()).hexdigest()


def get_filtered_queue(session, time=None, rank_limit=None, location=None,
                       altitude_limit=20, hourangle_limit=6, limit_number=None):
    """Get the currently valid queue, and filter against visibility and rank.

    The queue is defined as all pending pointings with a valid date range. We then filter
    out tiles that are not visible, order by rank and limit by rank if requested.

    Note that rank limits are designed to not return low rank tiles, but *only* if higher
    ranked tiles exist. If all tiles are below the rank limit, all tiles are returned.

    Parameters
    ----------
    session : `sqlalchemy.Session.session`
        a session object - see `load_session` or `open_session` for details

    time : `~astropy.time.Time`
        If given, the time to fetch the queue for.
        Defaults to the current time.

    rank_limit: int or None, default None
        Only return pointings with rank greater than given value.
        If all tiles are below the rank limit, all tiles are returned.

    location : `~astropy.coordinates.EarthLocation`
        Location of observatory. If provided, only visible pointings are returned.

    altitude_limit : float
        If filtering by visibility, only pointings which ever rise above this altitude are returned

    hourangle_limit : float
        If filtering by visibility, only pointings with hour angles closer to transit than this
        limit are returned.

    limit_number : int or None
        If not None, limit the number of results.

    Returns
    -------
    current_pointing : `Pointing`
        Pointing being observed now
    pending_pointings : `Pointing`
        all current pending Pointings, after filtering

    Examples
    --------
    An example using `open_session`. The items retrieved won't function outside the
    scope of the with block, since they rely on a session to access the underlying
    database.

    >>> with open_session() as session:
    >>>     current_job, pending_jobs = get_filtered_queue(session, limit_results=100)
    >>>     njobs = len(pending_jobs)
    >>>     exposure_list = current_job.exposure_sets
    >>> current_job
    DetachedInstanceError: Instance <Pointing at 0x10a00ac50> is not bound to a Session;
    attribute refresh operation cannot proceed

    An example using `load_session`. This will return a session object. With this method
    you don't have to worry so about about scope, but you do have to be careful about
    trapping errors and commiting any changes you make to the database. See the help for
    `load_session` for more details.

    >>> session = load_session()
    >>> current_job, pending_jobs = get_filtered_queue(session, limit_results=100)

    """
    queue = session.query(Pointing)

    # only get pending pointings
    queue = queue.filter(Pointing.status == 'pending')

    # only get pointings with exposure sets, we've seen some odd ones without any
    queue = queue.filter(Pointing.exposure_sets)

    if time is None:
        time = Time.now()
    now = time.iso

    # are we after the start time?
    queue = queue.filter(Pointing.start_time < now)
    # are we before the stop time (if any)?
    queue = queue.filter(or_(Pointing.stop_time > now, Pointing.stop_time == None))  # noqa: E711

    # now limit by RA and Dec
    if location is not None:
        # local sidereal time, units of degrees
        lst = time.sidereal_time('mean', location.lon)
        lo_lim = Longitude(lst - hourangle_limit * u.hourangle).deg
        up_lim = Longitude(lst + hourangle_limit * u.hourangle).deg
        if up_lim > lo_lim:
            queue = queue.filter(Pointing.ra < up_lim,
                                 Pointing.ra > lo_lim)
        else:
            queue = queue.filter(or_(Pointing.ra < up_lim,
                                     Pointing.ra > lo_lim))

        # is latitude ever greater than limit?
        lat = location.lat.deg
        queue = queue.filter(Pointing.dec > lat - 90 + altitude_limit,
                             Pointing.dec < lat + 90 - altitude_limit)

    queue = queue.order_by('rank')

    if rank_limit:
        number_passing_rank = queue.filter(Pointing.rank < rank_limit).count()

        # if there are results which pass the rank limit, filter by rank
        if number_passing_rank > 0:
            queue = queue.filter(Pointing.rank < rank_limit)

    if limit_number:
        queue = queue.limit(limit_number)

    # actually get the queue
    pending_pointings = queue.all()

    # find the current pointing too, with a separate query
    current_pointing = session.query(Pointing).filter(Pointing.status == 'running').one_or_none()

    return current_pointing, pending_pointings


def get_queue(session, time=None):
    """Get the currently valid queue.

    The queue is defined as all pending pointings with a valid date range.

    Parameters
    ----------
    session : `sqlalchemy.Session.session`
        a session object - see `load_session` or `open_session` for details

    time : `~astropy.time.Time`
        If given, the time to fetch the queue for.
        Defaults to the current time.

    Returns
    -------
    current_pointing : `Pointing`
        Pointing being observed now
    pending_pointings : `Pointing`
        all current pending Pointings

    """
    # call get_filtered_queue with no filtering
    current_pointing, pending_queue = get_filtered_queue(session, time)
    return current_pointing, pending_queue


def get_pointings(session, pointing_ids=None, status=None):
    """Get pointings, filtered by ID or status.

    Parameters
    ----------
    session : `sqlalchemy.Session.session`
        a session object - see `load_session` or `open_session` for details
    pointing_ids : int or list
        supply a pointing ID or list of IDs to filter by id
    status : string or list
        supply a status or list of statuses to filter by status

    Returns
    -------
    pointings : list
        a list of all matching Pointings

    """
    query = session.query(Pointing)
    if pointing_ids is not None:
        query = query.filter(Pointing.db_id.in_(list(pointing_ids)))
    if status is not None:
        if isinstance(status, str):
            status = [status]
        query = query.filter(Pointing.status.in_(status))
    pointings = query.all()
    return pointings


def get_pointing_by_id(session, pointing_id):
    """Get a single Pointing, filtered by ID.

    Parameters
    ----------
    session : `sqlalchemy.Session.session`
        a session object - see `load_session` or `open_session` for details
    pointing_id : int
        the id number of the pointing

    Returns
    -------
    pointing : `Pointing`
        the matching Pointing

    Raises
    ------
    ValueError : if no matching Pointing is found in the database

    """
    pointing = session.query(Pointing).filter(Pointing.db_id == pointing_id).one_or_none()
    if not pointing:
        raise ValueError('No matching Pointing found')
    return pointing


def get_mpointings(session, mpointing_ids=None, status=None):
    """Get mpointings, filtered by ID or status.

    Parameters
    ----------
    session : `sqlalchemy.Session.session`
        a session object - see `load_session` or `open_session` for details
    mpointing_ids : int or list
        supply a ID or list of IDs to filter results
    status : string or list
        supply a status or list of statuses to filter by status

    Returns
    -------
    mpointings : list
        a list of all matching Mpointings

    """
    query = session.query(Mpointing)
    if mpointing_ids is not None:
        query = query.filter(Mpointing.db_id.in_(list(mpointing_ids)))
    if status is not None:
        if isinstance(status, str):
            status = [status]
        query = query.filter(Mpointing.status.in_(status))
    mpointings = query.all()
    return mpointings


def get_mpointing_by_id(session, mpointing_id):
    """Get a single Mpointing, filtered by ID.

    Parameters
    ----------
    session : `sqlalchemy.Session.session`
        a session object - see `load_session` or `open_session` for details
    mpointing_id : int
        the id number of the Mpointing

    Returns
    -------
    mpointing : `Mpointing`
        the matching Mpointing

    Raises
    ------
    ValueError : if no matching Mpointing is found in the database

    """
    mpointing = session.query(Mpointing).filter(Mpointing.db_id == mpointing_id).one_or_none()
    if not mpointing:
        raise ValueError('No matching Mpointing found')
    return mpointing


def get_exposure_set_by_id(session, expset_id):
    """Get a single Exposure Set, filtered by ID.

    Parameters
    ----------
    session : `sqlalchemy.Session.session`
        a session object - see `load_session` or `open_session` for details
    expset_id : int
        the id number of the Exposure Set

    Returns
    -------
    exposure_set : `ExposureSet`
        the matching Exposure Set

    Raises
    ------
    ValueError : if no matching Exposure Set is found in the database

    """
    exposure_set = session.query(ExposureSet).filter(ExposureSet.db_id == expset_id).one_or_none()
    if not exposure_set:
        raise ValueError('No matching Exposure Set found')
    return exposure_set


def get_current_grid(session):
    """Get the current (last-defined) Grid from the grids table."""
    # Get the final entry in the grids table, assuming that's the latest and therefore current one
    # Note there's no query.last(), so need to order by id decending then take the first.
    grid = session.query(Grid).order_by(Grid.db_id.desc()).first()
    if not grid:
        raise ValueError('No defined grids found')
    return grid


def get_grid_tile_by_name(session, grid_name, tile_name):
    """Get a tile in a grid from the name of the grid and tile.

    Parameters
    ----------
    session : `sqlalchemy.Session.session`
        a session object - see `load_session` or `open_session` for details
    grid_name : str
        the name of the grid
    tile_name : str
        the name of the tile

    Returns
    -------
    grid_tile : `GridTile`
        the matching Grid Tile

    Raises
    ------
    ValueError : if no matching Grid or Grid Tile is found in the database

    """
    grid = session.query(Grid).filter(Grid.name == grid_name).one_or_none()
    if not grid:
        raise ValueError('No matching GridTile found')
    grid_tile = session.query(GridTile).filter(GridTile.name == tile_name,
                                               GridTile.grid_id == grid.db_id,
                                               ).one_or_none()
    if not grid_tile:
        raise ValueError('No matching GridTile found')
    return grid_tile


def get_events(session, event_type=None, source=None):
    """Get all events, filtered by type or source.

    Parameters
    ----------
    session : `sqlalchemy.Session.session`
        a session object - see `load_session` or `open_session` for details
    event_type: str or list of str, optional
        the type of event, e.g. GW, GRB
    source : str or list of str,, optional
        the event's origin, e.g. LVC, Fermi, GAIA

    Returns
    -------
    events : list
        a list of all matching Events

    """
    query = session.query(Event)
    if event_type is not None:
        if isinstance(event_type, str):
            event_type = [event_type]
        query = query.filter(Event.event_type.in_(event_type))
    if source is not None:
        if isinstance(source, str):
            source = [source]
        query = query.filter(Event.source.in_(source))
    events = query.all()
    return events


def delete_event_pointings(session, event):
    """Set all the Mpointings and Pointings associated with the given Event to deleted.

    Note this doesn't physically delete the rows from the database tables,
    it just sets the statuses to 'deleted'.
    """
    for mpointing in event.mpointings:
        mpointing.status = 'deleted'
        if mpointing.pointings[-1].status == 'pending':
            mpointing.pointings[-1].status = 'deleted'


def get_expired_pointings(session, time=None):
    """Find all the pointings still pending whose valid period has expired.

    Parameters
    ----------
    session : `sqlalchemy.Session.session`
        a session object - see `load_session` or `open_session` for details
    time : `~astropy.time.Time`
        If given, the time to evaluate time period at.
        Defaults to the current time.

    Returns
    -------
    pointings : list
        a list of all matching Pointings

    """
    if time is None:
        now = Time.now().iso
    else:
        now = time.iso

    pointings = session.query(Pointing).filter(Pointing.status == 'pending',
                                               Pointing.stop_time < now).all()
    return pointings


def get_expired_mpointings(session, time=None):
    """Find all the mpointings still unscheduled whose valid period has expired.

    Parameters
    ----------
    session : `sqlalchemy.Session.session`
        a session object - see `load_session` or `open_session` for details
    time : `~astropy.time.Time`
        If given, the time to evaluate time period at.
        Defaults to the current time.

    Returns
    -------
    mpointings : list
        a list of all matching Mpointings

    """
    if time is None:
        now = Time.now().iso
    else:
        now = time.iso

    mpointings = session.query(Mpointing).filter(Mpointing.status == 'unscheduled',
                                                 Mpointing.stop_time < now).all()
    return mpointings


def insert_items(session, items):
    """Insert one or more items into the database.

    This routine is slow for large numbers of items,
    use one of the bulk_insert routines instead.

    Parameters
    ----------
    session : `sqlalchemy.Session.session`
        the session object
    items : `gtecs.obs.database.models.Base`
        A model instance, e.g `Pointing`

    """
    items = list(items)
    for i, item in enumerate(items):
        session.add(item)
        if i % 1000 == 0:
            session.flush()


def update_pointing_status(session, pointing_id, status):
    """Set the status of the given pointing.

    For large numbers of pointings use `bulk_update_pointing_status`.

    Parameters
    ----------
    session : `sqlalchemy.Session.session`
        the session object
    pointing_id : int
        the id number of the pointing
    status : string
        status to set pointing to

    """
    pointing = get_pointing_by_id(session, pointing_id)
    pointing.status = status


def bulk_update_status(session, items, status):
    """Set the status of a large number of Pointings or Mpointings.

    For large numbers of items this routine is most efficient.

    Parameters
    ----------
    session : `sqlalchemy.Session.session`
        the session object
    items : list
        a list of `Pointing`s or `Mpointing`s to update
    status : string
        status to set pointings to

    """
    # Make sure they're all the same type
    if not all(isinstance(item, type(items[0])) for item in items):
        raise ValueError('Items must be all the same type (`Pointing` or `Mpointing`)')

    mappings = [{'db_id': item.db_id, 'status': status} for item in items]
    session.bulk_update_mappings(type(items[0]), mappings)


def mark_completed(pointing_id):
    """Update the given pointing's status to 'completed'.

    A utility function that creates its own session.

    Parameters
    ----------
    pointing_id : int
        the id number of the pointing

    """
    with open_session() as s:
        update_pointing_status(s, pointing_id, 'completed')


def mark_aborted(pointing_id):
    """Update the given pointing's status to 'aborted'.

    A utility function that creates its own session.

    Parameters
    ----------
    pointing_id : int
        the id number of the pointing

    """
    with open_session() as s:
        update_pointing_status(s, pointing_id, 'aborted')


def mark_interrupted(pointing_id):
    """Update the given pointing's status to 'interrupted'.

    A utility function that creates its own session.

    Parameters
    ----------
    pointing_id : int
        the id number of the pointing

    """
    with open_session() as s:
        update_pointing_status(s, pointing_id, 'interrupted')


def mark_running(pointing_id):
    """Update the given pointing's status to 'running'.

    A utility function that creates its own session.

    Parameters
    ----------
    pointing_id : int
        the id number of the pointing

    """
    with open_session() as s:
        update_pointing_status(s, pointing_id, 'running')
