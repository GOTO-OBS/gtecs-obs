"""Utility functions for using the database.

Functions that return database ORM model classes (e.g. `get_user` -> `User`,
`get_pointing_by_id` -> `Pointing` etc) require a session to be passed as their first argument.
The returned classes are then linked to this session, so any changes will have to be committed
manually or when the session is closed.

Functions that either act on the database and return nothing (e.g. `mark_pointing_running`)
or that query the database and return information instead of classes (e.g. `validate_user` -> bool,
`get_pointing_info` -> dict) do not require a session (they create their own internally).
"""

import hashlib
import json
from contextlib import contextmanager

from astropy import units as u
from astropy.coordinates import Longitude
from astropy.time import Time

from gtecs.common.database import get_session as get_session_common

from sqlalchemy import or_

from .models import ExposureSet, Grid, GridTile, Pointing, Site, Target, Telescope, User
from .. import params


__all__ = ['get_session', 'session_manager',
           'get_user', 'validate_user',
           'get_pointings', 'get_pointing_by_id', 'get_pending_pointings', 'get_current_pointing',
           'mark_pointing_running', 'mark_pointing_completed', 'mark_pointing_interrupted',
           'mark_pointing_confirmed', 'mark_pointing_failed',
           'get_pointing_info',
           'get_targets', 'get_target_by_id',
           'get_exposure_set_by_id',
           'get_telescope_by_id', 'get_telescope_info', 'get_site_by_id', 'get_site_info',
           'get_current_grid', 'get_grid_tile_by_name',
           'insert_items',
           ]


def get_session(user=None, password=None, host=None, echo=None, pool_pre_ping=None):
    """Create a database connection session.

    All arguments are passed to `gtecs.common.database.get_session()`,
    with the defaults taken from the module parameters.

    Note it is generally better to use the session_manager() context manager,
    which will automatically commit or rollback changes when done.
    """
    # This means the user doesn't need to worry about the params, but can overwrite if needed.
    if user is None:
        user = params.DATABASE_USER
    if password is None:
        password = params.DATABASE_PASSWORD
    if host is None:
        host = params.DATABASE_HOST
    if echo is None:
        echo = params.DATABASE_ECHO
    if pool_pre_ping is None:
        pool_pre_ping = params.DATABASE_PRE_PING
    session = get_session_common(
        user=user,
        password=password,
        host=host,
        echo=echo,
        pool_pre_ping=pool_pre_ping,
    )
    return session


@contextmanager
def session_manager(**kwargs):
    """Create a session context manager connection to the database."""
    session = get_session(**kwargs)
    try:
        yield session
        session.commit()
    except Exception:
        session.rollback()
        raise
    finally:
        session.close()


def get_user(session, username):
    """Return the `User` for a given username.

    Parameters
    ----------
    session : `sqlalchemy.Session.session`
        a session object
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


def validate_user(username, password):
    """Check that the user with the given name exists and their password is correct.

    Parameters
    ----------
    username : string
        short name of user
    password : string
        plain text password (stored in DB using a sha512 hash for security)

    Returns
    -------
    passed : bool
        True if the user exists and the password is correct

    """
    with session_manager() as session:
        user = get_user(session, username)
        return user.password_hash == hashlib.sha512(password.encode()).hexdigest()


def get_pointings(session, pointing_ids=None, status=None, time=None):
    """Get pointings, filtered by ID or status.

    Parameters
    ----------
    session : `sqlalchemy.Session.session`
        a session object
    pointing_ids : int or list
        supply a pointing ID or list of IDs to filter by id
    status : string or list
        supply a status or list of statuses to filter by status

    time : `~astropy.time.Time` or None
        If given, the time to check the status at.
        Defaults to Time.now().

    Returns
    -------
    pointings : list of `Pointing`
        a list of all matching Pointings

    """
    if time is None:
        time = Time.now()
    query = session.query(Pointing)
    if pointing_ids is not None:
        query = query.filter(Pointing.db_id.in_(list(pointing_ids)))
    if status is not None:
        if isinstance(status, str):
            status = [status]
        query = query.filter(Pointing.status_at_time(time).in_(status))
    pointings = query.all()
    return pointings


def get_pointing_by_id(session, pointing_id):
    """Get a single Pointing, filtered by ID.

    Parameters
    ----------
    session : `sqlalchemy.Session.session`
        a session object
    pointing_id : int
        the id number of the pointing

    Returns
    -------
    pointing : `Pointing`
        the matching Pointing

    """
    pointing = session.query(Pointing).filter(Pointing.db_id == pointing_id).one_or_none()
    if not pointing:
        raise ValueError('No matching Pointing found')
    return pointing


def get_pending_pointings(session, time=None, location=None, altitude_limit=20, hourangle_limit=6,
                          rank_limit=None, limit_number=None):
    """Get the queue of pending Pointings for the given telescope.

    The queue is defined as all pending pointings with a valid date range. We then filter
    out pointings that are not visible from the telescope site, as well as ordering and limiting
    by rank if requested.

    Note that rank limits are designed to not return low rank pointings, but *only* if higher
    ranked pointings exist. If all pointings are below the rank limit, all pointings are returned.

    Parameters
    ----------
    session : `sqlalchemy.Session.session`
        a session object

    time : `~astropy.time.Time`
        If given, the time to fetch the queue at.
        Defaults to Time.now().
    location : `~astropy.coordinates.EarthLocation`
        If given, the site location to fetch the queue for.
        Defaults to None, meaning the queue will not be filtered for visibility.
    altitude_limit : float
        If a location is given, only pointings which ever rise above this altitude are returned.
        Defaults to 20 degrees.
    hourangle_limit : float
        If a location is given, only pointings with hour angles closer to transit than this
        limit are returned.
        Defaults to 6 hours.
    rank_limit: int
        If given, only return pointings with rank greater than given value.
        If all pointings are below the rank limit, all pointings are returned.
        Defaults to None, meaning no rank filtering.
    limit_number : int
        If given, limit the number of results.
        Defaults to None, meaning no limit of the number of pointings to return.

    Returns
    -------
    pointings : list of `Pointing`
        all pending Pointings, after filtering

    """
    if time is None:
        time = Time.now()

    query = session.query(Pointing)

    # only get pending pointings (will account for start/stop times)
    query = query.filter(Pointing.status_at_time(time) == 'pending')

    if location is not None:
        # limit by hour angle
        if hourangle_limit is not None:
            lst = time.sidereal_time('mean', location.lon)
            lo_lim = Longitude(lst - hourangle_limit * u.hourangle).deg
            up_lim = Longitude(lst + hourangle_limit * u.hourangle).deg
            if up_lim > lo_lim:
                query = query.filter(Pointing.ra < up_lim,
                                     Pointing.ra > lo_lim)
            else:
                query = query.filter(or_(Pointing.ra < up_lim,
                                         Pointing.ra > lo_lim))

        # limit by altitude
        if altitude_limit is not None:
            lat = location.lat.deg
            query = query.filter(Pointing.dec > lat - 90 + altitude_limit,
                                 Pointing.dec < lat + 90 - altitude_limit)

    query = query.order_by('rank')
    if rank_limit:
        # if there are results which pass the rank limit, filter by rank
        number_passing_rank = query.filter(Pointing.rank < rank_limit).count()
        if number_passing_rank > 0:
            query = query.filter(Pointing.rank < rank_limit)

    if limit_number:
        query = query.limit(limit_number)

    # make the query and return
    pointings = query.all()
    return pointings


def get_current_pointing(session, telescope_id, time=None):
    """Fetch the current Pointing from the database for the given telescope.

    Parameters
    ----------
    session : `sqlalchemy.Session.session`
        a session object
    telescope_id : int
        The telescope to fetch the pointing for.

    time : `~astropy.time.Time` or None
        If given, the time to check the status at.
        Defaults to Time.now().

    Returns
    -------
    pointings : `Pointing` or None
        a single Pointing, or None if the telescope is idle.

    """
    if time is None:
        time = Time.now()
    query = session.query(Pointing)
    query = query.filter(Pointing.status_at_time(time) == 'running')
    query = query.filter(Pointing.telescope_id == telescope_id)
    pointing = query.one_or_none()  # Shouldn't ever be more than one running per telescope
    return pointing


def mark_pointing_running(pointing_id, telescope_id, time=None):
    """Update the given pointing's status to 'running'.

    Parameters
    ----------
    pointing_id : int
        the id number of the pointing
    telescope_id : int
        the id number of the telescope the pointing is running on

    time : str, `datetime.datetime`, `astropy.time.Time`, optional
        If given, the time to mark the Pointing as interrupted at.

    """
    with session_manager() as session:
        pointing = get_pointing_by_id(session, pointing_id)
        telescope = get_telescope_by_id(session, telescope_id)
        pointing.mark_running(telescope, time=time)
        session.commit()


def mark_pointing_completed(pointing_id, schedule_next=True, time=None):
    """Update the given pointing's status to 'completed'.

    Parameters
    ----------
    pointing_id : int
        the id number of the pointing

    schedule_next : bool, optional
        If True, automatically create the next Pointing for the Target
        (unless it has no more to do).
        default = True
    time : str, `datetime.datetime`, `astropy.time.Time`, optional
        If given, the time to mark the Pointing as completed at.

    """
    with session_manager() as session:
        pointing = get_pointing_by_id(session, pointing_id)
        pointing.mark_finished(completed=True, time=time)
        session.commit()

        if schedule_next:
            next_pointing = pointing.target.get_next_pointing(time=time)
            if next_pointing is not None:
                session.add(next_pointing)
                session.commit()


def mark_pointing_interrupted(pointing_id, schedule_next=True, delay=None, time=None):
    """Update the given pointing's status to 'interrupted'.

    Parameters
    ----------
    pointing_id : int
        the id number of the pointing

    schedule_next : bool, default=True
        If True, automatically create the next Pointing for the Target
        (unless it has no more to do).
    delay : float, default=None
        Time (in seconds) to delay the validity of the next Pointing from the time
        the previous Pointing (i.e. the one we're marking) finished.
        If None, the new Pointing takes the same start_time as the interrupted one.
        Only relevant if schedule_next=True.
    time : str, `datetime.datetime`, `astropy.time.Time`, optional
        If given, the time to mark the Pointing as interrupted at.

    """
    with session_manager() as session:
        pointing = get_pointing_by_id(session, pointing_id)
        pointing.mark_finished(completed=False, time=time)
        session.commit()

        if schedule_next:
            next_pointing = pointing.target.get_next_pointing(time=time, reschedule_delay=delay)
            if next_pointing is not None:
                session.add(next_pointing)
                session.commit()


def mark_pointing_confirmed(pointing_id, time=None):
    """Validate the given pointing's data as confirmed to be good quality.

    Parameters
    ----------
    pointing_id : int
        the id number of the pointing

    time : str, `datetime.datetime`, `astropy.time.Time`, optional
        If given, the time to mark the Pointing as interrupted at.

    """
    with session_manager() as session:
        pointing = get_pointing_by_id(session, pointing_id)
        pointing.mark_validated(good=True, time=time)
        session.commit()


def mark_pointing_failed(pointing_id, schedule_next=True, delay=None, time=None):
    """Validate the given pointing's data quality as poor, and needs re-observing.

    Parameters
    ----------
    pointing_id : int
        the id number of the pointing

    schedule_next : bool, default=True
        If True, automatically create the next Pointing for the Target
        (unless it has no more to do).
    delay : float, default=None
        Time (in seconds) to delay the validity of the next Pointing from the time
        the previous Pointing (i.e. the one we're marking) finished.
        If None, the new Pointing takes the same start_time as the failed one.
        Only relevant if schedule_next=True.
    time : str, `datetime.datetime`, `astropy.time.Time`, optional
        If given, the time to mark the Pointing as failed at.

    """
    with session_manager() as session:
        pointing = get_pointing_by_id(session, pointing_id)
        pointing.mark_validated(good=False, time=time)
        session.commit()

        if schedule_next:
            # There might be an upcoming pointing already, we'll need to delete it first
            # Although we don't want to interrupt if it's running
            for p in pointing.target.pointings:
                if p.status_at_time(time) in ['upcoming', 'pending']:
                    p.mark_deleted(time)
            session.commit()

            next_pointing = pointing.target.get_next_pointing(time=time, reschedule_delay=delay)
            if next_pointing is not None:
                session.add(next_pointing)
                session.commit()


def get_pointing_info(pointing_id):
    """Get a dictionary of info for the given Pointing.

    This should contain all the information needed for an image FITS header.

    Note the dict needs to be JSON-serializable, so we can't include any datetimes or other classes.
    """
    pointing_info = {}

    with session_manager() as session:
        # Get the Pointing from the database
        pointing = get_pointing_by_id(session, pointing_id)

        # Get Pointing info
        pointing_info['id'] = pointing.db_id
        # pointing_info['status'] = pointing.status  # Don't include anything time-dependent
        pointing_info['rank'] = pointing.rank
        pointing_info['start_time'] = pointing.start_time.isoformat(sep=' ')
        t = pointing.stop_time.isoformat(sep=' ') if pointing.stop_time is not None else None
        pointing_info['stop_time'] = t

        # Get Target info
        target = pointing.target
        pointing_info['target_id'] = target.db_id
        pointing_info['name'] = target.name
        pointing_info['ra'] = target.ra
        pointing_info['dec'] = target.dec
        pointing_info['start_rank'] = target.rank
        pointing_info['weight'] = target.weight
        pointing_info['is_template'] = target.is_template
        pointing_info['num_completed'] = target.num_completed
        pointing_info['target_start_time'] = target.start_time.isoformat(sep=' ')
        t = target.stop_time.isoformat(sep=' ') if target.stop_time is not None else None
        pointing_info['target_stop_time'] = t

        # Get Strategy info
        strategy = pointing.strategy
        pointing_info['strategy_id'] = strategy.db_id
        pointing_info['infinite'] = strategy.infinite
        pointing_info['min_time'] = strategy.min_time
        pointing_info['too'] = strategy.too
        pointing_info['requires_template'] = strategy.requires_template
        pointing_info['min_alt'] = strategy.min_alt
        pointing_info['max_sunalt'] = strategy.max_sunalt
        pointing_info['max_moon'] = strategy.max_moon
        pointing_info['min_moonsep'] = strategy.min_moonsep
        pointing_info['tel_mask'] = strategy.tel_mask

        # Get TimeBlock info
        time_block = pointing.time_block
        pointing_info['time_block_id'] = time_block.db_id
        pointing_info['block_num'] = time_block.block_num
        pointing_info['wait_time'] = time_block.wait_time
        pointing_info['valid_time'] = time_block.valid_time
        pointing_info['rank_change'] = time_block.rank_change

        # Get ExposureSet info
        pointing_info['exposure_sets'] = []
        for exposure_set in pointing.exposure_sets:
            expset_info = {}
            expset_info['id'] = exposure_set.db_id
            expset_info['num_exp'] = exposure_set.num_exp
            expset_info['exptime'] = exposure_set.exptime
            expset_info['filt'] = exposure_set.filt
            expset_info['binning'] = exposure_set.binning
            expset_info['dithering'] = exposure_set.dithering
            expset_info['ut_mask'] = exposure_set.ut_mask
            pointing_info['exposure_sets'].append(expset_info)

        # Get User info
        user = target.user
        pointing_info['user_id'] = user.db_id
        pointing_info['user_name'] = user.username
        pointing_info['user_fullname'] = user.full_name

        # Get Grid info
        grid_tile = pointing.grid_tile
        if grid_tile is not None:
            pointing_info['grid_id'] = grid_tile.grid.db_id
            pointing_info['grid_name'] = grid_tile.grid.name
            pointing_info['tile_id'] = grid_tile.db_id
            pointing_info['tile_name'] = grid_tile.name
        else:
            pointing_info['grid_id'] = None
            pointing_info['grid_name'] = None
            pointing_info['tile_id'] = None
            pointing_info['tile_name'] = None

        # Get Survey info
        survey = pointing.survey
        if survey is not None:
            pointing_info['survey_id'] = survey.db_id
            pointing_info['survey_name'] = survey.name
        else:
            pointing_info['survey_id'] = None
            pointing_info['survey_name'] = None

        # Try to get details from the alertdb
        # Because we import it in __init__, the backrefs should work
        try:
            # Get Notice info
            notice = pointing.target.notice
            pointing_info['notice_id'] = notice.db_id
            pointing_info['notice_ivorn'] = notice.ivorn
            pointing_info['notice_time'] = notice.received.isoformat(sep=' ')
            # Get Event info
            event = notice.event
            pointing_info['event_id'] = event.db_id
            pointing_info['event_name'] = event.name
            pointing_info['event_type'] = event.type
            pointing_info['event_origin'] = event.origin
            pointing_info['event_time'] = event.time.isoformat(sep=' ')
        except Exception:
            pointing_info['notice_id'] = None
            pointing_info['notice_ivorn'] = None
            pointing_info['notice_time'] = None
            pointing_info['event_id'] = None
            pointing_info['event_name'] = None
            pointing_info['event_type'] = None
            pointing_info['event_origin'] = None
            pointing_info['event_time'] = None

    # Check the dict is serializable, so we can send it though the web server API
    if json.loads(json.dumps(pointing_info)) != pointing_info:
        raise ValueError('Pointing info is not JSON-serializable')

    return pointing_info


def get_targets(session, target_ids=None, status=None, time=None):
    """Get targets, filtered by ID or status.

    Parameters
    ----------
    session : `sqlalchemy.Session.session`
        a session object
    target_ids : int or list
        supply a ID or list of IDs to filter results
    status : string or list
        supply a status or list of statuses to filter by status

    time : `~astropy.time.Time` or None
        If given, the time to check the status at.
        Defaults to Time.now().

    Returns
    -------
    targets : list of `Target`
        a list of all matching Targets

    """
    if time is None:
        time = Time.now()
    query = session.query(Target)
    if target_ids is not None:
        query = query.filter(Target.db_id.in_(list(target_ids)))
    if status is not None:
        if isinstance(status, str):
            status = [status]
        query = query.filter(Target.status_at_time(time).in_(status))
    targets = query.all()
    return targets


def get_target_by_id(session, target_id):
    """Get a single Target, filtered by ID.

    Parameters
    ----------
    session : `sqlalchemy.Session.session`
        a session object
    target_id : int
        the id number of the Target

    Returns
    -------
    target : `Target`
        the matching Target

    """
    target = session.query(Target).filter(Target.db_id == target_id).one_or_none()
    if not target:
        raise ValueError('No matching Target found')
    return target


def get_exposure_set_by_id(session, expset_id):
    """Get a single Exposure Set, filtered by ID.

    Parameters
    ----------
    session : `sqlalchemy.Session.session`
        a session object
    expset_id : int
        the id number of the Exposure Set

    Returns
    -------
    exposure_set : `ExposureSet`
        the matching Exposure Set

    """
    exposure_set = session.query(ExposureSet).filter(ExposureSet.db_id == expset_id).one_or_none()
    if not exposure_set:
        raise ValueError('No matching Exposure Set found')
    return exposure_set


def get_telescope_by_id(session, telescope_id):
    """Get a single Telescope, filtered by ID.

    Parameters
    ----------
    session : `sqlalchemy.Session.session`
        a session object
    telescope_id : int
        the id number of the telescope

    Returns
    -------
    telescope : `Telescope`
        the matching Telescope

    """
    telescope = session.query(Telescope).filter(Telescope.db_id == telescope_id).one_or_none()
    if not telescope:
        raise ValueError('No matching Telescope found')
    return telescope


def get_telescope_info():
    """Get the key info on all the telescopes defined in the database."""
    data = {}
    with session_manager() as session:
        telescopes = session.query(Telescope).all()
        for telescope in telescopes:
            tel_data = {}
            tel_data['name'] = telescope.name
            tel_data['status'] = telescope.status
            if telescope.current_pointing is not None:
                tel_data['current_pointing'] = telescope.current_pointing.db_id
            else:
                tel_data['current_pointing'] = None
            tel_data['horizon'] = telescope.get_horizon()
            tel_data['location'] = telescope.site.location
            data[telescope.db_id] = tel_data
    return data


def get_site_by_id(session, site_id):
    """Get a single Site, filtered by ID.

    Parameters
    ----------
    session : `sqlalchemy.Session.session`
        a session object
    site_id : int
        the id number of the site

    Returns
    -------
    site : `Site`
        the matching site

    """
    site = session.query(Site).filter(Site.db_id == site_id).one_or_none()
    if not site:
        raise ValueError('No matching Site found')
    return site


def get_site_info(site_id=None):
    """Get the key info on one or all of the sites and telescopes defined in the database."""
    data = {}
    with session_manager() as session:
        sites = session.query(Site).all()
        for site in sites:
            site_data = {}
            site_data['name'] = site.name
            site_data['location'] = site.location
            site_data['telescopes'] = {}
            for telescope in site.telescopes:
                tel_data = {}
                tel_data['name'] = telescope.name
                tel_data['status'] = telescope.status
                if telescope.current_pointing is not None:
                    tel_data['current_pointing'] = telescope.current_pointing.db_id
                else:
                    tel_data['current_pointing'] = None
                tel_data['horizon'] = telescope.get_horizon()
                site_data['telescopes'][telescope.db_id] = tel_data
            data[site.db_id] = site_data
    if site_id is not None:
        return data[site_id]
    return data


def get_current_grid(session):
    """Get the current (last-defined) Grid from the grids table."""
    # Get the final entry in the grids table, assuming that's the latest and therefore current one
    # Note there's no query.last(), so need to order by id descending then take the first.
    grid = session.query(Grid).order_by(Grid.db_id.desc()).first()
    if not grid:
        raise ValueError('No defined grids found')
    return grid


def get_grid_tile_by_name(session, grid_name, tile_name):
    """Get a tile in a grid from the name of the grid and tile.

    Parameters
    ----------
    session : `sqlalchemy.Session.session`
        a session object
    grid_name : str
        the name of the grid
    tile_name : str
        the name of the tile

    Returns
    -------
    grid_tile : `GridTile`
        the matching Grid Tile

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
