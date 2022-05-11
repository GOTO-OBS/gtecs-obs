"""Utility functions for using the database.

Functions that return database ORM model classes (e.g. `get_user` -> `User`,
`get_pointing_by_id` -> `Pointing` etc) require a session to be passed as their first argument.
The returned classes are then linked to this session, so any changes will have to be committed
manually (if using `load_session`) or when the session is closed (if using `open_session`).

Functions that either act on the database and return nothing (e.g. `mark_pointing_running`)
or that query the database and return infomation instead of classes (e.g. `validate_user` -> bool,
`get_pointing_info` -> dict) do not require a session (they create their own internally).
"""

import hashlib

from astropy import units as u
from astropy.coordinates import Longitude
from astropy.time import Time

from sqlalchemy import or_

from .management import open_session
from .models import Event, ExposureSet, Grid, GridTile, Pointing, Target, Telescope, User


__all__ = ['get_user', 'validate_user',
           'get_pointings', 'get_pointing_by_id', 'get_pending_pointings', 'get_current_pointing',
           'mark_pointing_running', 'mark_pointing_completed', 'mark_pointing_interrupted',
           'mark_pointing_confirmed', 'mark_pointing_failed',
           'get_pointing_info',
           'get_targets', 'get_target_by_id',
           'get_exposure_set_by_id',
           'get_telescope_by_id', 'get_telescope_info',
           'get_current_grid', 'get_grid_tile_by_name',
           'get_events',
           'insert_items',
           ]


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
    with open_session() as session:
        user = get_user(session, username)
        return user.password_hash == hashlib.sha512(password.encode()).hexdigest()


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
    pointings : list of `Pointing`
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

    """
    pointing = session.query(Pointing).filter(Pointing.db_id == pointing_id).one_or_none()
    if not pointing:
        raise ValueError('No matching Pointing found')
    return pointing


def get_pending_pointings(session, telescope_id, time=None, rank_limit=None,
                          altitude_limit=20, hourangle_limit=6, limit_number=None):
    """Get the queue of pending Pointings for the given telescope.

    The queue is defined as all pending pointings with a valid date range. We then filter
    out pointings that are not visible from the telescope site, as well as ordering and limiting
    by rank if requested.

    Note that rank limits are designed to not return low rank pointings, but *only* if higher
    ranked pointings exist. If all pointings are below the rank limit, all pointings are returned.

    Parameters
    ----------
    session : `sqlalchemy.Session.session`
        a session object - see `load_session` or `open_session` for details
    telescope_id : int
        The telescope to fetch the queue for.

    time : `~astropy.time.Time`
        If given, the time to fetch the queue at.
        Defaults to Time.now().
    rank_limit: int or None, default None
        Only return pointings with rank greater than given value.
        If all pointings are below the rank limit, all pointings are returned.
    altitude_limit : float
        If filtering by visibility, only pointings which ever rise above this altitude are returned
    hourangle_limit : float
        If filtering by visibility, only pointings with hour angles closer to transit than this
        limit are returned.
    limit_number : int or None
        If not None, limit the number of results.

    Returns
    -------
    pointings : list of `Pointing`
        all pending Pointings, after filtering

    Examples
    --------
    An example using `open_session`. The items retrieved won't function outside the
    scope of the with block, since they rely on a session to access the underlying
    database.

    >>> with open_session() as session:
    >>>     pointings = get_pending_pointings(session, telescope_id=1, limit_results=100)
    >>>     n_pointings = len(pointings)
    >>>     exposure_list = pointings[0].exposure_sets
    >>> pointings[0]
    DetachedInstanceError: Instance <Pointing at 0x10a00ac50> is not bound to a Session;
    attribute refresh operation cannot proceed

    An example using `load_session`. This will return a session object. With this method
    you don't have to worry so about about scope, but you do have to be careful about
    trapping errors and commiting any changes you make to the database. See the help for
    `load_session` for more details.

    >>> session = load_session()
    >>> pointings = get_pending_pointings(session, telescope_id=1, limit_results=100)

    """
    if time is None:
        time = Time.now()

    query = session.query(Pointing)

    # only get pending pointings (will account for start/stop times)
    query = query.filter(Pointing.status_at_time(time) == 'pending')

    # limit by telescope ID
    query = query.filter(Pointing.tel_is_valid(telescope_id))

    # get the telescope site location
    if altitude_limit is not None or hourangle_limit is not None:
        telescope = get_telescope_by_id(session, telescope_id)
        location = telescope.site.location

    # limit by RA and Dec
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

    # is latitude ever greater than limit?
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
        a session object - see `load_session` or `open_session` for details
    telescope_id : int
        The telescope to fetch the pointing for.

    time : `~astropy.time.Time` or None
        If given, the time to check the status at (useful for simulations).
        Defaults to Time.now().

    Returns
    -------
    pointings : `Pointing` or None
        a single Pointing, or None if the telescope is idle.

    """
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
    with open_session() as session:
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
    with open_session() as session:
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
    with open_session() as session:
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
    with open_session() as session:
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
    with open_session() as session:
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

    This should contain all the infomation needed for an image FITS header.
    """
    pointing_info = {}

    with open_session() as session:
        # Get the Pointing from the database
        pointing = get_pointing_by_id(session, pointing_id)

        # Get Pointing info
        pointing_info['id'] = pointing.db_id
        # pointing_info['status'] = pointing.status  # Don't include anything time-dependent
        pointing_info['rank'] = pointing.rank
        pointing_info['start_time'] = pointing.start_time
        pointing_info['stop_time'] = pointing.stop_time

        # Get Target info
        target = pointing.target
        pointing_info['target_id'] = target.db_id
        pointing_info['name'] = target.name
        pointing_info['ra'] = target.ra
        pointing_info['dec'] = target.dec
        pointing_info['start_rank'] = target.rank
        pointing_info['rank_decay'] = target.rank_decay
        pointing_info['weight'] = target.weight
        pointing_info['is_template'] = target.is_template
        pointing_info['num_completed'] = target.num_completed
        pointing_info['target_start_time'] = target.start_time
        pointing_info['target_stop_time'] = target.stop_time

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

        # Get ExposureSet info
        pointing_info['exposure_sets'] = []
        for exposure_set in pointing.exposure_sets:
            expset_info = {}
            expset_info['id'] = exposure_set.db_id
            expset_info['num_exp'] = exposure_set.num_exp
            expset_info['exptime'] = exposure_set.exptime
            expset_info['filt'] = exposure_set.filt
            expset_info['binning'] = exposure_set.binning
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
            pointing_info['skymap'] = survey.skymap
        else:
            pointing_info['survey_id'] = None
            pointing_info['survey_name'] = None
            pointing_info['skymap'] = None

        # Get Event info
        event = pointing.event
        if event is not None:
            pointing_info['event_id'] = event.db_id
            pointing_info['event_name'] = event.name
            pointing_info['event_source'] = event.source
            pointing_info['event_type'] = event.type
            pointing_info['event_time'] = event.time
        else:
            pointing_info['event_id'] = None
            pointing_info['event_name'] = None
            pointing_info['event_source'] = None
            pointing_info['event_type'] = None
            pointing_info['event_time'] = None

    return pointing_info


def get_targets(session, target_ids=None, status=None):
    """Get targets, filtered by ID or status.

    Parameters
    ----------
    session : `sqlalchemy.Session.session`
        a session object - see `load_session` or `open_session` for details
    target_ids : int or list
        supply a ID or list of IDs to filter results
    status : string or list
        supply a status or list of statuses to filter by status

    Returns
    -------
    targets : list of `Target`
        a list of all matching Targets

    """
    query = session.query(Target)
    if target_ids is not None:
        query = query.filter(Target.db_id.in_(list(target_ids)))
    if status is not None:
        if isinstance(status, str):
            status = [status]
        query = query.filter(Target.status.in_(status))
    targets = query.all()
    return targets


def get_target_by_id(session, target_id):
    """Get a single Target, filtered by ID.

    Parameters
    ----------
    session : `sqlalchemy.Session.session`
        a session object - see `load_session` or `open_session` for details
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
        a session object - see `load_session` or `open_session` for details
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
        a session object - see `load_session` or `open_session` for details
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
    with open_session() as session:
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
    events : list of `Event`
        a list of all matching Events

    """
    query = session.query(Event)
    if event_type is not None:
        if isinstance(event_type, str):
            event_type = [event_type]
        query = query.filter(Event.type.in_(event_type))
    if source is not None:
        if isinstance(source, str):
            source = [source]
        query = query.filter(Event.source.in_(source))
    events = query.all()
    return events


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
