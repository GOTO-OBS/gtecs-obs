"""Slack messaging tools."""

import datetime
import os
from collections import Counter

from astropy.coordinates import EarthLocation
from astropy.time import Time

from gtecs.common.slack import send_message

from . import database as db
from . import params


def send_slack_msg(text, channel=None, *args, **kwargs):
    """Send a message to Slack.

    Parameters
    ----------
    text : string
        The message text.
    channel : string, optional
        The channel to post the message to.
        If None, defaults to `gtecs.obs.params.SLACK_DEFAULT_CHANNEL`.

    Other parameters are passed to `gtecs.common.slack.send_slack_msg`.

    """
    if channel is None:
        channel = params.SLACK_DEFAULT_CHANNEL

    if params.ENABLE_SLACK:
        # Use the common function
        send_message(text, channel, params.SLACK_BOT_TOKEN, *args, **kwargs)
    else:
        print('Slack Message:', text)


def send_database_report(slack_channel=None):
    """Send a Slack message containing the pending pointings in the database."""
    attachments = []
    with db.open_session() as session:
        pointings = session.query(db.Pointing).filter(db.Pointing.status == 'pending').all()
        msg = '*There are {} pending pointings in the database*'.format(len(pointings))

        # Pending pointings that are associated with a survey
        surveys = [pointing.survey for pointing in pointings
                   if pointing.survey is not None]
        if len(surveys) > 0:
            # Print number of pointings and surveys
            survey_counter = Counter(surveys)
            title = '{} pointing{} from {} sky survey{}:'.format(
                len(surveys), 's' if len(surveys) != 1 else '',
                len(survey_counter), 's' if len(survey_counter) != 1 else '')
            # Print info for all surveys
            text = '\n'
            for survey, count in survey_counter.most_common():
                text += '- `{}`'.format(survey.name)
                text += ' (_'
                text += '{} pointing{}'.format(count, 's' if count != 1 else '')
                survey_pointings = [pointing for pointing in pointings
                                    if pointing.survey == survey]
                ranks = sorted({p.rank for p in survey_pointings})
                text += ', rank={}{}'.format(ranks[0], '+' if len(ranks) > 1 else '')
                text += '_)\n'
        else:
            title = '0 pointings from surveys'
            text = ''
        attach = {'fallback': title,
                  'text': title + text,
                  }
        attachments.append(attach)

        # Other pending pointings
        objects = [pointing.object_name for pointing in pointings
                   if pointing.survey is None]
        if len(objects) > 0:
            # Print number of pointings and objects
            objects_counter = Counter(objects)
            title = '{} non-survey pointing{} of {} object{}:'.format(
                len(objects), 's' if len(objects) != 1 else '',
                len(objects_counter), 's' if len(objects_counter) != 1 else '')
            # Print info for all objects
            text = '\n'
            for obj, count in objects_counter.most_common():
                text += '- `{}`'.format(obj)
                text += ' (_'
                text += '{} pointing{}'.format(count, 's' if count != 1 else '')
                object_pointings = [pointing for pointing in pointings
                                    if pointing.object_name == obj]
                ranks = sorted({p.rank for p in object_pointings})
                text += ', rank={}{}'.format(ranks[0], '+' if len(ranks) > 1 else '')
                text += '_)\n'
        else:
            title += '0 non-survey pointings'
            text = ''
        attach = {'fallback': title,
                  'text': title + text,
                  }
        attachments.append(attach)

    send_slack_msg(msg, attachments=attachments, channel=slack_channel)


def send_observation_report(date=None, location=None, alt_limit=30, sun_limit=-12,
                            slack_channel=None):
    """Send a Slack message containing last night's observation plots."""
    if date is None:
        now = datetime.datetime.utcnow()
        if now.hour < 12:
            now = now - datetime.timedelta(days=1)
        date = now.strftime('%Y-%m-%d')
    if location is None:
        location = EarthLocation.of_site('lapalma')

    plot_direc = os.path.join(params.FILE_PATH, 'plots')
    if not os.path.exists(plot_direc):
        os.mkdir(plot_direc)

    # Get the dates for the start and end of the night just finished
    midday_yesterday = datetime.datetime.strptime(date + ' 12:00:00', '%Y-%m-%d %H:%M:%S')
    midday_today = midday_yesterday + datetime.timedelta(days=1)

    with db.open_session() as session:
        # Get the current grid from the database and create a SkyGrid
        db_grid = db.get_current_grid(session)
        grid = db_grid.skygrid

        # Use Astroplan to get all the tiles that would have been visible last night
        visible_tiles = grid.get_visible_tiles(location,
                                               time_range=(Time(midday_yesterday),
                                                           Time(midday_today)),
                                               alt_limit=alt_limit,
                                               sun_limit=sun_limit,
                                               )
        notvisible_tiles = [tile for tile in grid.tilenames if tile not in visible_tiles]

        # Get all (on-grid) pointings observed last night
        pointings = session.query(db.Pointing).filter(
            db.Pointing.status == 'completed',
            db.Pointing.grid == db_grid,
            db.Pointing.stopped_time > midday_yesterday,
            db.Pointing.stopped_time < midday_today,
        ).all()

        all_surveys = []
        all_tiles = []

        # First find the all-sky survey tiles completed last night
        allsky_survey = db_grid.surveys[0]
        tiles = [pointing.grid_tile.name for pointing in pointings
                 if pointing.survey == allsky_survey]
        all_surveys.append(allsky_survey.name)
        all_tiles.append(tiles)

        # Then find the tiles of any other survey pointings completed last night
        surveys = [pointing.survey for pointing in pointings
                   if pointing.survey != allsky_survey and pointing.survey is not None]
        if len(surveys) > 0:
            survey_counter = Counter(surveys)
            for survey, _ in survey_counter.most_common():
                tiles = [pointing.grid_tile.name for pointing in pointings
                         if pointing.survey == survey]
                all_surveys.append(survey.name)
                all_tiles.append(tiles)

        # Get the object names for other non-survey pointings
        # Here we count them as small, single-tile surveys
        objects = [pointing.object_name for pointing in pointings
                   if pointing.survey is None]
        if len(objects) > 0:
            object_counter = Counter(objects)
            for obj, _ in object_counter.most_common():
                tiles = [pointing.grid_tile.name for pointing in pointings
                         if pointing.object_name == obj]
                all_surveys.append(obj)
                all_tiles.append(tiles)

    # Remove the empty all-sky survey list if we didn't observe any
    n_obs = sum(len(tiles) for tiles in all_tiles)
    n_obs_allsky = len(all_tiles[0])
    if n_obs_allsky == 0:
        all_surveys = all_surveys[1:]
        all_tiles = all_tiles[1:]

    # Make a plot of last night's observations (assuming we observed anything)
    if n_obs > 0:
        msg = 'Last night coverage plot'

        # Create plot
        title = 'GOTO observations for\nnight beginning {}'.format(date)
        filepath = os.path.join(plot_direc, '{}_observed.png'.format(date))
        grid.plot(filename=filepath,
                  color={tilename: '0.5' for tilename in notvisible_tiles},
                  highlight=all_tiles,
                  highlight_label=all_surveys,
                  alpha=0.5,
                  title=title)

        # Send message to Slack with the plot attached
        send_slack_msg(msg, filepath=filepath, channel=slack_channel)

        # Create plot of all-sky survey coverage (assuming we observed any new ones)
        if n_obs_allsky > 0:
            msg = 'All-sky survey coverage plot'

            with db.open_session() as session:
                # Get the current survey from the database
                db_grid = db.get_current_grid(session)
                db_survey = db_grid.surveys[0]

                # Get all completed all-sky survey pointings since it started
                query = session.query(db.Pointing).filter(
                    db.Pointing.status == 'completed',
                    db.Pointing.survey == db_survey,
                    db.Pointing.stopped_time < midday_today,
                )
                survey_pointings = query.all()

                # Count tiles
                counter = Counter([p.grid_tile.name for p in survey_pointings])
                count_dict = dict(counter)

                # Get start date of the survey
                startdate = min(p.stopped_time for p in survey_pointings)

            # Create plot
            title = 'GOTO all-sky survey coverage\n'
            title += 'from {} to {}'.format(startdate.strftime('%Y-%m-%d'), date)
            filepath = os.path.join(plot_direc, '{}_survey.png'.format(date))
            grid.plot(filename=filepath,
                      color=count_dict,
                      discrete_colorbar=True,
                      highlight=all_tiles[0],
                      highlight_color='red',
                      highlight_label='observed last night',
                      alpha=0.5,
                      title=title)

            # Send message to Slack with the plot attached
            send_slack_msg(msg, filepath=filepath, channel=slack_channel)
        else:
            send_slack_msg('No all-sky survey tiles were observed last night',
                           channel=slack_channel)
    else:
        send_slack_msg('No tiles were observed last night', channel=slack_channel)

    return n_obs, n_obs_allsky
