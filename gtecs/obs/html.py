"""HTML page functions for the scheduler."""

import os

from astroplan.moon import moon_illumination

from astropy import units as u
from astropy.coordinates import AltAz, ICRS, get_moon, get_sun
from astropy.time import Time

from . import database as db
from . import params


# setup
# common strings
html_size2 = '<font size=2 color=black face=\"Courier New\">\n'
html_size5 = '<font size=5 color=black face=\"Courier New\">\n'
popup_str = ('<div class=\"apple_overlay\" id=\"overlay\">' +
             '<div class=\"contentWrap\"></div>' + '</div>')

# set debug level
debug = 1

# interval between automatic refreshes
html_refresh = 15
html_refresh_string = ('<meta http-equiv=\"refresh\" content=\"' +
                       str(html_refresh) + '\">\n')


def import_queue_file():
    """Import the queue file."""
    import json
    lines = []
    with open(os.path.join(params.QUEUE_PATH, 'queue_info')) as f:
        for line in f.readlines():
            lines.append(line)

    time = Time(json.loads(lines[0]))
    all_constraint_names = json.loads(lines[1])
    # remaining lines are pointings
    pointing_list = []
    for line in lines[2:]:
        db_id, altaz_start, altaz_end, constraints = json.loads(line)
        db_id = int(db_id)
        constraint_names, valid_arr = list(zip(*constraints))
        valid_bools = [bool(x) for x in valid_arr]
        pointing_info = [db_id, altaz_start, altaz_end,
                         list(constraint_names), valid_bools]
        pointing_list.append(pointing_info)
    return time, all_constraint_names, pointing_list


def write_flag_file(pointing, time, all_constraint_names, pointing_info):
    """Write flag file for a given pointing."""
    db_id, altaz_start, altaz_end, constraint_names, valid_arr = pointing_info
    flag_filename = os.path.join(params.HTML_PATH, 'ID_{}_flags.html'.format(db_id))

    with open(flag_filename, 'w') as f:
        f.write('<html><body>\n')
        f.write(html_size2)
        f.write('ID_{:.0f}<br>\n'.format(db_id))

        time.format = 'iso'
        time.precision = 0
        f.write(str(time) + '<br>\n')

        for name, valid in zip(constraint_names, valid_arr):
            if valid:
                f.write('<font color=darkgreen>')
            else:
                f.write('<font color=red>')
            f.write(name + ' = ' + str(valid))
            f.write('</font><br>\n')

        for name in all_constraint_names:
            if name not in constraint_names:
                f.write('<font color=grey>')
                f.write(name + ' = N/A')
                f.write('</font><br>\n')

        f.write('<br>\n')

        if debug == 1:
            f.write('Debug info:<br>\n')
            f.write('now = ' + str(time) + '<br>\n')
            start = Time(pointing.start_time)
            start.format = 'iso'
            start.precision = 0
            f.write('start_time = ' + str(start) + '<br>\n')
            if pointing.stop_time:
                stop = Time(pointing.stop_time)
                stop.format = 'iso'
                stop.precision = 0
            else:
                stop = 'None'
            f.write('stop_time = ' + str(stop) + '<br>\n')

            target = ICRS(pointing.ra * u.deg, pointing.dec * u.deg)
            ra = target.ra.to_string(sep=':', precision=2, unit=u.hour)
            dec = target.dec.to_string(sep=':', precision=2)
            f.write('ra = ' + ra + '<br>\n')
            f.write('dec = ' + dec + '<br>\n')

            alt_start, az_start = altaz_start
            f.write('alt_start = {:.2f}<br>\n'.format(alt_start))
            f.write('az_start = {:.2f}<br>\n'.format(az_start))

            alt_end, az_end = altaz_end
            f.write('alt_end = {:.2f}<br>\n'.format(alt_end))
            f.write('az_end = {:.2f}<br>\n'.format(az_end))
        f.write("</body></html>")


def write_exp_file(db_id, exposure_sets):
    """Write exposure files for a pointing."""
    exp_filename = os.path.join(params.HTML_PATH, 'ID_{}_exp.html'.format(db_id))

    # unlike the flags, exposure info doesn't change
    # so don't re-write the files if they're already there!
    if not os.path.exists(exp_filename):
        with open(exp_filename, 'w') as f:
            f.write('<html><body><table cellpadding=5 cellspacing=5>\n')
            f.write('<tr>' +
                    '<td><b> </b></td>' +
                    '<td><b>UTs</b></td>' +
                    '<td><b>NumExp</b></td>' +
                    '<td><b>Exptime</b></td>' +
                    '<td><b>Filter</b></td>' +
                    '<td><b>Bin Factor</b></td>' +
                    '<td><b>Type</b></td>' +
                    '<td><b>RaOff</b></td>' +
                    '<td><b>DecOff</b></td>' +
                    '</tr>\n')
            i = 1
            for exposure_set in exposure_sets:
                f.write('<tr>' +
                        '<td align=right><b>' + str(i) + '</b></td>' +
                        '<td> ' + str(exposure_set.ut_mask) + ' </td>' +
                        '<td> ' + str(exposure_set.num_exp) + ' </td>' +
                        '<td> ' + str(exposure_set.exptime) + ' </td>' +
                        '<td> ' + str(exposure_set.filt) + ' </td>' +
                        '<td> ' + str(exposure_set.binning) + ' </td>' +
                        '</tr>\n')
                i += 1
            f.write("</table></body></html>")


def write_queue_page(observer):
    """Write the GOTO queue page."""
    # load any needed infomation saved by the scheduler
    time, all_constraint_names, pointing_list = import_queue_file()

    queue_filename = os.path.join(params.HTML_PATH, 'queue.html')
    with open(queue_filename, 'w') as f:
        f.write('<html><head>\n')
        f.write('<script src=\"jquery.tools.min.js\"></script>' +
                '<link rel=\"stylesheet\" type=\"text/css\" ' +
                'href=\"overlay-apple.css\"/>\n')
        f.write('<style>#overlay {background-image:url(transparent.png); ' +
                'color:#efefef;height:450px;}' +
                'div.contentWrap {height:441px;overflow-y:auto;}</style>\n')
        f.write(html_refresh_string)
        f.write('<style type=\"text/css\"><!-- ' +
                'td{font-family: arial; font-size: 10pt; color: black; ' +
                'white-space: nowrap;} tr{background-color: white;} --->' +
                '</style>\n')
        f.write('<title>GOTO queue</title></head>\n')
        f.write('<body bgcolor=white vlink=black alink=black link=black>\n')
        f.write('<center><br><br><img src=goto-logo.jpg width=50%><br><br>\n')

        f.write(html_size5)
        f.write('queue')
        f.write(html_size2)
        f.write('<br><br>')

        f.write('last updated:\n')
        time.format = 'iso'
        time.precision = 0
        f.write(str(time))
        f.write('<br>\n')

        f.write('LST: ')
        import warnings
        warnings.simplefilter("ignore", UnicodeWarning)
        lst = time.sidereal_time('mean', longitude=observer.location.lon)
        f.write(lst.to_string(sep=':', precision=2))

        altaz_frame = AltAz(obstime=time, location=observer.location)

        f.write('  SunAlt: ')
        sun = get_sun(time)
        sun_alt = sun.transform_to(altaz_frame).alt.degree
        f.write('{:.1f} deg'.format(sun_alt))

        f.write('  MoonAlt: ')
        moon = get_moon(time)
        moon_alt = moon.transform_to(altaz_frame).alt.degree
        f.write('{:.1f} deg'.format(moon_alt))

        f.write('  MoonPhase: ')
        moon_ill = moon_illumination(time)
        if 0 <= moon_ill < 0.25:
            moonstring = "D"
        elif 0.25 <= moon_ill < 0.65:
            moonstring = "G"
        elif 0.65 <= moon_ill <= 1.00:
            moonstring = "B"
        f.write('{:.2f} [{}]'.format(moon_ill, moonstring))

        f.write('<br><br>' +
                '<table border=1 bordercolor=#00639C cellspacing=\"0\"' +
                ' cellpadding=\"10\"><tr bgcolor=white>\n')
        f.write('<td><b> Target name </b></td>' +
                '<td><b> RA (j2000) </b></td>' +
                '<td><b> Dec (j2000) </b></td>' +
                '<td><b> Priority </b></td>' +
                '<td><b> MinTime </b></td>' +
                '<td><b> MinAlt </b></td>' +
                '<td><b> MaxSun </b></td>' +
                '<td><b> MaxMoon </b></td>' +
                '<td><b> MinMoonSep </b></td>' +
                '<td><b> User </b></td>' +
                '<td><b> Start UTC </b></td>' +
                '<td><b> Stop UTC </b></td>' +
                '<td><b> Filename/Flags </b></td>' +
                '</tr>\n')

        for pointing_info in pointing_list:
            db_id = pointing_info[0]

            # find database info
            session = db.load_session()
            pointing = db.get_pointing_by_id(session, db_id)
            if pointing.exposure_sets is not None:
                exposure_sets = pointing.exposure_sets
            username = pointing.user.username
            session.close()

            # create the small pointing files
            write_flag_file(pointing, time, all_constraint_names, pointing_info)
            flag_link = 'ID_{}_flags.html'.format(db_id)
            flag_str = ('<a href=' + flag_link + ' rel=\"#overlay\">' +
                        'ID_' + str(db_id) + '</a>' + popup_str)

            write_exp_file(db_id, exposure_sets)
            exp_link = 'ID_{}_exp.html'.format(db_id)
            exp_str = ('<a href=' + exp_link + ' rel=\"#overlay\">' +
                       str(pointing.name) + '</a>' + popup_str)

            # find ra/dec
            target = ICRS(pointing.ra * u.deg, pointing.dec * u.deg)
            ra = target.ra.to_string(sep=':', precision=2, unit=u.hour)
            dec = target.dec.to_string(sep=':', precision=2)

            # write table
            f.write("<tr bgcolor=white>\n")
            f.write('<td>' + exp_str + '</td>' +
                    '<td>' + ra + '</td>' +
                    '<td>' + dec + '</td>' +
                    '<td>' + str(pointing.min_time) + '</td>' +
                    '<td>' + str(pointing.min_alt) + '</td>' +
                    '<td>' + str(pointing.max_sunalt) + '</td>' +
                    '<td>' + str(pointing.max_moon) + '</td>' +
                    '<td>' + str(pointing.min_moonsep) + '</td>' +
                    '<td>' + str(username) + '</td>' +
                    '<td>' + str(pointing.start_time) + '</td>' +
                    '<td>' + str(pointing.stop_time) + '</td>' +
                    '<td>' + flag_str + '</td>' +
                    '</tr>\n')
            f.write("</tr>\n")
        f.write("</table>\n")
        f.write('<br><br></center>\n')
        f.write('<script>$(function() ' +
                '{$(\"a[rel]\").overlay({mask: \'black\',effect: \'apple\',' +
                'onBeforeLoad: function() ' +
                '{var wrap = this.getOverlay().find(\".contentWrap\");' +
                'wrap.load(this.getTrigger().attr(\"href\"));}});});' +
                '</script>\n')

        f.write('</body></html>\n')
