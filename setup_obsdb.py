#!/usr/bin/env python3
"""Create a blank observing database."""

import argparse

from gtecs.obs import database as db


def run(overwrite=False, add_defaults=False, verbose=False):
    """Create and fill the database."""
    # Create a blank database
    print('Creating blank database...')
    db.create_database(overwrite, verbose)

    # Create the schema from the base
    print('Filling database schema...')
    db.fill_database(verbose)

    if add_defaults:
        # Create the default user (NB no password)
        print('Creating default user...')
        user = db.User('default', '', 'Default user')
        with db.open_session() as session:
            session.add(user)

        # Add in default GOTO sites and telescopes
        print('Creating entires for GOTO telescopes...')
        # TODO: Add 8-UT grid
        site1 = db.Site.from_name('LaPalma')
        site2 = db.Site.from_name('SSO')
        goto1 = db.Telescope(name='GOTO-1', site=site1)
        goto2 = db.Telescope(name='GOTO-2', site=site1)
        goto3 = db.Telescope(name='GOTO-3', site=site2)
        goto4 = db.Telescope(name='GOTO-4', site=site2)
        with db.open_session() as session:
            db.insert_items(session, [site1, site2, goto1, goto2, goto3, goto4])
            session.commit()

    print('Done')


if __name__ == '__main__':
    description = """Fill the GOTO observation database with survey tiles."""

    parser = argparse.ArgumentParser(description=description,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-v', '--verbose', action='store_true', default=False,
                        help='Print SQL statements?')
    parser.add_argument('-o', '--overwrite', action='store_true', default=False,
                        help='Overwrite an existing database [WARNING: WILL LOSE DATA]?')
    parser.add_argument('-d', '--add-defaults', action='store_true', default=False,
                        help='Add default (GOTO) telescope definitions?')

    args = parser.parse_args()

    run(args.overwrite, args.add_defaults, args.verbose)
