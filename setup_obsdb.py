"""Create a blank observing database."""

import argparse

from gtecs.obs import database as db


def run(overwrite=False, verbose=False):
    """Create and fill the database."""
    # Create a blank database
    print('Creating blank database...')
    db.create_database(overwrite, verbose)

    # Create the schema from the base
    print('Filling database schema...')
    db.fill_database(verbose)

    # Create the default user (NB no password)
    print('Creating default user...')
    user = db.User('goto', '', 'Default user')
    with db.open_session() as session:
        session.add(user)

    print('Done')


if __name__ == '__main__':
    description = """Fill the GOTO observation database with survey tiles."""

    parser = argparse.ArgumentParser(description=description,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-v', '--verbose', action='store_true', default=False,
                        help='Print SQL statements?')
    parser.add_argument('-o', '--overwrite', action='store_true', default=False,
                        help='Overwrite an existing database [WARNING: WILL LOSE DATA]?')

    args = parser.parse_args()

    run(args.overwrite, args.verbose)
