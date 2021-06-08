# G-TeCS obs package

**G-TeCS** (*gee-teks*) is the GOTO Telescope Control System.

This package (`gtecs-obs`) contains functions for organising and scheduling observations.

Note this module is Python3 only and has been developed for Linux, otherwise use at your own risk.

## Requirements

This package requires several Python modules, which should be included during installation.

This package doesn't require any other G-TeCS packages to be installed, but it itself is a requirement of several of them.

This package requires the following packages created for GOTO:

- [GOTO-tile](https://github.com/GOTO-OBS/goto-tile)

## Installation

Once you've downloaded or cloned the repository, in the base directory run:

    pip3 install . --user

You should then be able to import the module from within Python.

Several scripts from the `scripts` folder should also be added to your path.

### Setting up the database

Make sure MariaDB is installed (recommended over MySQL).

Start a session as root:

    mysql -u root -p

Setup the `goto_obs` database by sourcing from the schema included in this repository:

    mysql> source schema.sql

### Configuration

The module will look for a file named `.obsdb.conf` either in the user's home directory or any path specified by the `OBSDB_CONF` environment variable. An example file is included in the base directory of this repository.

When installing ObsDB, copy the included `.obsdb.conf` file to one of the above locations, and change the database parameters to specify how to connect to your local database.

### Testing

After setting up the database and installing the module, you can test it works correctly using the `test_database.py` script in the `/tests/` directory.

Note this script assumes you start from a completely clean, empty database.

## Usage

TODO
