# G-TeCS obs package

**G-TeCS** (*gee-teks*) is the GOTO Telescope Control System.

This package (`gtecs-obs`) contains functions for organising and scheduling observations.

Note this module is Python3 only and has been developed for Linux, otherwise use at your own risk.

## Requirements

This package requires several Python modules, which should be included during installation.

This package requires the following G-TeCS packages to function fully:

- [gtecs-common](https://github.com/GOTO-OBS/gtecs-common)

It also requires the following packages created for GOTO:

- [GOTO-tile](https://github.com/GOTO-OBS/goto-tile)

## Installation

Once you've downloaded or cloned the repository, in the base directory run:

    pip3 install . --user

You should then be able to import the module from within Python.

Several scripts from the `scripts` folder should also be added to your path, in particular the `scheduler` script which will monitor the database and determine which pointings to observe (see **Usage** below).

### Setting up the database

You'll need to make sure PostgreSQL is installed and configured before setting up the database. The config file contains parameters for the user and password to use when interacting with the database, make sure you create this user with rights first (e.g. `sudo -u postgres createuser -edP gtecs`) (a useful hint from https://stackoverflow.com/a/26735105:  edit `/etc/postgresql/12/main/pg_hba.conf` to change `local all all peer` to ` local all all md5` to remove the annoying user bits, then restart `sudo service postgresql restart`).

Then you can create the database with `createdb -O gtecs gtecs`, or just log into PostgreSQL as the `gtecs` user and run `CREATE DATABASE gtecs;`.

Database migrations are handled using Alembic, which should be installed as part of the package. To create the initial database schema, run:

    alembic upgrade head

Once the database is created run the `setup_obsdb` script with the `-d/--add-defaults` option to add in definitions for the basic GOTO telescopes, otherwise the database will be entirely empty.

### Configuration

The module will look for a file named `.obs.conf` either in the user's home directory, the `gtecs` subdirectory, or a path specified by the `GTECS_CONF` environment variable. An example file is included in the base directory of this repository.

### Testing

After setting up the database and installing the module, you can test it works correctly using the `test_database.py` script in the `/tests/` directory.

Note this script assumes you start from a completely clean, empty database.

## Usage

To run the scheduler after the package is installed the script can be started with

    scheduler start

For a list of arguments use `scheduler help`.

To have the scheduler run continuously as a system service edit the `scheduler.service` file to add your username (the default is `goto`, note you need to edit both the `User=` and `ExecStart=` entries), then copy it to the right location with

    sudo cp scheduler.service /etc/systemd/system/

Then enable and start the service with

    sudo systemctl enable scheduler
    sudo systemctl start scheduler

You can use `systemctl status scheduler` to see the status of the service, or `scheduler info` to see what the scheduler currently thinks is the highest-priority pointing. To view the full queue, use the `db_queue` command.
