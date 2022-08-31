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

Several scripts from the `scripts` folder should also be added to your path, in particular the `scheduler` script which will monitor the database and determine which pointings to observe (see **Usage** below).

### Setting up the database

Make sure MariaDB is installed and configured (recommended over MySQL).

Then run the `setup_obsdb` script to create a blank database with all the required tables. Note you can use the `-d/--add-defaults` option to add in definitions for the basic GOTO telescopes, otherwise the database will be entirely empty.

### Configuration

The module will look for a file named `.obs.conf` either in the user's home directory, the `gtecs` subdirectory, or a path specified by the `GTECS_CONF` environment variable. An example file is included in the base directory of this repository.

After installing this package copy this sample config file to one of the above locations, and change the file path parameters to specify where you want the package to save files.

Once that has been done run the `initialise.py` script to create the expected directory structure at that location.

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
