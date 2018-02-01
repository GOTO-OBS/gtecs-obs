# GOTO Observation Database (ObsDB)

This repository contains the MySQL schema for the GOTO pointings database and a Python module *obsdb* which has various SQLAlchemy objects and functions to interact with the database.

This module was initially part of the G-TeCS control system repository but was split out on 31 Jan 2018.


Setting up the database
-----------------------

Make sure MySQL is installed (or a compatible alternative, like MariaDB).

Start a session as root:
> `$ mysql -u root -p`

Setup the `goto_obs` database by sourcing from the schema included in this repository:
> `mysql> source schema.sql`


Installation of the Python module
---------------------------------

Once you've downloaded or cloned the repository, in the `obsdb` directory:

> `$ python setup.py install`

Or alternatively:

> `$ pip install .`

You should then be able to import the API using `import obsdb` within Python.

Several scripts from the `obsdb/scripts` folder should also be added to your path.


Configuration
-------------
Configuration of the `obsdb` module is done using a config file and the Python module [configobj](http://configobj.readthedocs.io/en/latest/).
When running, the module will look for a file named *.obsdb.conf* either in the current directory, the user's home directory or any path specified by the *OBSDB_CONF* environment variable.

An example *.obsdb.conf* file is present in the base directory of this repository.

If no such file is found, `obsdb` will use the default config, as shown in the *data* directory.
Users can override any of these default settings using the *.obsdb.conf* file.


Testing
-------

After setting up the database and installing the module, you can test it works correctly using the included `test_database.py` script in the base directory.

> `$ python test_database.py`

Note this script assumes you start from a completely clean, empty database.


Martin Dyer
Last update: 1 Feb 2018
