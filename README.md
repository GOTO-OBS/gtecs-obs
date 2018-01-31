# GOTO Obs Database

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

Martin Dyer
Last update: 31 Jan 2018
