"""ObsDB module parameters."""

import os
import sys

import configobj

import pkg_resources

import validate

from .version import __version__


# Load configspec file for default configuration
if os.path.exists('obsdb/data/configspec.ini'):
    # We are running in install dir, during installation
    CONFIGSPEC_FILE = 'obsdb/data/configspec.ini'
else:
    # We are being imported, find pkg_resources
    CONFIGSPEC_FILE = pkg_resources.resource_filename('obsdb', 'data/configspec.ini')

# Try to find .obsdb.conf file, look in the home directory and
# anywhere specified by OBSDB_CONF environment variable
paths = [os.path.expanduser("~")]
if "OBSDB_CONF" in os.environ:
    OBSDB_CONF_PATH = os.environ["OBSDB_CONF"]
    paths.append(OBSDB_CONF_PATH)
else:
    OBSDB_CONF_PATH = None

# Load the config file as a ConfigObj
config = configobj.ConfigObj({}, configspec=CONFIGSPEC_FILE)
CONFIG_FILE_PATH = None
for loc in paths:
    try:
        with open(os.path.join(loc, ".obsdb.conf")) as source:
            config = configobj.ConfigObj(source, configspec=CONFIGSPEC_FILE)
            CONFIG_FILE_PATH = loc
    except IOError:
        pass

# Validate ConfigObj, filling defaults from configspec if missing from config file
validator = validate.Validator()
result = config.validate(validator)
if result is not True:
    print('Config file validation failed')
    print([k for k in result if not result[k]])
    sys.exit(1)

############################################################
# Module parameters
VERSION = __version__

# Database parameters
DATABASE_USER = config['DATABASE_USER']
DATABASE_PASSWORD = config['DATABASE_PASSWORD']
DATABASE_HOST = config['DATABASE_HOST']
DATABASE_NAME = config['DATABASE_NAME']

DATABASE_LOCATION = '{}:{}@{}/{}'.format(DATABASE_USER, DATABASE_PASSWORD,
                                         DATABASE_HOST, DATABASE_NAME)

DATABASE_ECHO = bool(config['DATABASE_ECHO'])
DATABASE_PRE_PING = bool(config['DATABASE_PRE_PING'])
