"""Package parameters."""

import os
import sys

import configobj

from importlib.metadata import version
try:
    import importlib.resources as pkg_resources
except ImportError:
    # Python < 3.7
    import importlib_resources as pkg_resources  # type: ignore

import validate


# Load configspec file for default configuration
CONFIGSPEC = pkg_resources.read_text('gtecs.obs.data', 'configspec.ini').split('\n')
config = configobj.ConfigObj({}, configspec=CONFIGSPEC)

# Try to find the config file, look in the home directory and
# anywhere specified by GTECS_CONF environment variable
CONFIG_FILE = '.obs.conf'
home = os.path.expanduser('~')
paths = [home, os.path.join(home, 'gtecs'), os.path.join(home, '.gtecs')]
if 'GTECS_CONF' in os.environ:
    paths.append(os.environ['GTECS_CONF'])

# Load the config file as a ConfigObj
CONFIG_FILE_PATH = None
for loc in paths:
    try:
        with open(os.path.join(loc, CONFIG_FILE)) as source:
            config = configobj.ConfigObj(source, configspec=CONFIGSPEC)
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
VERSION = version('gtecs-obs')

# Database parameters
DATABASE_USER = config['DATABASE_USER']
DATABASE_PASSWORD = config['DATABASE_PASSWORD']
DATABASE_HOST = config['DATABASE_HOST']
DATABASE_NAME = config['DATABASE_NAME']

DATABASE_LOCATION = '{}:{}@{}/{}'.format(DATABASE_USER, DATABASE_PASSWORD,
                                         DATABASE_HOST, DATABASE_NAME)

DATABASE_ECHO = bool(config['DATABASE_ECHO'])
DATABASE_PRE_PING = bool(config['DATABASE_PRE_PING'])
