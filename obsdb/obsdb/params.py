import os
import sys
import socket
import pkg_resources
import configobj
import validate


# get a default spec for config file, either from local path, or installed path
if os.path.exists('obsdb/data/configspec.ini'):
    # we are running in install dir, during installation
    configspec_file = 'obsdb/data/configspec.ini'
else:
    # we are being imported, find pkg_resources
    configspec_file = pkg_resources.resource_filename('obsdb', 'data/configspec.ini')

# try and load config file
# look in current dir, home directory and anywhere specified by OBSDB_CONF environment variable
paths = [os.curdir, os.path.expanduser("~")]
if "OBSDB_CONF" in os.environ:
    paths.append(os.environ["OBSDB_CONF"])

# now load config file
config = configobj.ConfigObj({}, configspec=configspec_file)
for loc in paths:
    try:
        with open(os.path.join(loc, ".obsdb.conf")) as source:
            config = configobj.ConfigObj(source, configspec=configspec_file)
    except IOError as e:
        pass

# validate ConfigObj, filling defaults from configspec if missing from config file
validator = validate.Validator()
result = config.validate(validator)
if result != True:
    print('Config file validation failed')
    sys.exit(1)


# Database parameters
DATABASE_USER = config['DATABASE_USER']
DATABASE_PASSWORD = config['DATABASE_PASSWORD']
DATABASE_HOST = config['DATABASE_HOST']
DATABASE_NAME = config['DATABASE_NAME']

DATABASE_LOCATION = '{}:{}@{}/{}'.format(DATABASE_USER, DATABASE_PASSWORD,
                                         DATABASE_HOST, DATABASE_NAME)

DATABASE_ECHO = bool(config['DATABASE_ECHO'])
DATABASE_PRE_PING = bool(config['DATABASE_PRE_PING'])
