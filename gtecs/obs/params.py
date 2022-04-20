"""Package parameters."""

import os

from gtecs.common.package import load_config, get_package_version


############################################################
# Load and validate config file
config, CONFIG_SPEC, CONFIG_FILE = load_config('obs', '.obs.conf')

############################################################
# Module parameters
VERSION = get_package_version('obs')

# File locations
FILE_PATH = config['FILE_PATH']
if FILE_PATH in ['path_not_set', '/path/goes/here/']:
    raise ValueError('FILE_PATH not set, check config file ({})'.format(CONFIG_FILE))
QUEUE_PATH = os.path.join(FILE_PATH, 'queue')
HTML_PATH = os.path.join(FILE_PATH, 'html')
LOG_PATH = os.path.join(FILE_PATH, 'logs')

############################################################
# Logging parameters
FILE_LOGGING = config['FILE_LOGGING']
STDOUT_LOGGING = config['STDOUT_LOGGING']

############################################################
# Database parameters
DATABASE_USER = config['DATABASE_USER']
DATABASE_PASSWORD = config['DATABASE_PASSWORD']
DATABASE_HOST = config['DATABASE_HOST']
DATABASE_URL = '{}:{}@{}'.format(DATABASE_USER, DATABASE_PASSWORD, DATABASE_HOST)
DATABASE_NAME = config['DATABASE_NAME']
DATABASE_ECHO = bool(config['DATABASE_ECHO'])
DATABASE_PRE_PING = bool(config['DATABASE_PRE_PING'])

############################################################
# Scheduling parameters
DEFAULT_MIN_ALT = config['DEFAULT_MIN_ALT']
DEFAULT_MAX_SUNALT = config['DEFAULT_MAX_SUNALT']
DEFAULT_MAX_MOON = config['DEFAULT_MAX_MOON']
DEFAULT_MIN_MOONSEP = config['DEFAULT_MIN_MOONSEP']

# these should be function arguments, so e.g. the simulations could change them
WEIGHTING_WEIGHT = config['WEIGHTING_WEIGHT']
AIRMASS_WEIGHT = config['AIRMASS_WEIGHT']
TTS_WEIGHT = config['TTS_WEIGHT']

HARD_ALT_LIM = config['HARD_ALT_LIM']
HARD_HA_LIM = config['HARD_HA_LIM']
READOUT_TIME = config['READOUT_TIME']
TEMPLATE_REQUIREMENT = config['TEMPLATE_REQUIREMENT']

############################################################
# Scheduler parameters
PYRO_HOST = config['PYRO_HOST']
PYRO_PORT = config['PYRO_PORT']
PYRO_URI = 'PYRO:scheduler@{}:{}'.format(PYRO_HOST, PYRO_PORT)
PYRO_TIMEOUT = config['PYRO_TIMEOUT']

SCHEDULER_CHECK_PERIOD = config['SCHEDULER_CHECK_PERIOD']
WRITE_QUEUE_FILE = config['WRITE_QUEUE_FILE']
WRITE_QUEUE_PAGE = config['WRITE_QUEUE_PAGE']

############################################################
# Slack bot parameters
ENABLE_SLACK = config['ENABLE_SLACK']
SLACK_BOT_TOKEN = config['SLACK_BOT_TOKEN']
SLACK_DEFAULT_CHANNEL = config['SLACK_DEFAULT_CHANNEL']
