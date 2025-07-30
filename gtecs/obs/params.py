"""Package parameters."""

from gtecs.common import config as pkg_config
from gtecs.common.package import get_package_version, load_config
from gtecs.common.system import get_local_ip


############################################################
# Load and validate config file
config, CONFIG_SPEC, CONFIG_FILE = load_config('obs', '.obs.conf')

############################################################
# Module parameters
VERSION = get_package_version('obs')

# General parameters
LOCAL_HOST = get_local_ip()

# File locations
FILE_PATH = pkg_config.CONFIG_PATH / 'obs'

############################################################
# Database parameters
DATABASE_USER = config['DATABASE_USER']
DATABASE_PASSWORD = config['DATABASE_PASSWORD']
DATABASE_HOST = config['DATABASE_HOST']
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
if PYRO_HOST == 'localhost':
    PYRO_HOST = LOCAL_HOST
PYRO_PORT = config['PYRO_PORT']
PYRO_URI = 'PYRO:scheduler@{}:{}'.format(PYRO_HOST, PYRO_PORT)
PYRO_TIMEOUT = config['PYRO_TIMEOUT']

SERVER_HOST = config['SERVER_HOST']
if SERVER_HOST == 'localhost':
    SERVER_HOST = LOCAL_HOST
SERVER_PORT = config['SERVER_PORT']
SERVER_PATH = config['SERVER_PATH']
SERVER_API_KEY = config['SERVER_API_KEY']

SCHEDULER_CHECK_PERIOD = config['SCHEDULER_CHECK_PERIOD']
WRITE_QUEUE_FILE = config['WRITE_QUEUE_FILE']
WRITE_QUEUE_PAGE = config['WRITE_QUEUE_PAGE']
SCHEDULER_SKIP_DAYTIME = config['SCHEDULER_SKIP_DAYTIME']
SCHEDULER_SUNALT_LIMIT = config['SCHEDULER_SUNALT_LIMIT']

############################################################
# Slack bot parameters
ENABLE_SLACK = config['ENABLE_SLACK']
SLACK_BOT_TOKEN = config['SLACK_BOT_TOKEN']
SLACK_DEFAULT_CHANNEL = config['SLACK_DEFAULT_CHANNEL']
