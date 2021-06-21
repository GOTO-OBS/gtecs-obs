"""Package parameters."""

from gtecs.common.package import load_config, get_package_version


############################################################
# Load and validate config file
config, CONFIG_SPEC, CONFIG_FILE = load_config('obs', '.obs.conf')

############################################################
# Module parameters
VERSION = get_package_version('obs')

############################################################
# Database parameters
DATABASE_USER = config['DATABASE_USER']
DATABASE_PASSWORD = config['DATABASE_PASSWORD']
DATABASE_HOST = config['DATABASE_HOST']
DATABASE_NAME = config['DATABASE_NAME']

DATABASE_LOCATION = '{}:{}@{}/{}'.format(DATABASE_USER, DATABASE_PASSWORD,
                                         DATABASE_HOST, DATABASE_NAME)

DATABASE_ECHO = bool(config['DATABASE_ECHO'])
DATABASE_PRE_PING = bool(config['DATABASE_PRE_PING'])
