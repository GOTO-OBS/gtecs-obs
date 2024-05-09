"""Observation database functions and ORM."""

from .models import *  # noqa: F401, F403
from .random import *  # noqa: F401, F403
from .utils import *  # noqa: F401, F403

try:
    # Try to use the alert database if it is available
    # This will add in the back-references to Events and Notices
    # We also import the specific table classes, so we can query them easily
    from gtecs.alert import database as alert_db  # noqa: F401
    from gtecs.alert.database import Event, Notice  # noqa: F401
except Exception:
    pass
