"""Observing database functions and ORM."""

from .models import *  # noqa: F401, F403
from .random import *  # noqa: F401, F403
from .utils import *  # noqa: F401, F403

try:
    # Try to use the alert database if it is available
    # This will add in the back-references to Events and Notices
    from gtecs.alert import database as alert_db
except Exception:
    pass
