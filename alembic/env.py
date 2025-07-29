import re
from logging.config import fileConfig

from alembic_utils.pg_function import PGFunction
from alembic_utils.pg_trigger import PGTrigger
from alembic_utils.replaceable_entity import register_entities
from sqlalchemy import MetaData, create_engine

from alembic import context
from gtecs.obs import params
from gtecs.obs.database.models import Base, functions, triggers

# This is the Alembic Config object, which provides
# access to the values within the .ini file in use.
config = context.config
if config.config_file_name is not None:
    fileConfig(config.config_file_name)

# Filter Base.metadata to only include tables in the 'obs' schema
metadata = MetaData(schema="obs")
for table in Base.metadata.tables.values():
    if table.schema == "obs":
        table.tometadata(metadata)
target_metadata = metadata

# Automatically create and register alembic-utils entities for functions and triggers
entities = []
for func_info in functions.values():
    pg_function = PGFunction(
        schema=func_info["schema"],
        signature=func_info["signature"],
        definition=func_info["definition"],
    )
    entities.append(pg_function)
for trigger_info in triggers.values():
    pg_trigger = PGTrigger(
        schema=trigger_info["schema"],
        signature=trigger_info["signature"],
        on_entity=trigger_info["on_entity"],
        is_constraint=False,
        definition=trigger_info["definition"],
    )
    entities.append(pg_trigger)
register_entities(entities)


def get_url() -> str:
    """Get the database URL from the package config."""
    url = "postgresql://{}:{}@{}/gtecs".format(
        params.DATABASE_USER, params.DATABASE_PASSWORD, params.DATABASE_HOST
    )
    return url


def include_name(name, type_, parent_names):
    """Include only specific object types in autogenerate."""
    if type_ == "schema":
        return name in ["obs"]
    elif type_ == "grant_table":
        return False
    else:
        return True


def run_migrations_offline() -> None:
    """Run migrations in 'offline' mode.

    This configures the context with just a URL
    and not an Engine, though an Engine is acceptable
    here as well.  By skipping the Engine creation
    we don't even need a DBAPI to be available.

    Calls to context.execute() here emit the given string to the
    script output.

    """
    url = get_url()
    context.configure(
        url=url,
        target_metadata=target_metadata,
        literal_binds=True,
        dialect_opts={"paramstyle": "named"},
        include_schemas=True,
        include_name=include_name,
        version_table="alembic_version_obs",  # Use separate version table for obs schema
    )

    with context.begin_transaction():
        context.run_migrations()


def run_migrations_online() -> None:
    """Run migrations in 'online' mode.

    In this scenario we need to create an Engine
    and associate a connection with the context.

    """
    connectable = create_engine(get_url())

    with connectable.connect() as connection:
        context.configure(
            connection=connection,
            target_metadata=target_metadata,
            include_schemas=True,
            include_name=include_name,
            version_table="alembic_version_obs",  # Use separate version table for obs schema
        )

        with context.begin_transaction():
            context.run_migrations()


if context.is_offline_mode():
    run_migrations_offline()
else:
    run_migrations_online()
