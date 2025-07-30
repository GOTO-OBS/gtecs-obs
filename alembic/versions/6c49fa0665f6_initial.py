"""initial

Revision ID: 6c49fa0665f6
Revises:
Create Date: 2025-07-29 10:27:49.647050

"""

from typing import Sequence, Union

import sqlalchemy as sa
from alembic_utils.pg_function import PGFunction
from alembic_utils.pg_trigger import PGTrigger
from sqlalchemy import text as sql_text

from alembic import op

# revision identifiers, used by Alembic.
revision: str = "6c49fa0665f6"
down_revision: Union[str, None] = None
branch_labels: Union[str, Sequence[str], None] = None
depends_on: Union[str, Sequence[str], None] = None


def upgrade() -> None:
    """Upgrade schema."""
    # Create obs schema if it doesn't exist
    op.execute("CREATE SCHEMA IF NOT EXISTS obs")

    # Create tables
    # Note the order of table creation is important due to foreign key constraints.

    # Grids table
    op.create_table(
        "grids",
        sa.Column("id", sa.Integer(), nullable=False),
        sa.Column("name", sa.String(length=255), nullable=False),
        sa.Column("ra_fov", sa.Float(), nullable=False),
        sa.Column("dec_fov", sa.Float(), nullable=False),
        sa.Column("ra_overlap", sa.Float(), nullable=False),
        sa.Column("dec_overlap", sa.Float(), nullable=False),
        sa.Column("algorithm", sa.String(length=255), nullable=False),
        sa.Column("ts", sa.DateTime(), server_default=sa.text("now()"), nullable=False),
        sa.PrimaryKeyConstraint("id"),
        schema="obs",
    )
    op.create_index(
        op.f("ix_obs_grids_name"), "grids", ["name"], unique=True, schema="obs"
    )

    # Sites table
    op.create_table(
        "sites",
        sa.Column("id", sa.Integer(), nullable=False),
        sa.Column("name", sa.String(length=255), nullable=False),
        sa.Column("latitude", sa.Float(), nullable=False),
        sa.Column("longitude", sa.Float(), nullable=False),
        sa.Column("height", sa.Float(), nullable=False),
        sa.Column("ts", sa.DateTime(), server_default=sa.text("now()"), nullable=False),
        sa.PrimaryKeyConstraint("id"),
        schema="obs",
    )
    op.create_index(
        op.f("ix_obs_sites_name"), "sites", ["name"], unique=True, schema="obs"
    )

    # Surveys table
    op.create_table(
        "surveys",
        sa.Column("id", sa.Integer(), nullable=False),
        sa.Column("name", sa.String(length=255), nullable=False),
        sa.Column("ts", sa.DateTime(), server_default=sa.text("now()"), nullable=False),
        sa.PrimaryKeyConstraint("id"),
        schema="obs",
    )
    op.create_index(
        op.f("ix_obs_surveys_name"), "surveys", ["name"], unique=False, schema="obs"
    )

    # Users table
    op.create_table(
        "users",
        sa.Column("id", sa.Integer(), nullable=False),
        sa.Column("username", sa.String(length=255), nullable=False),
        sa.Column("password", sa.String(length=255), nullable=False),
        sa.Column("full_name", sa.Text(), nullable=False),
        sa.Column("ts", sa.DateTime(), server_default=sa.text("now()"), nullable=False),
        sa.PrimaryKeyConstraint("id"),
        sa.UniqueConstraint("username"),
        schema="obs",
    )

    # GridTiles table
    op.create_table(
        "grid_tiles",
        sa.Column("id", sa.Integer(), nullable=False),
        sa.Column("name", sa.String(length=255), nullable=False),
        sa.Column("ra", sa.Float(), nullable=False),
        sa.Column("dec", sa.Float(), nullable=False),
        sa.Column("grid_id", sa.Integer(), nullable=False),
        sa.Column("ts", sa.DateTime(), server_default=sa.text("now()"), nullable=False),
        sa.ForeignKeyConstraint(
            ["grid_id"],
            ["obs.grids.id"],
        ),
        sa.PrimaryKeyConstraint("id"),
        schema="obs",
    )
    op.create_index(
        op.f("ix_obs_grid_tiles_grid_id"),
        "grid_tiles",
        ["grid_id"],
        unique=False,
        schema="obs",
    )
    op.create_index(
        op.f("ix_obs_grid_tiles_name"),
        "grid_tiles",
        ["name"],
        unique=False,
        schema="obs",
    )

    # Telescopes table
    op.create_table(
        "telescopes",
        sa.Column("id", sa.Integer(), nullable=False),
        sa.Column("name", sa.String(length=255), nullable=False),
        sa.Column("horizon", sa.Text(), nullable=True),
        sa.Column("site_id", sa.Integer(), nullable=False),
        sa.Column("grid_id", sa.Integer(), nullable=True),
        sa.Column("ts", sa.DateTime(), server_default=sa.text("now()"), nullable=False),
        sa.ForeignKeyConstraint(
            ["grid_id"],
            ["obs.grids.id"],
        ),
        sa.ForeignKeyConstraint(
            ["site_id"],
            ["obs.sites.id"],
        ),
        sa.PrimaryKeyConstraint("id"),
        schema="obs",
    )
    op.create_index(
        op.f("ix_obs_telescopes_grid_id"),
        "telescopes",
        ["grid_id"],
        unique=False,
        schema="obs",
    )
    op.create_index(
        op.f("ix_obs_telescopes_name"),
        "telescopes",
        ["name"],
        unique=True,
        schema="obs",
    )
    op.create_index(
        op.f("ix_obs_telescopes_site_id"),
        "telescopes",
        ["site_id"],
        unique=False,
        schema="obs",
    )

    # Targets table
    op.create_table(
        "targets",
        sa.Column("id", sa.Integer(), nullable=False),
        sa.Column("name", sa.Text(), nullable=False),
        sa.Column("ra", sa.Float(), nullable=False),
        sa.Column("dec", sa.Float(), nullable=False),
        sa.Column("rank", sa.Integer(), nullable=True),
        sa.Column("weight", sa.Float(), nullable=False),
        sa.Column("is_template", sa.Boolean(), nullable=False),
        sa.Column(
            "start_time", sa.DateTime(), server_default=sa.text("now()"), nullable=False
        ),
        sa.Column("stop_time", sa.DateTime(), nullable=True),
        sa.Column(
            "creation_time",
            sa.DateTime(),
            server_default=sa.text("now()"),
            nullable=False,
        ),
        sa.Column("deleted_time", sa.DateTime(), nullable=True),
        sa.Column("user_id", sa.Integer(), nullable=False),
        sa.Column("grid_tile_id", sa.Integer(), nullable=True),
        sa.Column("survey_id", sa.Integer(), nullable=True),
        sa.Column("ts", sa.DateTime(), server_default=sa.text("now()"), nullable=False),
        sa.ForeignKeyConstraint(
            ["grid_tile_id"],
            ["obs.grid_tiles.id"],
        ),
        sa.ForeignKeyConstraint(
            ["survey_id"],
            ["obs.surveys.id"],
        ),
        sa.ForeignKeyConstraint(
            ["user_id"],
            ["obs.users.id"],
        ),
        sa.PrimaryKeyConstraint("id"),
        schema="obs",
    )
    op.create_index(
        op.f("ix_obs_targets_grid_tile_id"),
        "targets",
        ["grid_tile_id"],
        unique=False,
        schema="obs",
    )
    op.create_index(
        op.f("ix_obs_targets_start_time"),
        "targets",
        ["start_time"],
        unique=False,
        schema="obs",
    )
    op.create_index(
        op.f("ix_obs_targets_stop_time"),
        "targets",
        ["stop_time"],
        unique=False,
        schema="obs",
    )
    op.create_index(
        op.f("ix_obs_targets_survey_id"),
        "targets",
        ["survey_id"],
        unique=False,
        schema="obs",
    )
    op.create_index(
        op.f("ix_obs_targets_user_id"),
        "targets",
        ["user_id"],
        unique=False,
        schema="obs",
    )

    # ExposureSets table
    op.create_table(
        "exposure_sets",
        sa.Column("id", sa.Integer(), nullable=False),
        sa.Column("num_exp", sa.Integer(), nullable=False),
        sa.Column("exptime", sa.Float(), nullable=False),
        sa.Column("filter", sa.String(length=1), nullable=False),
        sa.Column("binning", sa.Integer(), nullable=False),
        sa.Column("dithering", sa.Integer(), nullable=False),
        sa.Column("ut_mask", sa.Integer(), nullable=True),
        sa.Column("target_id", sa.Integer(), nullable=True),
        sa.Column("ts", sa.DateTime(), server_default=sa.text("now()"), nullable=False),
        sa.ForeignKeyConstraint(
            ["target_id"],
            ["obs.targets.id"],
        ),
        sa.PrimaryKeyConstraint("id"),
        schema="obs",
    )
    op.create_index(
        op.f("ix_obs_exposure_sets_target_id"),
        "exposure_sets",
        ["target_id"],
        unique=False,
        schema="obs",
    )

    # Strategies table
    op.create_table(
        "strategies",
        sa.Column("id", sa.Integer(), nullable=False),
        sa.Column("num_todo", sa.Integer(), nullable=False),
        sa.Column("stop_time", sa.DateTime(), nullable=True),
        sa.Column("min_time", sa.Float(), nullable=True),
        sa.Column("too", sa.Boolean(), nullable=False),
        sa.Column("requires_template", sa.Boolean(), nullable=False),
        sa.Column("min_alt", sa.Float(), nullable=False),
        sa.Column("max_sunalt", sa.Float(), nullable=False),
        sa.Column("max_moon", sa.String(length=1), nullable=False),
        sa.Column("min_moonsep", sa.Float(), nullable=False),
        sa.Column("tel_mask", sa.Integer(), nullable=True),
        sa.Column("target_id", sa.Integer(), nullable=True),
        sa.Column("ts", sa.DateTime(), server_default=sa.text("now()"), nullable=False),
        sa.ForeignKeyConstraint(
            ["target_id"],
            ["obs.targets.id"],
        ),
        sa.PrimaryKeyConstraint("id"),
        schema="obs",
    )
    op.create_index(
        op.f("ix_obs_strategies_target_id"),
        "strategies",
        ["target_id"],
        unique=False,
        schema="obs",
    )

    # TimeBlocks table
    op.create_table(
        "time_blocks",
        sa.Column("id", sa.Integer(), nullable=False),
        sa.Column("block_num", sa.Integer(), nullable=False),
        sa.Column("wait_time", sa.Float(), nullable=False),
        sa.Column("valid_time", sa.Float(), nullable=True),
        sa.Column("rank_change", sa.Integer(), nullable=True),
        sa.Column("strategy_id", sa.Integer(), nullable=False),
        sa.Column("ts", sa.DateTime(), server_default=sa.text("now()"), nullable=False),
        sa.ForeignKeyConstraint(
            ["strategy_id"],
            ["obs.strategies.id"],
        ),
        sa.PrimaryKeyConstraint("id"),
        schema="obs",
    )
    op.create_index(
        op.f("ix_obs_time_blocks_block_num"),
        "time_blocks",
        ["block_num"],
        unique=False,
        schema="obs",
    )
    op.create_index(
        op.f("ix_obs_time_blocks_strategy_id"),
        "time_blocks",
        ["strategy_id"],
        unique=False,
        schema="obs",
    )

    # Pointings table
    op.create_table(
        "pointings",
        sa.Column("id", sa.Integer(), nullable=False),
        sa.Column("rank", sa.Integer(), nullable=True),
        sa.Column(
            "start_time", sa.DateTime(), server_default=sa.text("now()"), nullable=False
        ),
        sa.Column("stop_time", sa.DateTime(), nullable=True),
        sa.Column(
            "creation_time",
            sa.DateTime(),
            server_default=sa.text("now()"),
            nullable=False,
        ),
        sa.Column("running_time", sa.DateTime(), nullable=True),
        sa.Column("completed", sa.Boolean(), nullable=False),
        sa.Column("finished_time", sa.DateTime(), nullable=True),
        sa.Column("valid", sa.Boolean(), nullable=True),
        sa.Column("validated_time", sa.DateTime(), nullable=True),
        sa.Column("target_id", sa.Integer(), nullable=False),
        sa.Column("time_block_id", sa.Integer(), nullable=False),
        sa.Column("strategy_id", sa.Integer(), nullable=False),
        sa.Column("telescope_id", sa.Integer(), nullable=True),
        sa.Column("ts", sa.DateTime(), server_default=sa.text("now()"), nullable=False),
        sa.ForeignKeyConstraint(
            ["strategy_id"],
            ["obs.strategies.id"],
        ),
        sa.ForeignKeyConstraint(
            ["target_id"],
            ["obs.targets.id"],
        ),
        sa.ForeignKeyConstraint(
            ["telescope_id"],
            ["obs.telescopes.id"],
        ),
        sa.ForeignKeyConstraint(
            ["time_block_id"],
            ["obs.time_blocks.id"],
        ),
        sa.PrimaryKeyConstraint("id"),
        schema="obs",
    )
    op.create_index(
        op.f("ix_obs_pointings_start_time"),
        "pointings",
        ["start_time"],
        unique=False,
        schema="obs",
    )
    op.create_index(
        op.f("ix_obs_pointings_stop_time"),
        "pointings",
        ["stop_time"],
        unique=False,
        schema="obs",
    )
    op.create_index(
        op.f("ix_obs_pointings_strategy_id"),
        "pointings",
        ["strategy_id"],
        unique=False,
        schema="obs",
    )
    op.create_index(
        op.f("ix_obs_pointings_target_id"),
        "pointings",
        ["target_id"],
        unique=False,
        schema="obs",
    )
    op.create_index(
        op.f("ix_obs_pointings_telescope_id"),
        "pointings",
        ["telescope_id"],
        unique=False,
        schema="obs",
    )
    op.create_index(
        op.f("ix_obs_pointings_time_block_id"),
        "pointings",
        ["time_block_id"],
        unique=False,
        schema="obs",
    )

    # Extra functions, triggers etc
    # Create a function to update the 'ts' column on update
    obs_update_ts = PGFunction(
        schema="obs",
        signature="update_ts()",
        definition="RETURNS TRIGGER\nLANGUAGE plpgsql AS\n$function$\nBEGIN\n    NEW.ts := now();\n    RETURN NEW;\nEND\n$function$",
    )
    op.create_entity(obs_update_ts)

    # Create timestamp update triggers for all tables with 'ts' column
    tables_with_ts = [
        "grids",
        "sites",
        "surveys",
        "users",
        "grid_tiles",
        "telescopes",
        "targets",
        "exposure_sets",
        "strategies",
        "time_blocks",
        "pointings",
    ]
    for table in tables_with_ts:
        update_ts_trigger = PGTrigger(
            schema="obs",
            signature=f"trig_update_ts_{table}",
            on_entity=f"obs.{table}",
            is_constraint=False,
            definition=f"BEFORE UPDATE ON obs.{table} FOR EACH ROW EXECUTE FUNCTION obs.update_ts()",
        )
        op.create_entity(update_ts_trigger)


def downgrade() -> None:
    """Downgrade schema."""
    # Drop the triggers and functions first
    tables_with_ts = [
        "grids",
        "sites",
        "surveys",
        "users",
        "grid_tiles",
        "telescopes",
        "targets",
        "exposure_sets",
        "strategies",
        "time_blocks",
        "pointings",
    ]
    for table in reversed(tables_with_ts):
        op.drop_entity(
            PGTrigger(
                schema="obs",
                signature=f"trig_update_ts_{table}",
                on_entity=f"obs.{table}",
                is_constraint=False,
                definition="",  # definition is not needed for drop
            )
        )
    op.drop_entity(
        PGFunction(
            schema="obs",
            signature="update_ts()",
            definition="",  # definition is not needed for drop
        )
    )

    # Drop tables in reverse order of creation

    # Pointings table
    op.drop_index(
        op.f("ix_obs_pointings_time_block_id"), table_name="pointings", schema="obs"
    )
    op.drop_index(
        op.f("ix_obs_pointings_telescope_id"), table_name="pointings", schema="obs"
    )
    op.drop_index(
        op.f("ix_obs_pointings_target_id"), table_name="pointings", schema="obs"
    )
    op.drop_index(
        op.f("ix_obs_pointings_strategy_id"), table_name="pointings", schema="obs"
    )
    op.drop_index(
        op.f("ix_obs_pointings_stop_time"), table_name="pointings", schema="obs"
    )
    op.drop_index(
        op.f("ix_obs_pointings_start_time"), table_name="pointings", schema="obs"
    )
    op.drop_table("pointings", schema="obs")

    # TimeBlocks table
    op.drop_index(
        op.f("ix_obs_time_blocks_strategy_id"), table_name="time_blocks", schema="obs"
    )
    op.drop_index(
        op.f("ix_obs_time_blocks_block_num"), table_name="time_blocks", schema="obs"
    )
    op.drop_table("time_blocks", schema="obs")

    # Strategies table
    op.drop_index(
        op.f("ix_obs_strategies_target_id"), table_name="strategies", schema="obs"
    )
    op.drop_table("strategies", schema="obs")

    # ExposureSets table
    op.drop_index(
        op.f("ix_obs_exposure_sets_target_id"), table_name="exposure_sets", schema="obs"
    )
    op.drop_table("exposure_sets", schema="obs")

    # Targets table
    op.drop_index(op.f("ix_obs_targets_user_id"), table_name="targets", schema="obs")
    op.drop_index(op.f("ix_obs_targets_survey_id"), table_name="targets", schema="obs")
    op.drop_index(op.f("ix_obs_targets_stop_time"), table_name="targets", schema="obs")
    op.drop_index(op.f("ix_obs_targets_start_time"), table_name="targets", schema="obs")
    op.drop_index(
        op.f("ix_obs_targets_grid_tile_id"), table_name="targets", schema="obs"
    )
    op.drop_table("targets", schema="obs")

    # Telescopes table
    op.drop_index(
        op.f("ix_obs_telescopes_site_id"), table_name="telescopes", schema="obs"
    )
    op.drop_index(op.f("ix_obs_telescopes_name"), table_name="telescopes", schema="obs")
    op.drop_index(
        op.f("ix_obs_telescopes_grid_id"), table_name="telescopes", schema="obs"
    )
    op.drop_table("telescopes", schema="obs")

    # GridTiles table
    op.drop_index(op.f("ix_obs_grid_tiles_name"), table_name="grid_tiles", schema="obs")
    op.drop_index(
        op.f("ix_obs_grid_tiles_grid_id"), table_name="grid_tiles", schema="obs"
    )
    op.drop_table("grid_tiles", schema="obs")

    # Users table
    op.drop_table("users", schema="obs")

    # Surveys table
    op.drop_index(op.f("ix_obs_surveys_name"), table_name="surveys", schema="obs")
    op.drop_table("surveys", schema="obs")

    # Sites table
    op.drop_index(op.f("ix_obs_sites_name"), table_name="sites", schema="obs")
    op.drop_table("sites", schema="obs")

    # Grids table
    op.drop_index(op.f("ix_obs_grids_name"), table_name="grids", schema="obs")
    op.drop_table("grids", schema="obs")

    # Finally, drop the obs schema
    op.execute("DROP SCHEMA obs CASCADE")
