-- GOTO ObsDB Database Schema
-- run this script in MariaDB using
-- "SOURCE path/to/script/schema.sql;"
-- NB this will drop the current database if it exists, so be careful!

-- -----------------------------------------------------
-- https://dev.mysql.com/doc/workbench/en/workbench-faq.html
SET @OLD_UNIQUE_CHECKS=@@UNIQUE_CHECKS, UNIQUE_CHECKS=0;
SET @OLD_FOREIGN_KEY_CHECKS=@@FOREIGN_KEY_CHECKS, FOREIGN_KEY_CHECKS=0;
SET @OLD_SQL_MODE=@@SQL_MODE, SQL_MODE='TRADITIONAL,ALLOW_INVALID_DATES';

-- -----------------------------------------------------
-- Create database
DROP DATABASE IF EXISTS `goto_obs`;
CREATE DATABASE `goto_obs` DEFAULT CHARACTER SET utf8 ;
USE `goto_obs` ;

-- -----------------------------------------------------
-- Create tables

-- Users table
DROP TABLE IF EXISTS `users` ;
CREATE TABLE IF NOT EXISTS `users` (
  -- primary key
  `id` INT NOT NULL AUTO_INCREMENT PRIMARY KEY,
  -- columns
  `username` VARCHAR(255) NOT NULL UNIQUE,
  `password` VARCHAR(255) NOT NULL,
  `full_name` TEXT NOT NULL
  -- indexes
  -- foreign keys
  )
  ENGINE = InnoDB
  DEFAULT CHARACTER SET = latin1;

-- Pointings table
DROP TABLE IF EXISTS `pointings` ;
CREATE TABLE `pointings` (
  -- primary key
  `id` INT NOT NULL AUTO_INCREMENT PRIMARY KEY,
  -- columns
  `status` ENUM('pending', 'running', 'completed', 'aborted', 'interrupted', 'expired', 'deleted') NOT NULL DEFAULT 'pending',
  `object` TEXT NOT NULL COMMENT 'object name',
  `ra` FLOAT NOT NULL COMMENT 'in decimal degrees',
  `decl` FLOAT NOT NULL COMMENT 'in decimal degrees',
  `rank` INT(1) UNSIGNED NOT NULL,
  `min_alt` FLOAT NOT NULL,
  `max_sunalt` FLOAT NOT NULL DEFAULT -15 COMMENT 'degrees',
  `min_time` FLOAT NOT NULL,
  `max_moon` CHAR(1) NOT NULL,
  `min_moonsep` FLOAT NOT NULL DEFAULT 30 COMMENT 'degrees',
  `too` BOOLEAN NOT NULL DEFAULT FALSE,
  `start_time` DATETIME NOT NULL,
  `stop_time` DATETIME NULL COMMENT 'If NULL then the pointing will never expire, and will remain until observed',
  `started_time` DATETIME NULL,
  `stopped_time` DATETIME NULL,
  `ts` TIMESTAMP(3) NOT NULL DEFAULT CURRENT_TIMESTAMP(3) ON UPDATE CURRENT_TIMESTAMP(3),
  -- indexes
  INDEX `status_idx` (`status`),
  INDEX `start_time_idx` (`start_time`),
  INDEX `stop_time_idx` (`stop_time`),
  -- foreign keys
  `user_id` INT NOT NULL,
  `mpointing_id` INT NULL,
  `time_block_id` INT NULL,
  `grid_tile_id` INT NULL,
  `survey_tile_id` INT NULL,
  `event_id` INT NULL,
  FOREIGN KEY (`user_id`) REFERENCES `users` (`id`),
  FOREIGN KEY (`mpointing_id`) REFERENCES `mpointings` (`id`),
  FOREIGN KEY (`time_block_id`) REFERENCES `time_blocks` (`id`),
  FOREIGN KEY (`grid_tile_id`) REFERENCES `grid_tiles` (`id`),
  FOREIGN KEY (`survey_tile_id`) REFERENCES `survey_tiles` (`id`),
  FOREIGN KEY (`event_id`) REFERENCES `events` (`id`)
  )
  ENGINE = InnoDB
  DEFAULT CHARACTER SET = latin1;

-- Exposure sets table
DROP TABLE IF EXISTS `exposure_sets` ;
CREATE TABLE `exposure_sets` (
  -- primary key
  `id` INT NOT NULL AUTO_INCREMENT PRIMARY KEY,
  -- columns
  `num_exp` INT UNSIGNED NOT NULL,
  `exptime` FLOAT NOT NULL,
  `filter` CHAR(1) NOT NULL,
  `binning` INT UNSIGNED NOT NULL,
  `imgtype` ENUM('SCIENCE', 'FOCUS', 'DARK', 'BIAS', 'FLAT', 'STD') NOT NULL,
  `ut_mask` INT NULL COMMENT 'bit mask to allocate to individual UTs, NULL means send to all',
  `ra_offset` FLOAT NOT NULL DEFAULT 0.0 COMMENT 'arcsecs',
  `dec_offset` FLOAT NOT NULL DEFAULT 0.0 COMMENT 'arcsec',
  -- indexes
  -- foreign keys
  `pointing_id` INT NULL,
  `mpointing_id` INT NULL,
  FOREIGN KEY (`pointing_id`) REFERENCES `pointings` (`id`),
  FOREIGN KEY (`mpointing_id`) REFERENCES `mpointings` (`id`)
  )
  ENGINE = InnoDB
  DEFAULT CHARACTER SET = latin1;

-- Mpointings table
DROP TABLE IF EXISTS `mpointings` ;
CREATE TABLE `mpointings` (
  -- primary key
  `id` INT NOT NULL AUTO_INCREMENT PRIMARY KEY,
  -- columns
  `status` ENUM('unscheduled', 'scheduled', 'completed', 'aborted', 'expired', 'deleted') NOT NULL DEFAULT 'unscheduled',
  `object` TEXT NOT NULL,
  `ra` FLOAT NOT NULL COMMENT 'decimal degrees',
  `decl` FLOAT NOT NULL COMMENT 'decimal degrees',
  `current_rank` INT NOT NULL,
  `initial_rank` INT NOT NULL,
  `num_todo` INT NOT NULL,
  `num_completed` INT NOT NULL DEFAULT 0,
  `infinite` BOOLEAN NOT NULL DEFAULT FALSE,
  `min_alt` FLOAT NOT NULL,
  `max_sunalt` FLOAT NOT NULL DEFAULT -15,
  `min_time` FLOAT NOT NULL,
  `max_moon` CHAR(1) NOT NULL,
  `min_moonsep` FLOAT NOT NULL DEFAULT 30 COMMENT 'degrees',
  `too` BOOLEAN NOT NULL DEFAULT FALSE,
  `start_time` DATETIME NOT NULL DEFAULT CURRENT_TIMESTAMP COMMENT 'only works on mysql later than 5.6.5',
  `stop_time` DATETIME NULL COMMENT 'If Null then the Mpointing will continue until it is complete',
  -- indexes
  INDEX `status_idx` (`status`),
  INDEX `start_time_idx` (`start_time`),
  INDEX `stop_time_idx` (`stop_time`),
  -- foreign keys
  `user_id` INT NOT NULL,
  `grid_tile_id` INT NULL,
  `survey_tile_id` INT NULL,
  `event_id` INT NULL,
  FOREIGN KEY (`user_id`) REFERENCES `users` (`id`),
  FOREIGN KEY (`grid_tile_id`) REFERENCES `grid_tiles` (`id`),
  FOREIGN KEY (`survey_tile_id`) REFERENCES `survey_tiles` (`id`),
  FOREIGN KEY (`event_id`) REFERENCES `events` (`id`)
  )
  ENGINE = InnoDB;

-- Time blocks table
DROP TABLE IF EXISTS `time_blocks` ;
CREATE TABLE `time_blocks` (
  -- primary key
  `id` INT NOT NULL AUTO_INCREMENT PRIMARY KEY,
  -- columns
  `block_num` INT NOT NULL,
  `valid_time` FLOAT NOT NULL COMMENT 'how long after the start_time the pointing should be valid for in minutes',
  `wait_time` FLOAT NOT NULL COMMENT 'time to wait after this pointing before scheduling the next',
  `current` BOOLEAN NOT NULL DEFAULT FALSE,
  -- indexes
  INDEX `block_num_idx` (`block_num`),
  INDEX `current_idx` (`current`),
  -- foreign keys
  `mpointing_id` INT NOT NULL,
  FOREIGN KEY (`mpointing_id`) REFERENCES `mpointings` (`id`)
  )
  ENGINE = InnoDB;

-- Grids table
DROP TABLE IF EXISTS `grids` ;
CREATE TABLE `grids` (
  -- primary key
  `id` INT NOT NULL AUTO_INCREMENT PRIMARY KEY,
  -- columns
  `name` VARCHAR(255) NOT NULL UNIQUE,
  `ra_fov` FLOAT NOT NULL COMMENT 'decimal degrees',
  `dec_fov` FLOAT NOT NULL COMMENT 'decimal degrees',
  `ra_overlap` FLOAT NOT NULL COMMENT 'fraction (0-1)',
  `dec_overlap` FLOAT NOT NULL COMMENT 'fraction (0-1)',
  `algorithm` VARCHAR(255) NOT NULL,
  -- indexes
  INDEX `name_idx` (`name`)
  -- foreign keys
  )
  ENGINE = InnoDB;

-- Grid tiles table
DROP TABLE IF EXISTS `grid_tiles` ;
CREATE TABLE `grid_tiles` (
  -- primary key
  `id` INT NOT NULL AUTO_INCREMENT PRIMARY KEY,
  -- columns
  `name` VARCHAR(255) NOT NULL,
  `ra` FLOAT NOT NULL COMMENT 'decimal degrees',
  `decl` FLOAT NOT NULL COMMENT 'decimal degrees',
  -- indexes
  INDEX `name_idx` (`name`),
  -- foreign keys
  `grid_id` INT NOT NULL,
  FOREIGN KEY (`grid_id`) REFERENCES `grids` (`id`)
  )
  ENGINE = InnoDB;

-- Surveys table
DROP TABLE IF EXISTS `surveys` ;
CREATE TABLE `surveys` (
  -- primary key
  `id` INT NOT NULL AUTO_INCREMENT PRIMARY KEY,
  -- columns
  `name` VARCHAR(255) NOT NULL,
  -- indexes
  INDEX `name_idx` (`name`),
  -- foreign keys
  `grid_id` INT NOT NULL,
  `event_id` INT NULL,
  FOREIGN KEY (`grid_id`) REFERENCES `grids` (`id`),
  FOREIGN KEY (`event_id`) REFERENCES `events` (`id`)
  )
  ENGINE = InnoDB;

-- Survey tiles table
DROP TABLE IF EXISTS `survey_tiles` ;
CREATE TABLE `survey_tiles` (
  -- primary key
  `id` INT NOT NULL AUTO_INCREMENT PRIMARY KEY,
  -- column
  `current_weight` FLOAT NOT NULL COMMENT '0-1',
  `initial_weight` FLOAT NOT NULL COMMENT '0-1',
  -- indexes
  -- foreign keys
  `survey_id` INT NOT NULL,
  `grid_tile_id` INT NOT NULL,
  FOREIGN KEY (`survey_id`) REFERENCES `surveys` (`id`),
  FOREIGN KEY (`grid_tile_id`) REFERENCES `grid_tiles` (`id`)
  )
  ENGINE = InnoDB;

-- Events table
DROP TABLE IF EXISTS `events` ;
CREATE TABLE `events` (
  -- primary key
  `id` INT NOT NULL AUTO_INCREMENT PRIMARY KEY,
  -- columns
  `name` VARCHAR(255) NOT NULL,
  `ivorn` VARCHAR(255) NOT NULL UNIQUE,
  `source` VARCHAR(255) NOT NULL COMMENT 'eg LVC, Swift',
  `type` VARCHAR(255) NOT NULL COMMENT 'eg GW, GRB',
  `time` DATETIME NULL COMMENT 'Event time',
  `skymap` VARCHAR(255) NULL COMMENT 'Local path to skymap file',
  -- indexes
  INDEX `name_idx` (`name`),
  INDEX `type_idx` (`type`)
  -- foreign keys
  )
  ENGINE = InnoDB;

-- Image logs table
DROP TABLE IF EXISTS `image_logs` ;
CREATE TABLE `image_logs` (
  -- primary key
  `id` INT NOT NULL AUTO_INCREMENT PRIMARY KEY,
  -- columns
  `filename` VARCHAR(30) NOT NULL UNIQUE COMMENT 'full FITS file name, including extension',
  `run_number` INT NOT NULL,
  `ut` INT NOT NULL,
  `ut_mask` INT NOT NULL,
  `start_time` DATETIME NOT NULL,
  `write_time` DATETIME NOT NULL,
  `set_position` INT NOT NULL DEFAULT 1,
  `set_total` INT NOT NULL DEFAULT 1,
  -- indexes
  INDEX `run_number_idx` (`run_number`),
  INDEX `write_time_idx` (`write_time`),
  -- foreign keys
  `exposure_set_id` INT NULL DEFAULT NULL,
  `pointing_id` INT NULL DEFAULT NULL,
  `mpointing_id` INT NULL DEFAULT NULL,
  FOREIGN KEY (`exposure_set_id`) REFERENCES `exposure_sets` (`id`),
  FOREIGN KEY (`pointing_id`) REFERENCES `pointings` (`id`),
  FOREIGN KEY (`mpointing_id`) REFERENCES `mpointings` (`id`)
  )
  ENGINE = InnoDB
  DEFAULT CHARACTER SET = latin1;


-- -----------------------------------------------------
-- Create user
SET SQL_MODE = '';
GRANT USAGE ON *.* TO goto;
DROP USER goto;
SET SQL_MODE='TRADITIONAL,ALLOW_INVALID_DATES';
CREATE USER 'goto' IDENTIFIED BY 'gotoobs';
GRANT ALL ON `goto_obs`.* TO 'goto';

-- -----------------------------------------------------
SET SQL_MODE=@OLD_SQL_MODE;
SET FOREIGN_KEY_CHECKS=@OLD_FOREIGN_KEY_CHECKS;
SET UNIQUE_CHECKS=@OLD_UNIQUE_CHECKS;

-- -----------------------------------------------------
-- Create triggers
DELIMITER $$

-- Mpointings triggers
-- Before insert
DROP TRIGGER IF EXISTS `mpointings_BEFORE_INSERT` $$
CREATE DEFINER = CURRENT_USER TRIGGER `mpointings_BEFORE_INSERT`
  BEFORE INSERT ON `mpointings` FOR EACH ROW
  BEGIN
    -- Select RA and Dec from linked grid tile if not given
    IF ((NEW.`grid_tile_id` is not NULL) and (NEW.`ra` is NULL) and (NEW.`decl` is NULL)) THEN
      SET NEW.`ra` = (SELECT `ra` FROM `grid_tiles` WHERE NEW.`grid_tile_id` = `grid_tiles`.`id`);
      SET NEW.`decl` = (SELECT `decl` FROM `grid_tiles` WHERE NEW.`grid_tile_id` = `grid_tiles`.`id`);
    END IF;
    -- Set current_rank from initial_rank
    IF ((NEW.`current_rank` is NULL) and (NEW.`initial_rank` is not NULL)) THEN
      SET NEW.`current_rank` = NEW.`initial_rank`;
    END IF;
  END$$

-- Before update
DROP TRIGGER IF EXISTS `mpointings_BEFORE_UPDATE` $$
CREATE DEFINER = CURRENT_USER TRIGGER `mpointings_BEFORE_UPDATE`
  BEFORE UPDATE ON `mpointings` FOR EACH ROW
  BEGIN
    -- Set status to completed if finished
    -- NB infinite Mpointings can never be completed
    IF (NEW.`num_completed` = NEW.`num_todo` and NEW.`infinite` = 0) THEN
      SET NEW.`status` = 'completed';
    END IF;
  END$$


-- Pointings triggers
-- Before insert
DROP TRIGGER IF EXISTS `pointings_BEFORE_INSERT` $$
CREATE DEFINER = CURRENT_USER TRIGGER `pointings_BEFORE_INSERT`
  BEFORE INSERT ON `pointings` FOR EACH ROW
  BEGIN
    -- Select RA and Dec from linked grid tile if not given
    IF ((NEW.`grid_tile_id` is not NULL) and (NEW.`ra` is NULL) and (NEW.`decl` is NULL)) THEN
      SET NEW.`ra` = (SELECT `ra` FROM `grid_tiles` WHERE NEW.`grid_tile_id` = `grid_tiles`.`id`);
      SET NEW.`decl` = (SELECT `decl` FROM `grid_tiles` WHERE NEW.`grid_tile_id` = `grid_tiles`.`id`);
    END IF;
  END$$

-- After insert
DROP TRIGGER IF EXISTS `pointings_AFTER_INSERT` $$
CREATE DEFINER = CURRENT_USER TRIGGER `pointings_AFTER_INSERT`
  AFTER INSERT ON `pointings` FOR EACH ROW
  BEGIN
    -- Mark all other blocks for this Mpointing as current=False,
    -- and the one for this Pointing as current=True
    IF (NEW.`time_block_id` is not NULL) AND (NEW.status = 'pending') THEN
      UPDATE `time_blocks` SET `current` = FALSE WHERE (NEW.`mpointing_id` = `time_blocks`.`mpointing_id`);
      UPDATE `time_blocks` SET `current` = TRUE WHERE (NEW.`time_block_id` = `time_blocks`.`id`);
    END IF;
    -- Mark any linked Mpointing as scheduled
    IF (NEW.`mpointing_id` is not NULL) AND (NEW.status = 'pending') THEN
      UPDATE `mpointings` SET `status` = 'scheduled' WHERE (NEW.`mpointing_id` = `mpointings`.`id`);
    END IF;
  END$$

-- Before update
DROP TRIGGER IF EXISTS `pointings_BEFORE_UPDATE` $$
CREATE DEFINER = CURRENT_USER TRIGGER `pointings_BEFORE_UPDATE`
  BEFORE UPDATE ON `pointings` FOR EACH ROW
  BEGIN
    -- Store time when set running
    IF (OLD.`status` != 'running' AND NEW.`status` = 'running') THEN
      SET NEW.`started_time` = UTC_TIMESTAMP();
    END IF;
    -- Store time when finished (completed, aborted, interrupted, expired...)
    IF (OLD.`status` IN ('pending', 'running') AND NEW.`status` NOT IN ('pending', 'running')) THEN
      SET NEW.`stopped_time` = UTC_TIMESTAMP();
    END IF;
  END$$

-- After update
DROP TRIGGER IF EXISTS `pointings_AFTER_UPDATE` $$
CREATE DEFINER = CURRENT_USER TRIGGER `pointings_AFTER_UPDATE`
  AFTER UPDATE ON `pointings` FOR EACH ROW
  BEGIN
    DECLARE isinfinite INT;
    -- Only trigger if the timestamp has changed
    IF (NEW.`ts` <> OLD.`ts`) THEN
      -- Mark Mpointing as unscheduled when the pointing is finished somehow
      -- NB only if the Mpointing is scheduled
      IF NEW.`status` NOT IN ('pending', 'running') THEN
        UPDATE `mpointings` SET `status` = 'unscheduled' WHERE (`mpointings`.`id` = NEW.`mpointing_id` and `mpointings`.`status` = 'scheduled');
      END IF;
      -- If the pointing was completed...
      IF (OLD.`status` != 'completed' AND NEW.`status` = 'completed') THEN
        -- Increase the Mpointing's completed count
        UPDATE `mpointings` SET `num_completed` = `num_completed` + 1 WHERE (NEW.`mpointing_id` = `mpointings`.`id`);
        -- Add 10 to the current_rank
        -- NB only if the Mpointing is not infinite
        SELECT `infinite` INTO isinfinite FROM `mpointings` WHERE (NEW.`mpointing_id` = `mpointings`.`id`);
        IF isinfinite = 0 THEN
          UPDATE `mpointings` SET `current_rank` = `current_rank` + 10 WHERE (NEW.`mpointing_id` = `mpointings`.`id`);
        END IF;
      END IF;
    END IF;
  END$$

-- Survey tiles triggers
-- Before insert
DROP TRIGGER IF EXISTS `survey_tiles_BEFORE_INSERT` $$
CREATE DEFINER = CURRENT_USER TRIGGER `survey_tiles_BEFORE_INSERT`
  BEFORE INSERT ON `survey_tiles` FOR EACH ROW
  BEGIN
    -- Set current_weight from initial_weight
    IF ((NEW.`current_weight` is NULL) and (NEW.`initial_weight` is not NULL)) THEN
      SET NEW.`current_weight` = NEW.`initial_weight`;
    END IF;
  END$$

DELIMITER ;
