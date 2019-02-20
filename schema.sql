-- GOTO ObsDB Database Schema
-- run this script in MariaDB using
-- "SOURCE path/to/script/goto-obsdb/schema.sql;"
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

-- Events table
DROP TABLE IF EXISTS `events` ;
CREATE TABLE `events` (
  `eventID` INT NOT NULL AUTO_INCREMENT PRIMARY KEY,
  `name` VARCHAR(255) NOT NULL,
  `source` VARCHAR(255) NOT NULL COMMENT 'LIGO, SWIFT etc.',
  `ivo` VARCHAR(255) NOT NULL UNIQUE,
  `skymap` VARCHAR(255) NULL
  )
  ENGINE = InnoDB;

-- Users table
DROP TABLE IF EXISTS `users` ;
CREATE TABLE IF NOT EXISTS `users` (
  `userKey` INT(11) NOT NULL AUTO_INCREMENT PRIMARY KEY,
  `user_name` VARCHAR(255) NOT NULL UNIQUE,
  `password` VARCHAR(255) NOT NULL,
  `fullName` TEXT NOT NULL
  )
  ENGINE = InnoDB
  DEFAULT CHARACTER SET = latin1;

-- Surveys table
DROP TABLE IF EXISTS `surveys` ;
CREATE TABLE `surveys` (
  `surveyID` INT NOT NULL AUTO_INCREMENT PRIMARY KEY,
  `name` VARCHAR(255) NOT NULL
  )
  ENGINE = InnoDB;

-- Survey tiles table
DROP TABLE IF EXISTS `survey_tiles` ;
CREATE TABLE `survey_tiles` (
  `tileID` INT NOT NULL AUTO_INCREMENT PRIMARY KEY,
  `name` VARCHAR(255) NOT NULL,
  `ra` FLOAT NOT NULL COMMENT 'decimal degrees',
  `decl` FLOAT NOT NULL COMMENT 'decimal degrees',
  `surveys_surveyID` INT NOT NULL,
  INDEX `fk_survey_tiles_surveys1_idx` (`surveys_surveyID`),
  CONSTRAINT `fk_survey_tiles_surveys1`
    FOREIGN KEY (`surveys_surveyID`)
    REFERENCES `surveys` (`surveyID`)
  )
  ENGINE = InnoDB;

-- Event tiles table
DROP TABLE IF EXISTS `event_tiles` ;
CREATE TABLE`event_tiles` (
  `tileID` INT NOT NULL AUTO_INCREMENT PRIMARY KEY,
  `ra` FLOAT NOT NULL,
  `decl` FLOAT NOT NULL,
  `probability` FLOAT NOT NULL,
  `unobserved_probability` FLOAT NOT NULL,
  `events_eventID` INT NOT NULL,
  `survey_tiles_tileID` INT NULL,
  INDEX `fk_event_tiles_events1_idx` (`events_eventID`),
  INDEX `fk_event_tiles_survey_tiles1_idx` (`survey_tiles_tileID`),
  CONSTRAINT `fk_event_tiles_events1`
    FOREIGN KEY (`events_eventID`)
    REFERENCES `events` (`eventID`),
  CONSTRAINT `fk_event_tiles_survey_tiles1`
    FOREIGN KEY (`survey_tiles_tileID`)
    REFERENCES `survey_tiles` (`tileID`)
  )
  ENGINE = InnoDB;

-- Mpointings table
DROP TABLE IF EXISTS `mpointings` ;
CREATE TABLE `mpointings` (
  `mpointingID` INT NOT NULL AUTO_INCREMENT PRIMARY KEY,
  `status` ENUM('unscheduled', 'scheduled', 'completed', 'aborted', 'expired', 'deleted') NOT NULL DEFAULT 'unscheduled',
  `object` TEXT NOT NULL,
  `ra` FLOAT NOT NULL COMMENT 'decimal degrees',
  `decl` FLOAT NOT NULL COMMENT 'decimal degrees',
  `rank` INT(11) NOT NULL,
  `start_rank` INT(11) NOT NULL,
  `minAlt` FLOAT NOT NULL,
  `maxSunAlt` FLOAT NOT NULL DEFAULT -15,
  `minTime` FLOAT NOT NULL,
  `maxMoon` CHAR(1) NOT NULL,
  `minMoonSep` FLOAT NOT NULL DEFAULT 30 COMMENT 'degrees',
  `ToO` TINYINT(1) NOT NULL DEFAULT 0,
  `infinite` TINYINT(1) NOT NULL DEFAULT 0,
  `startUTC` DATETIME NOT NULL DEFAULT CURRENT_TIMESTAMP COMMENT 'only works on mysql later than 5.6.5',
  `stopUTC` DATETIME NULL COMMENT 'If Null then the Mpointing will continue until it is complete',
  `num_todo` INT NOT NULL,
  `num_completed` INT NOT NULL DEFAULT 0,
  `users_userKey` INT(11) NOT NULL,
  `surveys_surveyID` INT NULL,
  `survey_tiles_tileID` INT NULL,
  `events_eventID` INT NULL,
  `event_tiles_tileID` INT NULL,
  INDEX `fk_mpointing_events1_idx` (`events_eventID`),
  INDEX `fk_mpointing_users1_idx` (`users_userKey`),
  INDEX `fk_mpointings_survey_tiles1_idx` (`survey_tiles_tileID`),
  INDEX `fk_mpointings_event_tiles1_idx` (`event_tiles_tileID`),
  INDEX `fk_mpointings_surveys1_idx` (`surveys_surveyID`),
  INDEX `status_idx` (`status`),
  INDEX `startUTC_idx` (`startUTC`),
  INDEX `stopUTC_idx` (`stopUTC`),
  CONSTRAINT `fk_mpointing_events1`
    FOREIGN KEY (`events_eventID`)
    REFERENCES `events` (`eventID`),
  CONSTRAINT `fk_mpointing_users1`
    FOREIGN KEY (`users_userKey`)
    REFERENCES `users` (`userKey`),
  CONSTRAINT `fk_mpointings_survey_tiles1`
    FOREIGN KEY (`survey_tiles_tileID`)
    REFERENCES `survey_tiles` (`tileID`),
  CONSTRAINT `fk_mpointings_event_tiles1`
    FOREIGN KEY (`event_tiles_tileID`)
    REFERENCES `event_tiles` (`tileID`),
  CONSTRAINT `fk_mpointings_surveys1`
    FOREIGN KEY (`surveys_surveyID`)
    REFERENCES `surveys` (`surveyID`)
  )
  ENGINE = InnoDB;

-- Observing blocks table
DROP TABLE IF EXISTS `observing_blocks` ;
CREATE TABLE `observing_blocks` (
  `blockID` INT NOT NULL AUTO_INCREMENT PRIMARY KEY,
  `blockNum` INT NOT NULL,
  `current` INT NOT NULL DEFAULT 0,
  `valid_time` FLOAT NOT NULL COMMENT 'how long after the startUTC the pointing should be valid for in minutes',
  `wait_time` FLOAT NOT NULL COMMENT 'time to wait after this pointing before scheduling the next',
  `mpointings_mpointingID` INT NOT NULL,
  INDEX `fk_observing_blocks_mpointing1_idx` (`mpointings_mpointingID`),
  INDEX `blockNum_idx` (`blockNum`),
  INDEX `current_idx` (`current`),
  CONSTRAINT `fk_observing_blocks_mpointing1`
    FOREIGN KEY (`mpointings_mpointingID`)
    REFERENCES `mpointings` (`mpointingID`)
  )
  ENGINE = InnoDB;

-- Pointings table
DROP TABLE IF EXISTS `pointings` ;
CREATE TABLE `pointings` (
  `pointingID` INT(24) NOT NULL AUTO_INCREMENT PRIMARY KEY,
  `status` ENUM('pending', 'running', 'completed', 'aborted', 'interrupted', 'expired', 'deleted') NOT NULL DEFAULT 'pending',
  `object` TEXT NOT NULL COMMENT 'object name',
  `ra` FLOAT NOT NULL COMMENT 'in decimal degrees',
  `decl` FLOAT NOT NULL COMMENT 'in decimal degrees',
  `rank` INT(1) UNSIGNED NOT NULL,
  `minAlt` FLOAT NOT NULL,
  `maxSunAlt` FLOAT NOT NULL DEFAULT -15 COMMENT 'degrees  ',
  `minTime` FLOAT NOT NULL,
  `maxMoon` CHAR(1) NOT NULL,
  `minMoonSep` FLOAT NOT NULL DEFAULT 30 COMMENT 'degrees',
  `ToO` TINYINT(1) UNSIGNED NOT NULL,
  `startUTC` DATETIME NOT NULL,
  `stopUTC` DATETIME NULL COMMENT 'If Null then the pointing will never expire, and will remain until observed',
  `startedUTC` DATETIME NULL,
  `stoppedUTC` DATETIME NULL,
  `ts` TIMESTAMP(3) NOT NULL DEFAULT CURRENT_TIMESTAMP(3) ON UPDATE CURRENT_TIMESTAMP(3),
  `users_userKey` INT(11) NOT NULL,
  `mpointings_mpointingID` INT NULL,
  `observing_blocks_blockID` INT NULL,
  `surveys_surveyID` INT NULL,
  `survey_tiles_tileID` INT NULL,
  `events_eventID` INT NULL,
  `event_tiles_tileID` INT NULL,
  INDEX `fk_pointings_events1_idx` (`events_eventID`),
  INDEX `fk_pointings_users1_idx` (`users_userKey`),
  INDEX `fk_pointings_observing_blocks1_idx` (`observing_blocks_blockID`),
  INDEX `fk_pointings_mpointings1_idx` (`mpointings_mpointingID`),
  INDEX `status_idx` (`status`),
  INDEX `startUTC_idx` (`startUTC`),
  INDEX `stopUTC_idx` (`stopUTC`),
  INDEX `fk_pointings_event_tiles1_idx` (`event_tiles_tileID`),
  INDEX `fk_pointings_survey_tiles1_idx` (`survey_tiles_tileID`),
  INDEX `fk_pointings_surveys1_idx` (`surveys_surveyID`),
  CONSTRAINT `fk_pointings_events1`
    FOREIGN KEY (`events_eventID`)
    REFERENCES `events` (`eventID`),
  CONSTRAINT `fk_pointings_users1`
    FOREIGN KEY (`users_userKey`)
    REFERENCES `users` (`userKey`),
  CONSTRAINT `fk_pointings_observing_blocks1`
    FOREIGN KEY (`observing_blocks_blockID`)
    REFERENCES `observing_blocks` (`blockID`),
  CONSTRAINT `fk_pointings_mpointings1`
    FOREIGN KEY (`mpointings_mpointingID`)
    REFERENCES `mpointings` (`mpointingID`),
  CONSTRAINT `fk_pointings_event_tiles1`
    FOREIGN KEY (`event_tiles_tileID`)
    REFERENCES `event_tiles` (`tileID`),
  CONSTRAINT `fk_pointings_survey_tiles1`
    FOREIGN KEY (`survey_tiles_tileID`)
    REFERENCES `survey_tiles` (`tileID`),
  CONSTRAINT `fk_pointings_surveys1`
    FOREIGN KEY (`surveys_surveyID`)
    REFERENCES `surveys` (`surveyID`)
  )
  ENGINE = InnoDB
  DEFAULT CHARACTER SET = latin1;

-- Exposure sets table
DROP TABLE IF EXISTS `exposure_sets` ;
CREATE TABLE `exposure_sets` (
  `expID` INT(24) NOT NULL AUTO_INCREMENT PRIMARY KEY,
  `utMask` INT NULL COMMENT 'bit mask to allocate to individual UTs. NULL means send to all',
  `typeFlag` ENUM('SCIENCE', 'FOCUS', 'DARK', 'BIAS', 'FLAT', 'STD') NOT NULL,
  `filter` CHAR(2) NOT NULL,
  `exptime` FLOAT NOT NULL,
  `binning` INT(11) UNSIGNED NOT NULL,
  `numexp` INT(11) UNSIGNED NOT NULL,
  `raoff` FLOAT NOT NULL DEFAULT 0.0 COMMENT 'RA offset (arcsecs)',
  `decoff` FLOAT NOT NULL DEFAULT 0.0 COMMENT 'dec offset (arcsec)',
  `pointings_pointingID` INT NULL,
  `mpointings_mpointingID` INT NULL,
  INDEX `fk_exposures_pointings1_idx` (`pointings_pointingID`),
  INDEX `fk_exposures_mpointing1_idx` (`mpointings_mpointingID`),
  CONSTRAINT `fk_exposures_pointings1`
    FOREIGN KEY (`pointings_pointingID`)
    REFERENCES `pointings` (`pointingID`),
  CONSTRAINT `fk_exposures_mpointing1`
    FOREIGN KEY (`mpointings_mpointingID`)
    REFERENCES `mpointings` (`mpointingID`)
  )
  ENGINE = InnoDB
  DEFAULT CHARACTER SET = latin1;

-- Image logs table
DROP TABLE IF EXISTS `image_logs` ;
CREATE TABLE `image_logs` (
  `logID` INT(24) NOT NULL AUTO_INCREMENT PRIMARY KEY,
  `filename` VARCHAR(30) NOT NULL UNIQUE COMMENT 'full FITS file name, including extension',
  `runNumber` INT NOT NULL,
  `ut` INT NOT NULL,
  `utMask` INT NOT NULL,
  `startUTC` DATETIME NOT NULL,
  `writeUTC` DATETIME NOT NULL,
  `set_position` INT NOT NULL DEFAULT 1,
  `set_total` INT NOT NULL DEFAULT 1,
  `exposure_sets_expID` INT NULL DEFAULT NULL,
  `pointings_pointingID` INT NULL DEFAULT NULL,
  `mpointings_mpointingID` INT NULL DEFAULT NULL,
  INDEX `fk_image_logs_exposure_sets1_idx` (`exposure_sets_expID`),
  INDEX `fk_image_logs_pointings1_idx` (`pointings_pointingID`),
  INDEX `fk_image_logs_mpointings1_idx` (`mpointings_mpointingID`),
  INDEX `runNumber_idx` (`runNumber`),
  INDEX `writeUTC_idx` (`writeUTC`),
  CONSTRAINT `fk_image_logs_exposure_sets1`
    FOREIGN KEY (`exposure_sets_expID`)
    REFERENCES `exposure_sets` (`expID`),
  CONSTRAINT `fk_image_logs_pointings1`
    FOREIGN KEY (`pointings_pointingID`)
    REFERENCES `pointings` (`pointingID`),
  CONSTRAINT `fk_image_logs_mpointings1`
    FOREIGN KEY (`mpointings_mpointingID`)
    REFERENCES `mpointings` (`mpointingID`)
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

-- Event tiles triggers
-- Before insert
DROP TRIGGER IF EXISTS `event_tiles_BEFORE_INSERT` $$
CREATE DEFINER = CURRENT_USER TRIGGER `event_tiles_BEFORE_INSERT` BEFORE INSERT ON `event_tiles` FOR EACH ROW
  BEGIN
    -- Select RA and Dec from linked survey tile if not given
    IF ((NEW.survey_tiles_tileID is not NULL) and (NEW.ra is NULL) and (NEW.decl is NULL)) THEN
      SET NEW.ra = (SELECT ra FROM `survey_tiles` WHERE NEW.survey_tiles_tileID = `survey_tiles`.`tileID`);
      SET NEW.decl = (SELECT decl FROM `survey_tiles` WHERE NEW.survey_tiles_tileID = `survey_tiles`.`tileID`);
    END IF;
    -- Set unobserved_probability from probability
    IF ((NEW.unobserved_probability is NULL) and (NEW.probability is not NULL)) THEN
      SET NEW.unobserved_probability = NEW.probability;
    END IF;
  END$$


-- Mpointings triggers
-- Before insert
DROP TRIGGER IF EXISTS `mpointings_BEFORE_INSERT` $$
CREATE DEFINER = CURRENT_USER TRIGGER `mpointings_BEFORE_INSERT` BEFORE INSERT ON `mpointings` FOR EACH ROW
  BEGIN
    -- Select RA and Dec from linked survey tile if not given
    IF ((NEW.survey_tiles_tileID is not NULL) and (NEW.ra is NULL) and (NEW.decl is NULL)) THEN
      SET NEW.ra = (SELECT ra FROM `survey_tiles` WHERE NEW.survey_tiles_tileID = `survey_tiles`.`tileID`);
      SET NEW.decl = (SELECT decl FROM `survey_tiles` WHERE NEW.survey_tiles_tileID = `survey_tiles`.`tileID`);
   END IF;
  END$$


-- Before update
DROP TRIGGER IF EXISTS `mpointings_BEFORE_UPDATE` $$
CREATE DEFINER = CURRENT_USER TRIGGER `mpointings_BEFORE_UPDATE` BEFORE UPDATE ON `mpointings` FOR EACH ROW
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
CREATE DEFINER = CURRENT_USER TRIGGER `pointings_BEFORE_INSERT` BEFORE INSERT ON `pointings` FOR EACH ROW
  BEGIN
    -- Select RA and Dec from linked survey tile if not given
    IF ((NEW.survey_tiles_tileID is not NULL) and (NEW.ra is NULL) and (NEW.decl is NULL)) THEN
      SET NEW.ra = (SELECT ra FROM `survey_tiles` WHERE NEW.survey_tiles_tileID = `survey_tiles`.`tileID`);
      SET NEW.decl = (SELECT decl FROM `survey_tiles` WHERE NEW.survey_tiles_tileID = `survey_tiles`.`tileID`);
    END IF;
  END$$


-- After insert
DROP TRIGGER IF EXISTS `pointings_AFTER_INSERT` $$
CREATE DEFINER = CURRENT_USER TRIGGER `pointings_AFTER_INSERT` AFTER INSERT ON `pointings` FOR EACH ROW
  BEGIN
    -- Mark all other blocks for this Mpointing as current=False,
    -- and the one for this Pointing as current=True
    IF (NEW.`observing_blocks_blockID` is not NULL) AND (NEW.status = 'pending') THEN
      UPDATE `observing_blocks` SET `current` = 0 WHERE (NEW.`mpointings_mpointingID` = `observing_blocks`.`mpointings_mpointingID`);
      UPDATE `observing_blocks` SET `current` = 1 WHERE (NEW.`observing_blocks_blockID` = `observing_blocks`.`blockID`);
    END IF;
    -- Mark any linked Mpointing as scheduled
    IF (NEW.`mpointings_mpointingID` is not NULL) AND (NEW.status = 'pending') THEN
      UPDATE `mpointings` SET `status` = 'scheduled' WHERE (NEW.`mpointings_mpointingID` = `mpointings`.`mpointingID`);
    END IF;
  END$$


-- Before update
DROP TRIGGER IF EXISTS `pointings_BEFORE_UPDATE` $$
CREATE DEFINER = CURRENT_USER TRIGGER `pointings_BEFORE_UPDATE` BEFORE UPDATE ON `pointings` FOR EACH ROW
  BEGIN
    -- Store time when set running
    IF (OLD.`status` != 'running' AND NEW.`status` = 'running') THEN
      SET NEW.`startedUTC` = UTC_TIMESTAMP();
    END IF;
    -- Store time when finished (completed, aborted, interrupted, expired...)
    IF (OLD.`status` IN ('pending', 'running') AND NEW.`status` NOT IN ('pending', 'running')) THEN
      SET NEW.`stoppedUTC` = UTC_TIMESTAMP();
    END IF;
  END$$

-- After update
DROP TRIGGER IF EXISTS `pointings_AFTER_UPDATE` $$
CREATE DEFINER = CURRENT_USER TRIGGER `pointings_AFTER_UPDATE` AFTER UPDATE ON `pointings` FOR EACH ROW
  BEGIN
    DECLARE isinfinite INT;
    -- Only trigger if the timestamp has changed
    IF (NEW.`ts` <> OLD.`ts`) THEN
      -- Mark Mpointing as unscheduled when the pointing is finished somehow
      -- NB only if the Mpointing is scheduled
      IF NEW.`status` NOT IN ('pending', 'running') THEN
        UPDATE `mpointings` SET `status` = 'unscheduled' WHERE (`mpointings`.`mpointingID` = NEW.`mpointings_mpointingID` and `mpointings`.`status` = 'scheduled');
      END IF;
      -- If the pointing was completed...
      IF NEW.`status` = 'completed' THEN
        -- Increase the Mpointing's completed count
        UPDATE `mpointings` SET `num_completed` = `num_completed` + 1 WHERE (NEW.`mpointings_mpointingID` = `mpointings`.`mpointingID`);
        -- Add 10 to the rank
        -- NB only if the Mpointing is not infinite
        SELECT `infinite` INTO isinfinite FROM `mpointings` WHERE (NEW.`mpointings_mpointingID` = `mpointings`.`mpointingID`);
        IF isinfinite = 0 THEN
          UPDATE `mpointings` SET `rank` = `rank` + 10 WHERE (NEW.`mpointings_mpointingID` = `mpointings`.`mpointingID`);
        END IF;
      END IF;
    END IF;
  END$$


DELIMITER ;
