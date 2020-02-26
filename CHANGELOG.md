# Changelog
All notable changes to this project will be documented in this file.
The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]
### Added
### Changed
### Deprecated
### Removed
### Fixed

## [1.1.0] – 2020-02-27
### Added
* Function to parse input, `utils.parse_input()`
### Changed
* ESObject casts input data to float32 to decrease memory load. This behaviour is controlled via the `dtype` argument.
* Updated README.
* Updated demo.
* Updated various docs.
### Deprecated
### Removed
### Fixed
* Duplicate cell id's cause crash.
* Mapping.

## [1.0.0] - 2019-09-07
First stable release for computing cell type expression specificity.
