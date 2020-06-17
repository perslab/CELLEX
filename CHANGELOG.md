# Changelog
All notable changes to this project will be documented in this file.
The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]
### Added
* `utils.mapping.human_symbol_to_human_ens()`
### Changed
* anova no longer uses multiprocessing. This should decrease memory usage significantly.
* es_mu and es_sd now have index column name: "gene".
* mapping module naming and docs
### Deprecated
### Removed
### Fixed

## [1.1.1] – 2020-03-09
### Fixed
* FileNotFoundError for mapping functions.

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
