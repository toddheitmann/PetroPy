# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project generally adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).

## [0.1.7] - 2019-

### Added

- Use TkAgg backend in matplotlib for consistent performance
- Download 2017 and 2018 logs in [kgs_download](https://toddheitmann.github.io/PetroPy/function/download.html#petropy.kgs_download) function
- Facies and pay flag input for statistics export.

### Fixed

- [Issue 1](https://github.com/toddheitmann/PetroPy/issues/1) calling newer methods from old versions of required packages. Versions added to requirements.txt.
- Issue returning empty depth array when trying to find next_formation_depth
- Issue not removing rows from calculations when only PE is null causing a crash

## [0.1.6] - 2018-02-15

### Fixed

- Issue with show method for LogViewer class using different back ends of matplotlib

## [0.1.5] - 2017-10-24

### Added

- Gradients for filled curves
- Inventory las files from given folders
- Wolfcamp examples
- autodetect_encoding to LASFile requiring cchardet on install

### Fixed

- Downloading las data from KGS and University Lands
- Electrofacies calculations

## [0.1.4] - 2017-10-04

### Fixed

- Graphs error with window location

## [0.1.3] - 2017-10-04

### Fixed

- Data included in distribution

## [0.1.2] - 2017-10-04

### Fixed

- Wheel source distribution

## [0.1.0] - 2017-10-04

### Added

- Log object
  - Subclass of [lasio](https://github.com/kinverarity1/lasio) LASFile
  - Calculates petrophysical properties
  - Exports statistics
- LogViewer object
  - displaying log data
  - interacting with log data
    - manually manipulate and redraw curves
    - bulk shift curve data
- Electrofacies function for digital rock types
- Dataset function to read packaged data
- Download functions to download public data
  - ul_lands_download() from Texas University Lands
  - kgs_download() from Kansas Geologic Society
- package added to PyPI registry
