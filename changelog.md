# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).

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
  - ul_lands_download() from Texas Unversity Lands
  - kgs_download() from Kansas Geologic Society
- package added to PyPI registry
