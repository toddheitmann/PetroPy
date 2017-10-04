<div align="center">
  <img src="https://github.com/toddheitmann/PetroPy/blob/master/petropy_logo.png"><br>
</div>

# PetroPy

A petrophysics with python package allowing scientific python computing of conventional and unconventional formation evaluation. Reads las files using [lasio](https://github.com/kinverarity1/lasio). Includes a petrophysical workflow and a log viewer based on XML templates.

<div align="center">
  <img src="https://github.com/toddheitmann/PetroPy/blob/master/university_6-17_no1.png"><br>
</div>

## Requirements

- [lasio](https://github.com/kinverarity1/lasio)
- [numpy](http://www.numpy.org)
- [scipy](https://www.scipy.org)
- [pandas](http://pandas.pydata.org)
- [matplotlib](http://matplotlib.org)
- [scikit-learn](http://scikit-learn.org/stable/)

## Currently Version 0.0

The basic workflow in the module is functional, but requires downloading source code and importing .py files.

## TO DO version 0.1

### Features
- [x] Read las files
- [x] Calculate fluid properties
- [x] Calculate multimineral, posority, saturation model
- [x] Calculate and export statistics
- [x] Curve edit manual redrawing
- [x] Curve edit bulk shift data
- [x] Electrofacies module
- [x] Replace log object by subclassing [lasio](https://github.com/kinverarity1/lasio) object
- [x] Add to pypi package registry

### Examples

- [x] Sphinx Documentation
- [ ] Wolfcamp example from University Lands
- [ ] Mississippi Limestone Example from KGS
- [ ] Dataframe manipulation for different petrophysical crossplots
- [ ] Template creation for log viewer

## TO DO version 0.2

### Features
- [ ] Histogram module for normalization

## TO DO version 0.3

### Extended Features with Geologic Software

- [ ] Create Geographix GGXLog object
  - [ ] Read Well with UWI input
  - [ ] Save calculations to Database

- [ ] Create Petra PetraLog object
  - [ ] Read Well with UWI input
  - [ ] Save calculations to Database
