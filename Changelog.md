# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).

## 1.0.1 - 2018-03-22 
### Added
- A-weighting filter to Filter
- JavaDoc comments to Utilities class
- a changelog

### Changed
- Power and Magnitude calculations in FFT to single sided form, divide by 2 and mirror the array to 
get back to 2 sided versions
- FormantExtractor now takes number of formants to extract as a constructor argument (before it was 
hardcoded to 3, now it will default to 4 if not specified)
- index.html was changed to add a few small tweaks to the web page

### Fixed 
- MFCC returning NaN's if filter bank was setup with too many filters, or if the minimum frequency 
is set too low.
- fixed an error in Filters recursive filter calculation

## 1.0.0 - 2018-03-09
### Added
- TactileWaves Package - the core toolbox dsp library and audio engine