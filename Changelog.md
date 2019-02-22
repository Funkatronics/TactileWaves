# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).

## 1.1.0 - 2019-02-22 
### Added
- Com package for easy Bluetooth communication on Android and Java
- Windowing and pre-emphasis filtering to FormantExtractor (increases formant accuracy)
- a2db methods to WaveFrame to convert linear values to different decibel scales
- new Biquad filters to the Filter class

### Changed
- Updated and corrected code comments and documentation
- Filter class was modified to reflect new filter types 
- Following the extensive performace testing that was performed in 2018, the FFT class has been 
further optimized resulting in a modest performace boost
- FFT class has also been updated for readability and ease of use. It's documenation and usage 
should be more clear now
- LPC class was rearranged a bit for clarity and consistency with other updates
- MFCC algorithm was modifed to match other common MFCC implementations
- YIN methods were converted to static becasue there is no need for instantiation
- Step 6 of the YIN PDA was added to the YIN class, which further refines the estimated pitch 
- Sort class was cleaned up and further optimized
- AndroidFormat method in WaveFormat was changed to properly reflect the default sampling rate on 
Android (44.1 kHz)
- Window classs was restructured to include an extensive collection of Windowing Fucntions
- Minor changes throughout the library to ensure compatibility and consistency with above changes

### Removed
- lin2db method from WaveFrame because it has been replaced by a2db methods

### Fixed 
- Cleaned up various sections of the code
- Corrected several typos in comments/docs
- Various bug fixes and performace tweaks across the library

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
- FormantExtractor now ony returns valid formants (frequency > 90 Hz, bandwidth < 400 Hz)
- index.html was changed to add a few small tweaks to the web page

### Fixed 
- MFCC returning NaN's if filter bank was setup with too many filters, or if the minimum frequency 
is set too low.
- fixed an error in Filters recursive filter calculation
- fixed an error in FFT's Power spectrum calculation - previously was dividing by N instead of N^2

## 1.0.0 - 2018-03-09
### Added
- TactileWaves Package - the core toolbox dsp library and audio engine