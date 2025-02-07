Description
=====
DIPFIT is an [EEGLAB](http://eeglab.org) plugin to perform inverse source localization.

A major obstacle to using EEG data to visualize macroscopic brain dynamics is the underdetermined nature of the inverse problem: Given an EEG scalp distribution of activity observed at given scalp electrodes, any number of brain source activity distributions can be found that would produce it. This is because there is any number of possible brain source area pairs or etc. that, jointly, add to the scalp data. Therefore, solving this EEG inverse problem requires making additional assumptions about the nature of the source distributions. A computationally tractable approach is to find some number of brain current dipoles (like vanishingly small batteries) whose summed projections to the scalp most nearly resemble the observed scalp distribution.

Documentation
====
For documentation see [https://eeglab.org/tutorials/09_source/DIPFIT.html](https://eeglab.org/tutorials/09_source/DIPFIT.html)

Version history
=====
v5.5
- Force unit to mm for Leadfield Matrix computation

v5.4
- Adding Neuromag coordinate format
- Better support to exclude channels
- Better support for STUDY and group-level Leadfield computation

v5.3
- Fixing ignoring channel when computing LeadField matrix
- Fixing setting up fiducial manually

v5.2
- Fixing issue with normalizing plot in dipplot

v5.1
- Adding 3-D plot for eLoreta

v5.0
- Adding MEG support
- Adding custom Fieldtrip model support
- Various optimizations

v4.4
- Conversion for MFF file electrode coordinates
- New Atlas for Loreta
- Compute leadfield for multiple datasets
- Add electrode O9 and O10 to template
- Allow processing STUDY with no event

v4.3
- Fix Leadfield matrix error for file not found and issue with Fieldtrip conversion

v4.2
- Reverting to using old file to avoid MNI coordinate conversion conflict

v4.1
- Updating corrupted CED file to prevent crash

v4.0
- Allowing computing and storing leadfield matrix

v3.7
- Better compatibility with compiled version of EEGLAB

v3.6
- Fix menu tag so the ERPSOURCE plugin can find the DIPFIT menu item

v3.5
- Use correct MRI for multifit, reorder menus

v3.4
- Adding possiblity to run DIPFIT on EEGLAB STUDY

v3.3
- Ensures backward compatiblity with old versions of EEGLAB

v3.2
- Fix bug with missing folders

v3.1
- Improved support for eeglab2fieldtrip.m

v3.0
- Adding support for eLoreta

v2.0
- Now uses Fieldtrip instead of its own functions

