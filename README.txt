Version History

FLYNN 3.0 .mat input (EEGLAB format
FLYNN 2.0 trial eeg text file input, multiple .mat output
FLYNN 1.0 average eeg text file input, single .mat output

Major/Minor Revision History

3.3.2
Bug fix in FFT baseline.

3.3.1
Fixed a bug for ALL anlyses of only one marker type.

3.3.0
FLYNN can now take continuous EEG data from Analyzer. Also fixed a bug where single participant summaries didn't display rejected trials properly.

3.2.2
Bug in artifact rejection (gradient criterion was only counting increases in voltage). 

3.2.1
Function now returns DISC.

3.2.0
FLYNN now takes as optional second argument a user-specified locs file. The variable chanlocs is compared to the user-specified locs file - chanlocs and data are reordered as needed.

3.1.5
Compare chanlocs of first participant with subsequent participants. Add info to DISC. Display warnings if chanlocs inconsistent.

3.1.4
Changed isRejected to isArtifact for clarity.

3.1.3
Allow for case where user does not provide a path name (default to FLYNNConfiguration.txt)

3.1.2
Added ALL as an analysis type (all trials)

3.1.1
Removed audio, since it won't work in Westgrid anyway

3.1.0
FLYNN is now a function that takes the config file as argument
WAV data now saved properly
doWav actually returns wavelet now
Baseline subtraction now uses repmat to be compatible with older versions of Matlab

3.0.1
Removed working directory (users will just add it to the config file, line one)
Removed tic/toc (users can add this if they wish)
Removed TODO list (will appear below)
Initialized all WAV variables
Use dsearchn to find indices (baseline, ERP, FFT, WAV)
Allow for no EEG baseline, 

TODO
- deal with missing config file
- restrict all ERPs to have the same timepoints?
- warn if about to overwrite exisiting .mat files?
- output a text summary of artifacts by participant and condition
- make artifact rejection a function
- deal with the case that there is only one condition (extra legend entries at end)