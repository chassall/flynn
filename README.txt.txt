

Version History

FLYNN 1.0 average eeg text file input, single .mat output
FLYNN 2.0 trial eeg text file input, multiple .mat output
FLYNN 3.0 .mat input (EEGLAB format

FLYNN 3.0 Revision History

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