SMA - Single Molecule Analysis
BDA 10/22/15 UNDER DEVELOPMENT
Updated 2/22/16 - developing for 2 channel data. Done, lightly tested on bead data. map_coords works; warping the image is either right or very close.


Dependencies - 
numpy, scipy, Matplotlib, PIL
https://github.com/ZhuangLab/storm-analysis
Also requires yattag - www.yattag.org
OpenCV -- http://opencv.org/


This code is run in several steps:

-[TO DO: For multi channel data: generate bead mapping]

-Find peaks, either in the beginning of the movie (typically for smFRET data) or throughout the entire movie (typically for ORBIT data). Output .pks or .pks3d file. Currently only .pks3d.

-Analyze those peaks: intensities, and optionally fits. Output either as .traces file (for smFRET data) or as .trdir directory of single trace files (for ORBIT data). Currently only .trdir.


Major pieces:
-ffpdax: find peaks in .dax movie. Can be set to automatically run apdax when done.
-apdax: analyze those peaks.
-single xml file specifies settings for both codes.
-multibatch - run multiple analyses in parallel. (Or, for better output, run each in a separate command line window)

Functional and tested for a very limited set of settings (for now). Works with: .pks3d, emchs=1 or 2, .trdir, ALEX4 = 0 or 1