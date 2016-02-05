SMA - Single Molecule Analysis
BDA 10/22/15 UNDER DEVELOPMENT
Updated 2/4/16 - developing for 2 channel data. IN PROGRESS. currently, can generate 2 channel .pks3d file. But, no mapping is applied - assumes half frame offset. Keep emchs =1 and it works fine for real data analysis.

Dependencies - See Storm analysis and DAOSTORM analysis requirements. You should be able to run STORM analysis using the sample data in the 3d_daostorm directory before trying to use this code.
https://github.com/ZhuangLab/storm-analysis
Also requires yattag - www.yattag.org
And Scikit-image http://scikit-image.org/, which in turn requires cython.

Maybe later I'll try to pare down to what is really needed from STORM analysis

This code is run in several steps:

-[TO DO: For multi channel data: generate bead mapping]

-Find peaks, either in the beginning of the movie (typically for smFRET data) or throughout the entire movie (typically for ORBIT data). Output .pks or .pks3d file. Currently only .pks3d.

-Analyze those peaks: intensities, and optionally fits. Output either as .traces file (for smFRET data) or as .trdir directory of single trace files (for ORBIT data). Currently only .trdir.


Major pieces:
-ffpdax: find peaks in .dax movie. Can be set to automatically run apdax when done.
-apdax: analyze those peaks.
-single xml file specifies settings for both codes.
-multibatch3 - run multiple analyses in parallel

Functional and tested for a very limited set of settings (for now). Works with: .pks3d, emchs=1, .trdir, ALEX4 = 0 (and probably ALEX4 = 1, but untested)