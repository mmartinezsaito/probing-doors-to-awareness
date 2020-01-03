# probing-doors-to-awareness
Utility code for the research article "Probing doors to awareness: choice set, visibility, and confidence."
Stimuli generation and display (Matlab with Psychtoolbox), plotting and data analysis (Matlab, R) code files. 

## Index
- Stimuli generation
  + /stimuli
    + Calibration utilities
      * Calibration\_Bitmap\_Black\_1280x1024\_wnumbers.bmp
      * E3calib\_det.m
      * E3calib\_dis1.m
      * E3calib\_dis2.m
    * E3\_4detec.m
    * E3\_4detec\_rps.mat
    * E3\_4detec\_rps\_pilot4b.mat
    * E3\_4discr1.m
    * E3\_4discr1\_rps.mat
    * E3\_4discr1\_rps\_pilot2b.mat
    * E3\_4discr2.m
    * E3\_4discr2\_rps.mat
    * E3\_4discr2\_rps\_pilot2b.mat
    * E3\_4.m
    * EEndScr.m
    * Gamma correction and contrast adjustment hardware specifications
      - AcerAF715\_gamma-corrected.mat
      - VideoSwitcher
- Main file
  * main.m: outlines file reading, preprocessing, and analysis stages 
- File reading routines
  * readRawData.m
  * readPupilData1.m
  * readPupilData2.m
- Preprocessing routines
  * calcSignalwiseVars.m
  * completeFtTable.m
  * confmvars.m
  * deblink\_smooth.m
  * poolIntoTables.m
  * get\_evwi.m
  * getfitfields.m
  * make\_facvec.m
  * makeTableG.m
- Preprocessed data
  * AUC\_all\_4.mat
  * AUCSDI\_all\_4.mat
  * pmcmc\_fdr\_adj.mat
- Data analysis routines
  * doingAnovas.m
  * evlk\_apd\_excldf2.m
  * evlk\_apd.m
  * fittingGlmes.m
  * fitMetadModel\.m: M: Maniscalco and Lau's (2012) function to estimate meta-d'
  * regana.Rmd
  * testingWaveformDifferences.m
- Plotting routines
  * plt\_evlk\_apd.m
  * plt\_evwist.m
  * plt\_raw.m
  * plt\_spectra.m
- License
  * LICENSE: MIT license
