lcoqa is a collection of software to test performance of APOGEE South
at LCO.

lcoqa has some executables which can be run for each MJD:

 * gcam_qa -t lco -m [mjd] - gathers guide camera exposures
 * exposures_qa -t lco -m [mjd] - gathers exposures from as1D files
 * signal_qa -t lco -m [mjd] - measures signal, S/N info from exposures
 * exposures_qa_plots -t lco -m [mjd] - make plots of exposure information
 * signal_qa_plots -t lco -m [mjd] - make plots of signal information

Data from lcoqa is stored in a directory set by $LCOQA_DATA. The basic
structure is:

* $LCOQA_DATA/
 - [telescope]/
  * apogee-summary.fits - list of all exposures
  * apogeeqa-summary-good.html - QA plots of all "good conditions" exposures
  * apogeeqa-summary-[mjd].html - 
  * [MJD]/
   - gcam-[MJD].fits - summary file of all processed guider images
   - exposures-[MJD].fits - summary file of all science exposures
   - signal-[MJD]-[expno].fits 
    * HDU1 - summary information on this exposure
    * HDU2 - information for each fiber used for summary information
  * plots/ - QA plots
