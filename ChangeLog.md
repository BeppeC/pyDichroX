v 4.0.1
-------
- Update of configuration file for Deimos beamline in order to work with changes in experimental data logfiles formatting

v 4.0.0
-------
- Added errorbands in energy scan analysis. Now it is possible to visualize error bands of energy scan spectra. Error bands are calcuated starting from standard deviation of averaged scans and propagating them through the analysis.
- Bug correction in hysteresis scans average. A bug was corrected which prevented to average multiple scans in hysteeresis on the fly analisis. 

v 3.1.2
-------
- Update in configuartion file for ESRF ID-12 @ ESRF Synchrotron (Grenoble, France) cube configuartion.

v 3.1.1
-------
- Added configurations for ESRF ID-12 @ ESRF Synchrotron (Grenoble, France) 17 T configuartion.
- Added configurations for ESRF ID-12 @ ESRF Synchrotron (Grenoble, France) cube configuartion.

v 3.1.0
-------
- Update of configuration file for ID-32 beamline @ ESRF Synchrotron (Grenoble, France). Data collecting from input file was implemented in a strongest way in order to work also with oldest file.

v 3.0.1
-------
- Minor bug correction in reporting log.

v 3.0.0
-------
- Added configuration and made adjustments in order to perform XMCD, XNLD, XNXD, XNCD and hystereis on the fly data analysis for ID-32 beamline @ ESRF Synchrotron (Grenoble, France)

v 2.1.3
-------
- For energy scans analysis in plot window for choosing edge and pre-edge energies also the preview of XAS averaged with removed baseline is shown, useful for low-dichroism signal data. Only for linear baseline data treatment.

v 2.1.2
-------
- Minor bug fixes

v 2.1.1
-------
- Minor bug fixes

v 2.1.0
-------
- New interactive window for select energy scan data to be average by means of check buttons.

v 2.0.0
-------
- Added new baseline extrapolation in addtion to linear one for the computation of edge jump.
Now it is possible to choose between linear and asymmetrically reweighted penalized least squares regression (2nd derivative constrained weighted regression) following the algorithm proposed by S.-J. Baek, A Park, Y.-J. Ahna and J. Choo (Analyst, 2015, 140, 250–257 DOI: 10.1039/c4an01061b).
- Graph window for the input of edge, pre-edge, post-edge energies and other parameters for baseline treatment is now interactive.

v 1.0.0
-------
- Added analysis for hysteresis point by point for Soleil.

v 0.4.1
-------
- Fixed some bugs on hysteresis analysis.

v 0.4.0
-------
- Added analysis for hystereis on the fly for Soleil.

v 0.3.1
-------
- Fixed a bug in dialogue for setting the edges energies.

v 0.3.0
-------
- Added configuration file for APE beamline at Elettra Synchrotron
  (Trieste, Italy).
- Now is possible to import scan data from a reference sample in order
  to normalize experimental data by reference. To activate this option
  just set True  the ref_norm attribute in configuration file.
- Fixed a bug in analysis selection.
- Fixed a bug which made crash the program if energy scan were provided
  with energy descending order.
- Fixed a bug which made crash the program if data scans have different
  lenght.

v 0.2.1
-------
- Fixed a bug in Boreas configuration file related to importing log
  information.

v 0.2.0
-------
- Added configurations for Boreas beamline at Alba Synchrotron 
  (Barcellona, Spain) to compute XMCD, XNCD, XNXD and XNLD.
- In order to correctly manage data coming from different facilities
  data log files are now mandatory for the functioning of the program.
- Data logs are now entirely managed by configuration files.
- Bug concerning overwriting of different datafile with the same name
  fixed.
