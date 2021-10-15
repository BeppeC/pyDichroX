v 0.4.0
-------
- Added hystereis on the fly for Soleil.

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
- Fixed a bug in Boreas configuration files related to importing log
  information.

v 0.2.0
-------
- Added configuration for Boreas beamline at Alba Synchrotron 
  (Barcellona, Spain) to compute XMCD, XNCD, XNXD and XNLD.
- In order to correctly manage data coming from different facilities
  data log files are now mandatory for the functioning of the program.
- Data logs are now entirely managed by configuration files.
- Bug concerning overwriting of different datafile with the same name
  fixed.