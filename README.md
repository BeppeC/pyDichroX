pyDichroX
---------

pyDichroX is a program written in Pyhton to analyse XAS data.

Starting from XAS scans and related logfiles, pyDichroX computes:
- from energy scans: XMCD, XNXD, XNCD, XNLD spectra;
- from magnetic field scans taken on the fly: XMCD hysteresis;
- from time scan acquisitions taken at different magnetic fields: XMCD hysteresis time-averaging the data for each field or reporting XMCD values for each time and field.

N.B. analysis routines of hysteresis data acquired with point by point methods does not discern different scan field branches. To avoid data from different branches to be averaged all together, analysis of different branches must be accomplished separately.

Results of the data analysis are saved to a csv file togheter with the final graphs.
Also a log file with information concerning experimental conditions (collected from scans' logfiles) and parameter choices durnig analysis is produced.

Making use of configuration files wich take care of specific file formats and layouts, pyDichroX can analyse data collected in different beamlines.

Currently pyDichroX analyses data from:
- DEIMOS Beamline @ Soleil Synchrotron, Paris (France) - XMCD, XNXD, XNCD, XNLD, Hystersis on the fly, Hysteresis point by point
- BOREAS Beamline @ Alba Synchrotron, Barcelona (Spain) - XMCD, XNXD, XNCD, XNLD
- APE Beamline @ Elettra Synchrotron, Trieste (Italy) - XNCD, XNLD


REQUIREMENTS
------------
- Python >= 3.9
- numpy
- scipy
- pandas
- matplotlib
- easygui


USAGE
-----
Just run pyDichroX.py.

Dialogues and interactions through the data analysis are provided by a sequence of pup-up windows, here listed: 

- Chose the analysis type:
A dialogue asks to choose which type of analysis perform between XMCD, XNXD, XNCD, XNLD, Hysteresis on the fly and the two Hysteresis point by point analysis.

N.B. analysis routines of hysteresis data acquired with point by point methods does not discern different scan field branches. To avoid data from different branches to be averaged all together, analysis of different branches must be accomplished separately.

- Edge file selection:
Select the edge-file where edge, preedge and postedge energies of the element of interests are listed. By default edge-files are located in 'edge files' folder.
If there's no edge-file present a wizard allows to create one.
Preedge and postedge energy are used during energy scans analisi to obtain edge jump. 

- Edge chose:
Select one of the edge present in the edge-file at which work. If the edge is not present in the list, the GUI allows to add new edges.

- Insert sample angle
Insert the angle of the X-Ray beam with respect the sample.

- Select XAS data files

- Choose to add another set of data
All the further set of data will be averaged together.

- Set interpolation points numbers
During analysis, XAS spectra will be interpolated using the number of points inserted. By default the number of points in the scans is setted.

- Chose XAS spectra to average
The graphs of the scans are shown. Scans selected by means of checking boxs are averaged and the average computed will be used for successive analysis. Currently only for energy scan analysis.

- If configuration file is setted in order to accept scans from a reference sample, repeat the same passages before for reference scans. Currently only for energy scan analysis.

- Energies - Only for energy scans (XMCD, XNXD, XNCD, XNLD)
A graph with analyzed data is shown togheter with the position of edge energy present in the edge-file and the edge energy found from experimental data.
The average of the XAS spectra is also shown with a baseline based on selected method: linear baseline based on pre and post-edge energies or ArpLS interpolated baseline (based on asymmetrically reweighted penalized least squares regression following the algorithm proposed by S.-J. Baek, A Park, Y.-J. Ahna and J. Choo Analyst, 2015, 140, 250–257 DOI: 10.1039/c4an01061b) where the smoothing parameter is required. Final data will be normalized based edge jump computed using both the pre-edge energy setted and chosed baseline.
Textboxes allow to change edge, pre-edge, post-edge energies and other parameters for baseline which is automoatically updated as the values are entered.

- Save data
In the end the final normalized data are shown and output file name is asked.
As output, a csv file with analized data, the final plots and a log file are
saved in the chosed directory.


TO DO
-----
- Add configurations in order to analyse data coming from other beamlines.


ACKNOWLEDGEMENTS
----------------
PyDichroX was written and is maintained by Giuseppe Cucinotta (giuseppe.cucinotta@protonmail.com).

This program was initially realized with the aim of fulfill the need of scientists at CETECS Laboratory at University of Florence of treating in an efficient and fast way the XAS data collected during their synchrotron experiments. In particular I acknowledge and I am gratefeul to:


dr. Lorenzo Poggini

dr. Giulia Serrano

dr. Michele Serri

Niccolò Giaconi

dr. Andrea Luigi Sorrentino

prof. Matteo Mannini


I also acknowledge the staff of DEIMOS Beamline (Soleil Synchrotron, Paris, France), BOREAS Beamline (Alba Synchrotron, Barcelona , Spain), APE Beamline (Elettra Synchrotron, Trieste, Italy).

The realization of this program would not have been possible without the collaboration, suggestions, discussion, hints and testing on the field they provided.


COPYRIGHTS
----------
Copyright (C) Giuseppe Cucinotta.

This Source Code Form is subject to the terms of the 
Mozilla Public License, v. 2.0. If a copy of the MPL was not distributed with
this file, You can obtain one at https://mozilla.org/MPL/2.0/.
