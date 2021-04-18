pyDichroX
---------

pyDichroX is a program written in Pyhton to analyse XAS data.

Starting with XAS scan pyDichroX computes XMCD, XNXD, XNCD, XNLD.

Currently pyDichroX analyse only data from Deimos beamline from Soleil
Synchrotron, Paris (France)

REQUIREMENTS
------------
Python >= 3
sys
os
importlib
numpy
scipy
pandas
matplotlib
easygui

USAGE
-----
Just run pyDichroX.py.

Dialogues and interactions through the data analysis are provided by a sequence
of pup-up windows, here listed: 

- Chose the analysis type:
A dialogue ask to choose which type of analysis perform between XMCD, XNXD,
XNCD, XNLD.

- Edge file selection:
Select the edge-file where edge and preedge energies of the element of interests
are listed. By default edge-files are located in 'edge files' folder.
If there's no edge-file present a wizard allows to create one.

- Edge chose:
Select one of the edge present in the edge-file at which work.

- Insert sample angle
Insert the angle of the X-Ray beam with respect the sample.

- Select XAS data files

- Choose to add another set of data
All the further set of data will be averaged together.

- Set interpolation points numbers
During analysis, XAS spectra will be interpolated using the number of points
inserted. By default the number of points in the scans is setted.

- Chose XAS spectra to average
The graphs of the scans are shown and then a dialogue ask which scans select
to be averaged and used for the analysis.

- Energies
A graph with analyzed data is shown togheter with the position of edge energy
present in the edge-file and the edge energy found from experimental data.
The average of the XAS spectra is also shown with a linear baseline based on pre
and post-edge energies. Final data will be normalized based edge jump computed
using both the pre-edge energy setted and linear baseline.
Another pop-up window allows to change edge, pre-edge and post-edge energies.

- Save data
In the end the final normalized data are shown and output file is asked.
As output, a csv file with analized data, the two final plots and a log file are
saved in the chosed directory.


TO DO
-----
At the moment I'm working on

- Add hysteresis analysis
- Add configurations in order to analyse data coming from other beamlines


COPYRIGHTS
----------
Copyright (C) Giuseppe Cucinotta.

This Source Code Form is subject to the terms of the 
Mozilla Public License, v. 2.0. If a copy of the MPL was not distributed with
this file, You can obtain one at https://mozilla.org/MPL/2.0/.
