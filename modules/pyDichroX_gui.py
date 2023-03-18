"""
pyDichroX_GUI.py

Class and methods to provide user interfaces for XMCD, XNCD, XNLD and
hysteresys data analysis.

Classes
--------
GUI : Provides methods to manage GUI.

Methods
-------
ask_continue()
    GUI dialogue to ask for continuing with other analysis.

ask_quit(title, mes)
    GUI dialogue to ask if quit or not the program.

no_config()
    GUI dialogue for no presence of configuration file.

set_config(cfg_list)
    GUI dialogue to select configuration file.

sel_edg_fls()
    GUI dialogue to select or create edge-list file.
"""

# Copyright (C) Giuseppe Cucinotta.
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

import sys
import easygui as eg
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.interpolate as itp
from matplotlib.widgets import TextBox, Button
from scipy.special import expit
from scipy import sparse
from scipy.sparse.linalg import spsolve

import modules.pyDichroX_escan_datatreat as esdt


class GUI:
    '''
    Provides methods to manage GUI.

    Attributes
    ----------
    type : dict
        Links each type of measurements (keys) to corresponding analysis
        type (values).

    sense : dict
        Links the type of detection (TEY or Fluorescence) to
        corresponding analysis.

    title : str
        Title of windows GUI

    analysis : str
        Name of analysis type.

    interactive : bool
        True for interactive mode false otherwise

    infile_ref : bool
        To select GUI message. True if input files are related to 
        reference sample, False otherwise.

    ------- Attributes for graphs in edges selection ------
    escale : array
        Energy scale, for graphs.

    y1int : UnivariateSpline.
        Interpolation of spectra to be plotted in edge choice.

    y2int : UnivariateSpline.
        Interpolation of spectra to be plotted in edge choice.

    st_e_vals : list
        Collect edge, pre-edge and post-edge energies dafult starting
        values.
        [0] : float, expected edge energy
        [1] : float, experimental edge energy
        [2] : float, pre-edge energy
        [3] : float, post-edge energy / float, smooth param for arpls
        [4] : int, pre-edge width
        [5] : bool, tag for energy scale recalibration
        [6] : array, positive spectra for arpls
        [7] : array, negative spectra for arpls
        [8] : pandas DataFrame, logdatatable for arpls

    e_vals : list
        Same as st_e_vals but collecting updated values

    bsl : 2d line
        Linear interpolation of baseline.

    intedg : 2d plot
        Interpolated value of baseline at edge energy.

    exel : line obj
        vline for edge energy.

    pewl : vspan obj
        highligt pre-edge range in plot.

    pel : line obj
        vline for pre-edge energy.

    psel : line obj
        vline for post-edge energy.

    bsliso : line obj
        baseline of average XAS data.

    isonbsl : line obj
        acvraged XAS baseline subtracted

    text_box_ed : TextBox obj
        TextBox widget used to set edge energy.

    text_box_pe : TextBox obj
        TextBox widget used to set pre-edge energy.

    text_box_pewd : TextBox obj
        TextBox widget used to set pre-edge width energy.

    text_box_pste : TextBox obj
        TextBox widget used to set post-edge energy.

    text_box_lam : TextBox obj
        TextBox widget used to set lambda smoothing parameter.

    bnrec : Button obj
        Button widget used to set recalibration energy scale choice.

    bsl_int : bool
        Select method to process baseline. If True ALS method is used,
        if False a linear approximation is adopted.


    Methods
    -------    
    chs_analysis()
        GUI dialogue to choose the data analysis to perform.

    chs_edge(edg_filename)
        Provides a GUI to set the edge and pre-edge energies for energy
        scan analysis.

    set_edges_lin(sel_edg, exper_edge, x, y1, y2, y1int, y2int)
        Allows user set the edge and pre-edge energies from experimental
        data in order to perform energy calibration and normalizations.

    arpls_bsl(y, lam)
        Create baseline for spectra using the asymmetrically reweighted
        penalized least squares regression.

    comp_bsl(pos, neg, log_dt)
        Create arpls baseline for positive and negative spectra and then
        compute baseline for xd and averaged XAS spectra.

    set_edge_energy(edg_e)
        Update edge energy line on graphs based on input in TextBox.

    set_pedge_energy(pedg_e)
        Update pre-edge energy line on graphs based on input in TextBox.

    set_psedge_energy(psedg_e):
        Update post-edge energy line on graphs based on input in
        TextBox.

    set_pew_energy(pew_e)
        Update pre-edge average range highlighted region based on energy
        width input in TextBox.

    set_lam_par(lam)
        Update baseline on graphs based on smoothing parametr input in
        TextBox.

    recal(event)
        Action associated to click on recalibrate button.

    reset_en(event)
        Reset TextBox and graphs to energy starting values.

    finish_but(event)
        Close the figure on pressing the button Finish.

    ask_angle()
        Provides a GUI to ask for experimental angle of the sample
        respect the beam.

    ask_T()
        GUI dialogue to ask for sample temperature.

    ask_H()
        GUI dialogue to ask for magnetic field.

    sel_ifls_fls(confobj)
        GUI dialogue to select data input files.

    add_set()
        GUI dialogue to add further sets of XMCD/XNCD/XNLD scans.

    in_dtst(confobj)
        GUI dialogue to select input files data sets.

    not_enough_fls(pol, preedge)
        GUI warning that not enough files for a given polarization have
        been supplied.

    not_enough_sp_br(self, scanobj, pol, preedge=False)
        GUI warning that not enough files for a given polarization have
        been supplied. For splitted branch hysteresis scans.

    ask_logfn(self)
        Gui dialogue to ask for datalog filename.

    no_log(fname, action='')
        Message box indicating a missing logfile.

    num_pnts(num_pts)
        GUI dialogue to set the number of points used for energy scan
        interpolations.

    num_times(log_dt)
        GUI dialogue to set the number of time steps to be plotted in
        hysteresis point by point splitted time analysis.

    chs_scns(choices)
        GUI dialogue to select the scans to be averaged from a list of
        labels.

    ask_bsl_interp()
        GUI dialogue to choose which method use to process baseline
        between linear approximation and peak screened ALS method.

    confirm_choice()
        Provides a GUI to let user confirms his choice and continue or
        make a different choice.

    outfile_name(default_nm)
        GUI dialogue to choose output file name.

    wrongpol(scn_num, polarisation)
        GUI dialogue to warn that a file with wrong polarisation has
        been passed.

    acq_times(pos, neg, fieldscn)
        Set start and end values to define the time window to select
        data to for hysteresis point by point analysis.

    ask_acq_times()
        Prompt for starting and ending acquisition times to be 
        considered in hysteresis point by point analysis.
    '''

    def __init__(self, confobj):
        '''
        At instantiation the title of windows, and the attributes type,
        interactive, log and sense are created.

        Parameters
        ----------
        confobj : configuration object.
        '''
        self.type = {'hyst': ['hyst_fly', 'hyst_t_aver', 'hyst_t_split'],
                     'xmcd': ['XMCD'],  # analysis values for XMCD
                     'xncd': ['XNCD'],  # analysis values for XNCD
                     'xnld': ['XNLD'],  # analysis values for XNLD
                     'xnxd': ['XNXD']}  # analysis values for XNXD

        self.title = 'pyDichroX'

        self.interactive = confobj.interactive
        self.sense = confobj.sense
        self.list_analysis = confobj.list_analysis

    def chs_analysis(self):
        '''
        GUI dialogue to choose the data analysis to perform. It also
        sets window's title and analysis attributes.

        Returns
        -------
        set title attribure (str) and analysis attribute (str).
        '''
        self.title = 'pyDichroX'

        msg = 'Select the analysis you want to perform.'
        choices = self.list_analysis

        self.title += ' {} analysis'.format(self.sense)

        while True:
            a_choice = eg.choicebox(msg=msg, title=self.title, choices=choices)

            if a_choice is None:
                ask_quit(self.title)
            else:
                break

        # Make easier id for hysteresis analysis
        if a_choice == 'XMCD Hysteresis on the fly':
            self.analysis = 'hyst_fly'
        elif a_choice == 'XMCD Hysteresis point by point - time average':
            self.analysis = 'hyst_t_aver'
        elif a_choice == 'XMCD Hysteresis point by point - time split':
            self.analysis = 'hyst_t_split'
        else:
            self.analysis = a_choice

        self.title = '{} {} analysis'.format(a_choice, self.sense)

    def chs_edge(self, edg_filename):
        '''
        GUI dialogue to set the edge and pre-edge energies for energy
        scan analysis.
        The edge and pre-edge energies can be choosen from a file, if
        present. Otherwise the file is created.
        If the needed edge is missing it can be added and saved to the
        file.

        Parameters
        ----------
        edg_filename : str
            File lisitng edge and pre-edge energies.

        Returns
        -------
        list (str) containing for the selcted edge:
        - [0] : name of the edge
        - [1] : edge energy
        - [2] : pre-edge energy
        - [3] : post-edge energy

        Notes
        -----
        The file edge_list.txt containing the edges' data must be a
        four-columns csv, with comma ',' separated values. First line
        contains the headers of the columns, namely 'Name',
        'Edge Energy', 'Pre-edge Energy' and 'Post-edge Energy'.
        '''
        try:
            edg_lst = pd.read_csv(edg_filename, sep=',')
        except:  # if the file is empty
            with open(edg_filename, 'w') as f:
                # Write file headers with column names
                f.write('Name,Edge Energy,Pre-edge Energy,Post-edge Energy')
            edg_lst = pd.read_csv(edg_filename, sep=',')

        msg = 'Choose the edge you want to use or insert a new one.'
        # Create list of edges shown to user
        edges = []
        for i in edg_lst['Name']:
            edges.append(i)
        edges.append('Add a new edge')

        chsn_edge = eg.choicebox(msg, self.title, edges)

        # Check that a choice has been made
        while True:
            if chsn_edge is None:
                ask_quit(self.title)
                chsn_edge = eg.choicebox(msg, self.title, edges)
            else:
                break

        # If user select Add, a new edge is added to the edge list
        if chsn_edge == 'Add a new edge':
            msg = 'Add a new edge.'
            field_nms = ['Name', 'Edge Energy', 'Pre-edge Energy',
                         'Post-edge Energy']
            field_vals = eg.multenterbox(msg, self.title, field_nms)

            # Check that enetered values are valid
            while True:
                errmsg = ''
                for val in field_vals:
                    if not val.strip():
                        errmsg += '\n"{}" is a required field.\n\n'.format(val)
                try:
                    float(field_vals[1])
                except:
                    errmsg += ('\nOnly numerical values are accepted for'
                               + ' Edge energy.\n\n')
                try:
                    float(field_vals[2])
                except:
                    errmsg += ('\nOnly numerical values are accepted for'
                               + ' Pre-Edge energy.\n\n')
                try:
                    float(field_vals[3])
                except:
                    errmsg += ('\nOnly numerical values are accepted for'
                               + ' Post-Edge energy.\n\n')

                if field_vals[0] in edg_lst['Name']:
                    errmsg += ('\nThere\'s already an Edge named {}.\n'.format(
                        field_vals[0]) + 'Please choose another name.')
                if not errmsg:
                    break

                field_vals = eg.multenterbox(msg + errmsg, self.title,
                                             field_nms)

            # Add the new edge to the list
            edg_lst.loc[len(edg_lst)] = field_vals
            edg_lst.to_csv(edg_filename, index=False, mode='w+')

            return field_vals
        else:
            # Select row in DataFrame
            # values.tolist() returns a list of selected rows, each of
            # them is on turn a list of the values in the rows.
            sel_edg = edg_lst[edg_lst['Name'] == chsn_edge]
            sel_edg = sel_edg.values.tolist()[0]
            return sel_edg

    def set_edges_lin(self, sel_edg, exper_edge, x, y1, y2, y1int, y2int):
        '''
        Allow setting the edge energy from experimental data, pre-edge
        energy, post-edge energy and energy range for pre-edge average.
        It also computes a linear interpolation considering pre-edge and
        post-edge values to take into account of baseline effects.

        Parameters
        ----------
        sel_edg : list
            Contains data of the edges considered
            sel_edg[0] expected edge energy
            sel_edg[1] pre-edge energy
            sel_edg[2] post-edge energy.

        exper_edge : float
            experimental edge energy.

        x : array
            energy scale of spectrum.

        y1 : array
            xd spectrum data.

        y2 : array
            xd_average spectrum data.

        y1int : UnivariateSpline
            interpolation of spectrum y1.

        y2int : UnivariateSpline
            interpolation of spectrum y2.

        Returns
        -------
        list :
        [new_vals[1], new_vals[2], new_vals[4], pe_int, new_vals[5]]
        - [0] exper_edge : experimental edge energy
        - [1] e_pe : pre-edge energy
        - [2] e_poste : post-edge energy
        - [3] pe_wdt : number of points representing the half-width
                energy range used for pre-edge average
        - [4] cal : bool.

        '''
        # Initializes number of points for half-width averaging interval
        # of pre-edge energy to 4.
        self.escale = x
        self.y2int = y2int
        self.y1int = y1int
        pe_wdt = 4
        # List with energy starting values
        self.st_e_vals = [sel_edg[0], exper_edge, sel_edg[1], sel_edg[2],
                          pe_wdt, False]
        # List with updated energy values
        self.e_vals = self.st_e_vals.copy()

        fig, ax1 = plt.subplots()
        fig.subplots_adjust(bottom=0.25)

        if self.infile_ref:
            add_title = "Normalized by reference data."
        else:
            add_title = ""
        fig.suptitle("Choose energies\n\n" + self.title + add_title)

        ax1.set_xlabel('E (eV)')
        ax1.set_ylabel(self.analysis + ' (a.u.)', color='black')
        ax1.tick_params(axis='y', labelcolor='black')
        ax1.plot(self.escale, y1, color='black')

        ax2 = ax1.twinx()  # second axes that shares the same x-axis
        ax2.spines['right'].set_color('pink')
        ax2.tick_params(axis='y', colors='pink')
        ax2.yaxis.label.set_color('pink')
        ax2.set_ylabel('Averaged XAS (a.u.)', color='pink')
        ax2.plot(self.escale, y2, color='pink')

        ax3 = ax1.twinx()  # 3rd axes that shares the same x-axis
        ax3.axis('off')

        # Refernce lines for energies
        # Made them as attributes to pass to textbox update function
        ax1.axvline(x=self.e_vals[0], color='blue', linestyle='dashed',
                    label='Expected edge')

        self.exel = ax1.axvline(x=self.e_vals[1], color='seagreen',
                                linestyle='dashed', label='Experimental edge')

        # Range for pre-edge averaging
        pe_idx = np.argmin(np.abs(self.escale - self.e_vals[2]))

        # Left endpoint range index
        lpe_idx = pe_idx - int(self.e_vals[4])
        # If average interval extends over energy range shrink it
        if lpe_idx < 0:
            lpe_e = self.escale[0]
            self.e_vals[4] = pe_idx
        else:
            lpe_e = self.escale[lpe_idx]

        # Right endpoint range index
        rpe_idx = pe_idx + int(self.e_vals[4])
        if rpe_idx >= len(self.escale):
            rpe_e = self.escale[-1]
            self.e_vals[4] = len(self.escale) - 1 - pe_idx
        else:
            rpe_e = self.escale[rpe_idx]

        self.pewl = ax1.axvspan(lpe_e, rpe_e, color='mistyrose')

        self.pel = ax1.axvline(x=self.e_vals[2], color='coral',
                               linestyle='dashed', label='Pre-edge energy')
        self.psel = ax1.axvline(x=self.e_vals[3], color='plum',
                                linestyle='dashed', label='Post-edge energy')
        # Compute linear interpolation of baseline considering pre-edge
        # and post-edge energies

        # Pre-edge and post-edge points on interpolated curve.
        x_int = [self.e_vals[2], self.e_vals[3]]
        y_int_iso = [self.y2int(self.e_vals[2]), self.y2int(self.e_vals[3])]
        y_int = [self.y1int(self.e_vals[2]), self.y1int(self.e_vals[3])]

        edg_int = esdt.lin_interpolate(x_int, y_int_iso, self.e_vals[1])

        self.bsl, = ax1.plot(x_int, y_int, color='indianred',
                             linestyle='dashdot', label='linear baseline')
        self.bsliso, = ax2.plot(x_int, y_int_iso, color='indianred',
                                linestyle='dashdot')
        self.intedg, = ax2.plot(self.e_vals[1], edg_int, marker='x',
                                color='indianred')
        # Averaged xd with lin baseline subtraction
        pst_idx = np.argmin(np.abs(self.escale - self.e_vals[3]))
        sub_bsl = self.escale[pe_idx:pst_idx+1]
        self.isonbsl, = ax3.plot(sub_bsl, self.y2int(sub_bsl) 
                                 - esdt.lin_interpolate(
                                    x_int, y_int_iso, sub_bsl),
                                 color='grey', linestyle='dashdot',
                                 label='Averagd XAS no baseline')

        # Position of text boxes
        boxedg = fig.add_axes([0.1, 0.09, 0.1, 0.05])
        boxpe = fig.add_axes([0.32, 0.09, 0.1, 0.05])
        boxpewdt = fig.add_axes([0.62, 0.09, 0.05, 0.05])
        boxpste = fig.add_axes([0.8, 0.09, 0.1, 0.05])
        # Update graphs based on input
        self.text_box_ed = TextBox(boxedg, 'Edge')
        self.text_box_ed.set_val(np.around(self.st_e_vals[1], decimals=2))
        self.text_box_ed.on_submit(self.set_edge_energy)
        self.text_box_pe = TextBox(boxpe, 'Pre-edge')
        self.text_box_pe.set_val(self.st_e_vals[2])
        self.text_box_pe.on_submit(self.set_pedge_energy)
        self.text_box_pewd = TextBox(boxpewdt, 'Pre-edge width\n# of points')
        self.text_box_pewd.set_val(self.st_e_vals[4])
        self.text_box_pewd.on_submit(self.set_pew_energy)
        self.text_box_pste = TextBox(boxpste, 'Post-edge')
        self.text_box_pste.set_val(self.st_e_vals[3])
        self.text_box_pste.on_submit(self.set_psedge_energy)

        # Recalibrate button
        axrec = fig.add_axes([0.1, 0.02, 0.12, 0.05])
        self.bnrec = Button(axrec, 'Recal OFF')
        self.bnrec.on_clicked(self.recal)
        # Reset button
        axreset = fig.add_axes([0.4, 0.02, 0.1, 0.05])
        bnreset = Button(axreset, 'Reset')
        bnreset.on_clicked(self.reset_en)
        # Finish button
        axfinish = fig.add_axes([0.7, 0.02, 0.1, 0.05])
        bnfinish = Button(axfinish, 'Finish')
        bnfinish.on_clicked(self.finish_but)

        ax1.legend()
        ax3.legend()
        plt.show()

        return [self.e_vals[1], self.e_vals[2], self.e_vals[3], self.e_vals[4],
                self.e_vals[5]]

    def set_edges_arpls(self, sel_edg, exper_edge, x, y1, y2, y1int, y2int,
                        pos, neg, log_dt):
        '''
        Allow setting the edge energy from experimental data, pre-edge
        energy, post-edge energy and energy range for pre-edge average.
        It also computes a linear interpolation considering pre-edge and
        post-edge values to take into account of baseline effects.

        Parameters
        ----------
        sel_edg : list
            Contains data of the edges considered
            sel_edg[0] expected edge energy
            sel_edg[1] pre-edge energy

        exper_edge : float
            experimental edge energy.

        x : array
            energy scale of spectrum.

        y1 : array
            xd spectrum data.

        y2 : array
            xd_average spectrum data.

        y1int : UnivariateSpline
            interpolation of spectrum y1.

        y2int : UnivariateSpline
            interpolation of spectrum y2.

        pos : 

        neg :

        Returns
        -------
        list :
        [new_vals[1], new_vals[2], new_vals[4], pe_int, new_vals[5]]
        - [0] exper_edge : experimental edge energy
        - [1] e_pe : pre-edge energy
        - [2] smooth parameter
        - [3] pe_wdt : number of points representing the half-width
                energy range used for pre-edge average
        - [4] cal : bool
        - [5] univariate spline object with interpolation of average xd
                baseline.

        '''
        # Initializes number of points for half-width averaging interval
        # of pre-edge energy to 4.
        self.escale = x
        self.y2int = y2int
        self.y1int = y1int
        # Starting values of pre-edge width and smoothing parameters
        pe_wdt = 4
        lam = 1E7
        # List with energy starting values
        self.st_e_vals = [sel_edg[0], exper_edge, sel_edg[1], lam, pe_wdt,
                          False, pos, neg, log_dt]
        # List with updated energy values
        self.e_vals = self.st_e_vals.copy()

        fig, ax1 = plt.subplots()
        fig.subplots_adjust(bottom=0.25)

        if self.infile_ref:
            add_title = "Normalized by reference data."
        else:
            add_title = ""
        fig.suptitle("Choose energies and baseline smooth parameter\n\n"
                     + self.title + add_title)

        ax1.set_xlabel('E (eV)')
        ax1.set_ylabel(self.analysis + ' (a.u.)', color='black')
        ax1.tick_params(axis='y', labelcolor='black')
        ax1.plot(self.escale, y1, color='black')

        ax2 = ax1.twinx()  # second axes that shares the same x-axis
        ax2.spines['right'].set_color('pink')
        ax2.tick_params(axis='y', colors='pink')
        ax2.yaxis.label.set_color('pink')
        ax2.set_ylabel('Averaged XAS (a.u.)', color='pink')
        ax2.plot(self.escale, y2, color='pink')

        # Refernce lines for energies
        # Made them as attributes to pass to textbox update function
        ax1.axvline(x=self.e_vals[0], color='blue', linestyle='dashed',
                    label='Expected edge')

        self.exel = ax1.axvline(x=self.e_vals[1], color='seagreen',
                                linestyle='dashed', label='Experimental edge')

        # Range for pre-edge averaging
        pe_idx = np.argmin(np.abs(self.escale - self.e_vals[2]))

        # Left endpoint range index
        lpe_idx = pe_idx - int(self.e_vals[4])
        # If average interval extends over energy range shrink it
        if lpe_idx < 0:
            lpe_e = self.escale[0]
            self.e_vals[4] = pe_idx
        else:
            lpe_e = self.escale[lpe_idx]

        # Right endpoint range index
        rpe_idx = pe_idx + int(self.e_vals[4])
        if rpe_idx >= len(self.escale):
            rpe_e = self.escale[-1]
            self.e_vals[4] = len(self.escale) - 1 - pe_idx
        else:
            rpe_e = self.escale[rpe_idx]

        self.pewl = ax1.axvspan(lpe_e, rpe_e, color='mistyrose')

        self.pel = ax1.axvline(x=self.e_vals[2], color='coral',
                               linestyle='dashed', label='Pre-edge energy')
        # Compute baselines
        pbsl, nbsl, xdbsl, isobsl = self.comp_bsl(pos, neg, log_dt)

        self.bsl, = ax1.plot(self.escale, xdbsl, color='indianred',
                             linestyle='dashdot', label='baseline')
        self.bsliso, = ax2.plot(self.escale, isobsl, color='indianred',
                                linestyle='dashdot')

        # Position of text boxes
        boxedg = fig.add_axes([0.1, 0.09, 0.1, 0.05])
        boxpe = fig.add_axes([0.32, 0.09, 0.1, 0.05])
        boxpewdt = fig.add_axes([0.62, 0.09, 0.05, 0.05])
        boxlam = fig.add_axes([0.85, 0.09, 0.07, 0.05])
        # Update graphs based on input
        self.text_box_ed = TextBox(boxedg, 'Edge')
        self.text_box_ed.set_val(np.around(self.st_e_vals[1], decimals=2))
        self.text_box_ed.on_submit(self.set_edge_energy)
        self.text_box_pe = TextBox(boxpe, 'Pre-edge')
        self.text_box_pe.set_val(self.st_e_vals[2])
        self.text_box_pe.on_submit(self.set_pedge_energy)
        self.text_box_pewd = TextBox(boxpewdt, 'Pre-edge width\n# of points')
        self.text_box_pewd.set_val(self.st_e_vals[4])
        self.text_box_pewd.on_submit(self.set_pew_energy)
        self.text_box_lam = TextBox(boxlam, 'Baseline\nsmoothing par')
        # For readability smoothing parameter is divided by 1E7
        self.text_box_lam.set_val(self.e_vals[3] / 1e7)
        self.text_box_lam.on_submit(self.set_lam_par)

        # Recalibrate button
        axrec = fig.add_axes([0.1, 0.02, 0.12, 0.05])
        self.bnrec = Button(axrec, 'Recal OFF')
        self.bnrec.on_clicked(self.recal)
        # Reset button
        axreset = fig.add_axes([0.4, 0.02, 0.1, 0.05])
        bnreset = Button(axreset, 'Reset')
        bnreset.on_clicked(self.reset_en)
        # Finish button
        axfinish = fig.add_axes([0.7, 0.02, 0.1, 0.05])
        bnfinish = Button(axfinish, 'Finish')
        bnfinish.on_clicked(self.finish_but)

        ax1.legend()
        plt.show()

        pos.bsl, neg.bsl, xdbsl, isobsl = self.comp_bsl(self.e_vals[6],
                                                        self.e_vals[7],
                                                        self.e_vals[8])

        avbsl = itp.UnivariateSpline(self.escale, isobsl, k=3, s=0)

        return [self.e_vals[1], self.e_vals[2], self.e_vals[3], self.e_vals[4],
                self.e_vals[5], avbsl]

    def arpls_bsl(self, y, lam):
        '''
        Create baseline for spectra using the .

        Parameters
        ----------
        y : array
            Data of spectrum

        lam : float
            Smooth parameter

        Return
        ------
        array, baseline extracted from y data with arpls method.

        Notes
        -----
        A logistic function is used to compute the weights to
        discriminate points to be considered in the baselins from points
        belonging to the peaks.
        "Logistic function gives nearly the same weight to the signal
        below or above a baseline when the difference between the signal
        and the baseline is smaller than the estimated noise mean. It
        gradually reduces the weight as the level of the signal
        increases. If a signal is in the 3sigma from the estimated noise
        mean which covers 99.7% of noise on Gaussian assumption, small
        weight is still given.
        Finally, zero weight is given when a signal is much higher than
        the baseline as it can be regarded as a part of the peak. In the
        extreme case the standard deviation is nearly zero, it becomes a
        shifted and reversed unit step function which smoothes and
        estimates a baseline while leaving the peak larger than noise
        mean untouched.
        '''
        # Converging paramter
        tol = 1E-6

        l = len(y)
        D = sparse.diags([1, -2, 1], [0, -1, -2], shape=(l, l-2))
        lamDD = lam * D.dot(D.transpose())
        # Starting weights all one
        w = np.ones(l)
        W = sparse.spdiags(w, 0, l, l)

        while True:
            W.setdiag(w)
            Z = W + lamDD
            z = spsolve(Z, w * y)
            # Difference between baseline and data - modified from
            # original version using absolute values in order to take
            # account of negative peaks
            d = np.array(np.abs(y) - np.abs(z))
            # Consider only negative differences
            dn = d[d < 0]
            m = np.mean(dn)
            s = np.std(dn)
            x = - 2 * (d - (2 * s - m)) / s
            # Weight logistic function
            wt = expit(x)
            # Check i converging condition is met
            if ((np.sqrt(w.dot(w)) - np.sqrt(wt.dot(wt))) < tol):
                break
            else:
                w = wt
        return z

    def comp_bsl(self, pos, neg, log_dt):
        '''
        Create arpls baseline for positive and negative spectra and then
        compute baseline for xd and averaged XAS spectra.

        Parameters
        ----------
        pos : array
            positive spectrum data.

        neg : array
            negative spectrum data.

        log_dt : pandas DataFrame
            DataFrame with log data.

        Return
        ------
        arrays with positive, negative, xd and averaged XAS baselines
        respectively
        '''
        # Compute interpolation of baseline for positive and negative
        # averaged spectra taken singularly
        pbsl = self.arpls_bsl(pos.aver, self.e_vals[3])
        nbsl = self.arpls_bsl(neg.aver, self.e_vals[3])
        # Baseline for xd
        xdbsl = nbsl - pbsl
        # Baseline for averaged XAS
        if self.analysis in self.type['xnld']:
            # If XNLD the angle must be considered for weighted mean
            # computation
            theta = log_dt['bm_angle']  # retrive angle from log table

            # Numerator and denominator terms of angle weight
            ang_w_n = 2 * (np.cos(theta))**2 - (np.sin(theta))**2
            ang_w_d = 3 * (np.cos(theta))**2

            bsliso = (pbsl + (ang_w_n * nbsl)) / ang_w_d
        else:
            bsliso = (pbsl + nbsl) / 2

        # Return univariate object for positive and negative baselines
        pbsl = itp.UnivariateSpline(self.escale, pbsl, k=3, s=0)
        nbsl = itp.UnivariateSpline(self.escale, nbsl, k=3, s=0)
        return pbsl, nbsl, xdbsl, bsliso

    def set_edge_energy(self, edg_e):
        '''
        Update edge energy line on graphs based on input in TextBox.

        Parameters
        ----------
        edg_e : str
            user input for edge energy. 
        '''
        # Check that a number is passed
        while True:
            try:
                edg_e = float(edg_e)
            except:
                # Leave previous value if the input is not a number
                edg_e = self.e_vals[1]
                msg = 'Enter a float number for edge energy.'
                title = self.title
                eg.msgbox(msg=msg, title=title)
            else:
                if (edg_e <= self.escale[0]) or (edg_e >= self.escale[-1]):
                    msg = 'Enter an edge value whithin the energy scale range.'
                    title = self.title
                    eg.msgbox(msg=msg, title=title)
                    # Leave previous value if the input is outside the
                    # energy range
                    edg_e = self.e_vals[1]
                else:
                    # Update edge value
                    self.e_vals[1] = edg_e
                    break
        # Update graph
        self.exel.set_xdata(self.e_vals[1])

        if not self.bsl_int:
            # Linear baseline
            # Pre-edge and post-edge points on interpolated curve.
            x_int = [self.e_vals[2], self.e_vals[3]]
            y_int_iso = [self.y2int(self.e_vals[2]),
                         self.y2int(self.e_vals[3])]
            y_int = [self.y1int(self.e_vals[2]), self.y1int(self.e_vals[3])]

            edg_int = esdt.lin_interpolate(x_int, y_int_iso, self.e_vals[1])

            self.intedg.set_xdata(self.e_vals[1])
            self.intedg.set_ydata(edg_int)

        plt.draw()

    def set_pedge_energy(self, pedg_e):
        '''
        Update pre-edge energy line on graphs based on input in TextBox.

        Parameters
        ----------
        pedg_e : str
            user input for edge energy. 
        '''
        # Check that a number is passed
        while True:
            try:
                pedg_e = float(pedg_e)
            except:
                # Leave previous value if the input is not a number
                pedg_e = self.e_vals[2]
                msg = 'Enter a float number for pre-edge energy.'
                title = self.title
                eg.msgbox(msg=msg, title=title)
            else:
                if (pedg_e <= self.escale[0]) or (pedg_e >= self.escale[-1]):
                    msg = ('Enter a pre-edge value whithin the energy scale'
                           + ' range.')
                    title = self.title
                    eg.msgbox(msg=msg, title=title)
                    # Leave previous value if the input is outside the
                    # energy range
                    pedg_e = self.e_vals[2]
                else:
                    # Update edge value
                    self.e_vals[2] = pedg_e
                    break
        # Update graph
        self.pel.set_xdata(self.e_vals[2])
        self.set_pew_energy(self.e_vals[4])

        if not self.bsl_int:
            # Linear baseline
            # Pre-edge and post-edge points on interpolated curve.
            x_int = [self.e_vals[2], self.e_vals[3]]
            y_int_iso = [self.y2int(self.e_vals[2]),
                         self.y2int(self.e_vals[3])]
            y_int = [self.y1int(self.e_vals[2]), self.y1int(self.e_vals[3])]

            edg_int = esdt.lin_interpolate(x_int, y_int_iso, self.e_vals[1])

            # Averaged xd with lin baseline subtraction
            pe_idx = np.argmin(np.abs(self.escale - self.e_vals[2]))
            pst_idx = np.argmin(np.abs(self.escale - self.e_vals[3]))
            sub_bsl = self.escale[pe_idx:pst_idx+1]

            self.bsl.set_xdata(x_int)
            self.bsl.set_ydata(y_int)
            self.bsliso.set_xdata(x_int)
            self.bsliso.set_ydata(y_int_iso)
            self.intedg.set_xdata(self.e_vals[1])
            self.intedg.set_ydata(edg_int)
            self.isonbsl.set_xdata(sub_bsl)
            self.isonbsl.set_ydata(
                self.y2int(sub_bsl)
                - esdt.lin_interpolate(x_int, y_int_iso, sub_bsl))

        plt.draw()

    def set_psedge_energy(self, psedg_e):
        '''
        Update post-edge energy line on graphs based on input in
        TextBox.

        Parameters
        ----------
        psedg_e : str
            user input for edge energy. 
        '''
        # Check that a number is passed
        while True:
            try:
                psedg_e = float(psedg_e)
            except:
                # Leave previous value if the input is not a number
                psedg_e = self.e_vals[3]
                msg = 'Enter a float number for post-edge energy.'
                title = self.title
                eg.msgbox(msg=msg, title=title)
            else:
                if (psedg_e <= self.escale[0]) or (psedg_e >= self.escale[-1]):
                    msg = ('Enter a post-edge value whithin the energy scale'
                           + ' range.')
                    title = self.title
                    eg.msgbox(msg=msg, title=title)
                    # Leave previous value if the input is outside the
                    # energy range
                    psedg_e = self.e_vals[3]
                else:
                    # Update edge value
                    self.e_vals[3] = psedg_e
                    break
        # Update graph
        self.psel.set_xdata(self.e_vals[3])

        # Pre-edge and post-edge points on interpolated curve.
        x_int = [self.e_vals[2], self.e_vals[3]]
        y_int_iso = [self.y2int(self.e_vals[2]), self.y2int(self.e_vals[3])]
        y_int = [self.y1int(self.e_vals[2]), self.y1int(self.e_vals[3])]

        edg_int = esdt.lin_interpolate(x_int, y_int_iso, self.e_vals[1])

        # Averaged xd with lin baseline subtraction
        pe_idx = np.argmin(np.abs(self.escale - self.e_vals[2]))
        pst_idx = np.argmin(np.abs(self.escale - self.e_vals[3]))
        sub_bsl = self.escale[pe_idx:pst_idx+1]

        self.bsl.set_xdata(x_int)
        self.bsl.set_ydata(y_int)
        self.bsliso.set_xdata(x_int)
        self.bsliso.set_ydata(y_int_iso)
        self.intedg.set_xdata(self.e_vals[1])
        self.intedg.set_ydata(edg_int)

        self.isonbsl.set_xdata(sub_bsl)
        self.isonbsl.set_ydata(
            self.y2int(sub_bsl)
            - esdt.lin_interpolate(x_int, y_int_iso, sub_bsl))

        plt.draw()

    def set_pew_energy(self, pew_e):
        '''
        Update pre-edge average range highlighted region based on energy
        width input in TextBox.

        Parameters
        ----------
        pew_e : str
            user input for pre-edge width. 
        '''
        # Check that a number is passed
        while True:
            try:
                pew_e = int(pew_e)
            except:
                # Leave previous value if the input is not a number
                pew_e = self.e_vals[4]
                msg = 'Enter an integer number for pre-edge energy width.'
                title = self.title
                eg.msgbox(msg=msg, title=title)
            else:
                # Range for pre-edge averaging
                pe_idx = np.argmin(np.abs(self.escale - self.e_vals[2]))

                # Left endpoint range index
                lpe_idx = pe_idx - pew_e
                if lpe_idx < 0:
                    lpe_e = self.escale[0]
                    self.e_vals[4] = pe_idx
                else:
                    lpe_e = self.escale[lpe_idx]
                    self.e_vals[4] = pew_e
                # Right endpoint range index
                rpe_idx = pe_idx + int(self.e_vals[4])
                if rpe_idx >= len(self.escale):
                    rpe_e = self.escale[-1]
                    self.e_vals[4] = min(self.e_vals[4], (len(self.escale)
                                                          - 1 - pe_idx))
                else:
                    rpe_e = self.escale[rpe_idx]
                    self.e_vals[4] = min(self.e_vals[4], pew_e)
                break
        # Update graph
        self.pewl.set_xy([[lpe_e, 0], [lpe_e, 1], [rpe_e, 1], [rpe_e, 0]])

        plt.draw()

    def set_lam_par(self, lam):
        '''
        Update baseline on graphs based on smoothing parametr input in
        TextBox.

        Parameters
        ----------
        lam : str
            user input for smoothing lambda parameter. 
        '''
        # Check that a number is passed
        while True:
            try:
                lam = float(lam)
            except:
                # Leave previous value if the input is not a number
                lam = self.e_vals[3] / 1e7
                msg = 'Enter a float number for smoothing parameter.'
                title = self.title
                eg.msgbox(msg=msg, title=title)
            else:
                # Update edge value
                self.e_vals[3] = lam * 1e7
                break
        # Update baselines
        pbsl, nbsl, xdbsl, isobsl = self.comp_bsl(self.e_vals[6],
                                                  self.e_vals[7],
                                                  self.e_vals[8])
        self.bsl.set_xdata(self.escale)
        self.bsl.set_ydata(xdbsl)
        self.bsliso.set_xdata(self.escale)
        self.bsliso.set_ydata(isobsl)

        plt.draw()

    def recal(self, event):
        '''
        Action associated to click on recalibrate button.
        When clicked change the label on button and set the value on
        recalibrate field in self.eval
        '''
        # Change the label of button together with recal value in e_vals
        if self.bnrec.label.get_text() == 'Recal OFF':
            self.bnrec.label.set_text('Recal ON')
            self.e_vals[5] = True
        else:
            self.bnrec.label.set_text('Recal OFF')
            self.e_vals[5] = False

        plt.draw()

    def reset_en(self, event):
        '''
        Reset TextBox and graphs to energy starting values.
        '''
        # Reset energy values to starting ones
        self.e_vals = self.st_e_vals.copy()

        self.text_box_ed.set_val(np.around(self.e_vals[1], decimals=2))
        self.text_box_pe.set_val(self.e_vals[2])
        self.text_box_pewd.set_val(self.e_vals[4])
        if self.bsl_int:
            self.text_box_lam.set_val(self.e_vals[3]/1E7)
        else:
            self.text_box_pste.set_val(self.e_vals[3])

        # Range for pre-edge averaging
        pe_idx = np.argmin(np.abs(self.escale - self.e_vals[2]))

        pew_e = self.e_vals[4]

        # Left endpoint range index
        lpe_idx = pe_idx - pew_e
        if lpe_idx < 0:
            lpe_e = self.escale[0]
            self.e_vals[4] = pe_idx
        else:
            lpe_e = self.escale[lpe_idx]
        # Right endpoint range index
        rpe_idx = pe_idx + int(self.e_vals[4])
        if rpe_idx >= len(self.escale):
            rpe_e = self.escale[-1]
            self.e_vals[4] = min(self.e_vals[4], (len(self.escale) - 1
                                                  - pe_idx))
        else:
            rpe_e = self.escale[rpe_idx]
            self.e_vals[4] = min(self.e_vals[4], pew_e)

        # Update graph
        self.exel.set_xdata(self.e_vals[1])
        self.pel.set_xdata(self.e_vals[2])

        if self.bsl_int:
            pbsl, nbsl, xdbsl, isobsl = self.comp_bsl(self.e_vals[6],
                                                      self.e_vals[7],
                                                      self.e_vals[8])
            self.bsl.set_xdata(self.escale)
            self.bsl.set_ydata(xdbsl)
            self.bsliso.set_xdata(escale)
            self.bsliso.set_ydata(isobsl)
        else:
            self.psel.set_xdata(self.e_vals[3])

            # Pre-edge and post-edge points on interpolated curve.
            x_int = [self.e_vals[2], self.e_vals[3]]
            y_int_iso = [self.y2int(self.e_vals[2]),
                         self.y2int(self.e_vals[3])]
            y_int = [self.y1int(self.e_vals[2]), self.y1int(self.e_vals[3])]

            edg_int = esdt.lin_interpolate(x_int, y_int_iso, self.e_vals[1])

            self.bsl.set_xdata(x_int)
            self.bsl.set_ydata(y_int)
            self.bsliso.set_xdata(x_int)
            self.bsliso.set_ydata(y_int_iso)
            self.intedg.set_xdata(self.e_vals[1])
            self.intedg.set_ydata(edg_int)

        self.pewl.set_xy([[lpe_e, 0], [lpe_e, 1], [rpe_e, 1], [rpe_e, 0]])

        plt.draw()

    def finish_but(self, event):
        '''
        Close the figure on pressing the button Finish.
        '''
        plt.close()

    def ask_angle(self):
        '''
        GUI dialogue to ask for experimental angle of the sample respect
        to the X-Ray beam.

        Returns
        -------
        float, the inserted angle value, for log purpose
        float, the complementary angle of the entered one in radiants.
        '''
        msg = 'Enter the rotation angle of the sample.'
        in_ang = 0

        new_angle = eg.enterbox(msg, self.title, in_ang)

        # Error messages for wrong inputs
        err_empty = '\nPlease provide an angle.\n'
        err_num = '\nOnly numbers are accepted.\n'

        while True:
            errmsg = ''
            if new_angle is None:  # If Cancel is pressed exit program
                ask_quit(self.title)
                errmsg = ' '
            if not new_angle:
                errmsg += err_empty
            else:
                try:
                    float(new_angle)
                except:
                    errmsg += err_num
            if not errmsg:
                break
            new_angle = eg.enterbox(msg + errmsg, self.title, in_ang)

        # Return input angle value (for log purpose) and complementary
        # angle (the angle between X-Ray beam and the sample sufrace
        # normal) in radians for XNLD computations
        return new_angle, np.deg2rad(90 - float(new_angle))

    def ask_T(self):
        '''
        GUI dialogue to ask for sample temperature.

        Returns
        -------
        float, the inserted sample temperature, for log purpose.
        '''
        msg = 'Enter the sample temperature.'

        new_T = eg.enterbox(msg, self.title)

        # Error messages for wrong inputs
        err_empty = '\nPlease provide a temperature.\n'
        err_num = '\nOnly numbers are accepted.\n'

        while True:
            errmsg = ''
            if new_T is None:  # If Cancel is pressed exit program
                ask_quit(self.title)
                errmsg = ' '
            if not new_T:
                errmsg += err_empty
            else:
                try:
                    float(new_T)
                except:
                    errmsg += err_num
            if not errmsg:
                break
            new_T = eg.enterbox(msg + errmsg, self.title)

        return float(new_T)

    def ask_H(self):
        '''
        GUI dialogue to ask for magnetic field.

        Returns
        -------
        float, the inserted magnetic, for log purpose.
        '''
        msg = 'Enter the magnetic field value.'

        new_H = eg.enterbox(msg, self.title)

        # Error messages for wrong inputs
        err_empty = '\nPlease provide a field value.\n'
        err_num = '\nOnly numbers are accepted.\n'

        while True:
            errmsg = ''
            if new_H is None:  # If Cancel is pressed exit program
                ask_quit(self.title)
                errmsg = ' '
            if not new_H:
                errmsg += err_empty
            else:
                try:
                    float(new_H)
                except:
                    errmsg += err_num
            if not errmsg:
                break
            new_H = eg.enterbox(msg + errmsg, self.title)

        return float(new_H)

    def sel_ifls(self, confobj):
        '''
        GUI dialogue to select input files.
        If a cumulative file including single data scans is selected
        the method for extract them and save them in a folder is called.

        Parameters
        ----------
        confobj : configuration object.

        Returns
        -------
        List (str) with the filenames to be opened.
        '''
        # Default file extension
        default = confobj.default_ext

        ask_data = True  # while loop controller
        while ask_data:
            # Select message
            if self.infile_ref:
                msg = ("Choose reference data files")
            else:
                msg = ("Choose data files")

            if self.analysis in self.type['hyst']:
                msg += "\nInclude both Edge and Pre-edge scans."

            # If confobj.filetypes there's the possibility of extract
            # data files from a cumulative file
            if confobj.filetypes:
                filetypes = confobj.filetypes

                f_nms = eg.fileopenbox(msg=msg, title=self.title,
                                       multiple=True, default=default,
                                       filetypes=filetypes)

                only_ex_ftyp = []  # only filetypes w/o *

                for ft in filetypes:
                    only_ex_ftyp.append(ft.lstrip('*'))
                # Break while loop if no file is passed or dialogue box
                # is closed
                if f_nms == None:
                    ask_data = False
                # If there are more than 1 file they must be all single
                # data files. If cumulative files containing single scan
                # data must be passed in order to extract them,
                # cumulative files must be passed one at time.
                elif len(f_nms) > 1:
                    # ctl var to check the presence of cumulative file
                    chk_ex = False
                    for fn in f_nms:
                        for ft in only_ex_ftyp:
                            # check cumulative file ext is present
                            if ft in fn[-6:]:
                                # if there's also 1 cumulative file brek
                                chk_ex = True
                                ft_fnd = ft  # ft found
                                break
                        if chk_ex == True:
                            # if also 1 cumulative file is present brek
                            only_ex = default.lstrip('*')
                            msg2 = ('Choose only {} '.format(only_ex)
                                    + 'files or one {}'.format(ft_fnd)
                                    + ' file at time')
                            eg.msgbox(msg2, self.title)
                            break
                    # If only scan data are presente break while loop
                    if chk_ex is False:
                        ask_data = False
                # If only 1 file is seleceted check if it's a cumulative
                # file. In this case extract data.
                else:
                    # ctl var to check the presence of cumulative file
                    chk_ex = False
                    fn = f_nms[0]
                    for ft in only_ex_ftyp:
                        # check cumulative file ext is present
                        if ft in fn[-6:]:
                            chk_ex = True
                            msg2 = ('Choose directory where to save'
                                    + ' extracted data.')
                            ex_dir = eg.diropenbox(msg2, self.title)
                            # If no directory or cancel is selected or
                            # dialogue box is closed brak for loop and
                            # continue with while loop
                            if ex_dir == None:
                                break
                            else:
                                ex_dir += '/'
                                confobj.data_ex(fn, ex_dir)
                            break
                    if chk_ex is False:
                        ask_data = False
            else:
                f_nms = eg.fileopenbox(msg=msg, title=self.title,
                                       multiple=True, default=default)
                # If cumulative file are not considered break while loop
                ask_data = False
        return f_nms

    def add_set(self):
        '''
        GUI dialogue to ask for adding further data set.

        Returns
        -------
        bool, True to add another set of data, False otherwise.
        '''
        msg = ("Dou you want to add another set of data?\n"
               + "(N.B. These data will be averaged with the set of data "
               + "already inserted)")

        add = eg.ynbox(msg, self.title)

        if add is None:
            return False  # If window is closed don't add further scans
        else:
            return add

    def in_dtst(self, confobj):
        '''
        GUI dialogue to select input files data sets.

        Parameters
        ----------
        confobj : configuration object.

        Returns
        -------
        nested list (str), each element of the list being a list with
        the filenames of a data set:
        [[dtset11, dtset12, ...], [dtset21, dtset22, ...], ...].

        In case of hysteresis each added dataset is represented by a
        tuple whose first element is a list including the filenames of
        edge scans and the second being a list including the filenames
        of pre-edge scans:
        [([dtset1edg1, dtset1edg2,...], [dtset1pedg1, dtset1pedg2,...]),
        ([dtset2edg1, dtset2edg2,...], [dtset2pedg1, dtset2pedg2,...]),
        ...].
        '''
        dataset = []
        # File selection
        while True:
            dataset.append(self.sel_ifls(confobj))
            if None in dataset:
                dataset.pop()
                ask_quit(self.title, 1)
            else:
                break
        # Ask for datalogfile
        confobj.scanlog_fname(self)

        # Additional data selection
        while True:
            if self.add_set():
                add = self.sel_ifls(confobj)

                if add is None:
                    break
                else:
                    dataset.append(add)
                # Ask for datalogfile
                confobj.scanlog_fname(self)
            else:
                break

        return dataset

    def not_enough_fls(self, pol, preedge=False):
        '''
        GUI warning that not enough files for a given polarization have
        been supplied.

        Parameters
        ----------
        pol : bool
            True for sigma+, CR and LH
            False for sigma-, CL and LV.

        preedge : bool
            True for pre-edge scan in hysteresis analysis
            False otherwise

        Returns
        -------
        bool, True to choose again input files, False to quit program.
        '''
        if self.analysis in self.type['xmcd']:
            if pol:
                data_type = 'sigma +'
            else:
                data_type = 'sigma -'
        elif self.analysis in self.type['xnld']:
            if pol:
                data_type = 'LH'
            else:
                data_type = 'LV'
        else:
            if pol:
                data_type = 'CR'
            else:
                data_type = 'CL'

        if preedge:
            data_type = 'pre-edge ' + data_type

        msg = ("You did not provide enough {} ".format(data_type)
               + "files to continue data analysis.\n\n"
               + "Do you want to choose files again or do you want to quit?")

        return eg.ccbox(
            msg=msg, title=self.title, choices=('Choose again', 'Quit'),
            cancel_choice='Quit')

    def not_enough_sp_br(self, scanobj, pol, preedge=False):
        '''
        GUI warning that not enough files for a given polarization have
        been supplied. For splitted branch hysteresis scans.

        Based on check variables in scanobj warn user which scans miss.

        Parameters
        ----------
        scanobj : HScan data object.

        pol : bool
            True for sigma+, CR and LH
            False for sigma-, CL and LV.

        preedge : bool
            True for pre-edge scan in hysteresis analysis
            False otherwise

        Returns
        -------
        bool, True to choose again input files, False to quit program.
        '''
        if pol:
            data_type = 'CR'
        else:
            data_type = 'CL'

        if preedge:
            data_type = 'pre-edge ' + data_type
            if scanobj.pe_up_hp is not True:
                data_type += ' scan up positive fields'
            if scanobj.pe_up_hn is not True:
                data_type += ' scan up negative fields'
            if scanobj.pe_down_hp is not True:
                data_type += ' scan down positive fields'
            if scanobj.pe_donw_hn is not True:
                data_type += ' scan down negative fields'
        else:
            if scanobj.up_hp is not True:
                data_type += ' scan up positive fields'
            if scanobj.up_hn is not True:
                data_type += ' scan up negative fields'
            if scanobj.down_hp is not True:
                data_type += ' scan down positive fields'
            if scanobj.donw_hn is not True:
                data_type += ' scan down negative fields'

        msg = ("You did not provide enough {} ".format(data_type)
               + "files to continue data analysis.\n\n"
               + "Do you want to choose files again or do you want to quit?")

        return eg.ccbox(
            msg=msg, title=self.title, choices=('Choose again', 'Quit'),
            cancel_choice='Quit')

    def ask_logfn(self):
        '''
        Gui dialogue to ask for datalog filename.

        Used in case a single datalog file is present for an entire
        dataset.

        Return
        ------
        str, datalog filename.
        '''
        msg = 'Choose datalog file associated with your set of data.'

        logfn = eg.fileopenbox(msg=msg, title=self.title, multiple=False)

        return logfn

    def no_log(self, fname, action=''):
        '''
        Message box indicating a missing logfile.

        Parameters
        ----------
        fname : str
            filename related to the missing logfile

        action : str
            describe what will happen due to logfile missing
        '''
        msg = "Logfile {} not found.\n{}".format(fname, action)

        eg.msgbox(msg, self.title)

    def num_pnts(self, num_pts):
        '''
        GUI dialogue to set the number of points used for scan
        interpolations.

        Parameters
        ----------
        num_pts : int
            Default number.

        Returns
        -------
        int, the number of points for scan interpolation.
        '''
        msg = ("Please enter the number of points you want to consider "
               + "for data interpolations.\n\nDefault value is given by the "
               + "average number of data points.")

        while True:
            num_points = eg.integerbox(msg, self.title, default=num_pts,
                                       upperbound=100000)
            if num_points is None:
                ask_quit(self.title)
            else:
                break

        return num_points

    def num_times(self, log_dt):
        '''
        GUI dialogue to set the number of time steps to be plotted in
        hysteresis point by point splitted time analysis.

        Parameters
        ----------
        log_dt : dict
            Collect data for logfile.

        Returns
        -------
        int, the number of time steps to be plotted.
        '''
        n_pt = log_dt['t_scale'][2]
        msg = ('Scans contains {} time steps.\n\n'.format(n_pt)
               + 'Insert number of times step each scan must be plotted')

        default = str(int(np.rint(n_pt * 10 / 100)))
        lwrbn = 1
        uprbn = n_pt

        while True:
            num_points = eg.integerbox(msg, self.title, default=default,
                                       lowerbound=lwrbn, upperbound=uprbn)
            if num_points is None:
                ask_quit(self.title)
            else:
                break

        return num_points

    def chs_scns(self, choices):
        '''
        GUI dialogue to select the scans to be averaged from a list of
        labels.

        Parameters
        ----------
        choices : list (str)
            Labels of scans.

        Returns
        -------
        list (str), the scan numbers of choiced scans.
        '''
        if self.infile_ref:
            msg = "Select the reference scans you want to average."
        else:
            msg = "Select the cans you want to average."

        poss_chs = choices
        choiced = eg.multchoicebox(msg, self.title, poss_chs)

        # Check that user has made a choice
        while True:
            errmsg = ""
            if choiced is None:  # If Cancel is pressed exit program
                ask_quit(self.title)
            if not choiced:
                errmsg += ("\nPlease, choose at least one scan")
            if not errmsg:
                break
            choiced = eg.multchoicebox(msg + errmsg, self.title, poss_chs)

        # Retrieves scan numbers from labels
        for i in range(len(choiced)):
            choiced[i] = choiced[i].split(maxsplit=1)[0]

        return choiced

    def ask_bsl_interp(self):
        '''
        GUI dialogue to choose which method use to process baseline
        between linear approximation and peak screened ALS method.

        Return
        ------
        Set the boolean attribute bsl_int.
        '''
        msg = ("Select the method used to process baselines in order "
               + "to extrapolate the edge jumps")
        choices = [
            ("Linear approximation using pre-edge and post-edge energies"),
            ("Baseline interpolation peak screened ALS method")]

        choice = eg.choicebox(msg, self.title, choices)

        while True:
            noerr = True
            if choice is None:  # If Cancel is pressed exit program
                ask_quit(self.title)
                noerr = False
            elif not choice:
                ask_quit(self.title)
                noerr = False
            if noerr:
                break
            choice = eg.choicebox(msg, self.title, choices)

        if choice == choices[0]:
            self.bsl_int = False
        else:
            self.bsl_int = True

    def confirm_choice(self):
        '''
        GUI dialogue to confirm a choice and continue or make a
        different choice.

        Returns
        -------
        eg.boolbox obj.
        '''
        msg = "Do you confirm your choice?"
        buttons = ["Yes, continue.", "No, make another choice."]

        # It returns in turn bool True for 'Yes' and False for 'No'
        return eg.boolbox(msg, self.title, buttons)

    def outfile_name(self, default_nm):
        '''
        GUI dialogue to choose output file name

        Parameters
        ----------
        default_nm : str
            default value

        Returns
        -------
        str, file name
        '''
        msg = ''
        fname = eg.filesavebox(msg=None, title=self.title, default=default_nm,
                               filetypes='*.dat')

        return fname

    def wrongpol(self, scn_num, polarisation):
        '''
        GUI dialogue to warn that a file with wrong polarisation has
        been passed.

        Parameters
        ----------
        scn_num : str
            scan number

        polarisation : str
            polarisation type.
        '''
        msg = ('Polarisation of scan number {} is not {}'.format(
            scn_num, polarisation) + '\n\nScan number {}'.format(scn_num)
            + ' will not be considered')

        eg.msgbox(msg=msg, title=self.title)

    def acq_times(self, pos, neg):
        '''
        Set start and end values to define the time window to select
        data to for hysteresis point by point analysis.

        Parameters
        ----------
        pos : ScanData obj
            ScanData obj with CR data.

        neg : ScanData obj
            ScanData obj with CL data.

        Return
        ------
        float, float : start (by default 0) and end (by default the
        smallest of max times of all scans) times for common time scale
        definition.
        '''
        while True:
            st_t, end_t = self.ask_acq_times()
            msg = ''

            t_min = pos.min_t.copy()
            t_min.extend(neg.min_t)

            t_max = pos.max_t.copy()
            t_max.extend(neg.max_t)

            #t_num = pos.num_t.copy()
            # t_num.extend(neg.num_t)

            # Check that start and end times are consistent with data
            if st_t < 0:
                st_t = np.around(np.average(t_min), 1)
            if st_t > np.amin(t_max):
                msg += (
                    'The start time inserted is outside the acquired time'
                    + 'range.\nPlease chose a value smaller than {} s'.format(
                    np.amin(t_max)))
            if end_t < 0:
                end_t = np.around(np.average(t_max), 1)
            if st_t >= end_t:
                msg += (
                    'Start time must be smaller than end time.\n'
                    + 'Please chose a value smaller than {} s'.format(end_t))
            if not msg:
                break
            eg.msgbox(msg, self.title)

        return st_t, end_t

    def ask_acq_times(self):
        '''
        Prompt for starting and ending acquisition times to be 
        considered in hysteresis point by point analysis.

        Return
        ------
        float, float
            starting and ending times inserted by user. If nothing is
            passed 0 and -1 are the default values passed by deafult for
            start and end times respectively.
        '''
        msg = ("Enter the START acquisition time you want to consider (sec)"
               + "\n\nIf nothing is passed 0 is taken as default.")
        st_t = eg.enterbox(msg, self.title)

        while True:
            errmsg = ''
            if st_t is None or not st_t:
                st_t = -1
            else:
                try:
                    st_t = float(st_t)
                except:
                    errmsg += (
                        '\n\nPlease insert only numbers for start time.')
                if st_t < 0:
                    errmsg += (
                        '\n\nInsert only positive numbers for start time.')
            if not errmsg:
                break
            st_t = eg.enterbox(msg + errmsg, self.title)

        msg = ("Enter the STOP acquisition time you want to consider (sec)"
               + "\n\nIf nothing is passed last acquired time is taken as"
               + " default.")
        end_t = eg.enterbox(msg, self.title)

        while True:
            errmsg = ''
            if end_t is None or not end_t:
                end_t = -1
            else:
                try:
                    end_t = float(end_t)
                except:
                    errmsg += ('\n\nPlease insert only numbers for end time.')
                if end_t < 0:
                    errmsg += (
                        '\n\nInsert only positive numbers for end time.')
            if not errmsg:
                break
            end_t = eg.enterbox(msg + errmsg, self.title)

        return st_t, end_t


def ask_continue():
    '''
    GUI dialogue to ask for continuing with other analysis.
    '''
    title = 'pyDichroX'
    msg = 'Do you want continue doing other analysis?'

    return eg.ynbox(msg=msg, title=title)


def ask_quit(title, mes=0):
    '''
    GUI dialogue to ask if quit or not the program.
    Quit program if Yes or windows is closed, otherwhise pass.

    Parameters
    ----------
    title : str
        pop-up window's title
    mes : int
        Identifiers for the message.
    '''
    id_msg = {0: '',
              1: 'You didn\'t choose any file\n\n',
              2: 'Edge scans are mandatory\n\n'}

    msg = id_msg[mes] + 'Do you really want to quit?'
    quit = eg.ynbox(msg=msg, title=title)

    if quit or quit is None:
        sys.exit(0)
    else:
        pass


def no_config():
    '''
    GUI dialogue for no presence of configuration file.
    Quit program.
    '''
    title = 'pyDichroX'
    msg = 'No configuration file found, pyDichroX will quit.'
    eg.msgbox(msg=msg, title=title)

    sys.exit(0)


def set_config(cfg_list):
    '''
    GUI dialogue to select configuration file.

    Parameters
    ----------
    cfg_list : list
        list of configuration files.

    Returns
    -------
    Configuration file chosen.
    '''
    title = 'pyDichroX'
    msg = 'Select the configuration file.'

    while True:
        cfg = eg.choicebox(msg=msg, title=title, choices=cfg_list)

        if cfg is None:
            ask_quit(title)
        else:
            break

    return cfg


def sel_edg_fls():
    '''
    GUI dialogue to select or create edge-list file.

    Returns
    -------
    str, the name of the edge-list file.
    '''
    title = 'pyDichroX'
    # Current directory and txt files are setted as default.
    msg = ("Choose to open an existen edge list file or create a new one.")
    choices = ("Open a file", "Create a new file")

    default = 'edge files/*.txt'

    while True:
        chs = eg.boolbox(msg, title, choices)
        if chs:
            msg2 = "Choose the edge list file"
            f_nm = eg.fileopenbox(msg2, title, default=default)
        elif chs is False:
            msg2 = ("A new edge list file will be created.\n" +
                    "Choose the directory and the filename.")
            f_nm = eg.filesavebox(msg2, title, default=default)
        elif chs is None:
            ask_quit(title, 1)
            continue

        if f_nm:
            break
        else:
            ask_quit(title, 1)
            continue

    return f_nm
