"""
pyDichroX_escan_datatreat.py

Classes and methods for XMCD, XNCD and XNLD data analysis.

Classes
--------
ScanData : Collect raw data for energy scan experiments and provides
    methods to average them.

EngyScan : Allow computation on energy scan data in order to obtain
    XNLD, XMCD, XNCD and XNXD spectra.

Methods
-------
e_scale(guiobj, pos, neg, log_dt, pos_ref, neg_ref)
    Create the energy scale used for data analysis.

lin_interpolate(x, y, x0):
    Linear interpolation of the value at x0.
"""

# Copyright (C) Giuseppe Cucinotta.
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as itp
import scipy.optimize as opt
from scipy import sparse
from scipy.special import expit
from scipy.sparse.linalg import spsolve
from matplotlib.widgets import CheckButtons, Button, TextBox
from matplotlib.path import Path


class ScanData:
    '''
    Collect raw data for energy scan experiments and provides methods to
    average them.

    Attributes
    ----------
    label : list (str)
        Labels for graphs.

    idx : list (str)
        Scan indexes

    raw_imp : pandas DataFrame
        Collect imported raw data.

    energy : array
        common energy scale for average of scans.

    ax : plt.subplot
        subplot where to show scans and their average during choosing
        scans procedure.

    lines : list (2d line objects)
        Collect lines object of raw scan for plotting.

    blines : list (2d line objects)
        Collect lines object of raw scan for plotting.

    weig_val : float
        energy value for weighting data scans
    
    weig_ln : line obj
        vline for weight energy.

    avg_ln : 2d line objects
        Lines object of average of scans.

    plab : list (str)
        Collect the list of labels plotted in graphs for choosing scans.

    checkbx : CheckButtons obj
        Widget for choosing plots.

    aver : array
        Data average of selected scans from raw_imp.
        If ScanData is a reference one, aver contain normalized data by
        reference.

    aver_std : array
        STD of selected scans from raw_imput

    dtype : str
        Identifies the data collected, used for graph labelling:
        sigma+, sigma- for XMCD
        CR, CL for XNCD
        H-, H+ for XNXD
        LH, LV fot XNLD

    chsn_scns : list (str)
        Labels of chosen scans for the analysis.

    pe_av : float
        Value of spectra at pre-edge energy. It is obtained from
        averaging data in a defined energy range centered at pre-edge
        energy.

    pe_av_er : float
        STD on values range taken for computation of pre-edge.

    pe_av_int : float
        Pre-edge value obtained from linear interpolation considering
        pre-edge and post-edge energies.

    bsl : Univariate spline object
        spline interpolation of ArpLS baseline 

    norm : array
        Averaged data normalized by value at pre-edge energy.

    norm_er : array
        Erron on norm.

    norm_int : array
        Averaged data normalized by interpolated pre-edge value.

    norm_int_er : array
        Error on norm_int.

    ej : float
        Edge-jump value.

    ej_er : float
        Error on edge-jump value.

    ej_norm : float
            edge-jump value normalized by value at pre-edge energy.

    ej_int : float
        edge-jump value computed with interpolated pre-edge value.

    ej_int_er : float
        error edge-jump value computed with interpolated pre-edge value.
        ATTENTION! No propagation is performed form interpolation, so
        currently error is the same coming for std at the edge of
        experimental data.

    ej_norm_int : float
        edge-jump computed and normalized by interpolated pre-edge
        value.

    ej_norm_int_er : float
        error on edge-jump computed and normalized by interpolated 
        pre-edge value.
        ATTENTION! No propagation is performed form interpolation, so
        currently error is the same coming for std at the edge of
        experimental data.

    Methods
    -------
    man_aver_e_scans(guiobj, enrg, weig_dt, log_dt)
        Manage the choice of scans to be averaged and return the average
        of selected scans.

    check_but(label)
        When check buttons are checked switch the visibility of
        corresponding line.

    set_weigen(weig_e)
        Update weight energy and line on graph based on user input.

    averbut(event)
        When average button is pressed calls aver_e_scans to compute
        average on selected scans.

    weigaverbut(event)
        When average button is pressed calls weig_aver_e_scans to weight
        and compute average on selected scans.

    aver_e_scans()
        Performe the average of data scans.

    weig_aver_e_scans()
        Weight data scans and average them.

    reset(event)
        Reset graph to starting conditions.

    finish_but(event)
        Close the figure on pressing the button Finish.

    edge_norm(guiobj, enrg, e_edge, e_pe, pe_rng, pe_int)
        Normalize energy scan data by value at pre-edge energy and
        compute edge-jump.
    '''

    def __init__(self):
        '''
        Initialize attributes label, idx, and raw_imp.
        '''
        self.label = []
        self.idx = []
        self.raw_imp = pd.DataFrame()

    def man_aver_e_scans(self, guiobj, enrg, weig_dt, log_dt):
        '''
        Manage the choice of scans to be averaged and return the average
        and standard deviation of selected scans.

        Parameters
        ----------
        guiobj : GUI object
            Provides GUI dialogs.

        enrg : array
            Energy values at which average is calculated.

        weig_dt : bool
            If True weight data before averaging them.

        log_dt : pandas DataFrame
            DataFrame with log data - needed to initialize norm energy
            with pre-edge energy

        Return
        ------
        Set class attributes:
        aver : array
            Average values of the chosen scans.

        aver_std : array
            STD values of the chosen scans.

        chsn_scns : list
            Labels of chosen scans for the analysis (for log purpose).
        '''
        self.energy = enrg
        self.chsn_scns = []
        self.aver = 0

        if guiobj.interactive:  # Interactive choose of scans
            fig, self.ax = plt.subplots(figsize=(10, 6))
            fig.subplots_adjust(right=0.75)

            self.ax.set_xlabel('E (eV)')
            self.ax.set_ylabel(self.dtype)

            if guiobj.infile_ref:
                fig.suptitle('Choose reference sample scans')
            else:
                fig.suptitle('Choose sample scans')
            # Initialize list which will contains line obj of scans
            # lines contain colored lines for choose
            # blines contain dark lines to be showed with average
            self.lines = []
            self.blines = []

            if weig_dt:
                # Initialize normalization. As starting value pre-edge
                # energy is considered
                self.weig_val = float(log_dt['PreEdge_en'])

                self.weig_ln = self.ax.axvline(x=self.weig_val, color='seagreen',
                                          label='Weight Energy')

                # Create new pd.DataFrame to store weighted scans and
                # initialize to raw_imp
                #self.raw_imp_w = self.raw_imp

                # Populate list with line objs
                for i in self.idx:
                    e_col = 'E' + i
                    # Show lines and not blines
                    self.lines.append(
                                    self.ax.plot(self.raw_imp[e_col],
                                    self.raw_imp[i],
                                    label=self.label[self.idx.index(i)])[0])
                # Initialize chsn_scs and average line with all scans and
                # set it invisible
                for line in self.lines:
                    if line.get_visible():
                        self.chsn_scns.append(line.get_label())
                self.weig_aver_e_scans()
            else:    
                # Populate list with line objs
                for i in self.idx:
                    e_col = 'E' + i
                    # Show lines and not blines
                    self.lines.append(
                                    self.ax.plot(self.raw_imp[e_col],
                                    self.raw_imp[i],
                                    label=self.label[self.idx.index(i)])[0])
                    self.blines.append(
                                        self.ax.plot(self.raw_imp[e_col],
                                            self.raw_imp[i], color='dimgrey',
                                            visible=False)[0])
                # Initialize chsn_scs and average line with all scans and
                # set it invisible
                for line in self.lines:
                    if line.get_visible():
                        self.chsn_scns.append(line.get_label())
                self.aver_e_scans()
                # Set a value just for log dt
                self.weig_val = '--'
         
            self.avg_ln, = self.ax.plot(
                self.energy, self.aver, color='red', lw=2, visible=False)
            self.avg_std_fl = self.ax.fill_between(
                self.energy, self.aver - self.aver_std,
                self.aver + self.aver_std, color='red', alpha=0.2,
                visible=False)
            # Create box for checkbutton
            chax = fig.add_axes([0.755, 0.32, 0.24, 0.55], facecolor='0.95')
            self.plab = [str(line.get_label()) for line in self.lines]
            visibility = [line.get_visible() for line in self.lines]
            self.checkbx = CheckButtons(chax, self.plab, visibility)
            # Customizations of checkbuttons
            rxy = []
            bxh = 0.05
            for r in self.checkbx.rectangles:
                r.set(height=bxh)
                r.set(width=bxh)
                rxy.append(r.get_xy())
            for i in range(len(rxy)):
                self.checkbx.lines[i][0].set_xdata(
                    [rxy[i][0], rxy[i][0] + bxh])
                self.checkbx.lines[i][0].set_ydata(
                    [rxy[i][1], rxy[i][1] + bxh])
                self.checkbx.lines[i][1].set_xdata(
                    [rxy[i][0] + bxh, rxy[i][0]])
                self.checkbx.lines[i][1].set_ydata(
                    [rxy[i][1], rxy[i][1] + bxh])

            for l in self.checkbx.labels:
                l.set(fontsize='medium')
                l.set_verticalalignment('center')
                l.set_horizontalalignment('left')

            self.checkbx.on_clicked(self.check_but)

            # Create box for average reset and finish buttons
            if weig_dt:
                wenbox = fig.add_axes([0.87, 0.24, 0.1, 0.06])
                text_wenbox = TextBox(wenbox, 'Weight Energy')
                text_wenbox.set_val(np.around(float(log_dt['PreEdge_en']),
                                              decimals=2))
                text_wenbox.on_submit(self.set_weigen)
                averbox = fig.add_axes([0.77, 0.14, 0.08, 0.06])
                bnaver = Button(averbox, 'Weight\n& Aver')
                bnaver.on_clicked(self.weigaverbut)
                rstbox = fig.add_axes([0.89, 0.14, 0.08, 0.06])
                finbox = fig.add_axes([0.82, 0.04, 0.12, 0.06])
            else:
                averbox = fig.add_axes([0.77, 0.2, 0.08, 0.08])
                bnaver = Button(averbox, 'Average')
                bnaver.on_clicked(self.averbut)
                rstbox = fig.add_axes([0.89, 0.2, 0.08, 0.08])
                finbox = fig.add_axes([0.82, 0.07, 0.12, 0.08])

            bnrst = Button(rstbox, 'Reset')
            bnrst.on_clicked(self.reset)
            bnfinish = Button(finbox, 'Finish')
            bnfinish.on_clicked(self.finish_but)

            self.ax.legend()
            plt.show()

            # If average is not pressed automatically compute average on
            # selected scans
            if self.chsn_scns == []:
                for line in self.lines:
                    if line.get_visible():
                        self.chsn_scns.append(line.get_label())
                if weig_dt:
                    self.weig_aver_e_scans()
                else:
                    self.aver_e_scans()
        else:
            # Not-interactive mode: all scans except 'Dummy Scans' are
            # evaluated
            for lbl in self.label:
                # Check it is not a 'Dummy Scan' and append
                # corresponding scan number in chosen scan list
                if not ('Dummy' in lbl):
                    self.chsn_scns.append(self.idx[self.label.index(lbl)])
            self.aver_e_scans()
            # Set a value just for log dt
            self.weig_val = '--'
        # collect weight energy for final log
        if self.dtype in ['sigma^+', 'CR', 'H-', 'LH']:
            label_w_dt = 'PWEn'
        else:
            label_w_dt = 'NWEn'
        log_dt[label_w_dt] = self.weig_val

    def check_but(self, label):
        '''
        When check buttons are checked switch the visibility of
        corresponding line.
        Also update self.chsn_scns with labels of visible scans.
        '''
        index = self.plab.index(label)
        self.lines[index].set_visible(not self.lines[index].get_visible())
        # Update chsn_scns
        self.chsn_scns = []
        for line in self.lines:
            if line.get_visible():
                self.chsn_scns.append(line.get_label())
        plt.draw()

    def set_weigen(self, weig_e):
        '''
        Update weight energy value and line on graph based on input in
        TextBox.

        Parameters
        ----------
        weig_e : str
            user input for weight energy. 
        '''
        # Check that a number is passed
        while True:
            try:
                weig_e = float(weig_e)
            except:
                # Leave previous value if the input is not a number
                weig_e = self.weig_val
                msg = 'Enter a float number for weight energy.'
                eg.msgbox(msg=msg, title='Error')
            else:
                if (weig_e < self.energy[0]) or (weig_e > self.energy[-1]):
                    msg = ('Enter a weight energy value whithin the energy' +
                           ' scale range.')
                    eg.msgbox(msg=msg, title='Error')
                    # Leave previous value if the input is outside the
                    # energy range
                    weig_e = self.weig_val
                else:
                    # Update weight energy value
                    self.weig_val = weig_e
                    break
        # Update graph
        self.weig_ln.set_xdata(self.weig_val)

    def averbut(self, event):
        '''
        When average button is pressed calls aver_e_scans to compute
        average on selected scans.
        Update self.chsn_scns and self.aver.
        '''
        # Initialize list of chosed scans
        self.chsn_scns = []
        # Set visible only chosen scans in blines and append to
        # chsn_scns
        for i in range(len(self.lines)):
            if self.lines[i].get_visible():
                self.lines[i].set(visible=False)
                self.chsn_scns.append(self.lines[i].get_label())
                self.blines[i].set(visible=True)

        self.aver_e_scans()

        # Update average and std lines and make them visible
        self.avg_ln.set_ydata(self.aver)
        self.avg_ln.set(visible=True)
        # update errorband
        averp = self.aver + self.aver_std
        avern = self.aver - self.aver_std
        verts = np.block([[self.energy, self.energy[::-1]],
                        [averp, avern[::-1]]]).T
        codes = np.full(len(verts), Path.LINETO)
        codes[0] = codes[len(self.aver)] = Path.MOVETO
        path = Path(verts, codes)
        self.avg_std_fl.set_paths([path.vertices])
        self.avg_std_fl.set(visible=True)

        plt.draw()

    def weigaverbut(self, event):
        '''
        When average button is pressed calls weig_aver_e_scans to weight
        and compute average on selected scans.
        Update self.chsn_scns and self.aver.
        '''
        # Initialize list of chosed scans
        self.chsn_scns = []
        
        # Set visible only chosen scans in blines and append to
        # chsn_scns
        for i in range(len(self.lines)):
            if self.lines[i].get_visible():
                self.lines[i].set(visible=False)
                self.chsn_scns.append(self.lines[i].get_label())
        # Empty self.blines: it will be filled with weighted scans by
        # self.weig_aver_e_scans
        self.blines.clear()
        self.weig_aver_e_scans()
        # Make visible blines
        for i in range(len(self.blines)):        
            self.blines[i].set(visible=True)
        # Update average and std lines and make them visible
        self.avg_ln.set_ydata(self.aver)
        self.avg_ln.set(visible=True)
        self.avg_ln.set_zorder(2.5)
        # update errorband
        averp = self.aver + self.aver_std
        avern = self.aver - self.aver_std
        verts = np.block([[self.energy, self.energy[::-1]],
                        [averp, avern[::-1]]]).T
        codes = np.full(len(verts), Path.LINETO)
        codes[0] = codes[len(self.aver)] = Path.MOVETO
        path = Path(verts, codes)
        self.avg_std_fl.set_paths([path.vertices])
        self.avg_std_fl.set(visible=True)
        self.avg_std_fl.set_zorder(2.45)

        plt.draw()

    def aver_e_scans(self):
        '''
        Perform the average of data scans. 
        
        Returns
        -------
        array, containing the average of data scans.

        Notes
        -----
        To compute the average the common energy scale self.enrgy is
        used.
        All passed scans are interpolated with a linear spline (k=1 and
        s=0 in itp.UnivariateSpline) and evaluated along the common
        energy scale.
        The interpolated data are eventually averaged.
        '''
        intrp = []

        for i in self.idx:
            if self.label[self.idx.index(i)] in self.chsn_scns:
                # chosen data
                x = self.raw_imp['E' + i][1:]
                y = self.raw_imp[i][1:]

                # Compute linear spline interpolation
                y_int = itp.UnivariateSpline(x, y, k=1, s=0)
                # Evaluate interpolation of scan data on common energy
                # scale and append to previous interpolations
                intrp.append(y_int(self.energy))

        # Average all inteprolated scans
        self.aver = np.average(intrp, axis=0)
        # Compute STD of the interpolated scans
        self.aver_std = np.std(intrp, axis=0)

    def weig_aver_e_scans(self):
        '''
        Weight data scans and average them.
        
        Returns
        -------
        array, containing the average of weighted data scans.

        Notes
        -----
        For each data scan the weight is computed as
        y(weig_val) / <y(weig_val)>
        where
        y(weig_val) is the interpolated value at weig_val energy 
        <y(weig_val)> is the average of all values at weig_val energy  of
        data scans
        For itnerpolation linear spline (k=1 and s=0 in 
        itp.UnivariateSpline) is employed.

        To compute the average the common energy scale self.energy is
        used.
        All passed scans are first weighted, then interpolated with a
        linear spline (k=1 and s=0 in itp.UnivariateSpline) and
        evaluated along the common energy scale.
        The interpolated data are eventually averaged.
        '''
        intrp_we = []
        intrp = []

        raw_imp_w = pd.DataFrame()

        for i in self.idx:
            if self.label[self.idx.index(i)] in self.chsn_scns:
                # chosen data
                x = self.raw_imp['E' + i][1:]
                y = self.raw_imp[i][1:]

                # Compute linear spline interpolation
                y_int = itp.UnivariateSpline(x, y, k=1, s=0)
                # Evaluate interpolation of scan data at self.weig_val
                int_weig_e = y_int(self.weig_val)
                # Divide scan data by int_weig_e and store them in
                # raw_imp_w
                raw_imp_w[i] = self.raw_imp[i] / int_weig_e
                # Append int_weig_e to intrp_we in order to average all
                # the values at self.weig_val
                intrp_we.append(int_weig_e)

                # Evaluate interpolation of scan data on common energy
                # scale, divide by int_weig_e it and append to previous
                # interpolations
                intrp.append(y_int(self.energy) / int_weig_e)
        # Average all inteprolated scan values at self.weig_val
        weig_e_av = np.average(intrp_we)
        # Populate self.bline with wirghted selected scans
        for i in self.idx:
            if self.label[self.idx.index(i)] in self.chsn_scns:
                e_col = 'E' + i
                # weighted chosen data
                wy = weig_e_av * raw_imp_w[i]
                self.blines.append(
                                    self.ax.plot(self.raw_imp[e_col],
                                            wy, color='dimgrey',
                                            visible=False)[0])
        # Average all inteprolated scans and multiply by weig_e_av in 
        # this way each scan has been weighted by int_weig_e/weig_e_av
        self.aver = weig_e_av * np.average(intrp, axis=0)
        # Compute STD of the interpolated scans
        self.aver_std = weig_e_av * np.std(intrp, axis=0)

    def reset(self, event):
        '''
        Reset graph for schoosing scans to starting conditions.
        Show all scans and set checked all buttons.
        '''
        # Clear graph
        self.avg_ln.set(visible=False)
        self.avg_std_fl.set(visible=False)

        stauts = self.checkbx.get_status()
        for i, stat in enumerate(stauts):
            if not stat:
                self.checkbx.set_active(i)
        # Show all spectra
        for i in range(len(self.lines)):
            self.lines[i].set(visible=True)
        for i in range(len(self.blines)):
            self.blines[i].set(visible=False)

        plt.draw()

    def finish_but(self, event):
        '''
        Close the figure on pressing the button Finish.
        '''
        plt.close()

    def edge_norm(self, guiobj, enrg, e_edge, e_pe, e_poste, pe_rng):
        '''
        Normalize energy scan data by the value at pre-edge energy.
        Also compute the  energy jump defined as the difference between
        the value at the edge and pre-edge energies respectively.

        This computations are implemented also considering baseline.
        If linear baseline is selected edge jump is computed considering
        the the height of data at edge energy from the stright line
        passing from pre-edge and post edge data.
        If asymmetrically reweighted penalized least squares baseline is
        selected the edge jump is calculated considering as the distance
        at edge energy between the averaged spectrum and baseline.

        Parameters
        ----------
        guiobj: GUI object
            Provides GUI dialogs.

        enrg : array
            Energy values of scan.

        e_edge : float
            Edge energy value.

        e_pe : float
            Pre-edge energy value.

        pe_rng : int
            Number of points constituting the semi-width of energy range
            centered at e_pe.

        pe_int : float
            Pre-edge value obtained from linear interpolation based on
            pre- and post-edge energies.

        Returns
        -------
        Set class attributes:
        pe_av : float
            value at pre-edge energy.

        norm : array
            self.aver scan normalized by value at pre-edge energy.

        norm_int : array
            Averaged data normalized by interpolated pre-edge value.

        ej : float
            edge-jump value.

        ej_norm : float
            edge-jump value normalized by value at pre-edge energy.

        ej_int : float
            edge-jump value computed with interpolated pre-edge value.

        ej_norm_int : float
            edge-jump computed and normalized by interpolated pre-edge
            value.

        Notes
        -----
        To reduce noise effects the value of scan at pre-edge energy is
        obtained computing an average over an energy range of width
        pe_rng and centered at e_pe pre-edge energy.
        The value of scan at edge energy is obtained by cubic spline
        interpolation of data (itp.UnivariateSpline with k=3 and s=0).
        '''
        # Index of the nearest element to pre-edge energy
        pe_idx = np.argmin((np.abs(enrg - e_pe)))
        # Left and right extremes of energy range for pre-edge average
        lpe_idx = int(pe_idx - pe_rng)
        rpe_idx = int(pe_idx + pe_rng + 1)

        # Average of values for computation of pre-edge
        self.pe_av = np.average(self.aver[lpe_idx: rpe_idx: 1])
        # STD on values range taken for computation of pre-edge
        self.pe_av_er = np.std(self.aver[lpe_idx: rpe_idx: 1])

        # Cubic spline interpolation of energy scan
        y_int = itp.UnivariateSpline(enrg, self.aver, k=3, s=0)
        # value at edge energy from interpolation
        y_edg = y_int(e_edge)

        # Take as error on value at the edge energy the interpolated
        # value of STD data
        # Cubic spline interpolation of STD energy scan
        y_int_er = itp.UnivariateSpline(enrg, self.aver_std, k=3, s=0)
        # value at edge energy from interpolation
        y_edg_err = y_int_er(e_edge)

        # Edge-jumps computations - no baseline
        self.ej = y_edg - self.pe_av
        self.ej_er = y_edg_err + self.pe_av_er
        self.ej_norm = self.ej / self.pe_av
        relerr = (self.ej_er / self.ej) + (self.pe_av_er / self.pe_av)
        self.ej_norm_er = self.ej_norm * relerr
        # Normalization by pre-edge value
        # It's a normalization not a physical quantity computation, so
        # no error propagation is computed
        self.norm = self.aver / self.pe_av
        self.norm_er = self.aver_std / self.pe_av

        # Edge-jumps computations - consider baseline
        if guiobj.bsl_int:
            # ArpLS baseline
            # Interpolation of pre-edge energy
            self.pe_av_int = self.bsl(e_edge)
        else:
            # Linear baseline
            # Interpolation of pre-edge energy
            x = [e_pe, e_poste]
            y = [y_int(e_pe), y_int(e_poste)]
            self.pe_av_int = lin_interpolate(x, y, e_edge)

        # Normalization by pre-edge value
        self.norm_int = self.aver / self.pe_av_int
        self.norm_int_er = self.aver_std / self.pe_av_int

        self.ej_int = y_edg - self.pe_av_int
        self.ej_int_er = y_edg_err
        self.ej_norm_int = self.ej_int / self.pe_av_int
        self.ej_norm_int_er = self.ej_int_er / self.pe_av_int


class EngyScan:
    '''
    Allow computation on energy scan data in order to obtain XNLD, XMCD,
    XNCD and XNXD spectra.

    Attributes
    ----------
    energy : array
        Common energy scale for positive and negative XAS scans data
        treatement.

    exper_edge : float
        Experimental edge energy.

    e_pe : float
        Pre-edge energy value.

    e_poste : float
        Post-edge energy value.

    pe_wdt : float
        Half-width of energy range for pre-edge average computation.

    pe_int : float
        Interpolated value of pre-edge considering pre- and post-edge
        energies using a linear approximation.

    avbsl : Univariate spline object
        ArpLS baseline of average xd.

    offest : float
        Offset given by the difference between expected and experimental
        edge energy.

    energycal : array
        Common energy scale calibrated considering the offset.

    ang_w_n : float
        Angle weight for XNLD computation.

    ang_w_d : float
        Angle weight for XNLD computation.

    xd : array
        X-Ray dichroism data.

    xd_err : array
        Experimental error of X-Ray dichroism.

    xd_aver : array
        Average of positive and negative XAS. In case of XNLD analysis a
        weighted average is considered.

    xd_aver_err : array
        Experimental error of xd_aver.

    xd_pc : array
        Percentage of X-Ray dichroism normalized by the average of
        positve and negative scans edge jumps respectively.

    xd_pc_er : array
        Error on xd_pc.

    xd_pc_int :array
        Same as xd_pc using data normalized with inerpolated edge-jump.

    xd_pc_int :array
        Error on xd_pc_int.

    pos_corr : array
        Normalized positive XAS spectrum weighted by the average of
        values at pre-edge energy of positive and negative spectra.
        For XNLD spectra the weight consider also the X-Ray beam's
        incidence angle.

    pos_corr_er : array
        Error on pos_corr.

    neg_corr : array
        Normalized negative XAS spectrum weighted by the average of
        values at pre-edge energy of positive and negative spectra.
        For XNLD spectra the weight consider also the X-Ray beam's
        incidence angle.

    neg_corr_er : array
        Error on neg_corr.

    pos_corr_int : array
        Same as pos_corr but using for the avereage of pre-edge values
        the ineterpolated values. For XNLD spectra the weight consider
        also the X-Ray beam's incidence angle.

    pos_corr_int_er : array
        Error on pos_corr_int.

    neg_corr_int : array
        Same as neg_corr but using for the avereage of pre-edge values
        the ineterpolated values. For XNLD spectra the weight consider
        also the X-Ray beam's incidence angle.

    neg_corr_int_er : array
        Error on neg_corr_int.

    xd_pc_av_ej : array
        Percentage  ofX-Ray dichroism normalized by edge jump of xd_aver
        scan.

    xd_pc_av_ej_er : array
        Error on xd_pc_av_ej 

    xd_pc_av_ej_int : array
        Same as xd_pc_av_ej but using for the edge-jump computation the 
        interpolated values of pre-edges.

    xd_pc_av_ej_int_er: array
        Error on xd_pc_av_ej_int

    Methods
    -------
    scan_average(guiobj, pos, neg, log_dt, pos_to_norm, neg_to_norm)
        Copmute averages of positive and negative scans.

    compt_mxd(guiobj, pos, neg, log_dt)
        Performs computation on energy scan data extracting XNLD, XMCD,
        XNXD or XNCD spctra.

    compt_xd(self, guiobj, pos, neg, log_dt)
        Compute not normalized X-Ray Dichroism.

    edges(guiobj, pos, neg, log_dt)
        Set values of edge and pre-edge energies.

    compt_xd_pc(guiobj, pos, neg)
        Compute the percentage of X-Ray Dichroism normalized for
        edge-jump.

    compt_pe_corr(guiobj, pos, neg)
        Calculate xd percentage normalizing with edge-jump obtained from
        weighted average of positive and negative spectra.

    comp_xd_pc_av_ej(self, log_dt)
        Compute percentage of X-Ray Dichroism normalized for edge-jump
        of xd_aver spectrum.
    '''

    def __init__(self, guiobj, e_scale, pos, neg, log_dt,
                 pos_to_norm=ScanData(), neg_to_norm=ScanData()):
        '''
        At instantiation scan_average and compt_xd method are called and
        energy attribute is setted, so upon its creation an EngyScan
        object collects all the informations about X-Ray dichroism.

        Parameters
        ----------
        guiobj : GUI obj
            Provides GUI dialogs.

        e_scale : array
            Common energy scale for data analysis.

        pos : ScanData obj
            Positive scans (CR for XMCD and XNCD, LH for XNLD).

        neg : ScanData obj
            Negative scnas (CL for XMCD and XNCD, LV for XNLD).       

        log_dt : dict
            Collect data for logfile.

        pos_to_norm : ScanData obj
            Positive scans (CR for XMCD and XNCD, LH for XNLD) data to
            be normalized by reference if present. If no reference data
            are present pass empty ScanData object.

        neg_to_norm : ScanData obj
            Negative scans (CR for XMCD and XNCD, LH for XNLD) data to
            be normalized by reference if present. If no reference data
            are present pass empty ScanData object.
        '''
        self.energy = e_scale

        if guiobj.interactive:
            # Ask for process baseline method
            guiobj.ask_bsl_interp()
        else:
            # For non interactive linear baseline is the default
            guiobj.bsl_int = False

        self.scan_average(guiobj, pos, neg, log_dt, pos_to_norm, neg_to_norm)

        self.compt_mxd(guiobj, pos, neg, log_dt)

    def scan_average(self, guiobj, pos, neg, log_dt, pos_to_norm, neg_to_norm):
        '''
        Copmute averages of positive and negative scans.
        If reference data are considered compute also the averages of
        reference scans.

        Parameters
        ----------
        guiobj : GUI obj
            Provides GUI dialogs.

        pos : ScanData obj
            Positive scans (CR for XMCD and XNCD, LH for XNLD).

        neg : ScanData obj
            Negative scnas (CL for XMCD and XNCD, LV for XNLD).     

        log_dt : dict
            Collect data for logfile.

        pos_to_norm : ScanData obj
            Positive scans (CR for XMCD and XNCD, LH for XNLD) data to
            be normalized by reference if present.

        neg_to_norm : ScanData obj
            Negative scans (CR for XMCD and XNCD, LH for XNLD) data to
            be normalized by reference if present.

        Return
        ------
        Instantiate aver attribute of positive and negative ScanData
        objects. In case referecnce data are passed aver attribute is
        the normalization ***_to_norm by reference data.

        Add keys to log_dt

        pos_chs : list (str) with positive chosen scan
        neg_chs : list (str) with negative chosen scan

        '''
        if guiobj.interactive:
            # only in interactive mode ask for normalization
            weig_dt = guiobj.ask_weig()
        else:
            weig_dt = False

        # If no reference data pos_to_norm and neg_to_norm are empty
        # ScanData object
        if pos_to_norm.raw_imp.empty or neg_to_norm.raw_imp.empty:
            # Computes averages of positive and negative polarization
            # scans
            guiobj.infile_ref = False

            pos.man_aver_e_scans(guiobj, self.energy, weig_dt, log_dt)
            neg.man_aver_e_scans(guiobj, self.energy, weig_dt, log_dt)
        else:
            # If present computes averages of positive and negative
            # polarization  of reference scans and normalize data for
            # them.
            # In this case pos and neg contains reference data and
            # pos_to_norm and neg_to_norm are data to be normalized.
            # They are supposed to be already analyzed so aver attribute
            # is already setted.
            guiobj.infile_ref = True
            pos.man_aver_e_scans(guiobj, self.energy, weig_dt, log_dt)
            neg.man_aver_e_scans(guiobj, self.energy, weig_dt, log_dt)

            # Just substitute pos.aver and neg.aver with normalized data
            pos.aver = pos_to_norm.aver / pos.aver
            neg.aver = neg_to_norm.aver / neg.aver

            guiobj.infile_ref = False

        # Add keys with chosen scans to log_dt
        log_dt['pos_chs'] = pos.chsn_scns
        log_dt['neg_chs'] = neg.chsn_scns

    def compt_mxd(self, guiobj, pos, neg, log_dt):
        '''
        Performs computation on energy scan data extracting XNLD, XMCD,
        XNXD or XNCD spctra.
        Given two ScanData object (one for positive and one for negative
        scans):
        - compute X-Ray Dichroism as negative - positive;
        - compute the arithmetical mean of positive and negative spectra
          for XMCD, XNCD, XNXD data and the weighted average for the
          angle for XNLD data (see Notes);
        - compute percentage of X-Ray Dichroism normalized by the
          average edge-jump for XMCD, XNCD, XNXD data while for XNLD an
          angle-weighted average is considered (see Notes);
        - compute percentage of X-Ray Dichroism normalized by edge-jump
          of the weighted average of positive and negative spectra.

        Parameters
        ----------
        guiobj : GUI obj
            Provides GUI dialogs.

        pos : ScanData obj
            Positive scans (CR for XMCD and XNCD, LH for XNLD).

        neg : ScanData obj
            Negative scnas (CL for XMCD and XNCD, LV for XNLD).

        log_dt : dict
            Collect data for logfile.

        Returns
        -------
        Add keys to log_dt:

        pos_ej : float, edge-jump value for postitive scans.
        pos_ej_int : float, edge-jump value for postitive scans,
            interpolated value of pre-edge energy is used.
        neg_ej : float, edge-jump value for postitive scans.
        neg_ej_int : float, edge-jump value for postitive scans,
            interpolated value of pre-edge energy is used.
        Notes
        -----
        In XNLD analsysis LH and LV spectra are weighted by the angle t
        between the sample surface and the X-Ray incident beam
        direction.

        - XNLD average = (LH + (2 * cos(t)^2 - sin(t)^2) * LV) / 
                         3 * cos(t)^2

        - XNLD (%) = 100 * 3 * cos(t)^2 * (LV/LV_pe - LH/LH_pe) /
                     (LH/LH_pe + (2 * cos(t)^2 - sin(t)^2) * LV/LV_pe)

        - XNLD norm by XNLD_aver edge-jump (%) = 
                    100 * ((LH_pe + (2 * cos(t)^2 - sin(t)^2) * LV_pe) / 
                    3 * cos(t)^2) * ((LV / LV_pe) - (LH / LH_pe)) /
                    (XNLD_aver_edge - XNLD_aver_pe))
        '''
        # Computes not normalized X-Ray Dichroism
        self.compt_xd(guiobj, pos, neg, log_dt)

        self.edges(guiobj, pos, neg, log_dt)

        # Normalize spectra
        pos.edge_norm(guiobj, self.energycal, self.exper_edge, self.e_pe,
                        self.e_poste, self.pe_wdt)
        neg.edge_norm(guiobj, self.energycal, self.exper_edge, self.e_pe,
                        self.e_poste, self.pe_wdt)

        log_dt['pos_ej'] = pos.ej
        log_dt['pos_ej_er'] = pos.ej_er
        log_dt['pos_ej_int'] = pos.ej_int
        log_dt['pos_ej_int_er'] = pos.ej_int_er
        log_dt['neg_ej'] = neg.ej
        log_dt['neg_ej_er'] = neg.ej_er
        log_dt['neg_ej_int'] = neg.ej_int
        log_dt['neg_ej_int_er'] = neg.ej_int_er

        # Compute percentage X-Ray Dichroism normalized for edge-jump
        self.compt_xd_pc(guiobj, pos, neg, log_dt)
        self.compt_pe_corr(guiobj, pos, neg)
        self.comp_xd_pc_av_ej(guiobj, log_dt)

    def compt_xd(self, guiobj, pos, neg, log_dt):
        '''
        Compute not normalized X-Ray Dichroism.

        Parameters
        ----------
        guiobj : GUI obj
            Provides GUI dialogs.

        pos_scan : ScanData obj
            Positive scans (CR for XMCD and XNCD, LH for XNLD).

        neg_scan : ScanData obj
            Negative scnas (CL for XMCD and XNCD, LV for XNLD).        

        log_dt : dict
            Collect data for logfile.

        Return
        ------
        Set attributes:

        xd : array
            X-Ray Dichroism as negative - positive scans.

        xd_aver : array
            Average of XAS spectra. Arithmetical mean of positive and
            negative scans for XMCD and XNCD analysis. Weigthed average
            by angle for XNLD analysis (see Notes).

        ang_w_n : float
            Angle weight for XNLD computation.

        ang_w_d : float
            Angle weight for XNLD computation.

        Notes
        -----
        In XNLD analsysis LH and LV spectra are weighted by the angle t
        between the sample surface and the X-Ray incident beam's
        direction.

        XNLD average = (LH + (2 * cos(t)^2 - sin(t)^2) * LV) /
                        3 * cos(t)^2
        '''
        self.xd = neg.aver - pos.aver
        self.xd_err = neg.aver_std + pos.aver_std
        if guiobj.analysis in guiobj.type['xnld']:
            # If XNLD the angle must be considered for weighted mean
            # computation
            theta = log_dt['bm_angle']  # retrive angle from log table

            # Numerator and denominator terms of angle weight
            self.ang_w_n = 2 * (np.cos(theta))**2 - (np.sin(theta))**2
            self.ang_w_d = 3 * (np.cos(theta))**2

            self.xd_aver = (pos.aver + (self.ang_w_n*neg.aver)) / self.ang_w_d
            self.xd_aver_err = ((pos.aver_std + (self.ang_w_n*neg.aver_std))
                                / self.ang_w_d)
        else:
            self.xd_aver = (pos.aver + neg.aver) / 2
            self.xd_aver_err = (pos.aver_std + neg.aver_std) / 2

    def edges(self, guiobj, pos, neg, log_dt):
        '''
        Set values of edge energy, pre-edge energy, pre-edge energy
        range (used to compute average of spectra at pre-edge energy)
        and post-edge energy.

        If in interactive mode a GUI is provided to set the experimental
        edge energy, pre-edge energy, post-edge energy and and
        half-width of interval adopted for pre-edge average.
        Otherwise for edge energy is considered the experimental
        computed one from minimization, for pre-edge end post-edge
        energies the ones provided from edge file are taken for good,
        and for the half-width of interval for pre-edge average 4 points
        is considered.

        Parameters
        ----------
        guiobj : GUI object
            Provides GUI dialogs.

        pos : ScanData obj
            Positive scans (CR for XMCD and XNCD, LH for XNLD).

        neg : ScanData obj
            Negative scnas (CL for XMCD and XNCD, LV for XNLD).

        log_dt : dict
            Collect data for logfile.

        Return
        ------
        Set class attributes

        exper_edge : float
            Experimental edge energy.

        e_pe : float
            Pre-edge energy value.

        e_poste : float
            Post-edge energy value.

        pe_wdt : float
            Half-width of energy range for pre-edge average computation.

        offset : float
            Offset given by the difference between expected and
            experimental edge energy.

        energycal : array
            Common energy scale calibrated considering the offset. In
            this way absorption edge of spectra falls at the expected
            energy.

        Add keys to log_dt

        exper_edge : experimental edge energy.
        setted_pedg : setted pre-edge energy.
        setted_postedg : setted post-edge energy.
        recal : bool, if True energy recalibration has been done.
        offset : energy offset adopted for energy recalibration.

        Notes
        -----
        Indication of experimental value of edge energy is provided
        finding the stationary points of experimental spectrum. This is
        obtained interpolating experimental data. Method used for
        interpolation is  scipy.interpolate.UnivariateSpline of order 3.
        '''
        # Retrive from log_dt edge, pre-edge energy and range for
        # pre-edge average

        # Check that pre-edge energy is included in the considered
        # energy range. If not, an energy value close to the nearest
        # range endpoint is considered
        if float(log_dt['PreEdge_en']) <= self.energy[5]:
            pe_e = self.energy[5]
        elif float(log_dt['PreEdge_en']) >= self.energy[-5]:
            pe_e = self.energy[-5]
        else:
            pe_e = float(log_dt['PreEdge_en'])

        sel_edg = [log_dt['Edge_en'], pe_e]

        if not guiobj.bsl_int:
            # Only for baseline linear interpolation post-edge is
            # Check that post edge energy is included in the considered
            # energy range
            if float(log_dt['PostEdge_en']) <= self.energy[0]:
                pste_e = self.energy[1]
            elif float(log_dt['PostEdge_en']) >= self.energy[-1]:
                pste_e = self.energy[-2]
            else:
                pste_e = float(log_dt['PostEdge_en'])
            sel_edg.append(pste_e)

        y = self.xd
        y_aver = self.xd_aver

        # Order 3 interpolation of data
        # opt.minimize_scalar serch for function minimum so negative
        # absolute value of interpolated data is considered.
        y_int_for_edge = itp.UnivariateSpline(self.energy, -abs(y), k=3, s=0)
        # Interpolation of average and xd spectra
        y_int_aver = itp.UnivariateSpline(self.energy, y_aver, k=3, s=0)
        y_int = itp.UnivariateSpline(self.energy, y, k=3, s=0)
        # Bounds for  minimum search - 5 eV window is considered
        u_bnd = log_dt['Edge_en'] + 2.5
        l_bnd = log_dt['Edge_en'] - 2.5
        min_y = opt.minimize_scalar(y_int_for_edge, bounds=(l_bnd, u_bnd),
                                    method='bounded')
        if guiobj.interactive:
            # x value which minimize y is passed as experimental edge
            # energy
            if guiobj.bsl_int:
                edgs = guiobj.set_edges_arpls(sel_edg, min_y.x, self.energy, y,
                                y_aver, y_int, y_int_aver, pos, neg, log_dt)
                self.exper_edge = edgs[0]
                self.e_pe = edgs[1]
                self.e_poste = None  # not needed in the following
                self.pe_wdt = int(edgs[3])
                recal = edgs[4]
                self.avbsl = edgs[5]
            else:
                # Linear interpolation of baseline
                edgs = guiobj.set_edges_lin(sel_edg, min_y.x, self.energy, y,
                                            y_aver, y_int, y_int_aver)
                self.exper_edge = edgs[0]
                self.e_pe = edgs[1]
                self.e_poste = edgs[2]
                self.pe_wdt = int(edgs[3])
                recal = edgs[4]
        else:
            # If not interactive for exper_edge the result of
            # minimization is considered
            self.exper_edge = float(min_y.x)
            self.e_pe = pe_e
            self.e_poste = pste_e
            self.pe_wdt = 4
            recal = False

        # New energy scale is created adding an offset given by the
        # difference  between expected and experimental edge energies if
        # recal is True, otherwise no offset is added
        if recal:
            self.offset = self.exper_edge - log_dt['Edge_en']
        else:
            self.offset = 0

        self.energycal = self.energy - self.offset

        # Log data related to edge and pre-edge energy setted
        log_dt['exper_edge'] = self.exper_edge
        log_dt['setted_pedg'] = self.e_pe
        if guiobj.bsl_int:
            log_dt['lambda'] = edgs[2]
        else:
            log_dt['setted_postedg'] = self.e_poste
        log_dt['recal'] = recal
        log_dt['offset'] = self.offset

    def compt_xd_pc(self, guiobj, pos, neg, log_dt):
        '''
        Compute the percentage of X-Ray Dichroism normalized for
        edge-jump.
        For XMCD and XNCD it is obtained as the difference between
        normalized by pre-edge negative and positve scans normalized by
        the arithmetical mean of the two respective edge jumps. For XNLD
        the average weighted by X-Ray beam's angle of incidence is
        considered (see Notes).

        Parameters
        ----------
        guiobj : GUI obj
            Provides GUI dialogs.

        pos : ScanData obj
            Positive scans (CR for XMCD and XNCD, LH for XNLD).

        neg : ScanData obj
            Negative scnas (CL for XMCD and XNCD, LV for XNLD).

        log_dt : dict
            Collect data for logfile.

        Return
        ------
        Set attributes
        xd_pc : array
            Percentage of X-Ray Dichroism normalized for edge-jump.
            For XMCD and XNCD it is obtained as the difference between
            normalized by pre-edge negative and positve scans normalized
            by the arithmetical mean of the two respective edge jumps.
            For XNLD the average weighted by X-Ray beam's angle of
            incidence is considered (see Notes).

        xd_pc_int : array
            Same as xd_pc using data normalized with inerpolated
            edge-jump.

        Add keys to logdt
        edge_pc : percentage of X-Ray Dichroism normalized for edge-jump
            at experimental edge energy
        edge_pc_er : error on edge_pc
        edge_pc_int : percentage of X-Ray Dichroism normalized with
            inerpolated edge-jump at experimental edge energy
        edge_pc_int_er : error on edge_pc_int 

        Notes
        -----
        In XNLD analsysis LH and LV spectra are weighted by the angle t
        between the sample surface and the X-Ray incident beam's
        direction.

         XNLD (%) = 100 * 3 * cos(t)^2 * (LV/LV_pe - LH/LH_pe) /
                     (LH/LH_pe + (2 * cos(t)^2 - sin(t)^2) * LV/LV_pe)
        '''
        # Compute mean
        if guiobj.analysis in guiobj.type['xnld']:
            # Percentage X-Ray dichroism normalized for edge jump
            self.xd_pc = (100 * self.ang_w_d * (neg.norm - pos.norm)
                          / (pos.ej_norm + self.ang_w_n * neg.ej_norm))
            # Error propagation
            d_rel = ((pos.ej_norm_er + self.ang_w_n * neg.ej_norm_er)
                    / (pos.ej_norm + self.ang_w_n * neg.ej_norm))
            n_rel = ((neg.norm_er + pos.norm_er) / (neg.norm - pos.norm))            

            self.xd_pc_int = (100 * self.ang_w_d * (neg.norm_int
                            - pos.norm_int) / (pos.ej_norm_int + self.ang_w_n
                            * neg.ej_norm_int))
            # Error propagation
            d_i_rel = ((pos.ej_norm_int_er + self.ang_w_n * neg.ej_norm_int_er)
                    / (pos.ej_norm_int + self.ang_w_n * neg.ej_norm_int))
            n_i_rel = ((neg.norm_int_er + pos.norm_int_er)
                    / (neg.norm_int - pos.norm_int))            
        else:
            # Percentage X-Ray dichroism normalized for edge jump
            self.xd_pc = 200 * ((neg.norm - pos.norm)
                                / (neg.ej_norm + pos.ej_norm))
            # Error propagation
            d_rel = ((neg.ej_norm_er + pos.ej_norm_er)
                    / (neg.ej_norm + pos.ej_norm))
            # neg/pos.norm_int_er = aver.std
            n_rel = ((neg.norm_er + pos.norm_er) / (neg.norm - pos.norm))

            self.xd_pc_int = 200 * ((neg.norm_int - pos.norm_int)
                                    /(neg.ej_norm_int + pos.ej_norm_int))
            # Error propagation
            d_i_rel = ((pos.ej_norm_int_er + neg.ej_norm_int_er)
                    / (neg.ej_norm_int + pos.ej_norm_int))
            n_i_rel = ((neg.norm_int_er + pos.norm_int_er)
                    / (neg.norm_int - pos.norm_int))

        xd_pc_errel = n_rel + d_rel
        self.xd_pc_er = self.xd_pc * xd_pc_errel
        xd_pc_int_errel = n_i_rel + d_i_rel
        self.xd_pc_int_er = self.xd_pc_int * xd_pc_int_errel

        # Compute dichroism percentage at edge energy
        in_xd_pc = itp.UnivariateSpline(self.energycal, self.xd_pc, k=1, s=0)
        in_xd_pc_er = itp.UnivariateSpline(self.energycal, self.xd_pc_er, k=1,
                                            s=0)
        log_dt['edge_pc'] = np.abs(in_xd_pc(self.exper_edge))
        log_dt['edge_pc_er'] = in_xd_pc_er(self.exper_edge)
        in_xd_pc_int = itp.UnivariateSpline(self.energycal, self.xd_pc_int,
                                            k=1, s=0)
        in_xd_pc_int_er = itp.UnivariateSpline(self.energycal,
                                                self.xd_pc_int_er, k=1, s=0)
        log_dt['edge_pc_int'] = np.abs(in_xd_pc_int(self.exper_edge))
        log_dt['edge_pc_int_er'] = in_xd_pc_int_er(self.exper_edge)


    def compt_pe_corr(self, guiobj, pos, neg):
        '''
        Calculate xd percentage normalizing with edge-jump obtained from
        weighted average of positive and negative spectra.

        Parameters
        ----------
        guiobj : GUI obj
            Provides GUI dialogs.

        pos : ScanData obj
            Positive scans (CR for XMCD and XNCD, LH for XNLD).

        neg : ScanData obj
            Negative scnas (CL for XMCD and XNCD, LV for XNLD).   

        Return
        ------
        Set attributes:

        pos_corr : array
            Average of positive spectra weighted by the ratio between
            the mean of positive pre-edge values and the mean of all
            pre-edge values.
            For XMCD and XNCD the mean of positive and negative
            pre-edges is arithmetical while for XNLD the average
            weighted by the angle of incidence- of X-Ray beam is
            considered (see Notes).

        pos_corr_int : array
            Same as pos_corr using interpolated pre-edge values.

        neg_corr : array
            Average of negative spectra weighted by the ratio between
            the mean of positive pre-edge values and the mean of all
            pre-edge values.
            For XMCD and XNCD the mean of positive and negative
            pre-edges is arithmetical while for XNLD the average
            weighted by angle of incidence if the X-Ray beam is
            considered (see Notes).

        neg_corr_int : array
            Same as neg_corr using interpolated pre-edge values.

        pos_corr_int : array
            Positive scan normalized by weighted average of interpolated
            edge-jump.

        neg_corr_int : array
            Negative scan normalized by weighted average of interpolated
            edge-jump.

        Notes
        -----
        In XNLD analsysis LH and LV spectra are weighted by the angle t
        between the sample surface and the X-Ray incident beam's
        direction.

        XNLD average = (LH + (2 * cos(t)^2 - sin(t)^2) * LV) /
                        3 * cos(t)^2
        '''
        # Compute mean
        if guiobj.analysis in guiobj.type['xnld']:
            # Angle weighted average of pre-edges values
            av_pe_av = (pos.pe_av + self.ang_w_n * neg.pe_av) / self.ang_w_d
            av_pe_av_er = ((pos.pe_av_er + self.ang_w_n * neg.pe_av_er)
                        / self.ang_w_d)
            av_pe_int = ((pos.pe_av_int + self.ang_w_n * neg.pe_av_int)
                         / self.ang_w_d)
        else:
            # Average of pre-edges values
            av_pe_av = (pos.pe_av + neg.pe_av) / 2
            av_pe_av_er = (pos.pe_av_er + neg.pe_av_er) / 2
            av_pe_int = (pos.pe_av_int + neg.pe_av_int) / 2

        self.pos_corr = pos.aver * av_pe_av / pos.pe_av
        self.pos_corr_er = (
            ((pos.aver_std / pos.aver) + (av_pe_av_er / av_pe_av)
            + (pos.pe_av_er / pos.pe_av)) * self.pos_corr)
        self.neg_corr = neg.aver * av_pe_av / neg.pe_av
        self.neg_corr_er = (
            ((neg.aver_std / neg.aver) + (av_pe_av_er / av_pe_av)
            + (neg.pe_av_er / neg.pe_av)) * self.neg_corr)

        self.pos_corr_int = pos.aver * av_pe_int / pos.pe_av_int
        self.pos_corr_int_er = pos.aver_std * av_pe_int / pos.pe_av_int
        self.neg_corr_int = neg.aver * av_pe_int / neg.pe_av_int
        self.neg_corr_int_er = neg.aver_std * av_pe_int / neg.pe_av_int

    def comp_xd_pc_av_ej(self, guiobj, log_dt):
        '''
        Compute percentage of X-Ray Dichroism normalized for edge-jump
        of xd_aver spectrum. It is computed employing the weighted
        averages pos_corr and neg_corr data. Edge-jump is calculated
        interpolating the values at edge and pre-edge energy of xd_aver
        spectrum with linear spline method.
        For XNLD the average weighted by the angle of incidence of the
        X-Ray beam is considered (see Notes).

        Parameters
        ----------
        guiobj : GUI obj
            Provides GUI dialogs.

        log_dt : dict
            Collect data for logfile.

        Returns
        -------
        Set class attributes:

        xd_pc_av_ej : array
            Percentage of X-Ray Dichroism normalized for edge-jump of
            xd_aver spectrum.

        xd_pc_av_ej_int : array
            Same as xd_pc_av_ej but using for the edge-jump computation
            the interpolated values of pre-edges.

        Add keys to log_dt
        xas_aver_ej : float, edge_jump value of xd average.
        xas_aver_ej_int : float, edge jump value of xd average
            considering interpolated value of pre-edge.

        Notes
        -----
        In XNLD analsysis LH and LV spectra are weighted by the angle t
        between the sample surface and the X-Ray incident beam's
        direction.

        XNLD norm by XNLD_aver edge-jump (%) = 
                    100 * ((LH_pe + (2 * cos(t)^2 - sin(t)^2) * LV_pe) / 
                    3 * cos(t)^2) * ((LV / LV_pe) - (LH / LH_pe)) /
                    (XNLD_aver_edge - XNLD_aver_pe))
        '''
        # Linear spline interpolation of xd_aver spectrum in order to
        # determine edge jump
        xd_aver_inter = itp.UnivariateSpline(self.energycal, self.xd_aver, k=1,
                                             s=0)
        # Ineterpolated values at edge and pre-edge energies
        edg_val_xd_aver = xd_aver_inter(self.exper_edge)
        pedg_val_xd_aver = xd_aver_inter(self.e_pe)

        # Linear spline interpolation of xd_aver_err array in order to
        # determine edge jump error
        xd_aver_inter_er = itp.UnivariateSpline(
                                    self.energycal, self.xd_aver_err, k=1, s=0)
        # Ineterpolated values of errors at edge and pre-edge energies
        edg_val_xd_aver_er = xd_aver_inter_er(self.exper_edge)
        pedg_val_xd_aver_er = xd_aver_inter_er(self.e_pe)

        if guiobj.bsl_int:
            # Compute value at egde of ArPLS baseline of average xd
            pedg_val_xd_aver_int = self.avbsl(self.exper_edge)
            pedg_val_xd_aver_int_er = 0
        else:
            # Compute linear interpolation of pre-edge for average xd
            pstedg_val_xd_aver = xd_aver_inter(self.e_poste)

            x = [self.e_pe, self.e_poste]
            y = [pedg_val_xd_aver, pstedg_val_xd_aver]
            pedg_val_xd_aver_int = lin_interpolate(x, y, self.exper_edge)

            # For linear interpolation error on edge and predge values
            # can be included in propagation for edge-jump computation.
            # Lines from pe+er:poste-er and pe-er:poste+er are
            # considered. Values at egde energy of both line is
            # interpolated then the greatest deviation from
            # pedg_val_xd_aver_int is taken as error estimation.

            # Ineterpolated values of errors at pre-edge and post-edge
            # energies            
            pstedg_val_xd_aver_er = xd_aver_inter_er(self.e_poste)
            y1 = [pedg_val_xd_aver + pedg_val_xd_aver_er,
                pstedg_val_xd_aver - pstedg_val_xd_aver_er]
            y2 = [pedg_val_xd_aver - pedg_val_xd_aver_er,
                pstedg_val_xd_aver + pstedg_val_xd_aver_er]
            y1_ext = lin_interpolate(x, y1, self.exper_edge)
            y2_ext = lin_interpolate(x, y2, self.exper_edge)
            pedg_val_xd_aver_int_er = np.amax(np.abs([
                                            pedg_val_xd_aver_int - y1_ext,
                                            pedg_val_xd_aver_int - y2_ext]))

        # xd percentage normalized for edge-jump of xd_aver
        self.xd_pc_av_ej = (100 * (self.neg_corr - self.pos_corr)
                            / (edg_val_xd_aver - pedg_val_xd_aver))
        d_er = ((edg_val_xd_aver_er + pedg_val_xd_aver_er)
                / (edg_val_xd_aver - pedg_val_xd_aver))
        n_er = (100 * (self.neg_corr_er + self.pos_corr_er)
                /(100 * (self.neg_corr - self.pos_corr)))
        self.xd_pc_av_ej_er = (n_er + d_er) * self.xd_pc_av_ej

        self.xd_pc_av_ej_int = (100 * (self.neg_corr_int - self.pos_corr_int)
                                / (edg_val_xd_aver - pedg_val_xd_aver_int))
        d_er = ((edg_val_xd_aver_er + pedg_val_xd_aver_int_er)
                / (edg_val_xd_aver - pedg_val_xd_aver_int))
        n_er = (100 * (self.neg_corr_int_er + self.pos_corr_int_er)
                / (self.neg_corr_int - self.pos_corr_int))
        self.xd_pc_av_ej_int_er = (n_er + d_er) * self.xd_pc_av_ej_int

        log_dt['xas_aver_ej'] = edg_val_xd_aver - pedg_val_xd_aver
        log_dt['xas_aver_ej_er'] = edg_val_xd_aver_er + pedg_val_xd_aver_er
        log_dt['xas_aver_ej_int'] = edg_val_xd_aver - pedg_val_xd_aver_int
        log_dt['xas_aver_ej_int_er'] = (
            edg_val_xd_aver_er + pedg_val_xd_aver_int_er)


def e_scale(guiobj, pos, neg, log_dt, pos_ref, neg_ref):
    '''
    Create the energy scale used for data analysis.
    Range is selected considering the intersection of the energy ranges
    of all the scans provided: the low-end is the highest of the minimum
    values and the high-end is the lowest of the maximum values between
    all the energy arrays.

    The number of points of the energy scale returned is by default the 
    average of the number of points of the energy scales of data scans.
    If in interactive mode GUI dialogues are provided to set the 
    number of points.

    Parameters
    ----------
    guiobj : GUI obj
        Provides GUI dialogs.

    pos : ScanData obj
        Positive scans (CR for XMCD and XNCD, LH for XNLD).

    neg : ScanData obj
        Negative scnas (CL for XMCD and XNCD, LV for XNLD).   

    log_dt : dict
        Collect data for logfile.

    pos_ref : ScanData obj
        Positive scans (CR for XMCD and XNCD, LH for XNLD) from
        reference.

    neg_ref : ScanData obj
        Negative scans (CR for XMCD and XNCD, LH for XNLD) from
        reference.

    Return
    ------
    array
        Common energy scale for positive and negative XAS scans data
        treatement (if reference scans are considered in confobj also
        their energy scale will be considered).

    Add e_scale key to log_dt with min, max and number of points of
    energy scale. 
    '''
    efirst_list = []
    elast_list = []
    e_len = []

    # Sort energy data arrays and count their elements
    for i in pos.idx:
        efirst_list.append(np.sort(pos.raw_imp['E' + i])[0])
        elast_list.append(np.sort(pos.raw_imp['E' + i])[-1])
        e_len.append(pos.raw_imp['E' + i].size)
    for i in neg.idx:
        efirst_list.append(np.sort(neg.raw_imp['E' + i])[0])
        elast_list.append(np.sort(neg.raw_imp['E' + i])[-1])
        e_len.append(neg.raw_imp['E' + i].size)
    # If reference data are present computes common energy scale for all
    # data
    if not (pos_ref.raw_imp.empty or neg_ref.raw_imp.empty):
        for i in pos_ref.idx:
            efirst_list.append(np.sort(pos_ref.raw_imp['E' + i])[0])
            elast_list.append(np.sort(pos_ref.raw_imp['E' + i])[-1])
            e_len.append(pos_ref.raw_imp['E' + i].size)
        for i in neg_ref.idx:
            efirst_list.append(np.sort(neg_ref.raw_imp['E' + i])[0])
            elast_list.append(np.sort(neg_ref.raw_imp['E' + i])[-1])
            e_len.append(neg_ref.raw_imp['E' + i].size)

    # Compute min, max and default langth of energy range
    e_min = np.around(np.amax(efirst_list), 1)
    e_max = np.around(np.amin(elast_list), 1)
    e_av_len = np.around(np.average(e_len), 0)

    # Set number of points of energy scale
    if guiobj.interactive:
        n_points = guiobj.num_pnts(e_av_len)
    else:
        n_points = int(e_av_len)

    e_step = (e_max - e_min) / (n_points - 1)

    log_dt['e_scale'] = [e_min, e_max, n_points, e_step]

    return np.linspace(e_min, e_max, n_points)


def lin_interpolate(x, y, x0):
    '''
    Linear interpolation of the value at x0
    A line is determined considering the couple of points passed in
    x and y.
    The value at x0 point is computed.

    Parameters
    ----------
    x : array
        x values of point 1 and point 2.

    y : array
        y values of point 1 and point 2.

    x0 : float
        point where linear interpolation is computed.

    Return
    ------
    float, interpolated value.
    '''
    lin_spl_int = itp.UnivariateSpline(x, y, k=1, s=0)
    return(lin_spl_int(x0))
