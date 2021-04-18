"""
pyDichroX_datatreat.py

Classes and methods for XMCD, XNCD, XNLD and hysteresys data analysis.

Classes
--------
ScanData : Collect raw data for energy scan experiments and provides methods to
    average them.

EngyScan : Allow computation on energy scan data extracting XNLD, XMCD or XNCD
    spctra.

Methods
-------
lin_interpolate(x, y, x0):
        Linear interpolation of the value at x0
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


class ScanData:
    '''
    Collect raw data for energy scan experiments and provides methods to
    average them.

    Attributes
    ----------
    label : list (str)
        Labels for graphs

    idx : list (str)
        Scan indexes

    raw_imp : pandas DataFrame
        Collect imported raw data

    aver : array
        Data average of selected scans from raw_imp

    dtype : str
        Identifies the data collected, used for graph labelling:
        sigma+, sigma- for XMCD
        CR, CL for XNCD
        H-, H+ for XNXD
        LH, LV fot XNLD

    chsd_scns : list (str)
        Labels of chosed scans for the analysis

    pe_av : float
        Value of spectra at pre-edge energy. It is obtained from averaging data
        in an defined energy range centered at pre-edge energy

    pe_av_int : float
        Pre-edge value obtained from linear interpolation considering pre-edge
        and post-edge energies

    norm : array
        Averaged data normalized by value at pre-edge energy

    norm_int : array
        Averaged data normalized by interpolated pre-edge value

    ej : float
        Edge-jump value

    ej_norm : float
            edge-jump value normalized by value at pre-edge energy

    ej_int : float
        edge-jump value computed with interpolated pre-edge value

    ej_norm_int : float
        edge-jump computed and normalized by interpolated pre-edge value

    Methods
    -------
    man_aver_e_scans(guiobj, enrg)
        Manage the choice of scans to be averaged and return the average of
        selected scans.

    aver_e_scans(enrg, chsd, interactive)
        Performe the average of data scans.

    edge_norm(enrg, e_edge, e_pe, pe_rng, pe_int)
        Normalize energy scan data by value at pre-edge energy and compute
        edge-jump.
    '''

    def __init__(self):
        '''
        Initialize attributes label, idx, and raw_imp.
        '''
        self.label = []
        self.idx = []
        self.raw_imp = pd.DataFrame()

    def man_aver_e_scans(self, guiobj, enrg):
        '''
        Manage the choice of scans to be averaged and return the average of
        selected scans.
        Works only with energy scans (i.e. XMCD, XNCD, XNXD and XNLD).

        Parameters
        ----------
        guiobj : GUI object
            Provides GUI dialogs.

        enrg : array
            Energy values at which average is calculated.

        Return
        ------
        Set class attributes:
        aver : array
            Average values of the choosed scans.

        chsd_scns : list
            Labels of chosed scans for the analysis (for log purpose)
        '''
        if guiobj.interactive:  # Interactive choose of scans
            # Loop until choice is confirmed
            isok = 0
            while not isok:
                # Plot all raw data
                plt.figure(1)
                for i in self.idx:
                    e_col = 'E' + i
                    plt.plot(self.raw_imp[e_col], self.raw_imp[i],
                             label=self.label[self.idx.index(i)])
                plt.xlabel('E (eV)')
                plt.ylabel(self.dtype)
                plt.legend()
                plt.show()

                # Dialogue to choose data to be averaged
                chsd = guiobj.chs_scns(self.label)

                # Compute average of choosed scans
                avgd = self.aver_e_scans(enrg, chsd, guiobj.interactive)

                # Ask for confirmation
                isok = guiobj.confirm_choice()
        else:
            # Not-interactive mode: all scans except 'Dummy Scans' are
            # evaluated.
            chsd = []

            for lbl in self.label:
                # Check it is not a 'Dummy Scan' and append corresponding scan
                # number in chosed scan list.
                if not ('Dummy' in lbl):
                    chsd.append(self.idx[self.label.index(lbl)])

            avgd = self.aver_e_scans(enrg, chsd, guiobj.interactive)

        self.aver = avgd
        self.chsd_scns = chsd  # Collect choosed scans for log purpose

    def aver_e_scans(self, enrg, chsd, interactive):
        '''
        Perform the average of data scans. 
        If interactive mode, data scans and their average are shown together
        in a plot. 

        Parameters
        ----------
        enrg : array
            Energy values at which average is calculated.

        chsd : list (str)
            Scan-numbers of scan to be averaged.

        interactive: bool
            If True chsd scans are plotted together with their average.

        Returns
        -------
        array, containing the average of data scans

        Notes
        -----
        To compute the average the common energy scale enrg is used.
        All passed scans are interpolated with a linear spline (k=1 and s=0 in
        itp.UnivariateSpline) and evaluated along the common energy scale.
        The interpolated data are eventually averaged.
        '''
        intrp = []

        if interactive:
            plt.figure(1)

        for scn in chsd:
            # Chosed data
            x = self.raw_imp['E' + scn][1:]
            y = self.raw_imp[scn][1:]

            if interactive:
                # Plot data
                plt.plot(x, y, color='black')

            # Compute linear spline interpolation
            y_int = itp.UnivariateSpline(x, y, k=1, s=0)
            # Evaluate interpolation of scan data on enrg energy scale and
            # append to previous interpolations
            intrp.append(y_int(enrg))

        # Average all inteprolated scans
        avgd = np.average(intrp, axis=0)

        if interactive:
            plt.plot(enrg, avgd, color='r', label='Average')
            plt.xlabel('E (eV)')
            plt.ylabel(self.dtype)
            plt.legend()
            plt.show()

        return avgd

    #def edge_norm(self, enrg, e_edge, e_pe, pe_rng, pe_int):
    def edge_norm(self, enrg, e_edge, e_pe, e_poste, pe_rng):
        '''
        Normalize energy scan data by the value at pre-edge energy.
        Also compute the  energy jump defined as the difference between
        the value at the edge and pre-edge energies respectively.
        This computations are implemented employing both a pre-edge value
        obtained from data (averaging in the region e_pe +/+ pe_rng) and an
        interpolated pre-edge value.

        Parameters
        ----------
        enrg : array
            Energy values of scan

        e_edge : float
            Edge energy value

        e_pe : float
            Pre-edge energy value

        pe_rng : int
            Number of points constituting the semi-width of energy range
            centered at e_pe.

        pe_int : float
            Pre-edge value obtained from linear interpolation based on pre- and
            post-edge energies.

        Returns
        -------
        Set class attributes:
        pe_av : float
            value at pre-edge energy

        norm : array
            self.aver scan normalized by value at pre-edge energy

        norm_int : array
            Averaged data normalized by interpolated pre-edge value

        ej : float
            edge-jump value

        ej_norm : float
            edge-jump value normalized by value at pre-edge energy

        ej_int : float
            edge-jump value computed with interpolated pre-edge value

        ej_norm_int : float
            edge-jump computed and normalized by interpolated pre-edge value

        Notes
        -----
        To reduce noise effects the value of scan at pre-edge energy is
        obtained computing an average over an energy range of width pe_rng and
        centered at e_pe pre-edge energy.
        The value of scan at edge energy is obtained by linear spline
        interpolation of data (itp.UnivariateSpline with k=1 and s=0).
        '''
        # Index of the nearest element to pre-edge energy
        pe_idx = np.argmin((np.abs(enrg - e_pe)))
        # Left and right extremes of energy range for pre-edge average
        lpe_idx = int(pe_idx - pe_rng)
        rpe_idx = int(pe_idx + pe_rng + 1)

        # Average of values for computation of pre-edge
        self.pe_av = np.average(self.aver[lpe_idx: rpe_idx: 1])

        # Cubic spline interpolation of energy scan
        y_int = itp.UnivariateSpline(enrg, self.aver, k=3, s=0)

        # Interpolation of pre-edge energy
        x = [e_pe, e_poste]
        y = [y_int(e_pe), y_int(e_poste)]
        self.pe_av_int = lin_interpolate(x, y, e_edge)

        # Normalization by pre-edge value
        self.norm = self.aver / self.pe_av
        self.norm_int = self.aver / self.pe_av_int

        
        # value at edge energy from interpolation
        y_edg = y_int(e_edge)

        # Edge-jumps computations
        self.ej = y_edg - self.pe_av
        self.ej_norm = self.ej / self.pe_av

        self.ej_int = y_edg - self.pe_av_int
        self.ej_norm_int = self.ej_int / self.pe_av_int


class EngyScan:
    '''
    Allow computation on energy scan data extracting XNLD, XMCD, XNCD and XNXD
    spectra.

    Attributes
    ----------
    energy : array
        Common energy scale for positive and negative XAS scans data treatement

    exper_edge : float
        Experimental edge energy

    e_pe : float
        Pre-edge energy value

    pe_wdt : float
        Half-width of energy range for pre-edge average computation

    pe_int : float
        Interpolated value of pre-edge considering pre- and post-edge energies
        using a linear approximation

    offest : float
        Offset given by the difference between expected and experimental edge 
        energy

    energycal : array
        Common energy scale calibrated considering the offset

    xd : array
        X-Ray dichroism data

    xd_aver : array
        Average of positive and negative XAS. In case of XNLD analysis a
        weighted average is considered

    xd_pc : array
        X-Ray dichroism normalized by the average of positve and negative scans
        edge jumps respectively. Obtained from positive and negative scans
        normalized by their values at pre-edge energy. The value is returned in
        percentage

    pos_corr : array
        Normalized positive XAS spectrum weighted by the average of values at
        pre-edge energy of positive and negative spectra. For XNLD spectra the
        weight consider also the X-Ray beam incidence angle.

    neg_corr : array
        Normalized negative XAS spectrum weighted by the average of values at
        pre-edge energy of positive and negative spectra. For XNLD spectra the
        weight consider also the X-Ray beam incidence angle.

    xd_pc_av_ej : array
        X-Ray dichroism normalized by edge jump of xd_aver scan.
        Obtained from positive and negative scans normalized by the weighted
        average of values at pre-edge energy. The value is returned in
        percentage.

    Methods
    -------
    e_scale(pos, neg, guiobj, log_dt)
        Creates the energy scale used for data analysis.

    edges(guiobj, log_dt)
        Set values of edge and pre-edge energies.

    compt_xd(pos, neg, guiobj, log_dt)

    '''

    def __init__(self, pos, neg, guiobj, log_dt):
        '''
        At instantiation e_scale and compt_xd method are called, so upon its
        creation an EngyScan object collects all the informations about X-Ray
        dichroism.

        Parameters
        ----------
        pos_scan : ScanData obj
            Positive scans (CR for XMCD and XNCD, LH for XNLD)

        neg_scan : ScanData obj
            Negative scnas (CL for XMCD and XNCD, LV for XNLD)

        guiobj : GUI obj
            Provides GUI dialogs.

        log_dt : dict
            Collect data for logfile
        '''

        self.e_scale(pos, neg, guiobj, log_dt)
        self.compt_xd(pos, neg, guiobj, log_dt)

    def e_scale(self, pos, neg, guiobj, log_dt):
        '''
        Create the energy scale used for data analysis.
        Range is selected considering the intersection of the energy ranges of
        all the provided scans: the low-end is the highest minimum value and 
        the high-end is the lowest maximum value of all the energy arrays.

        The number of points of the returned energy scale by default is the 
        average of the number of points of the energy scales of data scans.
        If in interactive mode GUI dialogues are provided to set the 
        number of points.

        Parameters
        ----------
        pos : ScanData object
            Positive scans data (CR for XMCD and XNCD, LH for XNLD).

        neg : ScanData object
            Negative scans data (CL for XMCD and XNCD, LV for XNLD).

        guiobj : GUI object
            Provides GUI dialogs.

        log_dt : dict
            Collect data for logfile

        Returns
        -------
        Set class attribute:

        energy : array
            Common energy scale for positive and negative XAS scans data
            treatement

        Add e_scale key to log_dt with min, max and number of points of energy
        scale. 
        '''
        e_lists = []
        e_len = []

        # Sort energy data arrays and count their elements
        for i in pos.idx:
            e_lists.append(np.sort(pos.raw_imp['E' + i]))
            e_len.append(pos.raw_imp['E' + i].size)
        for i in neg.idx:
            e_lists.append(np.sort(neg.raw_imp['E' + i]))
            e_len.append(neg.raw_imp['E' + i].size)

        # Compute min, max and default langth of energy range
        e_min = np.around(np.amax(e_lists, axis=0)[0], 1)
        e_max = np.around(np.amin(e_lists, axis=0)[-1], 1)
        e_av_len = np.around(np.average(e_len), 0)

        # Set number of points of energy scale
        if guiobj.interactive:
            n_points = guiobj.e_num_pnts(e_av_len)
        else:
            n_points = int(e_av_len)

        self.energy = np.linspace(e_min, e_max, n_points)

        log_dt['e_scale'] = [e_min, e_max, n_points]

    def edges(self, guiobj, log_dt):
        '''
        Set values of edge energy, pre-edge energy, pre-edge energy range (used
        to compute average of spectra at pre-edge energy) and post-edge energy.

        If in interactive mode a GUI is provided to set the experimental edge
        energy, pre-edge energy, post-edge energy and and half-width of
        interval adopted for pre-edge average.
        Otherwise for edge energy is considered the experimental computed one
        from minimization, for pre-edge end post-edge energies the ones provided
        from edge file are taken for good, and for the half-width of interval
        for pre-edge average 4 points is considered.

        Parameters
        ----------
        guiobj : GUI object
            Provides GUI dialogs.

        log_dt : dict
            Collect data for logfile

        Return
        ------
        Set class attributes

        exper_edge : float
            Experimental edge energy

        e_pe : float
            Pre-edge energy value

        e_poste : float
            Post-edge energy value

        pe_wdt : float
            Half-width of energy range for pre-edge average computation

        pe_int : float
            Interpolated value of pre-edge considering pre- and post-edge
            energies using a linear approximation

        offset : float
            Offset given by the difference between expected and experimental
            edge energy

        energycal : array
            Common energy scale calibrated considering the offset. In this way
            absorption edge of spectra falls at the expected energy.

        Add keys to log_dt

        exper_edge : experimental edge energy
        setted_pedg : setted pre-edge energy
        setted_postedg : setted post-edge energy
        pe_int : interpolated pre-edge
        recal : bool, if True energy recalibration has been done
        offset : energy offset adopted for energy recalibration

        Notes
        -----
        Indication of experimental value of edge energy is provided finding
        the stationary points of experimental spectrum. This is obtained
        interpolating experimental data. Method used for interpolation is 
        scipy.interpolate.UnivariateSpline of order 3.
        '''
        # Retrive from log_dt edge, pre-edge energy and range for pre-edge
        # average

        # Check that pre-edge energy is included in the considered energy
        # range. If not, an energy value close to the nearest range endpoint
        # is considered
        if float(log_dt['PreEdge_en']) <= self.energy[5]:
            pe_e = self.energy[5]
        elif float(log_dt['PreEdge_en']) >= self.energy[-5]:
            pe_e = self.energy[-5]
        else:
            pe_e = float(log_dt['PreEdge_en'])

        # Check that post edge energy is included in the considered energy
        # range
        if float(log_dt['PostEdge_en']) <= self.energy[0]:
            pste_e = self.energy[1]
        elif float(log_dt['PostEdge_en']) >= self.energy[-1]:
            pste_e = self.energy[-2]
        else:
            pste_e = float(log_dt['PostEdge_en'])

        sel_edg = [log_dt['Edge_en'], pe_e, pste_e]

        # Order 3 interpolation of data
        # opt.minimize_scalar serch for function minimum so negative absolute
        # value of interpolated data is considered.
        y_int_for_edge = itp.UnivariateSpline(self.energy, -abs(self.xd), k=3,
            s=0)
        # Bounds for  minimum search - 5 eV window is considered
        u_bnd = log_dt['Edge_en'] + 2.5
        l_bnd = log_dt['Edge_en'] - 2.5
        min_y = opt.minimize_scalar(y_int_for_edge, bounds=(l_bnd, u_bnd),
                                    method='bounded')
        y_int = itp.UnivariateSpline(self.energy, self.xd_aver, k=3, s=0)

        if guiobj.interactive:
            # x value which minimize y is passed as experimental edge energy
            edgs = guiobj.set_edges(sel_edg, min_y.x, self.energy, self.xd,
                                    self.xd_aver, y_int)
           
            self.exper_edge = edgs[0]
            self.e_pe = edgs[1]
            self.e_poste = edgs[2]
            self.pe_wdt = int(edgs[3])
            recal = edgs[4]
        else:
            # If not interactive for exper_edge the result of minimization is
            # considered
            self.exper_edge = float(min_y.x)
            self.e_pe = pe_e
            self.e_poste = pste_e
            self.pe_wdt = 4
            recal = False

        # New energy scale is created adding an offset given by the difference
        # between expected and experimental edge energies if recal is True,
        # otherwise no offset is added
        if recal:
            self.offset = self.exper_edge - log_dt['Edge_en']
        else:
            self.offset = 0

        self.energycal = self.energy - self.offset

        # Log data related to edge and pre-edge energy setted
        log_dt['exper_edge'] = self.exper_edge
        log_dt['setted_pedg'] = self.e_pe
        log_dt['setted_postedg'] = self.e_poste
        log_dt['recal'] = recal
        log_dt['offset'] = self.offset

    def compt_xd(self, pos, neg, guiobj, log_dt):
        '''
        Performs computation on energy scan data extracting XNLD, XMCD, XNXD or
        XNCD spctra.
        Given two ScanData object (one for positive and one for negative
        scans):
        - compute X-Ray Dichroism as negative - positive;
        - compute the arithmetical mean of positive and negative spectra for
          XMCD, XNCD, XNXD data and the weighted average for the angle for XNLD
          data (see Notes);
        - compute percentage of X-Ray Dichroism normalized by the average
          edge-jump for XMCD, XNCD, XNXD data while for XNLD an angle-weighted
          average is considered (see Notes);
        - compute percentage of X-Ray Dichroism normalized by edge-jump of the
          weighted average of positive and negative spectra.

        Parameters
        ----------
        pos : ScanData obj
            Positive scans data (CR for XMCD and XNCD, LH for XNLD).

        neg : ScanData obj
            Negative scans data (CL for XMCD and XNCD, LV for XNLD).

        guiobj : GUI object
            Provides GUI dialogs.

        log_dt : dict
            Collect data for logfile

        Returns
        -------
        Set class attributes:

        xd : array
            X-Ray Dichroism obtianed as negative - positive scans

        xd_aver : array
            average of XAS spectra. Arithmetical mean of positive and negative
            scans for XMCD and XNCD analysis. Weigthed average by angle for
            XNLD analysis (see Notes)

        xd_pc :array
            Percentage of X-Ray Dichroism normalized for edge-jump.
            For XMCD and XNCD it is obtained as the difference between
            normalized by pre-edge negative and positve scans normalized by the
            arithmetical mean of the two respective edge jumps. For XNLD the
            average weighted by angle of incidence X-Ray beam considered
            (see Notes)

        pos_corr : array
            Average of positive spectra weighted by the ratio between the mean
            of positive pre-edge values and the mean of all pre-edge values.
            For XMCD and XNCD the mean of positive and negative pre-edges is
            arithmetical while for XNLD the average weighted by angle of
            incidence X-Ray beam considered (see Notes)

        neg_corr : array
            Average of negative spectra weighted by the ratio between the mean
            of positive pre-edge values and the mean of all pre-edge values.
            For XMCD and XNCD the mean of positive and negative pre-edges is
            arithmetical while for XNLD the average weighted by angle of
            incidence X-Ray beam considered (see Notes)

        xd_pc_av_ej : array
            Percentage of X-Ray Dichroism normalized for edge-jump of xd_aver
            spectrum. It is computed employing the weighted averages pos_corr
            and neg_corr data. Edge-jump is calculated interpolating the values
            at edge and pre-edge energy of xd_aver spectrum with linear spline
            method.

        Add keys to log_dt

        pos_chs : list (str) with positive chosed scan
        neg_chs : list (str) with negative chosed scan

        Notes
        -----
        In XNLD analsysis LH and LV spectra are weighted by the angle t between
        the sample surface and the X-Ray incident beam direction.

        - XNLD average = (LH + (2 * cos(t)^2 - sin(t)^2) * LV) / 3 * cos(t)^2

        - XNLD (%) = 100 * 3 * cos(t)^2 * (LV/LV_pe - LH/LH_pe) /
                     (LH/LH_pe + (2 * cos(t)^2 - sin(t)^2) * LV/LV_pe)

        - XNLD norm by XNLD_aver edge-jump (%) = 
                         100 * ((LH_pe + (2 * cos(t)^2 - sin(t)^2) * LV_pe) / 
                         3 * cos(t)^2) * ((LV / LV_pe) - (LH / LH_pe)) /
                         (XNLD_aver_edge - XNLD_aver_pe))
        '''
        # Computes averages of positive and negative polarization scans.
        pos.man_aver_e_scans(guiobj, self.energy)
        neg.man_aver_e_scans(guiobj, self.energy)

        # Add keys with chosed scans to log_dt
        log_dt['pos_chs'] = pos.chsd_scns
        log_dt['neg_chs'] = neg.chsd_scns

        # Computes not normalized X-Ray Dichroism.
        self.xd = neg.aver - pos.aver

        if guiobj.case in guiobj.type['xnld']:
            # If XNLD the angle must be considered for weighted mean computation
            theta = log_dt['bm_angle']  # retrive angle from log table

            # Numerator and denominator terms of angle weight
            ang_w_n = 2 * (np.cos(theta))**2 - (np.sin(theta))**2
            ang_w_d = 3 * (np.cos(theta))**2

            self.xd_aver = (pos.aver + (ang_w_n * neg.aver)) / ang_w_d
        else:
            self.xd_aver = (pos.aver + neg.aver) / 2

        self.edges(guiobj, log_dt)

        # Normalize spectra
        pos.edge_norm(self.energycal, self.exper_edge, self.e_pe, self.e_poste,
                      self.pe_wdt)
        neg.edge_norm(self.energycal, self.exper_edge, self.e_pe, self.e_poste,
                      self.pe_wdt)


        log_dt['pos_ej'] = pos.ej
        log_dt['pos_ej_int'] = pos.ej_int
        log_dt['neg_ej'] = neg.ej
        log_dt['neg_ej_int'] = neg.ej_int

        # Compute mean
        if guiobj.case in guiobj.type['xnld']:
            # Percentage X-Ray dichroism normalized for edge jump
            self.xd_pc = (100 * ang_w_d * (neg.norm - pos.norm) /
                          (pos.ej_norm + ang_w_n * neg.ej_norm))
            self.xd_pc_int = (100 * ang_w_d * (neg.norm_int - pos.norm_int) /
                              (pos.ej_norm_int + ang_w_n * neg.ej_norm_int))

            # Angle weighted average of pre-edges values
            av_pe_av = (pos.pe_av + ang_w_n * neg.pe_av) / ang_w_d
            av_pe_int = (pos.pe_av_int + ang_w_n * neg.pe_av_int) / ang_w_d
        else:
            # Percentage X-Ray dichroism normalized for edge jump
            self.xd_pc = 200 * ((neg.norm - pos.norm) /
                                (neg.ej_norm + pos.ej_norm))
            self.xd_pc_int = 200 * ((neg.norm_int - pos.norm_int) /
                                    (neg.ej_norm_int + pos.ej_norm_int))

            # Average of pre-edges values
            av_pe_av = (pos.pe_av + neg.pe_av) / 2
            av_pe_int = (pos.pe_av_int + neg.pe_av_int) /2

        # Calculate xd percentage normalizing with edge-jump obtained from
        # weighted average of positive and negative spectra self.xd_aver
        self.pos_corr = pos.aver * av_pe_av / pos.pe_av
        self.neg_corr = neg.aver * av_pe_av / neg.pe_av

        self.pos_corr_int = pos.aver * av_pe_int / pos.pe_av_int
        self.neg_corr_int = neg.aver * av_pe_int / neg.pe_av_int

        # Linear spline interpolation of xd_aver spectrum in order to determine
        # edge jump
        xd_aver_inter = itp.UnivariateSpline(self.energycal, self.xd_aver, k=1,
                                            s=0)
        # Ineterpolated values at edge and pre-edge energies
        edg_val_xd_aver = xd_aver_inter(self.exper_edge)
        pedg_val_xd_aver = xd_aver_inter(self.e_pe)

        # Compute linear interpolation of pre-edge for average xd
        pstedg_val_xd_aver = xd_aver_inter(self.e_poste)

        x = [self.e_pe, self.e_poste]
        y = [pedg_val_xd_aver, pstedg_val_xd_aver]
        pedg_val_xd_aver_int = lin_interpolate(x, y, self.exper_edge)

        # xd percentage normalized for edge-jump of xd_aver
        self.xd_pc_av_ej = (100 * (self.neg_corr - self.pos_corr) /
                            (edg_val_xd_aver - pedg_val_xd_aver))
        self.xd_pc_av_ej_int = (100 * (self.neg_corr_int - self.pos_corr_int) /
                                (edg_val_xd_aver - pedg_val_xd_aver_int))

        log_dt['xas_aver_ej'] = edg_val_xd_aver - pedg_val_xd_aver
        log_dt['xas_aver_ej_int'] = edg_val_xd_aver - pedg_val_xd_aver_int

def lin_interpolate(x, y, x0):
        '''
        Linear interpolation of the value at x0
        A line is determined considering the couple of points passed in x and y.
        The value at x0 point is computed.

        Parameters
        ----------
        x : array
            x values of point 1 and point 2

        y : array
            y values of point 1 and point2

        x0 : float
            point where linear interpolation is computed

        Return
        ------
        float, interpolated value
        '''
        return y[0] + ((y[1] - y[0]) * (x0 - x[0]) / (x[1] - x[0]))