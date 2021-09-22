"""
pyDichroX_hscan_datatreat.py

Classes and methods for hysteresys data analysis.

Classes
--------
ScanData : Collect raw data for energy scan experiments and provides methods to
    average them.

EngyScan : Allow computation on energy scan data extracting XNLD, XMCD or XNCD
    spctra.

Methods
-------
e_scale(guiobj, pos, neg, log_dt, pos_ref, neg_ref)
    Create the energy scale used for data analysis.

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
    Collect raw data for magnetic field scan experiments and provides methods to
    average them.

    Attributes
    ----------
    label : list (str)
        Labels for graphs

    idx : list (str)
        Scan indexes

    raw_imp : pandas DataFrame
        Collect imported raw data

    up : pandas DataFrame

    down : pandas DataFrame
    
    pe_label : list (str)
        Labels for graphs

    pe_idx : list (str)
        Scan indexes

    pe_raw_imp : pandas DataFrame
        Collect imported raw data

    pe_up : pandas DataFrame

    pe_down : pandas DataFrame

    aver : array
        Data average of selected scans from raw_imp.
        If ScanData is a reference one, aver contain normalized data by
        reference.

    dtype : str
        Identifies the data collected, used for graph labelling:
        sigma+, sigma- for XMCD
        CR, CL for XNCD
        H-, H+ for XNXD
        LH, LV fot XNLD

    chsd_scns : list (str)
        Labels of chosed scans for the analysis

    norm : array
        Averaged data normalized by value at pre-edge energy

    Methods
    -------
    man_aver_e_scans(guiobj, enrg)
        Manage the choice of scans to be averaged and return the average of
        selected scans.

    aver_e_scans(enrg, chsd, guiobj)
        Performe the average of data scans.

    edge_norm(enrg, e_edge, e_pe, pe_rng, pe_int)
        Normalize energy scan data by value at pre-edge energy and compute
        edge-jump.
    '''

    def __init__(self):
        '''
        Initialize attributes label, idx, and raw_imp, 
        '''
        self.label = []
        self.idx = []
        self.raw_imp = pd.DataFrame()
        self.up = pd.DataFrame()
        self.down = pd.DataFrame()

        self.pe_label = []
        self.pe_idx = []
        self.pe_raw_imp = pd.DataFrame()
        self.pe_up = pd.DataFrame()
        self.pe_down = pd.DataFrame()

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
                if guiobj.infile_ref:
                    plt.title('Reference sample scans')
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
                avgd = self.aver_e_scans(enrg, chsd, guiobj)

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

            avgd = self.aver_e_scans(enrg, chsd, guiobj)

        self.aver = avgd
        self.chsd_scns = chsd  # Collect choosed scans for log purpose

    def aver_e_scans(self, enrg, chsd, guiobj):
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

        guiobj: GUI object
            Provides GUI dialogs.

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

        if guiobj.interactive:
            plt.figure(1)
            if guiobj.infile_ref:
                plt.title('Reference sample scans')

        for scn in chsd:
            # Chosed data
            x = self.raw_imp['E' + scn][1:]
            y = self.raw_imp[scn][1:]

            if guiobj.interactive:
                # Plot data
                plt.plot(x, y, color='black')

            # Compute linear spline interpolation
            y_int = itp.UnivariateSpline(x, y, k=1, s=0)
            # Evaluate interpolation of scan data on enrg energy scale and
            # append to previous interpolations
            intrp.append(y_int(enrg))

        # Average all inteprolated scans
        avgd = np.average(intrp, axis=0)

        if guiobj.interactive:
            plt.plot(enrg, avgd, color='r', label='Average')
            plt.xlabel('E (eV)')
            plt.ylabel(self.dtype)
            plt.legend()
            plt.show()

        return avgd

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

    e_poste : float
        Post-edge energy value

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

    ang_w_n : float
        Angle weight for XNLD computation

    ang_w_d : float
        Angle weight for XNLD computation

    xd : array
        X-Ray dichroism data

    xd_aver : array
        Average of positive and negative XAS. In case of XNLD analysis a
        weighted average is considered

    xd_pc : array
        Percentage of X-Ray dichroism normalized by the average of positve and
        negative scans edge jumps respectively. 

    xd_pc_int :array
        Same as xd_pc using data normalized with inerpolated edge-jump.

    pos_corr : array
        Normalized positive XAS spectrum weighted by the average of values at
        pre-edge energy of positive and negative spectra. For XNLD spectra the
        weight consider also the X-Ray beam incidence angle.

    neg_corr : array
        Normalized negative XAS spectrum weighted by the average of values at
        pre-edge energy of positive and negative spectra. For XNLD spectra the
        weight consider also the X-Ray beam incidence angle.

    pos_corr_int : array
        Same as pos_corr but using for the avereage of pre-edge values the
        ineterpolated values. For XNLD spectra the weight consider also the
        X-Ray beam incidence angle.

    neg_corr_int : array
        Same as neg_corr but using for the avereage of pre-edge values the
        ineterpolated values. For XNLD spectra the weight consider also the
        X-Ray beam incidence angle.

    xd_pc_av_ej : array
        Percentage  ofX-Ray dichroism normalized by edge jump of xd_aver scan.        

    xd_pc_av_ej_int : array
        Same as xd_pc_av_ej but using for the edge-jump computation the 
        interpolated values of pre-edges.

    Methods
    -------
    scan_average(guiobj, pos, neg, log_dt, pos_to_norm, neg_to_norm)
        Copmute averages of positive and negative scans.

    compt_mxd(guiobj, pos, neg, log_dt)
        Performs computation on energy scan data extracting XNLD, XMCD, XNXD or
        XNCD spctra.

    compt_xd(self, guiobj, pos, neg, log_dt)
        Compute not normalized X-Ray Dichroism.

    edges(guiobj, confobj, log_dt)
        Set values of edge and pre-edge energies.

    compt_xd_pc(guiobj, pos, neg)
        Compute the percentage of X-Ray Dichroism normalized for edge-jump.

    compt_pe_corr(guiobj, pos, neg)
        Calculate xd percentage normalizing with edge-jump obtained from
        weighted average of positive and negative spectra.

    comp_xd_pc_av_ej(self, log_dt)
        Compute percentage of X-Ray Dichroism normalized for edge-jump of
        xd_aver spectrum.
    '''

    def __init__(self, guiobj, confobj, e_scale, pos, neg, log_dt,
        pos_to_norm=ScanData(), neg_to_norm=ScanData()):
        '''
        At instantiation scan_average and compt_xd method are called and energy
        attribute is settee, so upon its creation an EngyScan object collects
        all the informations about X-Ray dichroism.

        Parameters
        ----------
        guiobj : GUI obj
            Provides GUI dialogs.

        confobj : Configuration object

        e_scale : array
            Common energy scale for data analysis

        pos : ScanData obj
            Positive scans (CR for XMCD and XNCD, LH for XNLD)

        neg : ScanData obj
            Negative scnas (CL for XMCD and XNCD, LV for XNLD)        

        log_dt : dict
            Collect data for logfile

        pos_to_norm : ScanData obj
            Positive scans (CR for XMCD and XNCD, LH for XNLD) data to be
            normalized by reference if present. If no reference data are present
            pass empty ScanData object.

        neg_to_norm : ScanData obj
            Negative scans (CR for XMCD and XNCD, LH for XNLD) data to be
            normalized by reference if present. If no reference data are present
            pass empty ScanData object.
        '''
        self.energy = e_scale

        self.scan_average(guiobj, pos, neg, log_dt, pos_to_norm, neg_to_norm)
        
        self.compt_mxd(guiobj, pos, neg, log_dt)

    def scan_average(self, guiobj, pos, neg, log_dt, pos_to_norm, neg_to_norm):
        '''
        Copmute averages of positive and negative scans.
        If reference data are considered compute also the averages of reference
        scans.

        Parameters
        ----------
        guiobj : GUI obj
            Provides GUI dialogs.

        pos_scan : ScanData obj
            Positive scans (CR for XMCD and XNCD, LH for XNLD)

        neg_scan : ScanData obj
            Negative scnas (CL for XMCD and XNCD, LV for XNLD)        

        log_dt : dict
            Collect data for logfile

        pos_to_norm : ScanData obj
            Positive scans (CR for XMCD and XNCD, LH for XNLD) data to be
            normalized by reference if present

        neg_to_norm : ScanData obj
            Negative scans (CR for XMCD and XNCD, LH for XNLD) data to be
            normalized by reference if present

        Return
        ------
        Instantiate aver attribute of positive and negative ScanData objects.
        In case referecnce data are passed aver attribute is the normalization
        ***_to_norm by reference data.

        Add keys to log_dt

        pos_chs : list (str) with positive chosed scan
        neg_chs : list (str) with negative chosed scan

        '''
        # If no reference data pos_to_norm and neg_to_norm are empty ScanData
        # object
        if pos_to_norm.raw_imp.empty or neg_to_norm.raw_imp.empty:
            # Computes averages of positive and negative polarization scans.
            guiobj.infile_ref = False

            pos.man_aver_e_scans(guiobj, self.energy)
            neg.man_aver_e_scans(guiobj, self.energy)            
        else:
            # If present computes averages of positive and negative polarization
            # of reference scans and normalize data for them.
            # In this case pos and neg contains reference data and pos_to_norm
            # neg_to_norm are data to be normalized. They are supposed to be
            # already analyzed so aver attribute is already setted
            guiobj.infile_ref = True
            pos.man_aver_e_scans(guiobj, self.energy)
            neg.man_aver_e_scans(guiobj, self.energy)

            # Just substitute pos.aver and neg.aver with normalized data
            pos.aver = pos_to_norm.aver / pos.aver
            neg.aver = neg_to_norm.aver / neg.aver

            guiobj.infile_ref = False

        # Add keys with chosed scans to log_dt
        log_dt['pos_chs'] = pos.chsd_scns
        log_dt['neg_chs'] = neg.chsd_scns

    def compt_mxd(self, guiobj, pos, neg, log_dt):
        '''
        Performs computation on energy scan data extracting XNLD, XMCD, XNXD or
        XNCD spctra.
        Given two ScanData object (one for positive and one for negative scans):
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
        guiobj : GUI obj
            Provides GUI dialogs.

        pos : ScanData obj
            Positive scans (CR for XMCD and XNCD, LH for XNLD)

        neg : ScanData obj
            Negative scnas (CL for XMCD and XNCD, LV for XNLD)        

        log_dt : dict
            Collect data for logfile

        Returns
        -------
        Add keys to log_dt

        pos_ej : float, edge-jump value for postitive scans
        pos_ej_int : float, edge-jump value for postitive scans, interpolated
            value of pre-edge energy is used
        neg_ej : float, edge-jump value for postitive scans
        neg_ej_int : float, edge-jump value for postitive scans, interpolated
            value of pre-edge energy is used
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
        # Computes not normalized X-Ray Dichroism.
        self.compt_xd(guiobj, pos, neg, log_dt)

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

        # Compute percentage X-Ray Dichroism normalized for edge-jump
        self.compt_xd_pc(guiobj, pos, neg)
        self.compt_pe_corr(guiobj, pos, neg)
        self.comp_xd_pc_av_ej(log_dt)

    def compt_xd(self, guiobj, pos, neg, log_dt):
        '''
        Compute not normalized X-Ray Dichroism.

        Parameters
        ----------
        guiobj : GUI obj
            Provides GUI dialogs

        pos_scan : ScanData obj
            Positive scans (CR for XMCD and XNCD, LH for XNLD)

        neg_scan : ScanData obj
            Negative scnas (CL for XMCD and XNCD, LV for XNLD)        

        log_dt : dict
            Collect data for logfile
        
        Return
        ------
        Set attributes:
        
        xd : array
            X-Ray Dichroism as negative - positive scans

        xd_aver : array
            Average of XAS spectra. Arithmetical mean of positive and
            negative scans for XMCD and XNCD analysis. Weigthed average by angle
            for XNLD analysis (see Notes)
        
        ang_w_n : float
            Angle weight for XNLD computation

        ang_w_d : float
            Angle weight for XNLD computation

        Notes
        -----
        In XNLD analsysis LH and LV spectra are weighted by the angle t between
        the sample surface and the X-Ray incident beam direction.

        XNLD average = (LH + (2 * cos(t)^2 - sin(t)^2) * LV) / 3 * cos(t)^2
        '''
        self.xd = neg.aver - pos.aver
        if guiobj.analysis in guiobj.type['xnld']:
            # If XNLD the angle must be considered for weighted mean computation
            theta = log_dt['bm_angle']  # retrive angle from log table

            # Numerator and denominator terms of angle weight
            self.ang_w_n = 2 * (np.cos(theta))**2 - (np.sin(theta))**2
            self.ang_w_d = 3 * (np.cos(theta))**2

            self.xd_aver = (pos.aver + (self.ang_w_n * neg.aver)) / self.ang_w_d
        else:
            self.xd_aver = (pos.aver + neg.aver) / 2

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

        y = self.xd
        y_aver = self.xd_aver

        # Order 3 interpolation of data
        # opt.minimize_scalar serch for function minimum so negative absolute
        # value of interpolated data is considered.
        y_int_for_edge = itp.UnivariateSpline(self.energy, -abs(y), k=3, s=0)
        # Bounds for  minimum search - 5 eV window is considered
        u_bnd = log_dt['Edge_en'] + 2.5
        l_bnd = log_dt['Edge_en'] - 2.5
        min_y = opt.minimize_scalar(y_int_for_edge, bounds=(l_bnd, u_bnd),
                                    method='bounded')
        y_int = itp.UnivariateSpline(self.energy, y_aver, k=3, s=0)

        if guiobj.interactive:
            # x value which minimize y is passed as experimental edge energy
            edgs = guiobj.set_edges(sel_edg, min_y.x, self.energy, y, y_aver,
                y_int)
           
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

    def compt_xd_pc(self, guiobj, pos, neg):
        '''
        Compute the percentage of X-Ray Dichroism normalized for edge-jump.
        For XMCD and XNCD it is obtained as the difference between
        normalized by pre-edge negative and positve scans normalized by the
        arithmetical mean of the two respective edge jumps. For XNLD the
        average weighted by angle of incidence X-Ray beam considered
        (see Notes)

        Parameters
        ----------
        guiobj : GUI obj
            Provides GUI dialogs.

        pos : ScanData obj
            Positive scans (CR for XMCD and XNCD, LH for XNLD)

        neg : ScanData obj
            Negative scnas (CL for XMCD and XNCD, LV for XNLD)        
        
        Return
        ------
        Set attributes
        xd_pc : array
            Percentage of X-Ray Dichroism normalized for edge-jump.
            For XMCD and XNCD it is obtained as the difference between
            normalized by pre-edge negative and positve scans normalized by the
            arithmetical mean of the two respective edge jumps. For XNLD the
            average weighted by angle of incidence X-Ray beam considered
            (see Notes)

        xd_pc_int : array
            Same as xd_pc using data normalized with inerpolated edge-jump.

        Notes
        -----
        In XNLD analsysis LH and LV spectra are weighted by the angle t between
        the sample surface and the X-Ray incident beam direction.

         XNLD (%) = 100 * 3 * cos(t)^2 * (LV/LV_pe - LH/LH_pe) /
                     (LH/LH_pe + (2 * cos(t)^2 - sin(t)^2) * LV/LV_pe)
        '''
        # Compute mean
        if guiobj.analysis in guiobj.type['xnld']:
            # Percentage X-Ray dichroism normalized for edge jump
            self.xd_pc = (100 * self.ang_w_d * (neg.norm - pos.norm) /
                (pos.ej_norm + self.ang_w_n * neg.ej_norm))
            self.xd_pc_int = (100 * self.ang_w_d * (neg.norm_int - 
                pos.norm_int) / (pos.ej_norm_int + self.ang_w_n * 
                neg.ej_norm_int))
        else:
            # Percentage X-Ray dichroism normalized for edge jump
            self.xd_pc = 200 * (neg.norm - pos.norm) / (neg.ej_norm + 
                pos.ej_norm)
            self.xd_pc_int = 200 * ((neg.norm_int - pos.norm_int) /
                (neg.ej_norm_int + pos.ej_norm_int))

    def compt_pe_corr(self, guiobj, pos, neg):
        '''
        Calculate xd percentage normalizing with edge-jump obtained from
        weighted average of positive and negative spectra

        Parameters
        ----------
        guiobj : GUI obj
            Provides GUI dialogs.

        pos : ScanData obj
            Positive scans (CR for XMCD and XNCD, LH for XNLD)

        neg : ScanData obj
            Negative scnas (CL for XMCD and XNCD, LV for XNLD)        
        
        Return
        ------
        Set attributes:

        pos_corr : array
            Average of positive spectra weighted by the ratio between the mean
            of positive pre-edge values and the mean of all pre-edge values.
            For XMCD and XNCD the mean of positive and negative pre-edges is
            arithmetical while for XNLD the average weighted by angle of
            incidence X-Ray beam considered (see Notes)

        pos_corr_int : array
            Same as pos_corr using interpolated pre-edge values

        neg_corr : array
            Average of negative spectra weighted by the ratio between the mean
            of positive pre-edge values and the mean of all pre-edge values.
            For XMCD and XNCD the mean of positive and negative pre-edges is
            arithmetical while for XNLD the average weighted by angle of
            incidence X-Ray beam considered (see Notes)

        neg_corr_int : array
            Same as neg_corr using interpolated pre-edge values
        
        pos_corr_int : array
            Positive scan normalized by weighted average of interpolated
            edge-jump
        
        neg_corr_int : array
            Negative scan normalized by weighted average of interpolated
            edge-jump

        Notes
        -----
        In XNLD analsysis LH and LV spectra are weighted by the angle t between
        the sample surface and the X-Ray incident beam direction.

        XNLD average = (LH + (2 * cos(t)^2 - sin(t)^2) * LV) / 3 * cos(t)^2
        '''
        # Compute mean
        if guiobj.analysis in guiobj.type['xnld']:
            # Angle weighted average of pre-edges values
            av_pe_av = (pos.pe_av + self.ang_w_n * neg.pe_av) / self.ang_w_d
            av_pe_int = ((pos.pe_av_int + self.ang_w_n * neg.pe_av_int) / 
                self.ang_w_d)
        else:
            # Average of pre-edges values
            av_pe_av = (pos.pe_av + neg.pe_av) / 2
            av_pe_int = (pos.pe_av_int + neg.pe_av_int) /2

        self.pos_corr = pos.aver * av_pe_av / pos.pe_av
        self.neg_corr = neg.aver * av_pe_av / neg.pe_av

        self.pos_corr_int = pos.aver * av_pe_int / pos.pe_av_int
        self.neg_corr_int = neg.aver * av_pe_int / neg.pe_av_int

    def comp_xd_pc_av_ej(self, log_dt):
        '''
        Compute percentage of X-Ray Dichroism normalized for edge-jump of
        xd_aver spectrum. It is computed employing the weighted averages
        pos_corr and neg_corr data. Edge-jump is calculated interpolating the
        values at edge and pre-edge energy of xd_aver spectrum with linear
        spline method.
        For XNLD the average weighted by angle of incidence X-Ray beam is
        considered (see Notes).

        Parameters
        ----------
        log_dt : dict
            Collect data for logfile

        Returns
        -------
        Set class attributes:

        xd_pc_av_ej : array
            Percentage of X-Ray Dichroism normalized for edge-jump of xd_aver
            spectrum.

        xd_pc_av_ej_int : array
            Same as xd_pc_av_ej but using for the edge-jump computation the 
            interpolated values of pre-edges.

        Add keys to log_dt
        xas_aver_ej : float, edge_jump value of xd average
        xas_aver_ej_int : float, edge jump value of xd average considering
            interpolated value of pre-edge

        Notes
        -----
        In XNLD analsysis LH and LV spectra are weighted by the angle t between
        the sample surface and the X-Ray incident beam direction.

        XNLD norm by XNLD_aver edge-jump (%) = 
                         100 * ((LH_pe + (2 * cos(t)^2 - sin(t)^2) * LV_pe) / 
                         3 * cos(t)^2) * ((LV / LV_pe) - (LH / LH_pe)) /
                         (XNLD_aver_edge - XNLD_aver_pe))
        '''
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


def e_scale(guiobj, pos, neg, log_dt, pos_ref, neg_ref):
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
    guiobj : GUI obj
        Provides GUI dialogs.

    pos_scan : ScanData obj
        Positive scans (CR for XMCD and XNCD, LH for XNLD)

    neg_scan : ScanData obj
        Negative scnas (CL for XMCD and XNCD, LV for XNLD)        

    log_dt : dict
        Collect data for logfile

    pos_ref : ScanData obj
        Positive scans (CR for XMCD and XNCD, LH for XNLD) from reference

    neg_ref : ScanData obj
        Negative scans (CR for XMCD and XNCD, LH for XNLD) from reference

    Return
    ------
    array
        Common energy scale for positive and negative XAS scans data
        treatement (if reference scans are considered in confobj also their
        energy scale will be considered)

    Add e_scale key to log_dt with min, max and number of points of energy
    scale. 
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
        n_points = guiobj.e_num_pnts(e_av_len)
    else:
        n_points = int(e_av_len)

    log_dt['e_scale'] = [e_min, e_max, n_points]

    return np.linspace(e_min, e_max, n_points)


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