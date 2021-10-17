"""
pyDichroX_hscan_datatreat.py

Classes and methods for hysteresys data analysis.

Classes
--------
ScanData : Collect raw data for magnetic field scan experiments and
    provides methods to average them.

FieldScan : Allow computation of XMCD hysteresis on magnetic field scan
    data.

Methods
-------
h_scale(guiobj, pos, neg, log_dt)
    Create the magnetic field scale used for data analysis.

h_num_points(h_arr)
    Count the number of different fields present in h_arr.

aver_duplicates(data, idx)
    Check if in the field column idx there are repeated values.
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
    Collect raw data for magnetic field scan experiments and provides
    methods to average them.

    Attributes
    ----------
    label : list (str)
        Labels for graphs - for scans collected at edge energy.

    idx : list (str)
        Scan indexes - for scans collected at edge energy.

    raw_imp : pandas DataFrame
        Collect imported raw data at edge energy.

    pre_edge : bool
        True if pre-edge data are present.

    up : pandas DataFrame
        Collect data related to magntic filed scan up branch.

    down : pandas DataFrame
        Collect data related to magntic filed scan down branch.

    up_chsn : list
        Collect labels of ub branch scans chosen for analysis.

    dw_chsn : list
        Collect labels of down branch scans chosen for analysis.

    up_aver : array
        Contains average of chosen up branch scans computed on common
        field scale.
    
    dw_aver : array
        Contains average of chosen down branch scans computed on common
        field scale.
    
    pe_label : list (str)
        Labels for graphs - for scans collected at pre-edge energy.

    pe_idx : list (str)
        Scan indexes - for scans collected at pre-edge energy.

    pe_raw_imp : pandas DataFrame
        Collect imported raw data at pre-edge energy.

    pe_up : pandas DataFrame.
        Collect data related to magntic filed scan up branch - for scans
        measured at pre-edge energy.

    pe_down : pandas DataFrame.
        Collect data related to magntic filed scan up branch - for scans
        measured at pre-edge energy.

    pe_up_chsn : list
        Collect labels of ub branch scans chosen for analysis.

    pe_dw_chsn : list
        Collect labels of down branch scans chosen for analysis.

    pe_up_aver : array
        Contains average of chosen up branch scans computed on common
        field scale.
    
    pe_dw_aver : array
        Contains average of chosen down branch scans computed on common
        field scale.

    dtype : str
        Identifies CR and CL data, used for graph labelling.

    Methods
    -------
    up_n_down()
        Separate input data in up branches and down branches.


    man_aver_e_scans(guiobj, enrg)
        Manage the choice of scans to be averaged and return the average of
        selected scans.

    plot_chs_avr(fields, guiobj, edge, up)
        Plot raw data, allow to choose which data will be used for
        analysis and average them.

    aver_h_scans(data, fields, chsn, guiobj, title)
        Perform the average of data scans.
    '''

    def __init__(self):
        '''
        Initialize attributes label, idx, raw_imp, up, down and pre-edge
        corresponding pe_label, pe_idx, pe_raw_imp, pe_up, pe_down.
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

    def up_n_down(self):
        '''
        Separate input data in up branches and down branches.

        Return
        ------
        Populate the attributes self.up and self.down with data from
        up and down branches respectively.
        If pre-edge data are present, also populate self.pe_up and 
        self.pe_down attributes.
        '''
        # Run through edge scans
        for i in self.idx:
            # Magnetic field column's name
            h_col = 'H' + i
            # Select data from scan i
            i_scan = self.raw_imp[[h_col, i]]

            if i_scan[h_col].iloc[0] > 0:
                # If the first element is > 0 => scan down
                # Look for the index where the minimum of magnetic
                # fields is to find where the up branch starts
                sep_idx = np.argmin(i_scan[h_col])

                # Seprate down branch from up branch.
                # For dw consider sep_idx + 1 in order to include the
                # element related to minimum field value
                dw = i_scan.iloc[0:sep_idx+1, :]
                up = i_scan.iloc[sep_idx:, :]
                self.down = pd.concat([self.down, dw], axis=1)
                self.up = pd.concat([self.up, up], axis=1)
            else:
                # If the first element is < 0 => scan up
                # Look for the index where the maximum of magnetic
                # fields is to find where the down branch starts
                sep_idx = np.argmax(i_scan[h_col])

                # Seprate up branch from down branch.
                # For up consider sep_idx + 1 in order to include the
                # element related to maximum field value
                up = i_scan.iloc[0:sep_idx+1, :]
                dw = i_scan.iloc[sep_idx:, :]
                self.down = pd.concat([self.down, dw], axis=1)
                self.up = pd.concat([self.up, up], axis=1)

        # Run through pre-edge scans if present
        if self.pre_edge:
            for i in self.pe_idx:
                # Magnetic field column's name
                h_col = 'H' + i
                # Select data from scan i
                i_scan = self.pe_raw_imp[[h_col, i]]

                if i_scan[h_col].iloc[0] > 0:
                    # If the first element is > 0 => scan down
                    # Look for the index where the minimum of magnetic
                    # fields is to find where the up branch starts
                    sep_idx = np.argmin(i_scan[h_col])

                    # Seprate down branch from up branch.
                    # For dw consider sep_idx + 1 in order to include
                    # the element related to minimum field value
                    dw = i_scan.iloc[0:sep_idx+1, :]
                    up = i_scan.iloc[sep_idx:, :]
                    self.pe_down = pd.concat([self.pe_down, dw], axis=1)
                    self.pe_up = pd.concat([self.pe_up, up], axis=1)
                else:
                    # If the first element is < 0 => scan up
                    # Look for the index where the maximum of magnetic
                    # fields is to find where the down branch starts
                    sep_idx = np.argmax(i_scan[h_col])

                    # Seprate up branch from down branch.
                    # For up consider sep_idx + 1 in order to include
                    # the element related to maximum field value
                    up = i_scan.iloc[0:sep_idx+1, :]
                    dw = i_scan.iloc[sep_idx:, :]
                    self.pe_down = pd.concat([self.pe_down, dw], axis=1)
                    self.pe_up = pd.concat([self.pe_up, up], axis=1)

    def man_aver_h_scans(self, guiobj, fields):
        '''
        Manage the choice of scans to be averaged and return the average
        of selected scans.

        Parameters
        ----------
        guiobj : GUI object
            Provides GUI dialogs.

        fields : array
            Magnetic field values at which average is calculated.

        Return
        ------
        Set class attributes:
        aver : array
            Average values of the chosen scans.

        chsn_scns : list
            Labels of chosen scans for the analysis (for log purpose)
        '''
        # Compute average of scans only if there are more than one scan
        # for each branch (each scan contain one branch up and one
        # branch down)
        if guiobj.interactive:  # Interactive choose of scans        
            # Choose edge scans.
            if len(self.idx) > 1:
                # Edge Up data
                # Loop until choice is confirmed
                isok = 0
                while not isok:
                    up_chsn, up_avgd = self.plot_chs_avr(fields, guiobj,
                                                        edge=True, up=True)
                    isok = guiobj.confirm_choice()

                # Edge Down data
                # Loop until choice is confirmed
                isok = 0
                while not isok:
                    dw_chsn, dw_avgd = self.plot_chs_avr(fields, guiobj,
                                                        edge=True, up=False)
                    isok = guiobj.confirm_choice()
            else:
                # If there is just one scan is useless the choose, just
                # interpolate it with the common field scale.
                up_chsn = self.label
                up_avgd = self.aver_h_scans(self.up, fields, up_chsn, guiobj,
                                            title='Edge - Up branch')
                dw_chsn = self.label
                dw_avgd = self.aver_h_scans(self.down, fields, dw_chsn, guiobj,
                                            title='Edge - Down branch')
            # Choose pre-edge scans
            if self.pre_edge:
                if len(self.pe_idx) > 1:
                    # pre-edge Up data
                    # Loop until choice is confirmed
                    isok = 0
                    while not isok:
                        pe_up_chsn, pe_up_avgd = self.plot_chs_avr(fields,
                                                guiobj, edge=False, up=True)
                        isok = guiobj.confirm_choice()

                    # pre-edge Down data
                    # Loop until choice is confirmed
                    isok = 0
                    while not isok:
                        pe_dw_chsn, pe_dw_avgd = self.plot_chs_avr(fields,
                                                guiobj, edge=False, up=False)
                        isok = guiobj.confirm_choice()
                else:  # There is just one pre-edge scan
                    # If there is just one scan is useless the choose, just
                    # interpolate it with the common field scale.
                    pe_up_chsn = self.pe_label
                    pe_up_avgd = self.aver_h_scans(self.pe_up, fields,
                            pe_up_chsn, guiobj, title='Pre-Edge - Up branch')
                    pe_dw_chsn = self.pe_label
                    pe_dw_avgd = self.aver_h_scans(self.pe_down, fields,
                            pe_dw_chsn, guiobj, title='Pre-Edge - Down branch')
            else:
                # If no pre-edge scans just return empty lists
                pe_up_chsn = []
                pe_up_avgd = []
                pe_dw_chsn = []
                pe_dw_avgd = []
        else: 
        # Not interactive - Consider and average all scans
            up_chsn = self.label
            up_avgd = self.aver_h_scans(self.up, fields, up_chsn, guiobj,
                                        title='Edge - Up branch')
            dw_chsn = self.label
            dw_avgd = self.aver_h_scans(self.down, fields, dw_chsn, guiobj,
                                        title='Edge - Down branch')
            if self.pre_edge:
                pe_up_chsn = self.pe_label
                pe_up_avgd = self.aver_h_scans(self.pe_up, fields, pe_up_chsn,
                                        guiobj, title='Pre-Edge - Up branch')
                pe_dw_chsn = self.pe_label
                pe_dw_avgd = self.aver_h_scans(self.pe_down, fields,
                            pe_dw_chsn, guiobj, title='Pre-Edge - Down branch')
            else:
                # If no pre-edge scans just return empty lists
                pe_up_chsn = []
                pe_up_avgd = []
                pe_dw_chsn = []
                pe_dw_avgd = []

        self.up_chsn = up_chsn
        self.up_aver = up_avgd
        self.dw_chsn = dw_chsn
        self.dw_aver = dw_avgd

        self.pe_up_chsn = pe_up_chsn
        self.pe_up_aver = pe_up_avgd
        self.pe_dw_chsn = pe_dw_chsn
        self.pe_dw_aver = pe_dw_avgd

    def plot_chs_avr(self, fields, guiobj, edge, up):
        '''
        Plot raw data, allow to choose which data will be used for
        analysis and average them.

        Parameters
        ----------
        fields : array
            Array with common magneti field scale, for interpolation and
            average of scans.

        guiobj : GUI object
            Provides GUI dialogs.

        edge : bool
            True for edge scans treatment
            False for pre-edge scan treatment.

        up : bool
            True for up branches treatment
            False for down branches treatment.

        Return
        ------
        list, array
        A list with the labels of chosen scans and an array with the
        average of the chosen scans.
        '''
        plt.figure(1)
        if edge:
            # Plot configurations for edge scans
            label = self.label
            idx = self.idx
            if up:
                # Plot configuration for up branches
                title = 'Edge ' + self.dtype + ' Up'
                plt.title(title)
                data = self.up
            else:
                # plot configurations for down branches
                title = 'Edge ' + self.dtype + ' Down'
                plt.title(title)
                data = self.down
        else:
            # Plot configurations for pre-edge scans
            label = self.pe_label
            idx = self.pe_idx
            if up:
                # Plot configurations for up branches
                title = 'Pre-Edge ' + self.dtype + ' Up'
                plt.title(title)
                data = self.pe_up
            else:
                # Plot configurations for down branches
                title = 'Pre-Edge ' + self.dtype + ' Down'
                plt.title(title)
                data = self.pe_down
        # Plot data
        for i in idx:
            h_col = 'H' + i
            plt.plot(data[h_col], data[i], label=label[idx.index(i)])

        plt.xlabel('H (T)')
        plt.ylabel(self.dtype)
        plt.legend()                
        plt.show()

        chsn = guiobj.chs_scns(label)

        avgd = aver_h_scans(data, fields, chsn, guiobj, title)

        return chsn, avgd

    def aver_h_scans(self, data, fields, chsn, guiobj, title):
        '''
        Perform the average of data scans. 
        If interactive mode, data scans and their average are shown
        together in a plot. 

        Parameters
        ----------
        data : pandas DataFrame
            scan data.

        fields : array
            Magnetic fields values at which average is calculated.

        chsn : list (str)
            Scan-numbers of scan to be averaged.

        guiobj: GUI object
            Provides GUI dialogs.

        title : str
            graph title.

        Returns
        -------
        array, containing the average of data scans.

        Notes
        -----
        To compute the average the common scale fields is used.
        All passed scans are interpolated with a linear spline
        (k=1 and s=0 in itp.UnivariateSpline) and evaluated along the
        common field scale.
        The interpolated data are eventually averaged.
        '''
        intrp = []

        if guiobj.interactive:
            plt.figure(1)
            plt.title(title)

        for scn in chsn:
            # Retrive scan number from label, remove PE- for pre-edge
            idx = scn.removeprefix('PE-')
            # chosen data

            # Univariate spline requires increasing x
            # Data are sorted by field
            ##sorted_data = data.sort_values(by=['H' + idx])
            
            ##x = sorted_data['H' + idx][1:].dropna()
            ##y = sorted_data[idx][1:].dropna()
            x, y = aver_duplicates(data, idx)

            if guiobj.interactive:
                # Plot data
                plt.plot(x, y, color='black')

            # Compute linear spline interpolation
            y_int = itp.UnivariateSpline(x, y, k=1, s=0)
            # Evaluate interpolation of field scan data on common field
            # scale and append to previous interpolations
            intrp.append(y_int(fields))

        # Average all inteprolated scans
        avgd = np.average(intrp, axis=0)

        if guiobj.interactive:
            plt.plot(fields, avgd, color='r', label='Average')
            plt.xlabel('H (T)')
            plt.ylabel(self.dtype)
            plt.legend()
            plt.show()

        return avgd


class FieldScan:
    '''
    Allow computation of XMCD hysteresis on magnetic field scan data.

    Attributes
    ----------
    pre_edge : bool
        True if pre-edge data are present.

    fields : array
        Common magnetic field scale for CR and CL scans data treatement.
        If present also scans collected at pre-edge energy are
        considered.

    cr_up : array
        Average of CR branch up scans.

    cr_down : array
        Average of CR branch down scans.

    cl_up : array
        Average of CL branch up scans.

    cl_down : array
        Average of CL branch down scans.
    
    cr_pe_up :array
        Average of CR branch up scans at pre-edge energy.
    
    cr_pe_down : array
        Average of CR branch down scans at pre-edge energy.

    cl_pe_up : array
        Average of CL branch up scans at pre-edge energy.
    
    cl_pe_down: array
        Average of CL branch down scans at pre-edge energy.

    edg_up : array
        XMCD for scan up branch considering only edge scans.

    edg_down 8: array
        XMCD for scan down branch considering only edge scans.

    edg_up_norm : array
        XMCD for scan up branch normalized to 1 considering only edge
        scans.

    edg_down_norm : array
        XMCD for scan down branch normalized to 1 considering only edge
        scans.

    up_w_pe : array
        XMCD for scan up branch with data normalized by pre-edge scans.

    dw_w_pe : array
        XMCD for scan down branch with data normalized by pre-edge
        scans.

    up_perc : array
        XMCD in percentage for scan up branch with data normalized by
        pre-edge scans.

    down_perc : array
        XMCD in percentage for scan down branch with data normalized by
        pre-edge scans.

    Methods
    -------
    scan_average(guiobj, pos, neg, log_dt)
        Separate CR and CL scans in up and down branches and average
        scans.

    compt_scanfield()
        Perform computation on magnetic field scan data in order to
        obtain XMCD hysteresis.
    '''

    def __init__(self, guiobj, confobj, h_scale, pos, neg, log_dt):
        '''
        At instantiation the field attribute is created with common
        magnetic field scale. CR and CL scans are separated in up and
        down branches and averaged. Finally XMCD is calculated.
        To do this scan_average and compt_scanfield method are called. 
        The FieldScan object created collects all the informations about X-Ray dichroism.

        Parameters
        ----------
        guiobj : GUI obj
            Provides GUI dialogs.

        confobj : Configuration object.

        h_scale : array
            Common magnetic field scale for data analysis.

        pos : ScanData obj
            Contains CR scans.

        neg : ScanData obj
            Contains CL scnas.

        log_dt : dict
            Collect data for logfile.
        '''
        # Initialize pre_edge False if pre-edges scans are present it
        # will be set True by scan_average
        self.pre_edge = False
        self.fields = h_scale

        self.scan_average(guiobj, pos, neg, log_dt)

        self.compt_scanfield()

    def scan_average(self, guiobj, pos, neg, log_dt):
        '''
        Separate CR and CL scans in up and down branches and average
        scans.
        Plot average results for edge and, if present, for pre-edge
        scans.

        Parameters
        ----------
        guiobj : GUI obj
            Provides GUI dialogs.

        pos : ScanData obj
            Contains CR scans.

        neg : ScanData obj
            Contains CL scnas.

        log_dt : dict
            Collect data for logfile.
        
        Return
        ------
        Instantiate attributes:

        cr_up : array
            Average of CR branch up scans.

        cr_down : array
            Average of CR branch down scans.

        cl_up : array
            Average of CL branch up scans.

        cl_down : array
            Average of CL branch down scans.
        
        cr_pe_up :array
            Average of CR branch up scans at pre-edge energy.
        
        cr_pe_down : array
            Average of CR branch down scans at pre-edge energy.

        cl_pe_up : array
            Average of CL branch up scans at pre-edge energy.
        
        cl_pe_down: array
            Average of CL branch down scans at pre-edge energy.

        Add keys to log_dt

        pos_up_chsn : list (str) with CR branch up chosen scans
        pos_dw_chsn : list (str) with CR branch down chosen scans
        pos_pe_up_chsn : list (str) with CR pre-edge branch up chosen
                        scans
        pos_pe_dw_chsn : list (str) with CR pre-edge branch down chosen
                        scans
        neg_up_chsn : list (str) with CL branch up chosen scans
        neg_dw_chsn : list (str) with CL branch down chosen scans
        neg_pe_up_chsn : list (str) with CL pre-edge branch up chosen
                        scans
        neg_pe_dw_chsn : list (str) with CL pre-edge branch down chosen
                        scans
        '''
        # Separate up and down branches and compute the  average of
        # scans
        pos.man_aver_h_scans(guiobj, self.fields)
        neg.man_aver_h_scans(guiobj, self.fields)

        # Fill log data with chosen scans
        log_dt['pos_up_chsn'] = pos.up_chsn
        log_dt['pos_dw_chsn'] = pos.dw_chsn
        log_dt['pos_pe_up_chsn'] = pos.pe_up_chsn
        log_dt['pos_pe_dw_chsn'] = pos.pe_dw_chsn
        log_dt['neg_up_chsn'] = neg.up_chsn
        log_dt['neg_dw_chsn'] = neg.dw_chsn
        log_dt['neg_pe_up_chsn'] = neg.pe_up_chsn
        log_dt['neg_pe_dw_chsn'] = neg.pe_dw_chsn

        self.cr_up = pos.up_aver
        self.cr_down = pos.dw_aver
        self.cl_up = neg.up_aver
        self.cl_down = neg.dw_aver

        plt.figure(1)
        plt.title('Up and Down branches')
        plt.subplot(221)
        plt.plot(self.fields, self.cr_up, label='CR Up')
        plt.ylabel('I (a.u.)')
        plt.xlabel('H (T)')
        plt.axhline(y=0, color='darkgray')
        plt.axvline(x=0, color='darkgray')
        plt.legend()

        plt.subplot(222)
        plt.plot(self.fields, self.cr_down, label='CR Down')
        plt.ylabel('I (a.u.)')
        plt.xlabel('H (T)')
        plt.axhline(y=0, color='darkgray')
        plt.axvline(x=0, color='darkgray')
        plt.legend()

        plt.subplot(223)
        plt.plot(self.fields, self.cl_up, label='CL Up')
        plt.ylabel('I (a.u.)')
        plt.xlabel('H (T)')
        plt.axhline(y=0, color='darkgray')
        plt.axvline(x=0, color='darkgray')
        plt.legend()

        plt.subplot(224)
        plt.plot(self.fields, self.cl_down, label='CL Down')
        plt.ylabel('I (a.u.)')
        plt.xlabel('H (T)')
        plt.axhline(y=0, color='darkgray')
        plt.axvline(x=0, color='darkgray')
        plt.legend()

        self.cr_pe_up = pos.pe_up_aver
        self.cr_pe_down = pos.pe_dw_aver
        self.cl_pe_up = neg.pe_up_aver
        self.cl_pe_down = neg.pe_dw_aver

        # Check if there are pre-edge scans.
        if (pos.pre_edge and neg.pre_edge):
            # Report the presence of pre-edge scans
            self.pre_edge = True

            plt.figure(2)
            plt.title('Up and Down pre-edge branches')
            plt.subplot(221)
            plt.plot(self.fields, self.cr_pe_up, label='CR pre-edge Up')
            plt.ylabel('I (a.u.)')
            plt.xlabel('H (T)')
            plt.axhline(y=0, color='darkgray')
            plt.axvline(x=0, color='darkgray')
            plt.legend()

            plt.subplot(222)
            plt.plot(self.fields, self.cr_pe_down, label='CR pre-edge Down')
            plt.ylabel('I (a.u.)')
            plt.xlabel('H (T)')
            plt.axhline(y=0, color='darkgray')
            plt.axvline(x=0, color='darkgray')
            plt.legend()

            plt.subplot(223)
            plt.plot(self.fields, self.cl_pe_up, label='CL pre-edge Up')
            plt.ylabel('I (a.u.)')
            plt.xlabel('H (T)')
            plt.axhline(y=0, color='darkgray')
            plt.axvline(x=0, color='darkgray')
            plt.legend()

            plt.subplot(224)
            plt.plot(self.fields, self.cl_pe_down, label='CL pre-edge Down')
            plt.ylabel('I (a.u.)')
            plt.xlabel('H (T)')
            plt.axhline(y=0, color='darkgray')
            plt.axvline(x=0, color='darkgray')
            plt.legend()
        plt.show()

    def compt_scanfield(self):
        '''
        Perform computation on magnetic field scan data in order to
        obtain XMCD hysteresis.
        Given two ScanData object (one for  CR and one for CL
        polarization):
        - compute X-Ray Dichroism as CL - CR for each scan field up and
          down branch;
        - compute normalized to 1 X-Ray dichroism for each scan field up
          and down branch;
        If pre-edge scans are present:
        - compute X-Ray Dichroism as CL - CR for each scan field up and
          down branch normalized by pre-edge data;
        - compute X-Ray dichroism in percentage normalized by pre-edge
            data for each scan field up and down branch.

        Returns
        -------
        Instantiate attributes:
        
        edg_up : array
            XMCD for scan up branch considering only edge scans.

        edg_down : array
            XMCD for scan down branch considering only edge scans.

        edg_up_norm : array
            XMCD for scan up branch normalized to 1 considering only
            edge scans.

        edg_down_norm : array
            XMCD for scan down branch normalized to 1 considering only
            edge scans.
        
        If pre-edge data are present:

        up_w_pe : array
            XMCD for scan up branch with data normalized by pre-edge
            scans.

        dw_w_pe : array
            XMCD for scan down branch with data normalized by pre-edge
            scans.

        up_perc : array
            XMCD in percentage for scan up branch with data normalized
            by pre-edge scans.

        down_perc : array
            XMCD in percentage for scan down branch with data normalized
            by pre-edge scans.
        '''
        # Compute XMCD for up and down branches @ edge energy
        self.edg_up = self.cl_up - self.cr_up
        self.edg_down = self.cl_down - self.cr_down

        # Normalize in [-1,1] edge XMCD up and down branches.
        # Normalization is performed considering values at high fields:
        # for each branch the average of 5 points at maximum and minimum
        # fields is considered then the normalization is computed
        # considering the greatest value of the two
        up_av_field1 = np.abs(np.average(self.edg_up[:5]))
        up_av_field2 = np.abs(np.average(self.edg_up[-5:-1]))

        dw_av_field1 = np.abs(np.average(self.edg_down[:5]))
        dw_av_field2 = np.abs(np.average(self.edg_down[-5:-1]))
        
        self.edg_up_norm = self.edg_up / np.maximum(up_av_field1, up_av_field2)
        self.edg_down_norm = self.edg_up / np.maximum(dw_av_field1,
                                                      dw_av_field2)

        if self.pre_edge:
            # Normalize branches by pre-edge data
            cl_up_over_pe = self.cl_up / self.cl_pe_up
            cr_up_over_pe = self.cr_up / self.cr_pe_up
            cl_dw_over_pe = self.cl_down / self.cl_pe_down
            cr_dw_over_pe = self.cr_down / self.cr_pe_down

            # XMCD for branches up and down considering data normalized
            # by pre-edge scans
            self.up_w_pe = cl_up_over_pe - cr_dw_over_pe
            self.dw_w_pe = cl_dw_over_pe - cr_dw_over_pe

            # XMCD in percentage for branches up and down considering
            # data normalized by pre-edge scans            
            self.up_perc = 200 * self.up_w_pe / (cl_up_over_pe +
                                                    cr_up_over_pe - 2)
            self.down_perc = 200 * self.dw_w_pe / (cl_dw_over_pe +
                                                    cr_dw_over_pe - 2)


def h_scale(guiobj, pos, neg, log_dt):
    '''
    Create the magnetic field scale used for data analysis.
    Range is selected considering the intersection of the field ranges
    of all the scans provided: the low-end is the highest of the minimum
    values and  the high-end is the lowest of the maximum values between
    all the magnetic field arrays.

    The number of points of the magnetic field scale returned is by
    default the average of the number of points of the field scales of
    data scans.
    If in interactive mode GUI dialogues are provided to set the 
    number of points.

    Parameters
    ----------
    guiobj : GUI obj
        Provides GUI dialogs.

    pos: ScanData obj
        Positive scans (CR for XMCD and XNCD, LH for XNLD).

    neg: ScanData obj
        Negative scnas (CL for XMCD and XNCD, LV for XNLD).      

    log_dt : dict
        Collect data for logfile.

    Return
    ------
    array
        Common magnetic field scale for positive and negative XAS scans
        data treatement

    Add h_scale key to log_dt with min, max and number of points of
    field scale. 
    '''
    h_maxs = []  # Maxima of magnetic field scans
    h_mins = []  # Minima of magnetic field scans
    h_len = []  # Number of different fields in each scan
    
    # Collect maxima and minima and count the fields for edge scans
    for i in pos.idx:
        h_maxs.append(np.amax(pos.raw_imp['H' + i]))
        h_mins.append(np.amin(pos.raw_imp['H' + i]))
        h_len.append(h_num_points(pos.raw_imp['H' + i].dropna()))
    for i in neg.idx:
        h_maxs.append(np.amax(neg.raw_imp['H' + i]))
        h_mins.append(np.amin(neg.raw_imp['H' + i]))
        h_len.append(h_num_points(neg.raw_imp['H' + i].dropna()))
    # If present pre-edge scans do the same for them
    if pos.pre_edge:
        for i in pos.pe_idx:
            h_maxs.append(np.amax(pos.pe_raw_imp['H' + i]))
            h_mins.append(np.amin(pos.pe_raw_imp['H' + i]))
            h_len.append(h_num_points(pos.pe_raw_imp['H' + i].dropna()))
    if neg.pre_edge:
        for i in neg.pe_idx:
            h_maxs.append(np.amax(neg.pe_raw_imp['H' + i]))
            h_mins.append(np.amin(neg.pe_raw_imp['H' + i]))
            h_len.append(h_num_points(neg.pe_raw_imp['H' + i].dropna()))
    
    # Compute min, max and default langth of energy range
    # Set decimal place to round the first OoM higher than tolerance in
    # h_num_points
    h_min = np.around(np.amax(h_mins), 3)
    h_max = np.around(np.amin(h_maxs), 3)
    h_av_len = np.around(np.average(h_len), 0)

    # Set number of points of energy scale
    if guiobj.interactive:
        n_points = guiobj.num_pnts(h_av_len)
    else:
        n_points = int(h_av_len)

    log_dt['h_scale'] = [h_min, h_max, n_points]

    return np.linspace(h_min,h_max, n_points)

def h_num_points(h_arr):
    '''
    Count the number of different fields present in h_arr.

    During on-fly scans the magnetic field speed variation is not
    uniform. In particular at the extremes of the field range (so at the
    beginning of the scan, at the end of the scan and each time the
    magntic field scan direction is reversed) there are different points
    at the same field value within the experimental uncertainty. This
    means that the number of acquisitions is not the same as the number
    of different magnetic fields measured. In order to obtain a reliable
    value for the number of different magnetic fileds measured during
    the scan a tolerance value is introduced and consecutive fields
    differing less than tolerance are considered the same and counted
    only once.

    Parameters
    ----------
    h_arr : array
        array with magnetic field values.

    Return
    ------
    int, the number of different magnetic fields present in h_arr.
    '''
    tolerance = 5e-4  # Field tolerance (5 Oe)

    n = 1  # Loop starts from element 1 of h_arr

    for i in range(1, len(h_arr)):
        if (h_arr[i] - h_arr[i-1]) < tolerance:
            continue
        else:
            n += 1

    return n

def aver_duplicates(data, idx):
    '''
    Check if in the field column idx there are repeated values. In case
    compute average of corresponding XAS values. This is required by the
    fact the univariate spline method employed to interpolate date on a
    common field scale requires strictly increasing x values.

    Parameters
    ----------
    data : Pandas DataFrame
        contains data.

    idx : str
        identifies scan to be analyse.

    Return
    ------
    array, field scale with unique values.
    array, XAS values corresponding to unique values field scale.
    '''

    # Remove NaNs and sort data by field column
    sorted_data = data.dropna().sort_values(by=['H' + idx])

    # Extract unique field values and corresonding indexes of first 
    # unique value occurence
    x, indexes = np.unique(sorted_data['H' + idx], return_index=True)

    y = []

    # From one index and the next there are repeated values of fields.
    # This for loop average the XAS value between consecutive indexes.
    for i in range(len(indices)):
        if i == len(indices) - 1:
            y.append(np.average(sorted_data[idx][indexes[i]:]))
        else:
            y.append(np.average(sorted_data[idx][indexes[i]:indexes[i+1]]))

    return x, y