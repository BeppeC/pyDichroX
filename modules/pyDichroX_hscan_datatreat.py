"""
pyDichroX_hscan_datatreat.py

Classes and methods for hysteresys data analysis.

Classes
--------
ScanData : Collect raw data for magnetic field scan experiments and
    provides methods to average them.

FieldScan : Allow computation of XMCD hysteresis on magnetic field scan
    data for on the fly collection experiments.

FieldPtScan: Allow computation of XMCD hysteresis on magnetic field scan
    data for point by point collection experiments.

Methods
-------
h_scale(guiobj, pos, neg, log_dt)
    Create the magnetic field scale used for data analysis.

h_scale_sp_bracnh(guiobj, pos, neg, log_dt)
    Create the magnetic field scales used for data analysis - for data
    acquisitions where hysteresis branch are splitted for positive and
    negative fields values.

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

    up_idx : list (str)
        Scan indexes - for up scans collected at edge energy - used only
        for splitted branch analysis to create common field scale.

    dw_idx : list (str)
        Scan indexes - for down scans collected at edge energy
        - used only for splitted branch analysis to create common field
        scale.

    pe_up_idx : list (str)
        Scan indexes - for up scans collected at pre-edge energy
        - used only for splitted branch analysis to create common field
        scale.

    pe_dw_idx : list (str)
        Scan indexes - for down scans collected at pre-edge energy
        - used only for splitted branch analysis to create common field
        scale.

    raw_imp : pandas DataFrame
        Collect imported raw data at edge energy.
        For hysteresis on the fly analysis it consist of a couple of
        columns for each scan ('HscanNum' and 'scanN') containing the
        field and data values respectively
        For hysteresis point by point it consists of just three columns
        'H', 'I' and 't' where field values, data and times of all scans
        are appended.

    pe_raw_imp : pandas DataFrame
        Collect imported raw data at pre-edge energy.
        For hysteresis on the fly analysis it consist of a couple of
        columns for each scan ('HscanNum' and 'scanN') containing the
        field and data values respectively
        For hysteresis point by point it consists of just three columns
        'H', 'I' and 't' where field values, data and times of all scans
        are appended.

    pre_edge : bool
        True if pre-edge data are present.
        Only for Hysteresis on the fly analisys.

    up_hp : bool
        Variable to check if up scans with postive fields are present
        Only for splitted branches.

    up_hn : bool
        Variable to check if up scans with negative fields are present
        Only for splitted branches.

    down_hp : bool
        Variable to check if down scans with postive fields are present
        Only for splitted branches.

    down_hn : bool
        Variable to check if down scans with negative fields are present
        Only for splitted branches.

    pe_up_hp : bool
        Variable to check if pe up scans with postive fields are present
        Only for splitted branches.

    pe_up_hn : bool
        Variable to check if pe up scans with negative fields are
        present - Only for splitted branches.

    pe_down_hp : bool
        Variable to check if pe down scans with postive fields are
        present - Only for splitted branches.

    pe_down_hn : bool
        Variable to check if pe down scans with negative fields are
        present - Only for splitted branches.

    up_label : list (str)
        Labels for graphs - for scans collected at edge energy up scans
        Only for splitted branch.

    dw_label : list (str)
        Labels for graphs - for scans collected at edge energy down
        scans - Only for splitted branch.

    pe_up_label : list (str)
        Labels for graphs - for scans collected at pre-edge energy up
        scans - Only for splitted branch.

    pe_dw_label : list (str)
        Labels for graphs - for scans collected at pre-edge energy down
        scans - Only for splitted branch.

    sel_data : pandas DataFrame
        Collect data at edge energy inside the chosen time window.
        Only for Hysteresis point by pont analysis.

    pe_sel_data : pandas DataFrame
        Collect data at pre-edge energy inside the chosen time window.
        Only for Hysteresis point by pont analysis.

    min_t : list
        Collect the minima of the timescales of data scans.
        Only for Hysteresis point by pont analysis.

    max_t : list
        Collect the maxima of the timescales of data scans.
        Only for Hysteresis point by pont analysis.

    aver : pandas DataFrame
        Collect averaged data at edge energy for time averaged analysis
        or data at edge energy interpolated on a common time scale for
        time splitted analysis.
        Only for Hysteresis point by point analysis.

    pe_aver : pandas DataFrame
        Collect averaged data at pre-edge energy for time averaged
        analysis or data at pre-edge energy interpolated on a common
        time scale for time splitted analysis.
        Only for Hysteresis point by point analysis.

    up : pandas DataFrame
        Collect data related to magntic filed scan up branch.
        Only for Hysteresis on the fly analisys.

    down : pandas DataFrame
        Collect data related to magntic filed scan down branch.
        Only for Hysteresis on the fly analisys.

    up_chsn : list
        Collect labels of ub branch scans chosen for analysis.
        Only for Hysteresis on the fly analisys.

    dw_chsn : list
        Collect labels of down branch scans chosen for analysis.
        Only for Hysteresis on the fly analisys.

    up_aver : array
        Contains average of chosen up branch scans computed on common
        field scale.
        Only for Hysteresis on the fly analisys.

    dw_aver : array
        Contains average of chosen down branch scans computed on common
        field scale.
        Only for Hysteresis on the fly analisys.

    pe_label : list (str)
        Labels for graphs - for scans collected at pre-edge energy.

    pe_idx : list (str)
        Scan indexes - for scans collected at pre-edge energy.

    pe_raw_imp : pandas DataFrame
        Collect imported raw data at pre-edge energy.
        For hysteresis on the fly analysis it consist of a couple of
        columns for each scan ('HscanNum' and 'scanN') containing the
        field and data values respectively
        For hysteresis point by point it consists of just three columns
        'H', 'I' and 't' where field values, data and times of all scans
        are appended.

    pe_up : pandas DataFrame.
        Collect data related to magntic filed scan up branch - for scans
        measured at pre-edge energy.
        Only for Hysteresis on the fly analisys.

    pe_down : pandas DataFrame.
        Collect data related to magntic filed scan up branch - for scans
        measured at pre-edge energy.
        Only for Hysteresis on the fly analisys.

    pe_up_chsn : list
        Collect labels of ub branch scans chosen for analysis.
        Only for Hysteresis on the fly analisys.

    pe_dw_chsn : list
        Collect labels of down branch scans chosen for analysis.
        Only for Hysteresis on the fly analisys.

    pe_up_aver : array
        Contains average of chosen up branch scans computed on common
        field scale.
        Only for Hysteresis on the fly analisys.

    pe_dw_aver : array
        Contains average of chosen down branch scans computed on common
        field scale.
        Only for Hysteresis on the fly analisys.

    dtype : str
        Identifies CR and CL data, used for graph labelling.

    Methods
    -------
    up_n_down()
        Separate input data in up branches and down branches.
        For hysteresis on the fly analysis.

    man_aver_e_scans(guiobj, confobj, fields)
        Manage the choice of scans to be averaged and return the average
        of selected scans.

    plot_chs_avr(fields, guiobj, confobj, edge, up)
        Plot raw data, allow to choose which data will be used for
        analysis and average them.

    aver_h_scans(data, fields, chsn, guiobj, confobj, title)
        Perform the average of data scans.
    '''

    def __init__(self, guiobj, confobj):
        '''
        Initialize attributes label, idx, raw_imp, up, down and pre-edge
        corresponding pe_label, pe_idx, pe_raw_imp, pe_up, pe_down.

        Prameters
        ---------
        guiobj : GUI object
            Provides GUI dialogs.

        confobj : Configurations object.
        '''
        self.label = []
        self.idx = []
        self.raw_imp = pd.DataFrame()

        self.pe_label = []
        self.pe_idx = []
        self.pe_raw_imp = pd.DataFrame()
        # Initialize False setted True by separate_hscans in io
        self.pre_edge = False

        # Only for hysteresis on the fly is currently cosidered the
        # possibility to seprate up from down branches.
        if guiobj.analysis == 'hyst_fly':
            self.up = pd.DataFrame()
            self.down = pd.DataFrame()
            self.pe_up = pd.DataFrame()
            self.pe_down = pd.DataFrame()

            if confobj.spl_brnch:
                self.up_hp = False
                self.up_hn = False
                self.down_hp = False
                self.down_hn = False
                self.pe_up_hp = False
                self.pe_up_hn = False
                self.pe_down_hp = False
                self.pe_down_hn = False

                self.up_idx = []
                self.dw_idx = []
                self.pe_up_idx = []
                self.pe_dw_idx = []

                self.up_label = []
                self.dw_label = []
                self.pe_up_label = []
                self.pe_dw_label = []
        else:
            self.min_t = []
            self.max_t = []

    def up_n_down(self):
        '''
        Separate input data in up branches and down branches.
        For Hysteresis on the fly analysis with complete branches
        acquisition.

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

    def man_aver_h_scans(self, guiobj, confobj, fields):
        '''
        Manage the choice of scans to be averaged and return the average
        of selected scans.

        Parameters
        ----------
        guiobj : GUI object
            Provides GUI dialogs.

        confobj : configuration object.

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
            # for splitted branch in order to chose there must be more
            # than 2 scans, one for positive and one for negative fields
            if confobj.spl_brnch:
                min_num_sc = 4
            else:
                min_num_sc = 1
            # Choose edge scans.
            if len(self.idx) > min_num_sc:
                # Edge Up data
                # Loop until choice is confirmed
                isok = 0
                while not isok:
                    up_chsn, up_avgd = self.plot_chs_avr(fields, guiobj,
                                                         confobj, edge=True,
                                                         up=True)
                    isok = guiobj.confirm_choice()

                # Edge Down data
                # Loop until choice is confirmed
                isok = 0
                while not isok:
                    dw_chsn, dw_avgd = self.plot_chs_avr(fields, guiobj,
                                                         confobj, edge=True,
                                                         up=False)
                    isok = guiobj.confirm_choice()
            else:
                # If there is just the minimun number of scans is
                # useless the choose, just interpolate it with the
                # common field scale.
                if confobj.spl_brnch:
                    up_chsn = self.up_label
                    dw_chsn = self.dw_label
                else:
                    up_chsn = self.label
                    dw_chsn = self.label
                up_avgd = self.aver_h_scans(self.up, fields, up_chsn, guiobj,
                                            confobj, title='Edge - Up branch')
                dw_avgd = self.aver_h_scans(self.down, fields, dw_chsn, guiobj,
                                        confobj, title='Edge - Down branch')
            # Choose pre-edge scans
            if self.pre_edge:
                if len(self.pe_idx) > min_num_sc:
                    # pre-edge Up data
                    # Loop until choice is confirmed
                    isok = 0
                    while not isok:
                        pe_up_chsn, pe_up_avgd = self.plot_chs_avr(
                            fields, guiobj, confobj, edge=False, up=True)
                        isok = guiobj.confirm_choice()

                    # pre-edge Down data
                    # Loop until choice is confirmed
                    isok = 0
                    while not isok:
                        pe_dw_chsn, pe_dw_avgd = self.plot_chs_avr(
                            fields, guiobj, confobj, edge=False, up=False)
                        isok = guiobj.confirm_choice()
                else:  # There is just one pre-edge scan
                   # If there is just the minimun number of scans is
                    # useless the choose, just interpolate it with the
                    # common field scale.
                    if confobj.spl_brnch:
                        pe_up_chsn = self.pe_up_label
                        pe_dw_chsn = self.pe_dw_label
                    else:
                        pe_up_chsn = self.pe_label
                        pe_dw_chsn = self.pe_label

                    pe_up_avgd = self.aver_h_scans(
                        self.pe_up, fields, pe_up_chsn, guiobj, confobj,
                        title='Pre-Edge - Up branch')
                    pe_dw_avgd = self.aver_h_scans(
                        self.pe_down, fields, pe_dw_chsn, guiobj, confobj,
                        title='Pre-Edge - Down branch')
            else:
                # If no pre-edge scans just return empty lists
                pe_up_chsn = []
                pe_up_avgd = []
                pe_dw_chsn = []
                pe_dw_avgd = []
        else:
            # Not interactive - Consider and average all scans
            if confobj.spl_brnch:
                up_chsn = self.up_label
                dw_chsn = self.dw_label
            else:
                up_chsn = self.label
                dw_chsn = self.label

            up_avgd = self.aver_h_scans(self.up, fields, up_chsn, guiobj,
                                        confobj, title='Edge - Up branch')
            dw_avgd = self.aver_h_scans(self.down, fields, dw_chsn, guiobj,
                                        confobj, title='Edge - Down branch')
            if self.pre_edge:
                if confobj.spl_brnch:
                    pe_up_chsn = self.pe_up_label
                    pe_dw_chsn = self.pe_dw_label
                else:
                    pe_up_chsn = self.pe_label
                    pe_dw_chsn = self.pe_label

                pe_up_avgd = self.aver_h_scans(
                    self.pe_up, fields, pe_up_chsn, guiobj, confobj,
                    title='Pre-Edge - Up branch')
                pe_dw_avgd = self.aver_h_scans(
                    self.pe_down, fields, pe_dw_chsn, guiobj, confobj,
                    title='Pre-Edge - Down branch')
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

    def plot_chs_avr(self, fields, guiobj, confobj, edge, up):
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

        confobj : configuration object.

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
            # if splitted branch up and down have different idx
            if confobj.spl_brnch:
                if up:
                    label = self.up_label
                    idx = self.up_idx
                else:
                    label = self.dw_label
                    idx = self.dw_idx
            else:
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
            # if splitted branch up and down have different idx
            if confobj.spl_brnch:
                if up:
                    label = self.pe_up_label
                    idx = self.pe_up_idx
                else:
                    label = self.pe_dw_label
                    idx = self.pe_dw_idx
            else:
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

        avgd = self.aver_h_scans(data, fields, chsn, guiobj, confobj, title)

        return chsn, avgd

    def aver_h_scans(self, data, fields, chsn, guiobj, confobj, title):
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

        confobj : configuration object.

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
        # If splitted branches interpolated variable consider one array
        # for positive field branch - index 0 - and one for negative
        # field - index 1 -
        if confobj.spl_brnch:
            intrp = [[], []]
        # If no splitted branches just one array for interpolated data
        else:
            intrp = []

        if guiobj.interactive:
            plt.figure(1)
            plt.title(title)

        for scn in chsn:
            # Retrive scan number from label, remove PE- for pre-edge
            idx = scn.removeprefix('PE-')

            # Univariate spline requires increasing x
            # x and y are data averaged in order to remove duplicates in
            # field values
            x, y = aver_duplicates(data, idx)

            if guiobj.interactive:
                # Plot data
                plt.plot(x, y, color='black')

            # Compute linear spline interpolation
            y_int = itp.UnivariateSpline(x, y, k=1, s=0)

            # Evaluate interpolation of field scan data on common field
            # scale and append to previous interpolations

            # if splitted branches append interpolated half branches
            # scan in appropriate array
            if confobj.spl_brnch:
                # positive field branc
                if np.amax(data['H' + idx]) > 0:
                    intrp[1].append(y_int(fields[1]))
                # negative field branch
                else:
                    intrp[0].append(y_int(fields[0]))
            # if no splitted branches just append interpolations
            else:
                intrp.append(y_int(fields))

        # Average all inteprolated scans
        if confobj.spl_brnch:
            avgd = np.concatenate((np.average(intrp[0], axis=0),
                                   np.average(intrp[1], axis=0)))
        else:
            avgd = np.average(intrp, axis=0)

        if guiobj.interactive:
            if confobj.spl_brnch:
                plt.plot(fields[2], avgd, color='r', label='Average')
            else:
                plt.plot(fields, avgd, color='r', label='Average')
            plt.xlabel('H (T)')
            plt.ylabel(self.dtype)
            plt.legend()
            plt.show()

        return avgd

    def aver_pt_scans(self, guiobj, time_scale):
        '''
        Perform the average of data scans. 
        If interactive mode, data scans and their average are shown
        together in a plot. 

        Parameters
        ----------
        guiobj: GUI object
            Provides GUI dialogs.

        time_scale : array
            Common time scale.

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
        if guiobj.analysis == 'hyst_t_aver':
            # Data at the same field collect in window time are averaged
            aver = pd.DataFrame(columns=['H', 'I', 'dI'])
            # Group data by field values and average them.
            for Hval, dat in self.sel_data.groupby('H'):
                # In case the baseline subtraction must be implemented
                # this is the place where to put it.
                # if guiobj.bsl_sub:
                # Subtract baseline and average data
                # pass
                # else:
                row = pd.Series([Hval, np.nanmean(dat['I'], dtype=np.float64),
                                 np.nanstd(dat['I'])], index=['H', 'I', 'dI'])
                #aver = aver.append(row, ignore_index=True)
                aver = pd.concat([aver, row.to_frame().T], ignore_index=True)
            # Sort data by H
            self.aver = aver.sort_values(by=['H'], ignore_index=True)
            if self.pre_edge:
                pe_aver = pd.DataFrame(columns=['H', 'I', 'dI'])
                # Group data by field values and average them.
                for Hval, dat in self.pe_sel_data.groupby('H'):
                    # In case the baseline subtraction must be implemented
                    # this is the place where to put it.
                    # if guiobj.bsl_sub:
                    # Subtract baseline and average data
                    # pass
                    # else:
                    row = pd.Series(
                        [Hval, np.nanmean(dat['I'], dtype=np.float64),
                        np.nanstd(dat['I'])], index=['H', 'I', 'dI'])
                    #pe_aver = pe_aver.append(row, ignore_index=True)
                    pe_aver = pd.concat([pe_aver, row.to_frame().T],
                                        ignore_index=True)
                # Sort data by H
                self.pe_aver = pe_aver.sort_values(by=['H'], ignore_index=True)
        else:
            aver = pd.DataFrame(columns=['H', 'I', 't'])
            # Group data by fields
            for Hval, dat in self.sel_data.groupby('H'):
                # Sort data by time
                srt = dat.sort_values(by=['t'])
                # Interpolate data on common time scale
                interp = itp.interp1d(srt['t'], srt['I'], kind='slinear',
                                      fill_value='extrapolate')
                # Append interpolated data
                temp = pd.DataFrame({'H': Hval, 'I': interp(time_scale),
                                     't': time_scale})
                #aver = aver.append(temp, ignore_index=True)
                aver = pd.concat([aver, temp], ignore_index=True)
            # Sort data by t and H
            self.aver = aver.sort_values(by=['t', 'H'], ignore_index=True)
            if self.pre_edge:
                pe_aver = pd.DataFrame(columns=['H', 'I', 't'])
                # Group data by fields
                for Hval, dat in self.pe_sel_data.groupby('H'):
                    # Sort data by time
                    srt = dat.sort_values(by=['t'])
                    # Interpolate data on common time scale
                    interp = itp.interp1d(srt['t'], srt['I'], kind='slinear',
                                          fill_value='extrapolate')
                    # Append interpolated data
                    temp = pd.DataFrame({'H': Hval, 'I': interp(time_scale),
                                         't': time_scale})
                    #pe_aver = pe_aver.append(temp, ignore_index=True)
                    pe_aver = pd.concat([pe_aver, temp], ignore_index=True)
                # Sort data by t and H
                self.pe_aver = pe_aver.sort_values(by=['t', 'H'],
                                                   ignore_index=True)


class FieldScan:
    '''
    Allow computation of XMCD hysteresis on magnetic field scan data.

    Attributes
    ----------
    pre_edge : bool
        True if pre-edge data are present.

    fields : array
        Common magnetic field scale for CR and CL scans data treatement.

    flds : array
        Final common magnetic filed scale for CR and CL scans.
        If no splitted branch is equal to field otherwise it is the
        union of positive and negative fields half-branches.

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
    scan_average(guiobj, confobj, pos, neg, log_dt)
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
        The FieldScan object created collects all the informations about
        X-Ray dichroism.

        Parameters
        ----------
        guiobj : GUI obj
            Provides GUI dialogs.

        confobj : configuration object.

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

        self.scan_average(guiobj, confobj, pos, neg, log_dt)

        self.compt_scanfield()

    def scan_average(self, guiobj, confobj, pos, neg, log_dt):
        '''
        Separate CR and CL scans in up and down branches and average
        scans.
        Plot average results for edge and, if present, for pre-edge
        scans.

        Parameters
        ----------
        guiobj : GUI obj
            Provides GUI dialogs.

        confobj : configuration object.

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

        flds : array
            Final common magnetic filed scale for CR and CL scans.
            If no splitted branch is equal to field otherwise it is the
            union of positive and negative fields half-branches.

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
        pos.man_aver_h_scans(guiobj, confobj, self.fields)
        neg.man_aver_h_scans(guiobj, confobj, self.fields)

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

        if confobj.spl_brnch:
            self.flds = self.fields[2]
        else:
            self.flds = self.fields

        plt.figure(1)
        plt.title('Up and Down branches')
        plt.subplot(221)
        plt.plot(self.flds, self.cr_up, label='CR Up')
        plt.ylabel('I (a.u.)')
        plt.xlabel('H (T)')
        plt.axhline(y=0, color='darkgray')
        plt.axvline(x=0, color='darkgray')
        plt.legend()

        plt.subplot(222)
        plt.plot(self.flds, self.cr_down, label='CR Down')
        plt.ylabel('I (a.u.)')
        plt.xlabel('H (T)')
        plt.axhline(y=0, color='darkgray')
        plt.axvline(x=0, color='darkgray')
        plt.legend()

        plt.subplot(223)
        plt.plot(self.flds, self.cl_up, label='CL Up')
        plt.ylabel('I (a.u.)')
        plt.xlabel('H (T)')
        plt.axhline(y=0, color='darkgray')
        plt.axvline(x=0, color='darkgray')
        plt.legend()

        plt.subplot(224)
        plt.plot(self.flds, self.cl_down, label='CL Down')
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
            plt.plot(self.flds, self.cr_pe_up, label='CR pre-edge Up')
            plt.ylabel('I (a.u.)')
            plt.xlabel('H (T)')
            plt.axhline(y=0, color='darkgray')
            plt.axvline(x=0, color='darkgray')
            plt.legend()

            plt.subplot(222)
            plt.plot(self.flds, self.cr_pe_down, label='CR pre-edge Down')
            plt.ylabel('I (a.u.)')
            plt.xlabel('H (T)')
            plt.axhline(y=0, color='darkgray')
            plt.axvline(x=0, color='darkgray')
            plt.legend()

            plt.subplot(223)
            plt.plot(self.flds, self.cl_pe_up, label='CL pre-edge Up')
            plt.ylabel('I (a.u.)')
            plt.xlabel('H (T)')
            plt.axhline(y=0, color='darkgray')
            plt.axvline(x=0, color='darkgray')
            plt.legend()

            plt.subplot(224)
            plt.plot(self.flds, self.cl_pe_down, label='CL pre-edge Down')
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
        self.edg_down_norm = self.edg_down / np.maximum(dw_av_field1,
                                                        dw_av_field2)

        if self.pre_edge:
            # Normalize branches by pre-edge data
            cl_up_over_pe = self.cl_up / self.cl_pe_up
            cr_up_over_pe = self.cr_up / self.cr_pe_up
            cl_dw_over_pe = self.cl_down / self.cl_pe_down
            cr_dw_over_pe = self.cr_down / self.cr_pe_down

            # XMCD for branches up and down considering data normalized
            # by pre-edge scans
            self.up_w_pe = cl_up_over_pe - cr_up_over_pe
            self.dw_w_pe = cl_dw_over_pe - cr_dw_over_pe

            # XMCD in percentage for branches up and down considering
            # data normalized by pre-edge scans
            self.up_perc = 200 * self.up_w_pe / (cl_up_over_pe
                                                 + cr_up_over_pe - 2)
            self.down_perc = 200 * self.dw_w_pe / (cl_dw_over_pe
                                                   + cr_dw_over_pe - 2)


class FieldPtScan:
    '''
    Allow computation of XMCD hysteresis on magnetic field scan data for
    point by point collection experiments.

    Attributes
    ----------
    pre_edge : bool
        True if pre-edge data are present.

    time_scale : array
        Common time scale for CR and CL scans data treatement.
        If present also scans collected at pre-edge energy are
        considered.

    xmcd : pandas DataFrame
        Collect data obtained from computation of XMCD.
        It consists of H, CR, CL, XMCD columns. If pre-edge scans are
        present there are also CRpe, CLpe and XMCD_norm. For time
        averaged analysis also error columns are present.

    subtitle : str
        Used for graph labelling and output file name. It includes the
        number of the first and last scans and the window time used for
        analysis.

    Methods
    -------
    time_scale(guiobj, pos, neg, log_dt)
        Create commom time scale for hysteresis point by point analysis.

    scan_average(guiobj, pos, neg, log_dt)
        Separate CR and CL scans in up and down branches and average
        scans.

    compt_scanfield()
        Perform computation on magnetic field scan data in order to
        obtain XMCD hysteresis.
    '''

    def __init__(self, guiobj, pos, neg, log_dt):
        '''
        At instantiation the time_scale attribute is created with common
        time scale. XMCD is calculated on averaged data in scanData 
        objects. 
        The FieldScan object created collects all the informations about
        X-Ray dichroism.

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
        '''
        if (pos.pre_edge and neg.pre_edge):
            self.pre_edge = True
        else:
            self.pre_edge = False

        self.time_scale(guiobj, pos, neg, log_dt)

        self.compt_pt_scan(guiobj, pos, neg)

    def time_scale(self, guiobj, pos, neg, log_dt):
        '''
        Create commom time scale for hysteresis point by point analysis.
        Start and end points of the scale are provided by
        guiobj.acq_times, the number of points of the scale is instead
        obtained as average of the number of points between start time
        and end time of all the scans.

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
        Set time_scale attribute with the common time scale.
        Add t_scale key to log_dt, list with starting, ending, number of
        points and stepsize used for time scale definition.

        '''
        len_t = []
        st_t, end_t = guiobj.acq_times(pos, neg)

        pos.sel_data = pos.raw_imp[(pos.raw_imp['t'] >= st_t)
                                   & (pos.raw_imp['t'] <= end_t)]
        neg.sel_data = neg.raw_imp[(neg.raw_imp['t'] >= st_t)
                                   & (neg.raw_imp['t'] <= end_t)]
        # Collect the number of time points in the time window for each
        # field
        for Hval, dat in pos.sel_data.groupby('H'):
            len_t.append(len(dat['t']))
        for Hval, dat in neg.sel_data.groupby('H'):
            len_t.append(len(dat['t']))

        if pos.pre_edge and neg.pre_edge:
            pos.pe_sel_data = pos.pe_raw_imp[(pos.pe_raw_imp['t'] >= st_t)
                                             & (pos.pe_raw_imp['t'] <= end_t)]
            neg.pe_sel_data = neg.pe_raw_imp[(neg.pe_raw_imp['t'] >= st_t)
                                             & (neg.pe_raw_imp['t'] <= end_t)]
            for Hval, dat in pos.pe_sel_data.groupby('H'):
                len_t.append(len(dat['t']))
            for Hval, dat in neg.pe_sel_data.groupby('H'):
                len_t.append(len(dat['t']))

        num_t = int(np.around(np.average(len_t), 0))
        step_t = (end_t - st_t) / (num_t - 1)

        log_dt['t_scale'] = [st_t, end_t, num_t, step_t]

        self.time_scale = np.linspace(st_t, end_t, num_t)

    def compt_pt_scan(self, guiobj, pos, neg):
        '''
        Compute XMCD for point by point hysteresis analysis.
        If time average is considered XMCD is computed for each field
        averaging all data collected by time scan. Standard deviation
        of the averages are considered as errors and propagated during
        computation.
        If split time is considered XMCD is computed for each fiedl and
        for each point at common time scale.

        Parameters
        ----------
        guiobj : GUI obj
            Provides GUI dialogs.

        pos : ScanData obj
            Contains CR scans.

        neg : ScanData obj
            Contains CL scnas.

        Return
        ------
        Instantiate attributes xmcd attribute with all xmcd computed
        values.

        '''
        # Average pos and neg scanData
        pos.aver_pt_scans(guiobj, self.time_scale)
        neg.aver_pt_scans(guiobj, self.time_scale)

        # First create xmcd DataFrame
        if guiobj.analysis == 'hyst_t_aver':
            # Rename columns to join them
            pos.aver.rename(columns={'H': 'H', 'I': 'CR', 'dI': 'dCR'},
                            inplace=True)
            neg.aver.rename(columns={'H': 'H', 'I': 'CL', 'dI': 'dCL'},
                            inplace=True)
            # Join pos and neg aver using H as reference
            xmcd = pos.aver.join(neg.aver.set_index('H'), on='H')

            if self.pre_edge:
                pos.pe_aver.rename(columns={'H': 'H', 'I': 'CRpe',
                                            'dI': 'dCRpe'}, inplace=True)
                neg.pe_aver.rename(columns={'H': 'H', 'I': 'CLpe',
                                            'dI': 'dCLpe'}, inplace=True)
                # Join pos and neg pe_aver using H as reference
                xmcd = xmcd.join(pos.pe_aver.set_index('H'), on='H')
                xmcd = xmcd.join(neg.pe_aver.set_index('H'), on='H')
        else:
            # For time split analysis group by H, sort by t and add the
            # columns
            pos.aver.rename(columns={'H': 'H', 'I': 'CR', 't': 't'},
                            inplace=True)
            pos_gb = pos.aver.groupby('H')
            neg_gb = neg.aver.groupby('H')

            if self.pre_edge:
                pospe_gb = pos.pe_aver.groupby('H')
                negpe_gb = neg.pe_aver.groupby('H')

                xmcd = pd.DataFrame(columns=['H', 'CR', 'CL', 'CRpe', 'CLpe',
                                             't'])
                for Hval, dat in pos_gb:
                    xmcd_p = pos.aver[pos.aver['H'] == Hval].sort_values(
                        by=['t'])
                    xmcd_p['CL'] = neg_gb.get_group(Hval).sort_values(
                        by=['t'])['I']
                    xmcd_p['CRpe'] = pospe_gb.get_group(Hval).sort_values(
                        by=['t'])['I']
                    xmcd_p['CLpe'] = negpe_gb.get_group(Hval).sort_values(
                        by=['t'])['I']
                    #xmcd = xmcd.append(xmcd_p, ignore_index=True)
                    xmcd = pd.concat([xmcd, xmcd_p], ignore_index=True)
            else:
                xmcd = pd.DataFrame(columns=['H', 'CR', 'CL', 't'])
                for Hval, dat in pos_gb:
                    xmcd_p = pos.aver[pos.aver['H'] == Hval].sort_values(
                        by=['t'])
                    xmcd_p['CL'] = neg_gb.get_group(Hval).sort_values(
                        by=['t'])['I']
                    #xmcd = xmcd.append(xmcd_p, ignore_index=True)
                    xmcd = pd.concat([xmcd, xmcd_p], ignore_index=True)
        # Compute XMCD
        xmcd['edge'] = xmcd['CL'] - xmcd['CR']

        if self.pre_edge:
            xmcd['pre-edge'] = xmcd['CLpe'] - xmcd['CRpe']

            crcrp = xmcd['CR'] / xmcd['CRpe']
            clclp = xmcd['CL'] / xmcd['CLpe']

            xmcd['no-norm'] = (clclp - crcrp)
            xmcd['norm'] = 200 * xmcd['no-norm'] / (crcrp + clclp - 2)

            if guiobj.analysis == 'hyst_t_aver':
                # For time average compute also error coming from
                # averaging
                xmcd['Dedge'] = xmcd['dCL'] + xmcd['dCR']
                xmcd['Dpre-edge'] = xmcd['dCLpe'] + xmcd['dCRpe']

                dcrcrp = np.abs(
                    (np.abs(xmcd['dCR'] / xmcd['CR']) + np.abs(xmcd['dCRpe']
                    / xmcd['CRpe'])) * crcrp)
                dclclp = np.abs(
                    (np.abs(xmcd['dCL'] / xmcd['CL']) + np.abs(xmcd['dCLpe']
                    / xmcd['CLpe'])) * clclp)
                denerr = (clclp + crcrp - 2) ** 2
                xmcd['Dnorm'] = 400 * (((np.abs(crcrp-1) / denerr) * dclclp)
                                       + ((np.abs(clclp-1) / denerr) * dcrcrp))

        self.xmcd = xmcd


def h_scale(guiobj, confobj, pos, neg, log_dt):
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

    confobj : configuration object.

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
        If splitted branches are provided return
        array[array_p, arrany_n, array_all]
        where array_p is the array with common scale for positive
        magnetic field branch, array_n contains the common scale
        for negative magnetic field branch and array_all the union of
        the two.

    Add h_scale key to log_dt with min, max and number of points of
    field scale. 
    '''
    h_maxs = []  # Maxima of magnetic field scans
    h_mins = []  # Minima of magnetic field scans
    h_len = []  # Number of different fields in each scan

    # splitted branch case
    if confobj.spl_brnch:
        zero_ps = []  # Near zero values for positive fields
        zero_ns = []  # Near zero values for negative fields

        # Collect maxima and minima and count the fields for edge scans
        # For splitted branches data are alreafy separated up from down
        for i in pos.up_idx:
            # Check from max h value which half-branch is
            maxh = np.amax(pos.up['H' + i])
            # if maxh < 0 negative field half-branch
            if maxh < 0:
                zero_ns.append(maxh)
                h_mins.append(np.amin(pos.up['H' + i]))
            # if maxh > 0 positive field half-branch
            else:
                h_maxs.append(maxh)
                zero_ps.append(np.amin(pos.up['H' + i]))
            h_len.append(h_num_points(pos.up['H' + i].dropna()))
        for i in pos.dw_idx:
            # Check from max h value which half-branch is
            maxh = np.amax(pos.down['H' + i])
            # if maxh < 0 negative field half-branch
            if maxh < 0:
                zero_ns.append(maxh)
                h_mins.append(np.amin(pos.down['H' + i]))
            # if maxh > 0 positive field half-branch
            else:
                h_maxs.append(maxh)
                zero_ps.append(np.amin(pos.down['H' + i]))
            h_len.append(h_num_points(pos.down['H' + i].dropna()))
        for i in neg.up_idx:
            # Check from max h value which half-branch is
            maxh = np.amax(neg.up['H' + i])
            # if maxh < 0 negative field half-branch
            if maxh < 0:
                zero_ns.append(maxh)
                h_mins.append(np.amin(neg.up['H' + i]))
            # if maxh > 0 positive field half-branch
            else:
                h_maxs.append(maxh)
                zero_ps.append(np.amin(neg.up['H' + i]))
            h_len.append(h_num_points(neg.up['H' + i].dropna()))
        for i in neg.dw_idx:
            # Check from max h value which half-branch is
            maxh = np.amax(neg.down['H' + i])
            # if maxh < 0 negative field half-branch
            if maxh < 0:
                zero_ns.append(maxh)
                h_mins.append(np.amin(neg.down['H' + i]))
            # if maxh > 0 positive field half-branch
            else:
                h_maxs.append(maxh)
                zero_ps.append(np.amin(neg.down['H' + i]))
            h_len.append(h_num_points(neg.down['H' + i].dropna()))
        # If present pre-edge scans do the same for them
        if pos.pre_edge:
            for i in pos.pe_up_idx:
                # Check from max h value which half-branch is
                maxh = np.amax(pos.pe_up['H' + i])
                # if maxh < 0 negative field half-branch
                if maxh < 0:
                    zero_ns.append(maxh)
                    h_mins.append(np.amin(pos.pe_up['H' + i]))
                # if maxh > 0 positive field half-branch
                else:
                    h_maxs.append(maxh)
                    zero_ps.append(np.amin(pos.pe_up['H' + i]))
                h_len.append(h_num_points(pos.pe_up['H' + i].dropna()))
            for i in pos.pe_dw_idx:
                # Check from max h value which half-branch is
                maxh = np.amax(pos.pe_down['H' + i])
                # if maxh < 0 negative field half-branch
                if maxh < 0:
                    zero_ns.append(maxh)
                    h_mins.append(np.amin(pos.pe_down['H' + i]))
                # if maxh > 0 positive field half-branch
                else:
                    h_maxs.append(maxh)
                    zero_ps.append(np.amin(pos.pe_down['H' + i]))
                h_len.append(h_num_points(pos.pe_down['H' + i].dropna()))
        if neg.pre_edge:
            for i in neg.pe_up_idx:
                # Check from max h value which half-branch is
                maxh = np.amax(neg.pe_up['H' + i])
                # if maxh < 0 negative field half-branch
                if maxh < 0:
                    zero_ns.append(maxh)
                    h_mins.append(np.amin(neg.pe_up['H' + i]))
                # if maxh > 0 positive field half-branch
                else:
                    h_maxs.append(maxh)
                    zero_ps.append(np.amin(neg.pe_up['H' + i]))
                h_len.append(h_num_points(neg.pe_up['H' + i].dropna()))
            for i in neg.pe_dw_idx:
                # Check from max h value which half-branch is
                maxh = np.amax(neg.pe_down['H' + i])
                # if maxh < 0 negative field half-branch
                if maxh < 0:
                    zero_ns.append(maxh)
                    h_mins.append(np.amin(neg.pe_down['H' + i]))
                # if maxh > 0 positive field half-branch
                else:
                    h_maxs.append(maxh)
                    zero_ps.append(np.amin(neg.pe_down['H' + i]))
                h_len.append(h_num_points(neg.pe_down['H' + i].dropna()))

        # Compute min, max, zero_p, zero_n and default length of energy
        # range
        # Set decimal place to round the first OoM higher than tolerance
        # in h_num_points
        h_min = np.around(np.amax(h_mins), 3)
        h_max = np.around(np.amin(h_maxs), 3)
        zero_p = np.around(np.amax(zero_ps), 3)
        zero_n = np.around(np.amin(zero_ns), 3)
        h_av_len = np.around(np.average(h_len), 0)

        # Set number of points of energy scale
        if guiobj.interactive:
            n_points = int(np.around((guiobj.num_pnts(2*h_av_len))/2))
        else:
            n_points = int(h_av_len)

        h_step = (h_max - h_min) / (2*n_points - 1)
        log_dt['h_scale'] = [h_min, h_max, n_points, h_step]

        p_branch = np.linspace(zero_p, h_max, n_points)
        n_branch = np.linspace(h_min, zero_n, n_points)
        all_branch = np.concatenate((n_branch, p_branch))
        return [n_branch, p_branch, all_branch]
    # no splitted branch case
    else:
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

        # Compute min, max and default length of energy range
        # Set decimal place to round the first OoM higher than tolerance
        # in h_num_points
        h_min = np.around(np.amax(h_mins), 3)
        h_max = np.around(np.amin(h_maxs), 3)
        h_av_len = np.around(np.average(h_len), 0)

        # Set number of points of energy scale
        if guiobj.interactive:
            n_points = guiobj.num_pnts(h_av_len)
        else:
            n_points = int(h_av_len)

        h_step = (h_max - h_min) / (n_points - 1)
        log_dt['h_scale'] = [h_min, h_max, n_points, h_step]

        return np.linspace(h_min, h_max, n_points)


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
    value for the number of different magnetic fields measured during
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
    for i in range(len(indexes)):
        if i == len(indexes) - 1:
            y.append(np.average(sorted_data[idx][indexes[i]:]))
        else:
            y.append(np.average(sorted_data[idx][indexes[i]:indexes[i+1]]))

    return x, y
