"""
pyDichroX_IO.py

Methods to manage input and output files for XMCD, XNCD, XNLD and
hysteresys data analysis.

Methods
--------
open_import_scan(guiobj, confobj)
    Open input files and import data for energy scan experiments.

e_scan_importer(guiobj, confobj, pos, neg, T, H ,log_dt)
    Manage energy scan import.

h_scan_importer(guiobj, confobj, pos, neg, T, H, log_dt)
    Manage magnetic field on-fly scan import.

dt_raw_import(confobj, data)
    Based on configuration file select data to import depending on
    sensing.

escan_split(confobj, data)
    Split energy scan containing multiple scans into single scans.

hscan_imp(confobj, data)
    Select data columns from imported data for hysteresis analysis.

set_scn_num(confobj, f_name, pos, neg)
    Associate an identifier for each scan.

separate_escans(guiobj, confobj, e_raw, dt_raw, scn_lbl, scn_num, ispol,
    lgrws, pos, neg)
    Separate scans and fill positive and negative ScanData objects with
    raw data and labels depending on analysis type.

separate_hscans(guiobj, confobj, h_raw, dt_raw, time_raw, scn_lbl,
    scn_num, ispol, lgrws, pos, neg)
    Separate magnetic field scans and fill positive and negative
    ScanData objects with raw data and labels depending on polarisation.

output_fls_escan(guiobj, pos, neg, scanobj)
    Organize output data and columns names for energy scans data
    analysis.

output_fls_hscan(scanobj)
    Organize output data and columns names for hysteresis on fly
    analysis.

output_fls_ptbypt(guiobj, scanobj)
    Organize output data and columns names for hysteresis on point by
    point analysis.

output_plot_escan(guiobj, pos, neg, scanobj, log_dt)
    Create a final figure with plot reporting averaged positive and
    negative XAS and percentage X-Ray dichroism.

output_plot_hscan(scanobj, log_dt)
    Create a final figure with plot reporting XMCD hysteresis on the fly
    analized.

output_plot_t_aver(scanobj, log_dt)
    Create a final figure with plot reporting XMCD hysteresis point by
    point analized with time averaging.

output_plot_t_split(guiobj, scanobj, log_dt)
    Create a final figure with plot reporting XMCD hysteresis point by
    point analized with time splitting.

save_data_escan(confobj, guiobj, pos, neg, scanobj, log_dt, pos_ref,
    neg_ref, scanobj_ref_norm, log_dt_ref)
    Save output data, logdata and final plots for energy scan analysis.

save_data_hscan(confobj, guiobj, scanobj, log_dt)
    Save output data, logdata and final plots for scan field on the fly
    analysis.

save_data_tscan(confobj, guiobj, scanobj, log_dt)
    Save output data, logdata and final plots for scan field point by
    point analysis.
"""

# Copyright (C) Giuseppe Cucinotta.
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

import sys
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import scipy.interpolate as itp
from scipy import stats
from matplotlib.widgets import Button

import modules.pyDichroX_escan_datatreat as esdt
import modules.pyDichroX_hscan_datatreat as hsdt


def open_import_scan(guiobj, confobj, edg_filename):
    '''
    Open input files import data and groups them:
    - for energy scan analysis basing on absorption coefficient
      for circularly polirized light on value of linear polarization for
      linearly polarized light
    - for field scan analysis basing on light polarization.

    Collect experimental data for compiling output logfile.

    For energy scan analysis if confobj.ref_norm data scans of a
    reference sample are prompted. 

    Parameters
    ----------
    guiobj : GUI object
        Provides GUI dialogs.

    confobj : Configurations object.

    edg_filename : str
            File lisitng edge and pre-edge energies.

    Returns
    -------
    ScanData object with positive scans (sigma+, CR, H- or LH scans,
        depending on analysis)

    ScanData object with negative scans (sigma-, CL, H+ or LV scans,
        depending on analysis)

    dict, contains log information:
    . angle : sample user setted angle
    . bm_angle : beam incidence angle
    . temp : mean of temperatures collected in log_tbl['t']
    . field : mean of field collected in log_tbl['field']
    . pos_chs : list (str) positive scans chosen for analysis - for
                energy scans
    . neg_chs : list (str) negative scans chosen for analysis - for
                energy scans
    . pos_up_chsn : list (str) with CR branch up chosen scans
    . pos_dw_chsn : list (str) with CR branch down chosen scans
    . pos_pe_up_chsn : list (str) with CR pre-edge branch up chosen
                       scans
    . pos_pe_dw_chsn : list (str) with CR pre-edge branch down chosen
                       scans
    . neg_up_chsn : list (str) with CL branch up chosen scans
    . neg_dw_chsn : list (str) with CL branch down chosen scans
    . neg_pe_up_chsn : list (str) with CL pre-edge branch up chosen
                       scans
    . neg_pe_dw_chsn : list (str) with CL pre-edge branch down chosen
                       scans    
    . Edge_name : name of the edge used for analysis
    . Edge_en : edge energy value tabulated in the edge list file
    . exper_edge : for energy scan analysis -> energy value used as edge
                   energy for analysis
                   for field scan analysis -> mean of edge energies in
                   log_tbl['edge_mon_en']
    . PreEdge_en : pre-edge energy value tabulated in the edge list file
    . setted_pedg : for energy scan analysis -> energy value chosen as
                    pre-edge energy
                    for field scan analysis -> mean of pre-edge energies
                    in log_tbl['pre_edge_mon_en']
    . PostEdge_en : post-edge energy value tabulated in the edge list
                    file
    . setted_postedg : energy value chooed as post-edge energy
    . recal : bool, True if energy scale has been recalibrated
    . offset : offset value added to energy scale to recalibrate it
    . e_scale : list containing start value, end value, number of points
                and stepsize used to construct the energy scale
    . h_scale : list containing start value, end value, number of points
                and stepsize used to construct the magnetic field scale
    . t_scale : list containing start value, end value, number of points
                and step size used to construct the time scale
    . pe_int : interpolated pre-edge value
    . pos_ej : edge-jump for positive scans
    . pos_ej_er : error on edge-jump for positive scans
    . pos_ej_int : edge-jump for positve scans computed using
                   interpolated pre-edge energy
    . pos_ej_int_er : error on edge-jump for positve scans computed
                      using interpolated pre-edge energy
    . neg_ej : edge-jump for negative scans
    . neg_ej_er : error on edge-jump for negative scans
    . neg_ej_int : edge-jump for negative scans computed using
                   interpolated pre-edge energy
    . neg_ej_int_er : error onedge-jump for negative scans
                      computed using interpolated pre-edge energy
    . xas_aver_ej : edge-jump computed from averaged XAS spectra
    . xas_aver_ej_er : error on xas_aver_ej
    . xas_aver_ej_int : edge-jump computed from averaged XAS spectra and
                        using interpolated pre-edge energy
    . xas_aver_ej_int_er : error on xas_aver_ej_int
    . edge_pc : percentage of X-Ray Dichroism normalized for edge-jump
                at experimental edge energy
    . edge_pc_er : error on edge_pc
    . edge_pc_int : percentage of X-Ray Dichroism normalized with
                    inerpolated edge-jump at experimental edge energy
    . edge_pc_int_er : error on edge_pc_int 
    . log_tbl / log_ref_tbl pandas DataFrame:
        - mon_en : monocromator energy
        - edge_mon_en : monocromator energy for edge scans - only field
                        scan analysis
        - pre_edge_mon_en : monocromator energy for pre-edge scans -
                        only field scan analysis
        - pol : polarisation id
        - field : magnetic field value
        - tb1 : sample temperature 1 (Deimos)
        - tb2 : sample temperature 2 (Deimos)
        - t : tb1, tb2 average (Deimos)
              temperature inserted by user (Ape)
              sample temperature (Boreas)
        - rz : sample rotation angle (Deimos)
        - tx : sample x position (Deimos)
        - tz : sample z position (Deimos)
        - x : sample x position (Ape)
        - y : sample y position (Ape)
        - z : sample z position (Ape) 
        - scan_num : list with scan numbers of input data
        - type : pol id for each scan

    If reference data are considered in confobj also ScanData objects
    with repsectively positive and negative scan data belonging to the
    reference sample are returned.
    '''
    # Set edge and pre-edge energies
    edge = guiobj.chs_edge(edg_filename)

    # dictionary collecting log data
    log_dt = {}

    log_dt['Edge_name'] = edge[0]
    log_dt['Edge_en'] = float(edge[1])
    log_dt['PreEdge_en'] = float(edge[2])
    log_dt['PostEdge_en'] = float(edge[3])

    # Insert experimental rotation angle of the sample
    angle, bm_angle = guiobj.ask_angle()

    # Ask for sample temperature if not provided by datalog file
    T = -100  # just initialize, if provided by log they won't be used
    H = -100
    if confobj.ask_for_T:
        T = guiobj.ask_T()
    # Ask for magnetic field if not provided by datalog file
    if confobj.ask_for_H:
        H = guiobj.ask_H()

    # Attribute to select message for GUI
    guiobj.infile_ref = False

    # Instantiation of objects for different polarization signs
    # pos will contain sigma+, CR, LH, H- polarization data
    # neg will contain sigma-, CL, LV, H+ polarization data
    # ref are referred to reference data. They will be filled only if
    # confobj.ref_norm is True
    #
    # If field scan data use hsdt data object, if energy scan data use
    # es.dt data object
    if guiobj.analysis in guiobj.type['hyst']:
        pos = hsdt.ScanData(guiobj, confobj)
        neg = hsdt.ScanData(guiobj, confobj)

        log_tbl = h_scan_importer(guiobj, confobj, pos, neg, T, H, log_dt)
        log_dt['exper_edge'] = np.round(log_tbl['edge_mon_en'].abs().mean(), 2)
        log_dt['setted_pedg'] = np.round(
            log_tbl['pre_edge_mon_en'].abs().mean(), 2)
        if confobj.spl_brnch:
            log_dt['field'] = np.round(log_tbl['field'].abs().max(), 1)
        else:
            log_dt['field'] = np.round(log_tbl['field'].abs().mean(), 1)

        # Just create empty data list for pos_ref and neg_ref currently
        # not considered for hysteresis
        pos_ref = []
        neg_ref = []
    else:
        pos = esdt.ScanData()
        neg = esdt.ScanData()

        log_tbl = e_scan_importer(guiobj, confobj, pos, neg, T, H, log_dt)

        log_dt['field'] = np.round(log_tbl['field'].abs().mean(), 1)

        pos_ref = esdt.ScanData()
        neg_ref = esdt.ScanData()

    # dictionary collecting log data
    log_dt['log_tbl'] = log_tbl
    log_dt['temp'] = np.round(log_tbl['t'].mean(), 1)
    log_dt['temp_er'] = np.round(log_tbl['t'].std(), 1)
    log_dt['angle'] = angle
    log_dt['bm_angle'] = bm_angle

    if (guiobj.analysis not in guiobj.type['hyst']) and confobj.ref_norm:
        # Currently no normalization by other data set is considered for
        # hysteresis analysis

        # Attribute to select message for GUI
        guiobj.infile_ref = True

        log_ref_tbl = e_scan_importer(guiobj, confobj, pos_ref, neg_ref, T, H,
                                      log_dt)
        log_dt['log_ref_tbl'] = log_ref_tbl

        guiobj.infile_ref = False

    return pos, neg, log_dt, pos_ref, neg_ref


def e_scan_importer(guiobj, confobj, pos, neg, T, H, log_dt):
    '''
    Manage energy scan import.

    Open input files, import data and group them based on absorption
    coefficient for circularly polirized light or based on value of
    linear polarization for linearly polarized light.

    Collect experimental data for compiling output lile.

    Parameters
    ----------
    guiobj : GUI object
        Provides GUI dialogs.

    confobj : Configurations object.

    pos : ScanData object
        with positive scans (sigma+, CR, H- or LH scans, depending on
        analysis).

    neg : ScanData object 
        with negative scans (sigma-, CL, H+ or LV scans, depending on
        analysis).

    T : float
        temperature value manually inserted if not provided by data
        logfiles.

    H : float
        field value manually inserted if not provided by data logfiles.

    log_dt : dict
        dictionary collecting log data.

    Returns
    -------
    Update and fill attributes of pos and neg ScanData objects.

    Update log_dt with log information (depending on beamline).
    '''
    while True:
        # log table with scan log details
        lgrws = {}
        log_tbl = pd.DataFrame()

        in_sets = guiobj.in_dtst(confobj)  # Import data sets

        for dtset in in_sets:
            chk_sgn_h = 1.  # To check a field change

            for file in dtset:
                f_name = os.path.basename(file)
                data = pd.read_csv(file, sep=confobj.sep,
                                   usecols=confobj.scn_cols(guiobj, f_name))
                try:
                    lgrws.update(confobj.log_scavenger(file))
                except:
                    # If logfile is not present analysis will fail so
                    # continue and skip data collecting for this scan
                    logfn = confobj.single_lognm(f_name)
                    guiobj.no_log(logfn, confobj.nologmess)
                    continue

                if confobj.ask_for_T:
                    lgrws['t'] = T
                if confobj.ask_for_H:
                    lgrws['field'] = H

                # Light polarisation
                pol = lgrws['pol']

                # Mean and sign of magnetic field
                h_sgn = np.sign(lgrws['field'])

                e_raw, dt_raw = escan_split(confobj, data)

                for i in range(len(e_raw)):
                    # Set scan number from filename
                    scn_num = set_scn_num(confobj, f_name, pos, neg)

                    # Check if magnetic field has changed: the first
                    # scan and the scans after the field has changed are
                    # usually rejected so a different label is provided
                    # for these scans.
                    if (dtset.index(file) == 0) and (i == 0):
                        scn_lbl = scn_num + ' Dummy scan'
                        lgrws['scn_num'] = scn_num + '(D)'
                        chk_sgn_h = h_sgn
                    elif (chk_sgn_h != h_sgn):
                        scn_lbl = scn_num + ' Dummy scan'
                        lgrws['scn_num'] = scn_num + '(D)'
                        chk_sgn_h = h_sgn
                    else:
                        scn_lbl = scn_num
                        lgrws['scn_num'] = scn_num

                    # Separate positive from negative scans
                    if guiobj.analysis in guiobj.type['xnld']:
                        try:
                            islv = confobj.lv_cond(pol)
                        except:
                            guiobj.wrongpol(scn_num, 'linear')
                            # continue if wrong file is found
                            continue

                        separate_escans(
                            guiobj, confobj, e_raw[i], dt_raw[i], scn_lbl,
                            scn_num, islv, lgrws, pos, neg)
                    else:
                        # if not XNLD check for circular polarisation
                        try:
                            iscr = confobj.cr_cond(pol)
                        except:
                            guiobj.wrongpol(scn_num, 'circular')
                            # continue if wrong file is found
                            continue

                        separate_escans(
                            guiobj, confobj, e_raw[i], dt_raw[i], scn_lbl,
                            scn_num, iscr, lgrws, pos, neg)

                    log_tbl = pd.concat(
                        [log_tbl, pd.DataFrame(lgrws, index=[0])],
                        ignore_index=True)

            # increase counter for cumulative scanlogs if present
            confobj.scanlog_cnt += 1
        # At least one scan file per polarization is needed
        if len(pos.idx) < 1:
            cont = guiobj.not_enough_fls(True)
            if (not cont) or (cont is None):
                sys.exit(0)
            else:
                continue
        if len(neg.idx) < 1:
            cont = guiobj.not_enough_fls(False)
            if (not cont) or (cont is None):
                sys.exit(0)
            else:
                continue

        break  # If no problem with number of files breaks the loop

    return log_tbl


def h_scan_importer(guiobj, confobj, pos, neg, T, H, log_dt):
    '''
    Manage magnetic field on-fly scan import.

    Open input files, import data, group them based on light
    polirization and separate edge from pre-edge scans.

    Collect experimental data for compiling output logfile.

    Parameters
    ----------
    guiobj : GUI object
        Provides GUI dialogs.

    confobj : Configurations object.

    pos : ScanData object
        with positive scans (sigma+, CR, H- or LH scans, depending on
        analysis).

    neg : ScanData object 
        with negative scans (sigma-, CL, H+ or LV scans, depending on
        analysis).

    T : float
        temperature value manually inserted if not provided by data
        logfiles.

    H : float
        field value manually inserted if not provided by data logfiles.

    log_dt : dict
        dictionary collecting log data.

    Returns
    -------
    Update and fill attributes of pos and neg ScanData objects.

    Update log_dt with log information (depending on beamline).
    '''
    while True:
        # log table with scan log details
        lgrws = {}
        log_tbl = pd.DataFrame()

        in_sets = guiobj.in_dtst(confobj)  # Import data sets

        for dtset in in_sets:

            for file in dtset:
                f_name = os.path.basename(file)
                data = pd.read_csv(file, sep=confobj.sep, na_values='nan',
                                   usecols=confobj.scn_cols(guiobj, f_name))
                try:
                    lgrws.update(confobj.log_scavenger(file))
                except:
                    logfn = confobj.single_lognm(f_name)
                    guiobj.no_log(logfn, confobj.nologmess)
                    continue

                if confobj.ask_for_T:
                    lgrws['t'] = T
                if confobj.ask_for_H:
                    lgrws['field'] = H

                # Light polarisation
                pol = lgrws['pol']

                # Import magnetic field scan
                time_raw, h_raw, dt_raw = hscan_imp(confobj, data)

                # Separate positive from negative scans
                try:
                    iscr = confobj.cr_cond(pol)
                except:
                    guiobj.wrongpol(scn_num, 'circular')
                    continue  # continue if wrong file is found

                # Set scan number from filename
                scn_num = set_scn_num(confobj, f_name, pos, neg)
                lgrws['scn_num'] = scn_num

                separate_hscans(guiobj, confobj, h_raw, dt_raw, time_raw,
                                scn_num, iscr, lgrws, log_dt, pos, neg)

                log_tbl = pd.concat([log_tbl, pd.DataFrame(lgrws, index=[0])],
                                    ignore_index=True)

            # increase counter for cumulative scanlogs if present
            confobj.scanlog_cnt += 1

        # if splitted braches check that at least on half brah per
        # polarisation is present
        if confobj.spl_brnch:
            pck = pos.up_hp and pos.up_hn and pos.down_hp and pos.down_nh
            if pck:
                break
            else:
                cont = guiobj.not_enough_sp_br(pos, True)
                if (not cont) or (cont is None):
                    sys.exit(0)
                else:
                    continue
            nck = neg.up_hp and neg.up_hn and neg.down_hp and neg.down_nh
            if nck:
                break
            else:
                cont = guiobj.not_enough_sp_br(neg, False)
                if (not cont) or (cont is None):
                    sys.exit(0)
                else:
                    continue
            # If no pre-edge files are present hysteresis analysis can
            # be done without pre-edge normalization.
            # But if pre-edge files are present they must be at least
            # one for each half branch.
            if (len(pos.pe_idx) > 0) and (len(pos.pe_idx) < 4):
                cont = guiobj.not_enough_fls(True, True)
                if (not cont) or (cont is None):
                    sys.exit(0)
                else:
                    continue
            elif len(pos.pe_idx) >= 4:
                p_peck = (pos.pe_up_hp and pos.pe_up_hn and
                          pos.pe_down_hp and pos.pe_down_nh)
                if p_peck:
                    break
                else:
                    cont = guiobj.not_enough_sp_br(pos, True, True)
                    if (not cont) or (cont is None):
                        sys.exit(0)
                    else:
                        continue
            if (len(neg.pe_idx) > 0) and (len(neg.pe_idx) < 4):
                cont = guiobj.not_enough_fls(False, True)
                if (not cont) or (cont is None):
                    sys.exit(0)
                else:
                    continue
            elif len(neg.pe_idx) >= 4:
                n_peck = (neg.pe_up_hp and neg.pe_up_hn and
                          neg.pe_down_hp and neg.pe_down_nh)
                if n_peck:
                    break
                else:
                    cont = guiobj.not_enough_sp_br(neg, False, True)
                    if (not cont) or (cont is None):
                        sys.exit(0)
                    else:
                        continue
        else:
            # At least one scan file per polarization is needed
            if len(pos.idx) < 1:
                cont = guiobj.not_enough_fls(True)
                if (not cont) or (cont is None):
                    sys.exit(0)
                else:
                    continue
            if len(neg.idx) < 1:
                cont = guiobj.not_enough_fls(False)
                if (not cont) or (cont is None):
                    sys.exit(0)
                else:
                    continue
            # If no pre-edge files are present hysteresis analysis can
            # be done without pre-edge normalization.
            # But if are present pre-edge files for one polarization
            # also at least one pre-edge file for the other polarization
            # is required.
            if (len(pos.pe_idx) < 1) and (len(neg.pe_idx) >= 1):
                cont = guiobj.not_enough_fls(True, True)
                if (not cont) or (cont is None):
                    sys.exit(0)
                else:
                    continue
            if (len(neg.pe_idx) < 1) and (len(pos.pe_idx) >= 1):
                cont = guiobj.not_enough_fls(False, True)
                if (not cont) or (cont is None):
                    sys.exit(0)
                else:
                    continue

        if guiobj.analysis != 'hyst_fly':
            # For point by point analysis if pre-edge scans are present
            # they must be in the same number as edge scans
            if (len(pos.pe_idx) > 0) and (len(pos.pe_idx) != len(pos.idx)):
                cont = guiobj.not_enough_fls(True, True)
                if (not cont) or (cont is None):
                    sys.exit(0)
                else:
                    continue
            if (len(neg.pe_idx) > 0) and (len(neg.pe_idx) != len(neg.idx)):
                cont = guiobj.not_enough_fls(False, True)
                if (not cont) or (cont is None):
                    sys.exit(0)
                else:
                    continue

        break  # If no problem with number of files breaks the loop

    if (confobj.spl_brnch == False) and (guiobj.analysis == 'hyst_fly'):
        # Up and down separation only for complete branches
        # For hysteresis point by point no branch recognition mechanism
        # is at the moment present.
        pos.up_n_down()
        neg.up_n_down()

    return log_tbl


def dt_raw_import(confobj, data):
    '''
    Based on configuration file select data to import depending on
    sensing.
    It also normalize data by i0 if normalization is not provided in
    datafile.

    Parameters
    ----------
    guiobj : GUI object
        Provides GUI dialogs.

    confobj : Configuration obj.

    data : Pandas DataFrame
        DataFrame containing imported data from input file.

    Return
    ------
    array with normalized data based on sensing reported in confobj.
    '''
    if confobj.norm_curr:  # Normalize data
        if (confobj.sense == 'TEY'):
            dt_raw = data[confobj.it]/data[confobj.i0]
        else:
            dt_raw = data[confobj.if1]/data[confobj.if0]
    else:  # Normalized data are provided
        if confobj.sense == 'TEY':
            dt_raw = data[confobj.iti0]
        else:
            dt_raw = data[confobj.ifi0]

    return dt_raw


def escan_split(confobj, data):
    '''
    Split energy scan containing multiple scans into single scans.

    Energy data trend in this case is supposed to have a sawtooth shape.
    To find when a scan ends and the next starts escan_split:
    . interpolate energy profile
    . computes the derivative of interpolation and its statistical mode
    . the points where the derivative is different from mode are the
      end-scan points.

    If input data contain just one scan, return the single scan.

    Parameters
    ----------
    confobj : Configuration Object.

    data : Pandas DataFrame
        DataFrame containing imported data from input file.

    Returns
    -------
    Set of two list. The first list contains pd.Series with energy
    vlaues of splitted scans. The second lisc contains pd.Series with
    data values of splitted scans.
    If only one scan is present the two list will contain just one
    element each.
    '''
    # Import energy and XAS normalized raw data
    energy_raw = abs(data[confobj.energy])
    data_raw = dt_raw_import(confobj, data)

    # indexes of e_raw, used as x for energy derivative
    x = np.array([i for i in range(len(energy_raw))])
    # e_raw linear interpolation
    e_raw_int = itp.UnivariateSpline(x, energy_raw, k=1, s=0)
    # derivative object of e_raw_int
    e_raw_der = e_raw_int.derivative()
    # compute first derivative on x and round it
    e_raw_der_round = np.around(e_raw_der(x), 0)
    # compute statistical mode
    der_mode = stats.mode(e_raw_der_round, keepdims=True).mode

    st_idx = 0  # Energy scan start index

    # List with separated scans
    e_raw = []
    dt_raw = []

    for i in range(len(e_raw_der_round)):
        if e_raw_der_round[i] != der_mode:  # Scan ends here
            stp_idx = i + 1
            e_raw.append(energy_raw.iloc[st_idx: stp_idx])
            dt_raw.append(data_raw.iloc[st_idx: stp_idx])

            st_idx = stp_idx
    # Append last scan, or the single scan in case no multiple scans are
    # present
    e_raw.append(energy_raw.iloc[st_idx:])
    dt_raw.append(data_raw.iloc[st_idx:])

    return e_raw, dt_raw


def hscan_imp(confobj, data):
    '''
    Select data columns from imported data for hysteresis analysis.

    Parameters
    ----------
    confobj : Configuration Object

    data : Pandas DataFrame
        DataFrame containing imported data from input file

    Returns
    -------
    Returns a tuple of pd.Series of three elements:
    - 1st. timestamps (needed for point by point hysteresis)
    - 2nd. field vlaues of the scan
    - 3rd. measured data values.
    '''
    # Import field and XAS normalized raw data
    h_raw = data[confobj.field]
    data_raw = dt_raw_import(confobj, data)
    time_raw = data[confobj.time]

    return time_raw, h_raw, data_raw


def set_scn_num(confobj, f_name, pos, neg):
    '''
    Associate an identifier for each scan.
    This identifier is usually the scan number associated to datafile.
    If the input consists of more than one dataset and different scans
    coming from different set have the same number identifier is
    modified in order to avoid overwiritng data problems during data
    import.

    Parameters
    ----------
    confobj : configuration obj.

    f_name : str
        scan name, used to exctract scan identifier.

    pos : ScanData obj
        contains positive scan data.

    neg : ScanData obj
        contains negative scan data.

    Returns
    -------
    str, identifier for the scan.
    '''
    scn_num = confobj.extract_num(f_name)
    scn_num_chk = scn_num

    suffx = 1  # suffix to be added if scan_id already used is found

    # if scan id is already used add suffix and check until no free id
    # is found
    while True:
        if scn_num_chk in pos.raw_imp.columns:
            scn_num_chk = scn_num + '_{}'.format(suffx)
            suffx += 1
        else:
            break

    while True:
        if scn_num_chk in neg.raw_imp.columns:
            scn_num_chk = scn_num + '_{}'.format(suffx)
            suffx += 1
        else:
            break

    return scn_num_chk


def separate_escans(guiobj, confobj, e_raw, dt_raw, scn_lbl, scn_num, ispol,
                    lgrws, pos, neg):
    '''
    Separate energy scans and fill positive and negative ScanData
    objects with raw data and labels depending on analysis type.

    Parameters
    ----------
    confobj : configuration obj.

    guiobj : GUI object
        Provides GUI dialogs.

    e_raw : Pandas Series
        Series containing energy values of the scan.

    dt_raw : Pandas Series
        Series containing imported data.

    scn_lbl : str
        Scan label.

    scn_num : str
        Scan number.

    ispol : bool
        True for CR and LV polarisations
        False for CL and LH polarisations.

    lgrws : dict
        Dictionary with log data collected from data logfile.

    pos, neg : ScanData objects
        Collect raw data from positive and negative scans.

    Return
    ------
    Fill pos and neg ScanData objects with raw data as well lgrws
    dictionary with scan labels.

    '''

    if e_raw.iloc[0] > e_raw.iloc[-1]:
        e_raw_to_imp = np.flip(e_raw.to_numpy())
        dt_raw_to_imp = np.flip(dt_raw.to_numpy())
    else:
        e_raw_to_imp = e_raw.to_numpy()
        dt_raw_to_imp = dt_raw.to_numpy()

    rawen = pd.DataFrame()
    rawdt = pd.DataFrame()

    rawen['E' + scn_num] = e_raw_to_imp
    rawdt[scn_num] = dt_raw_to_imp

    # For XNLD, data are divided and different labels assigned based on
    # polarization sign.
    if guiobj.analysis in guiobj.type['xnld']:

        scn_lbl += ' H = {:.2f} T'.format(lgrws['field'])

        if ispol:  # LV data
            neg.raw_imp = pd.concat([neg.raw_imp, rawen], axis=1)
            neg.raw_imp = pd.concat([neg.raw_imp, rawdt], axis=1)
            neg.label.append(scn_lbl)
            neg.idx.append(scn_num)
            neg.dtype = 'LV'
            lgrws['type'] = 'LV'
        else:  # LH data
            pos.raw_imp = pd.concat([pos.raw_imp, rawen], axis=1)
            pos.raw_imp = pd.concat([pos.raw_imp, rawdt], axis=1)
            pos.label.append(scn_lbl)
            pos.idx.append(scn_num)
            pos.dtype = 'LH'
            lgrws['type'] = 'LH'
    # For XMCD, data are divided and different labels assigned based on
    # sigma sign.
    elif guiobj.analysis in guiobj.type['xmcd']:
        h_sgn = np.sign(lgrws['field'])

        if ispol:
            scn_lbl += ' CR, H = {:.2f} T'.format(lgrws['field'])
            lgrws['type'] = 'CR'
        else:
            scn_lbl += ' CL, H = {:.2f} T'.format(lgrws['field'])
            lgrws['type'] = 'CL'

        sigma_sgn = h_sgn * confobj.phi_sgn

        if sigma_sgn > 0:  # sigma + data
            pos.raw_imp = pd.concat([pos.raw_imp, rawen], axis=1)
            pos.raw_imp = pd.concat([pos.raw_imp, rawdt], axis=1)
            pos.label.append(scn_lbl)
            pos.idx.append(scn_num)
            pos.dtype = 'sigma^+'
        else:  # sigma - data
            neg.raw_imp = pd.concat([neg.raw_imp, rawen], axis=1)
            neg.raw_imp = pd.concat([neg.raw_imp, rawdt], axis=1)
            neg.label.append(scn_lbl)
            neg.idx.append(scn_num)
            neg.dtype = 'sigma^-'
    # For XNCD, data are divided and different labels assigned based on
    # polarizaion sign.
    elif guiobj.analysis in guiobj.type['xncd']:
        scn_lbl += ' H = {:.2f} T'.format(lgrws['field'])

        if ispol:  # CR data
            pos.raw_imp = pd.concat([pos.raw_imp, rawen], axis=1)
            pos.raw_imp = pd.concat([pos.raw_imp, rawdt], axis=1)
            pos.label.append(scn_lbl)
            pos.idx.append(scn_num)
            pos.dtype = 'CR'
            lgrws['type'] = 'CR'
        else:  # CL data
            neg.raw_imp = pd.concat([neg.raw_imp, rawen], axis=1)
            neg.raw_imp = pd.concat([neg.raw_imp, rawdt], axis=1)
            neg.label.append(scn_lbl)
            neg.idx.append(scn_num)
            neg.dtype = 'CL'
            lgrws['type'] = 'CL'
    # For XNXD, data are divided and different labels assigned
    # based on magnetic field sign.
    elif guiobj.analysis in guiobj.type['xnxd']:
        h_sgn = np.sign(lgrws['field'])

        if ispol:
            scn_lbl += ' CR, H = {:.2f} T'.format(lgrws['field'])
            lgrws['type'] = 'CR'
        else:
            scn_lbl += ' CL, H = {:.2f} T'.format(lgrws['field'])
            lgrws['type'] = 'CL'

        if h_sgn < 0:  # H- data
            pos.raw_imp = pd.concat([pos.raw_imp, rawen], axis=1)
            pos.raw_imp = pd.concat([pos.raw_imp, rawdt], axis=1)
            pos.label.append(scn_lbl)
            pos.idx.append(scn_num)
            pos.dtype = 'H -'
        else:  # H+ data
            neg.raw_imp = pd.concat([neg.raw_imp, rawen], axis=1)
            neg.raw_imp = pd.concat([neg.raw_imp, rawdt], axis=1)
            neg.label.append(scn_lbl)
            neg.idx.append(scn_num)
            neg.dtype = 'H +'


def separate_hscans(guiobj, confobj, h_raw, dt_raw, time_raw, scn_num, ispol,
                    lgrws, log_dt, pos, neg):
    '''
    Separate magnetic field scans and fill positive and negative
    ScanData objects with raw data and labels depending on polarisation.

    Parameters
    ----------
    guiobj : GUI object
        Provides GUI dialogs.

    confobj : Configurations object.

    h_raw : Pandas Series
        Series containing field values of the scan.

    dt_raw : Pandas Series
        Series containing imported data.

    time_raw : Pandas Series
        Series containing timestamps of imported data.

    scn_num : str
        Scan number.

    ispol : bool
        True for CR
        False for CL.

    lgrws : dict
        Dictionary with log data collected from data logfile.

    log_dt : dict
        dictionary collecting log data.

    pos, neg : ScanData objects
        Collect raw data from positive and negative scans.

    Return
    ------
    Fill pos and neg ScanData objects with raw data as well as lgrws
    dictionary with scan labels.

    '''
    # Energy tolerance in eV. To discriminate between edge and pre-edge
    # scans
    e_tol = 2

    # label for graphs
    scn_lbl = scn_num

    h_raw_to_imp = h_raw.to_numpy()
    dt_raw_to_imp = dt_raw.to_numpy()

    rawdt = pd.DataFrame()

    # For hysteresis on fly scans are collected as separate columns
    # For hysteresis pont by point they are just appended on the same
    # columns
    if guiobj.analysis == 'hyst_fly':
        rawh = pd.DataFrame()
        rawdt[scn_num] = dt_raw_to_imp
        rawh['H' + scn_num] = h_raw_to_imp
    else:
        rawdt['H'] = np.around(h_raw_to_imp, decimals=3)
        rawdt['I'] = dt_raw_to_imp
        rawdt['t'] = time_raw.to_numpy()

        # List with min, max and number of times data
        # to be used for the realization of a common time scale for time
        # splitted point by point analysis
        min_t = np.nanmin(rawdt['t'])
        max_t = np.nanmax(rawdt['t'])

    if ispol:  # CR data
        pos.dtype = 'CR'
        # Edge scan considered
        if abs(log_dt['Edge_en'] - lgrws['mon_en']) < e_tol:
            lgrws['type'] = 'CR'
            lgrws['edge_mon_en'] = lgrws['mon_en']
            lgrws['pre_edge_mon_en'] = np.nan

            if guiobj.analysis == 'hyst_fly':
                # for splitted branch up and down division is made here
                if confobj.spl_brnch:
                    fst = h_raw_to_imp[0]
                    lst = h_raw_to_imp[-1]
                    # up branch
                    if (lst > fst):
                        pos.up = pd.concat([pos.up, rawdt], axis=1)
                        pos.up = pd.concat([pos.up, rawh], axis=1)
                        # positive field branch
                        if fst > 0:
                            pos.up_hp = True
                        # negative field branch
                        else:
                            pos.up_hn = True
                        pos.up_idx.append(scn_num)
                        pos.up_label.append(scn_lbl)
                    # down branch
                    else:
                        pos.down = pd.concat([pos.down, rawdt], axis=1)
                        pos.down = pd.concat([pos.down, rawh], axis=1)
                        # positive field branch
                        if fst > 0:
                            pos.down_hp = True
                        # negative field branch
                        else:
                            pos.down_nh = True
                        pos.dw_idx.append(scn_num)
                        pos.dw_label.append(scn_lbl)
                else:
                    pos.raw_imp = pd.concat([pos.raw_imp, rawdt], axis=1)
                    pos.raw_imp = pd.concat([pos.raw_imp, rawh], axis=1)
            else:
                #pos.raw_imp = pos.raw_imp.append(rawdt, ignore_index=True)
                pos.raw_imp = pd.concat([pos.raw_imp, rawdt],
                                        ignore_index=True)
                pos.min_t.append(min_t)
                pos.max_t.append(max_t)

            pos.label.append(scn_lbl)
            pos.idx.append(scn_num)
        else:
            # Pre-edge scan
            # Set True pre_edge attribute to indicate the presence of
            # pre-edge scans
            pos.pre_edge = True
            scn_lbl = 'PE-' + scn_lbl
            lgrws['type'] = 'PE-CR'
            lgrws['edge_mon_en'] = np.nan
            lgrws['pre_edge_mon_en'] = lgrws['mon_en']

            if guiobj.analysis == 'hyst_fly':
                # for splitted branch up and down division is made here
                if confobj.spl_brnch:
                    fst = h_raw_to_imp[0]
                    lst = h_raw_to_imp[-1]
                    # up branch
                    if (lst > fst):
                        pos.pe_up = pd.concat([pos.pe_up, rawdt], axis=1)
                        pos.pe_up = pd.concat([pos.pe_up, rawh], axis=1)
                        # positive field branch
                        if fst > 0:
                            pos.pe_up_hp = True
                        # negative field branch
                        else:
                            pos.pe_up_hn = True
                        pos.pe_up_idx.append(scn_num)
                        pos.pe_up_label.append(scn_lbl)
                    # down branch
                    else:
                        pos.pe_down = pd.concat([pos.pe_down, rawdt], axis=1)
                        pos.pe_down = pd.concat([pos.pe_down, rawh], axis=1)
                        # positive field branch
                        if fst > 0:
                            pos.pe_down_hp = True
                        # negative field branch
                        else:
                            pos.pe_down_nh = True
                        pos.pe_dw_idx.append(scn_num)
                        pos.pe_dw_label.append(scn_lbl)
                else:
                    pos.pe_raw_imp = pd.concat([pos.pe_raw_imp, rawh], axis=1)
                    pos.pe_raw_imp = pd.concat([pos.pe_raw_imp, rawdt], axis=1)
            else:
                #pos.pe_raw_imp = pos.pe_raw_imp.append(rawdt,
                #                                       ignore_index=True)
                pos.pe_raw_imp = pd.concat([pos.pe_raw_imp, rawdt],
                                            ignore_index=True)
                pos.min_t.append(min_t)
                pos.max_t.append(max_t)

            pos.pe_label.append(scn_lbl)
            pos.pe_idx.append(scn_num)
    else:  # CL data
        neg.dtype = 'CL'
        # Edge scan considered
        if abs(log_dt['Edge_en'] - lgrws['mon_en']) < e_tol:
            lgrws['type'] = 'CL'
            lgrws['edge_mon_en'] = lgrws['mon_en']
            lgrws['pre_edge_mon_en'] = np.nan

            if guiobj.analysis == 'hyst_fly':
                # for splitted branch up and down division is made here
                if confobj.spl_brnch:
                    fst = h_raw_to_imp[0]
                    lst = h_raw_to_imp[-1]
                    # up branch
                    if (lst > fst):
                        neg.up = pd.concat([neg.up, rawdt], axis=1)
                        neg.up = pd.concat([neg.up, rawh], axis=1)
                        # positive field branch
                        if fst > 0:
                            neg.up_hp = True
                        # negative field branch
                        else:
                            neg.up_hn = True
                        neg.up_idx.append(scn_num)
                        neg.up_label.append(scn_lbl)
                    # down branch
                    else:
                        neg.down = pd.concat([neg.down, rawdt], axis=1)
                        neg.down = pd.concat([neg.down, rawh], axis=1)
                        # positive field branch
                        if fst > 0:
                            neg.down_hp = True
                        # negative field branch
                        else:
                            neg.down_nh = True
                        neg.dw_idx.append(scn_num)
                        neg.dw_label.append(scn_lbl)
                else:
                    neg.raw_imp = pd.concat([neg.raw_imp, rawh], axis=1)
                    neg.raw_imp = pd.concat([neg.raw_imp, rawdt], axis=1)
            else:
                #neg.raw_imp = neg.raw_imp.append(rawdt, ignore_index=True)
                neg.raw_imp = pd.concat([neg.raw_imp, rawdt],
                                        ignore_index=True)
                neg.min_t.append(min_t)
                neg.max_t.append(max_t)
            neg.label.append(scn_lbl)
            neg.idx.append(scn_num)
        else:
            # Pre-edge scan
            # Set True pre_edge attribute to indicate the presence of
            # pre-edge scans
            neg.pre_edge = True
            scn_lbl = 'PE-' + scn_lbl
            lgrws['type'] = 'PE-CL'
            lgrws['edge_mon_en'] = np.nan
            lgrws['pre_edge_mon_en'] = lgrws['mon_en']

            if guiobj.analysis == 'hyst_fly':
                # for splitted branch up and down division is made here
                if confobj.spl_brnch:
                    fst = h_raw_to_imp[0]
                    lst = h_raw_to_imp[-1]
                    # up branch
                    if (lst > fst):
                        neg.pe_up = pd.concat([neg.pe_up, rawdt], axis=1)
                        neg.pe_up = pd.concat([neg.pe_up, rawh], axis=1)
                        # positive field branch
                        if fst > 0:
                            neg.pe_up_hp = True
                        # negative field branch
                        else:
                            neg.pe_up_hn = True
                        neg.pe_up_idx.append(scn_num)
                        neg.pe_up_label.append(scn_lbl)
                    # down branch
                    else:
                        neg.pe_down = pd.concat([neg.pe_down, rawdt], axis=1)
                        neg.pe_down = pd.concat([neg.pe_down, rawh], axis=1)
                        # positive field branch
                        if fst > 0:
                            neg.pe_down_hp = True
                        # negative field branch
                        else:
                            neg.pe_down_nh = True
                        neg.pe_dw_idx.append(scn_num)
                        neg.pe_dw_label.append(scn_lbl)
                else:
                    neg.pe_raw_imp = pd.concat([neg.pe_raw_imp, rawh], axis=1)
                    neg.pe_raw_imp = pd.concat([neg.pe_raw_imp, rawdt], axis=1)
            else:
                #neg.pe_raw_imp = neg.pe_raw_imp.append(rawdt,
                #                                       ignore_index=True)
                neg.pe_raw_imp = pd.concat([neg.pe_raw_imp, rawdt],
                                            ignore_index=True)
                neg.min_t.append(min_t)
                neg.max_t.append(max_t)

            neg.pe_label.append(scn_lbl)
            neg.pe_idx.append(scn_num)


def output_fls_escan(guiobj, pos, neg, scanobj):
    '''
    Organize output data and columns names for energy scans data
    analysis.    

    Parameters
    ----------
    guiobj : GUI object
        Provide GUI dialogues.

    pos : ScanData object
        ScanData for CR scans.

    neg : ScanData object
        ScanData for CL scans.

    scanobj : EngyScan object
        Contains results of X-Ray dichroism computations.

    Return
    ------
    Numpy array with data to be saved
    col_nms str with column names
    col_desc str with column description
    '''
    # Collect data
    out_data = np.stack((scanobj.energy,
                         scanobj.energycal,
                         pos.aver,
                         pos.aver_std,
                         neg.aver,
                         neg.aver_std,
                         scanobj.xd,
                         scanobj.xd_err,
                         scanobj.xd_aver,
                         scanobj.xd_aver_err,
                         pos.norm,
                         pos.norm_er,
                         neg.norm,
                         neg.norm_er,
                         scanobj.xd_pc,
                         scanobj.xd_pc_er,
                         scanobj.pos_corr,
                         scanobj.pos_corr_er,
                         scanobj.neg_corr,
                         scanobj.neg_corr_er,
                         scanobj.xd_pc_av_ej,
                         scanobj.xd_pc_int_er,
                         pos.norm_int,
                         pos.norm_int_er,
                         neg.norm_int,
                         neg.norm_int_er,
                         scanobj.xd_pc_int,
                         scanobj.xd_pc_int_er,
                         scanobj.xd_pc_av_ej_int,
                         scanobj.xd_pc_av_ej_int_er), axis=1)
    if guiobj.infile_ref:
        ref = ' norm by ref'
    else:
        ref = ''
    # Output file column names
    col_nms = 'E_int{} (eV),'.format(ref)
    col_nms += 'E{} (eV),'.format(ref)
    col_nms += '{}_av{} (a.u.),'.format(pos.dtype, ref)
    col_nms += '{}_av_er{} (a.u.),'.format(pos.dtype, ref)
    col_nms += '{}_av{} (a.u.),'.format(neg.dtype, ref)
    col_nms += '{}_av_er{} (a.u.),'.format(neg.dtype, ref)
    col_nms += '{}{} (a.u.),'.format(guiobj.analysis, ref)
    col_nms += '{}{}_er (a.u.),'.format(guiobj.analysis, ref)
    col_nms += 'XAS avgd{} (a.u.),'.format(ref)
    col_nms += 'XAS avgd{}_er (a.u.),'.format(ref)
    col_nms += '{}oPE{} (a.u.),'.format(pos.dtype, ref)
    col_nms += '{}oPE{}_er (a.u.),'.format(pos.dtype, ref)
    col_nms += '{}oPE{} (a.u.),'.format(neg.dtype, ref)
    col_nms += '{}oPE{}_er (a.u.),'.format(neg.dtype, ref)
    col_nms += '{}{} (%),'.format(guiobj.analysis, ref)
    col_nms += '{}{}_er (%),'.format(guiobj.analysis, ref)
    col_nms += '{} Corr{} (a.u.),'.format(pos.dtype, ref)
    col_nms += '{} Corr{}_er (a.u.),'.format(pos.dtype, ref)
    col_nms += '{} Corr{} (a.u.),'.format(neg.dtype, ref)
    col_nms += '{} Corr{}_er (a.u.),'.format(neg.dtype, ref)
    col_nms += '{} XAS avgd{} (%),'.format(guiobj.analysis, ref)
    col_nms += '{} XAS avgd{}_er (%),'.format(guiobj.analysis, ref)
    col_nms += '{}oPE_int{} (a.u.),'.format(pos.dtype, ref)
    col_nms += '{}oPE_int{}_er (a.u.),'.format(pos.dtype, ref)
    col_nms += '{}oPE_int{} (a.u.),'.format(neg.dtype, ref)
    col_nms += '{}oPE_int{}_er (a.u.),'.format(neg.dtype, ref)
    col_nms += '{}_int{} (%),'.format(guiobj.analysis, ref)
    col_nms += '{}_int{}_er (%),'.format(guiobj.analysis, ref)
    col_nms += '{} XAS avgd_int{} (%),'.format(guiobj.analysis, ref)
    col_nms += '{} XAS avgd_int_er{} (%)\n'.format(guiobj.analysis, ref)

    # Output file column descriptions
    col_desc = 'Interpolated E,'
    col_desc += 'Calibrated E - offset: {} eV,'.format(scanobj.offset)
    col_desc += 'Averaged {},'.format(pos.dtype)
    col_desc += 'STD of averaged {},'.format(pos.dtype)
    col_desc += 'Averaged {},'.format(neg.dtype)
    col_desc += 'STD of averaged {},'.format(neg.dtype)
    col_desc += '{}_av - {}_av,'.format(neg.dtype, pos.dtype)
    col_desc += 'Error of {}_av - {}_av,'.format(neg.dtype, pos.dtype)
    if guiobj.analysis in guiobj.type['xnld']:
        col_desc += '{}_av + ((2cos^2 - sin^2) {}_av)) / 3 cos^2,'.format(
            pos.dtype, neg.dtype)
    else:
        col_desc += '({}_av +  {}_av) / 2,'.format(neg.dtype, pos.dtype)
    col_desc += 'Error of XAS avgd{} (a.u.),'.format(ref)
    col_desc += '{}_av/PE_av,'.format(pos.dtype)
    col_desc += 'Error on {}_av/PE_av,'.format(pos.dtype)
    col_desc += '{}_av/PE_av,'.format(neg.dtype)
    col_desc += 'Error on {}_av/PE_av,'.format(neg.dtype)
    if guiobj.analysis in guiobj.type['xnld']:
        col_desc += ('Normalized {} 300 * cos^2 *'.format(guiobj.analysis)
                     + ' ({0}oPE - {1}oPE) / ({1}_EJNor'.format(
                        neg.dtype, pos.dtype)
                     + ' + (2cos^2 - sin^2) * {}_EJNor),'.format(neg.dtype))
    else:
        col_desc += ('Normalized {} 200 * '.format(guiobj.analysis)
                     + '({0}oPE - {1}oPE) / ({1}_EJNorm'.format(
                        neg.dtype, pos.dtype)
                     + ' + {}_EJNorm),'.format(neg.dtype))
    col_desc += 'Error on normalized {},'.format(guiobj.analysis)
    col_desc += '{0}_av * XASavgdPE_av / {0}PE_av,'.format(pos.dtype)
    col_desc += 'Error on {0}_av * XASavgdPE_av / {0}PE_av,'.format(pos.dtype)
    col_desc += '{0}_av * XASavgdPE_av / {0}PE_av,'.format(neg.dtype)
    col_desc += 'Error on {0}_av * XASavgdPE_av / {0}PE_av,'.format(neg.dtype)
    col_desc += 'XAS avgd edge-jump norm,'
    col_desc += 'Error on XAS avgd edge-jump norm,'
    col_desc += '{}_av/PE_int,'.format(pos.dtype)
    col_desc += 'Error on {}_av/PE_int,'.format(pos.dtype)
    col_desc += '{}_av/PE_int,'.format(neg.dtype)
    col_desc += 'Error on {}_av/PE_int,'.format(neg.dtype)
    if guiobj.analysis in guiobj.type['xnld']:
        col_desc += ('Normalized {} 300 * cos^2 *'.format(guiobj.analysis)
                     + ' ({0}oPE_int - {1}oPE_int) / ({1}_EJNor_int'.format(
                        neg.dtype, pos.dtype)
                     + ' + (2cos^2 - sin^2) * {}_EJNor_int),'.format(
                        neg.dtype))
    else:
        col_desc += ('Normalized {} 200 * '.format(guiobj.analysis)
                     + '({0}oPE_int - {1}oPE_int) / ({1}_EJNor_int'.format(
                        neg.dtype, pos.dtype)
                     + ' + {}_EJNor_int),'.format(neg.dtype))
    col_desc += 'Error on normalized {},'.format(guiobj.analysis)
    col_desc += '{0}_int{1} normalized by EJ of average {0},'.format(
        guiobj.analysis, ref)
    col_desc += 'Error on {0}_int{1} normalized by EJ of average {0}\n'.format(
        guiobj.analysis, ref)

    return out_data, col_nms, col_desc


def output_fls_hscan(scanobj):
    '''
    Organize output data and columns names for hysteresis on fly
    analysis.

    Parameters
    ----------
    scanobj : FieldScan object
        Contains results of X-Ray dichroism computations.

    Return
    ------
    Numpy array with data to be saved
    col_nms str with column names
    col_desc str with column description.
    '''
    # Collect data
    if scanobj.pre_edge:
        # Aggregate data
        out_data = np.stack((scanobj.flds, scanobj.cr_up, scanobj.cl_up,
                             scanobj.cr_down, scanobj.cl_down, scanobj.edg_up,
                             scanobj.edg_down, scanobj.edg_up_norm,
                             scanobj.edg_down_norm, scanobj.cr_pe_up,
                             scanobj.cl_pe_up, scanobj.cr_pe_down,
                             scanobj.cl_pe_down, scanobj.up_w_pe,
                             scanobj.dw_w_pe, scanobj.up_perc,
                             scanobj.down_perc), axis=1)
    else:
        out_data = np.stack((scanobj.flds, scanobj.cr_up, scanobj.cl_up,
                             scanobj.cr_down, scanobj.cl_down, scanobj.edg_up,
                             scanobj.edg_down, scanobj.edg_up_norm,
                             scanobj.edg_down_norm), axis=1)

    # Output file column names
    col_nms = 'H (T),'
    col_nms += 'CR up (a.u.),'
    col_nms += 'CL up (a.u.),'
    col_nms += 'CR down (a.u.),'
    col_nms += 'CL down (a.u.),'
    col_nms += 'Edge only Up (a.u.),'
    col_nms += 'Edge only Down (a.u.),'
    col_nms += 'Edge only Up - Norm (a.u.),'
    col_nms += 'Edge only Down - Norm (a.u.)\n'

    # Output file column description
    col_desc = 'Magnetic field,'
    col_desc += 'Interpolated & averaged CR,'
    col_desc += 'Interpolated & averaged CL,'
    col_desc += 'Interpolated & averaged CR,'
    col_desc += 'Interpolated & averaged CL,'
    col_desc += 'CL up - CR up - only edge,'
    col_desc += 'CL down - CR down - only edge,'
    col_desc += 'Normalized to 1,'
    col_desc += 'Normalized to 1\n'

    if scanobj.pre_edge:
        col_nms = col_nms.removesuffix('\n')
        col_nms += ',CR up pre-edge (a.u.),'
        col_nms += 'CL up pre-edge (a.u.),'
        col_nms += 'CR down pre-edge (a.u.),'
        col_nms += 'CL down pre-edge (a.u.),'
        col_nms += 'Up (a.u.),'
        col_nms += 'Down (a.u.),'
        col_nms += 'Up Norm (%),'
        col_nms += 'Down Norm (%)\n'

        col_desc = col_desc.removesuffix('\n')
        col_desc += ',Interpolated & averaged CR @ pre-edge,'
        col_desc += 'Interpolated & averaged CL @ pre-edge,'
        col_desc += 'Interpolated & averaged CR @ pre-edge,'
        col_desc += 'Interpolated & averaged CL @ pre-edge,'
        col_desc += 'CLup/CLupPE - CRup/CRupPE,'
        col_desc += 'CLdown/CLdownPE - CRdown/CRdownPE,'
        col_desc += '100 * 2 * Up/((CLup/CLupPE + CRup/CRupPE) - 2),'
        col_desc += '100 * 2 * Down/((CLdw/CLdwPE + CRdw/CRdwPE) - 2)\n'

    return out_data, col_nms, col_desc


def output_fls_ptbypt(guiobj, scanobj):
    '''
    Organize output data and columns names for hysteresis on point by
    point analysis.

    Parameters
    ----------
    guiobj : GUI object
        Provide GUI dialogues.

    scanobj : FieldScan object
        Contains results of X-Ray dichroism computations.

    Return
    ------
    Numpy array with data to be saved
    col_nms str with column names
    col_desc str with column description.
    '''
    if guiobj.analysis == 'hyst_t_aver':
        # Output column names
        col_nms = 'H (T),'
        col_nms += 'CR (a.u.),'
        col_nms += 'ErrCR (a.u.),'
        col_nms += 'CL (a.u.),'
        col_nms += 'ErrCL (a.u.),'
        col_nms += 'Edge only (a.u.),'
        col_nms += 'ErrEdge only (a.u.)\n'

        # Output file column description
        col_desc = 'Magnetic field,'
        col_desc += 'Averaged CR,'
        col_desc += 'Std of averaged CR,'
        col_desc += 'Averaged CL,'
        col_desc += 'Std of averaged CL,'
        col_desc += 'CL - CR - only edge,'
        col_desc += 'Std of CL - CR - only edge\n'

        if scanobj.pre_edge:
            # Order data
            scanobj.xmcd = scanobj.xmcd[['H', 'CR', 'dCR', 'CL', 'dCL', 'edge',
                                         'Dedge', 'CRpe', 'dCRpe', 'CLpe',
                                         'dCLpe', 'pre-edge', 'Dpre-edge',
                                         'no-norm', 'norm', 'Dnorm']]
            # Output column names
            col_nms = col_nms.removesuffix('\n')
            col_nms += ',CR_pe (a.u.),'
            col_nms += 'ErrCR_pe (a.u.),'
            col_nms += 'CL_pe (a.u.),'
            col_nms += 'ErrCL_pe (a.u.),'
            col_nms += 'PreEdge only (a.u.),'
            col_nms += 'ErrPreEdge only (a.u.),'
            col_nms += 'Hyst not norm (a.u.),'
            col_nms += 'Hyst norm (%),'
            col_nms += 'ErrHyst norm (%)\n'

            # Output column description
            col_desc = col_desc.removesuffix('\n')
            col_desc += ',Averaged CR pre-edge,'
            col_desc += 'Std of averaged CR pre-edge,'
            col_desc += 'Averaged CL pre-edge,'
            col_desc += 'Std of averaged CL pre-edge,'
            col_desc += 'CL_pe - CR_pe,'
            col_desc += 'Std of CL_pe - CR_pe'
            col_desc += 'CL/CL_pe - CR/CR_pe,'
            col_desc += '100 * 2 * Hyst_not_norm/((CL/CL_pe + CR/CR_PE) - 2),'
            col_desc += 'Std of Hyst norm\n'
        else:
            # Order data
            scanobj.xmcd = scanobj.xmcd[['H', 'CR', 'dCR', 'CL', 'dCL', 'edge',
                                         'Dedge']]
        scanobj.xmcd = scanobj.xmcd.sort_values(by=['H'])
    else:
        # Output column names
        col_nms = 'H (T),'
        col_nms += 't (s),'
        col_nms += 'CR (a.u.),'
        col_nms += 'CL (a.u.),'
        col_nms += 'Edge only (a.u.)\n'

        # Output file column description
        col_desc = 'Magnetic field,'
        col_desc += 'Time,'
        col_desc += 'CR,'
        col_desc += 'CL,'
        col_desc += 'CL - CR - only edge\n'

        if scanobj.pre_edge:
            # Order data
            scanobj.xmcd = scanobj.xmcd[['H', 't', 'CR', 'CL', 'edge', 'CRpe',
                                         'CLpe', 'pre-edge', 'no-norm', 'norm']]
            # Output column names
            col_nms = col_nms.removesuffix('\n')
            col_nms += ',CR_pe (a.u.),'
            col_nms += 'CL_pe (a.u.),'
            col_nms += 'PreEdge only (a.u.),'
            col_nms += 'Hyst not norm (a.u.),'
            col_nms += 'Hyst norm (%)\n'

            # Output column description
            col_desc = col_desc.removesuffix('\n')
            col_desc += ',CR pre-edge,'
            col_desc += 'CL pre-edge,'
            col_desc += 'CL_pe - CR_pe,'
            col_desc += 'CL/CL_pe - CR/CR_pe,'
            col_desc += '100 * 2 * Hyst_not_norm/((CL/CL_pe + CR/CR_PE) - 2)\n'
        else:
            # Order data
            scanobj.xmcd = scanobj.xmcd[['H', 't', 'CR', 'CL', 'edge']]
        scanobj.xmcd = scanobj.xmcd.sort_values(by=['t', 'H'])
    out_data = scanobj.xmcd.to_numpy()

    return out_data, col_nms, col_desc


def output_plot_escan(guiobj, pos, neg, scanobj, log_dt):
    '''
    Create a final figure with plot reporting averaged positive and
    negative XAS and percentage X-Ray dichroism.

    Parameters
    ----------
    guiobj : GUI object
        Provide GUI dialogues.

    pos : ScanData object
        ScanData for positive scans (CR, LH, H-).

    neg : ScanData object
        ScanData for negative scans (CL, LV, H+).

    scanobj : EngyScan object
        Contains results of X-Ray dichroism computations.

    log_dt : dict
        Collect data for logfile.

    Return
    ------
    Two pyplot figure objects one with data coming from computation of
    X-Ray dichroism using edge-jump and interpolated value of edge-jump
    respectively.
    '''
    if guiobj.infile_ref:
        ref = ' norm by ref'
    else:
        ref = ''

    class sh:
        # Class with functions to show/hide error bands
        def __init__(self, f):
            self.f = f

        def press(self, event):
            if self.f == 1:
                if f11.get_visible():
                    f11.set(visible=False)
                    f12.set(visible=False)
                    f21.set(visible=False)
                    plt.draw()
                else:
                    f11.set(visible=True)
                    f12.set(visible=True)
                    f21.set(visible=True)
                    plt.draw()
            else:
                if f31.get_visible():
                    f31.set(visible=False)
                    f32.set(visible=False)
                    f4.set(visible=False)
                    plt.draw()
                else:
                    f31.set(visible=True)
                    f32.set(visible=True)
                    f4.set(visible=True)
                    plt.draw()

    f1, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
    if guiobj.analysis in guiobj.type['xnld']:
        ax1.plot(scanobj.energycal, neg.aver, color='darkorange',
                 label=neg.dtype)
        f11 = ax1.fill_between(scanobj.energycal, neg.aver - neg.aver_std,
            neg.aver + neg.aver_std, color='darkorange', alpha=0.2,
            visible=True)
        ax1.plot(scanobj.energycal, pos.aver, color='darkviolet',
                 label=pos.dtype)
        f12 = ax1.fill_between(scanobj.energycal, pos.aver - pos.aver_std,
            pos.aver + pos.aver_std, color='darkviolet', alpha=0.2,
            visible=True)
    elif guiobj.analysis in guiobj.type['xmcd']:
        ax1.plot(scanobj.energycal, neg.aver,
                 color='blue', label=r'$\sigma^-$')
        f11 = ax1.fill_between(scanobj.energycal, neg.aver - neg.aver_std,
            neg.aver + neg.aver_std, color='blue', alpha=0.2, visible=True)
        ax1.plot(scanobj.energycal, pos.aver, color='red', label=r'$\sigma^+$')
        f12 = ax1.fill_between(scanobj.energycal, pos.aver - pos.aver_std,
            pos.aver + pos.aver_std, color='red', alpha=0.2, visible=True)
    else:
        ax1.plot(scanobj.energycal, neg.aver, color='blue', label=neg.dtype)
        f11 = ax1.fill_between(scanobj.energycal, neg.aver - neg.aver_std,
            neg.aver + neg.aver_std, color='blue', alpha=0.2, visible=True)
        ax1.plot(scanobj.energycal, pos.aver, color='red', label=pos.dtype)
        f12 = ax1.fill_between(scanobj.energycal, pos.aver - pos.aver_std,
            pos.aver + pos.aver_std, color='red', alpha=0.2, visible=True)
    ax1.set_ylabel('XAS{} (a.u.)'.format(ref))
    ax1.legend()
    if guiobj.analysis in guiobj.type['xnld']:
        ax2.plot(scanobj.energycal, scanobj.xd_pc, color='black',
                 label=guiobj.analysis)
        f21 = ax2.fill_between(scanobj.energycal,
            scanobj.xd_pc - scanobj.xd_pc_er, scanobj.xd_pc + scanobj.xd_pc_er,
            color='black', alpha=0.2, visible=True)
        ax2.axhline(y=0, color='darkgray')
    else:
        ax2.plot(scanobj.energycal, scanobj.xd_pc, color='green',
                 label=guiobj.analysis)
        f21 = ax2.fill_between(scanobj.energycal,
            scanobj.xd_pc - scanobj.xd_pc_er, scanobj.xd_pc + scanobj.xd_pc_er,
            color='green', alpha=0.2, visible=True)
        ax2.axhline(y=0, color='black')
    ax2.set_ylabel('{}{} (%)'.format(guiobj.analysis, ref))
    ax2.set_xlabel('E (eV)')
    ax2.legend()
    f1.suptitle('Edge : {:.2f}, PreEdge : {:.2f}, Edge-jump : {:.4f}'.format(
        log_dt['exper_edge'], log_dt['setted_pedg'], log_dt['xas_aver_ej'])
        + '+/-{:.4f}\n'.format(log_dt['xas_aver_ej_er'])
        + '{} @ edge energy : {:.2f} +/- {:.2f} %\n'.format(guiobj.analysis,
        log_dt['edge_pc'], log_dt['edge_pc_er'])
        + 'T = {} +/- {} K, H = {} T, '.format(log_dt['temp'],
        log_dt['temp_er'], log_dt['field'])
        + r'$\theta$ = {}°'.format(log_dt['angle']))

    # Button to show/hide errors
    callback = sh(1)
    axbutt1 = f1.add_axes([0.81, 0.0, 0.1, 0.065])
    butt1 = Button(axbutt1, label='S/H errs', color='whitesmoke')
    butt1.on_clicked(callback.press)

    f2, (ax3, ax4) = plt.subplots(2, 1, sharex=True)
    if guiobj.analysis in guiobj.type['xnld']:
        ax3.plot(scanobj.energycal, neg.aver, color='darkorange',
                 label=neg.dtype)
        f31 = ax3.fill_between(scanobj.energycal, neg.aver - neg.aver_std,
            neg.aver + neg.aver_std, color='darkorange', alpha=0.2,
            visible=True)
        ax3.plot(scanobj.energycal, pos.aver, color='darkviolet',
                 label=pos.dtype)
        f32 = ax3.fill_between(scanobj.energycal, pos.aver - pos.aver_std,
            pos.aver + pos.aver_std, color='darkviolet', alpha=0.2,
            visible=True)
    elif guiobj.analysis in guiobj.type['xmcd']:
        ax3.plot(scanobj.energycal, neg.aver,
                 color='blue', label=r'$\sigma^-$')
        f31 = ax3.fill_between(scanobj.energycal, neg.aver - neg.aver_std,
            neg.aver + neg.aver_std, color='blue', alpha=0.2, visible=True)
        ax3.plot(scanobj.energycal, pos.aver, color='red', label=r'$\sigma^+$')
        f32= ax3.fill_between(scanobj.energycal, pos.aver - pos.aver_std,
            pos.aver + pos.aver_std, color='red', alpha=0.2, visible=True)
    else:
        ax3.plot(scanobj.energycal, neg.aver, color='blue', label=neg.dtype)
        f31 = ax3.fill_between(scanobj.energycal, neg.aver - neg.aver_std,
            neg.aver + neg.aver_std, color='blue', alpha=0.2, visible=True)
        ax3.plot(scanobj.energycal, pos.aver, color='red', label=pos.dtype)
        f32 = ax3.fill_between(scanobj.energycal, pos.aver - pos.aver_std,
            pos.aver + pos.aver_std, color='red', alpha=0.2, visible=True)
    ax3.set_ylabel('XAS{} (a.u.)'.format(ref))
    ax3.legend()
    if guiobj.analysis in guiobj.type['xnld']:
        ax4.plot(scanobj.energycal, scanobj.xd_pc_int, color='black',
                 label=guiobj.analysis)
        f4 = ax4.fill_between(scanobj.energycal,
            scanobj.xd_pc_int - scanobj.xd_pc_int_er,
            scanobj.xd_pc_int + scanobj.xd_pc_int_er, color='black',
            alpha=0.2, visible=True)
        ax4.axhline(y=0, color='darkgray')
    else:
        ax4.plot(scanobj.energycal, scanobj.xd_pc_int, color='green',
                 label=guiobj.analysis)
        f4 = ax4.fill_between(scanobj.energycal,
            scanobj.xd_pc_int - scanobj.xd_pc_int_er,
            scanobj.xd_pc_int + scanobj.xd_pc_int_er, color='green',
            alpha=0.2, visible=True)
        ax4.axhline(y=0, color='black')
    ax4.set_ylabel('{}{}_int (%)'.format(guiobj.analysis, ref))
    ax4.set_xlabel('E (eV)')
    ax4.legend()
    if guiobj.bsl_int:
        f2.suptitle('Edge : {:.2f}, PreEdge : {:.2f},'.format(
            log_dt['exper_edge'], log_dt['setted_pedg'])
            + r' $\lambda$ : {:.2f}'.format(log_dt['lambda']/1E7)
            + 'Edge-jump : {:.4f}'.format(log_dt['xas_aver_ej_int'])
            + '+/-{:.4f}\n'.format(log_dt['xas_aver_ej_int_er'])
            + '{} @ edge energy : {:.2f} +/- {:.2f} %\n'.format(
            guiobj.analysis, log_dt['edge_pc_int'], log_dt['edge_pc_int_er'])
            + 'T = {} +/- {} K, H = {} T, '.format(log_dt['temp'],
            log_dt['temp_er'], log_dt['field'])
            + r'$\theta$ = {}°'.format(log_dt['angle']))
    else:
        f2.suptitle('Edge : {:.2f}, PreEdge : {:.2f},'.format(
            log_dt['exper_edge'], log_dt['setted_pedg'])
            + ' Post-Edge : {:.2f}\n'.format(log_dt['setted_postedg'])
            + ' Edge-jump : {:.4f}'.format(log_dt['xas_aver_ej_int'])
            + '+/-{:.4f}\n'.format(log_dt['xas_aver_ej_int_er'])
            + '{} @ edge energy : {:.2f} +/- {:.2f} %\n'.format(
            guiobj.analysis, log_dt['edge_pc_int'], log_dt['edge_pc_int_er'])
             + 'T = {} +/- {} K, H = {} T, '.format(log_dt['temp'],
            log_dt['temp_er'], log_dt['field'])
            + r'$\theta$ = {}°'.format(log_dt['angle']))

    # Button to show/hide errors
    callback2 = sh(2)
    axbutt2 = f2.add_axes([0.81, 0.0, 0.1, 0.065])
    butt2 = Button(axbutt2, label='S/H errs', color='whitesmoke')
    butt2.on_clicked(callback2.press)

    plt.show()

    return f1, f2


def output_plot_hscan(scanobj, log_dt):
    '''
    Create a final figure with plot reporting XMCD hysteresis on the fly
    analized.

    Parameters
    ----------
    scanobj : FieldScan object
        Contains results of X-Ray dichroism computations.

    log_dt : dict
        Collect data for logfile.

    Return
    ------
    A pyplot figure object with branch up and down of XMCD hysteresis
    normalized to 1.
    If pre-edge data are present also a second pyplot figure object with
    branch up and down of XMCD hysteresis normalized with pre-edge data
    is returned.
    '''
    # Limit y range to avoid spikes near 0
    limy = 1.3

    # Collect edge scans used - for titling graph
    # set removes duplicates and then sort
    edge_used = sorted(set(log_dt['pos_up_chsn'] + log_dt['pos_dw_chsn']
                           + log_dt['neg_up_chsn'] + log_dt['neg_dw_chsn']))
    f1, ax1 = plt.subplots()

    ax1.plot(scanobj.flds, scanobj.edg_up_norm, label='UP', color='crimson')
    ax1.plot(scanobj.flds, scanobj.edg_down_norm, label='DOWN', color='green')
    ax1.set_title('E = {} eV, T = {} +/- {} K, '.format(log_dt['exper_edge'],
        log_dt['temp'], log_dt['temp_er']) + r'$\theta$ = {}°'.format(
        log_dt['angle']))
    ax1.legend()
    ax1.set_ylim(limy, -limy)
    ax1.axhline(y=0, color='darkgray')
    ax1.axvline(x=0, color='darkgray')
    ax1.set_xlabel('H (T)')
    ax1.set_ylabel('XMCD (a.u.)')

    f1.suptitle('Edge only - Scan ' + ' '.join(edge_used))

    if scanobj.pre_edge:
        # Collect pre-edge scans used - for titling graph
        # set removes duplicates and then sort
        pre_edge_used = sorted(set(log_dt['pos_pe_up_chsn']
                                   + log_dt['pos_pe_dw_chsn']
                                   + log_dt['neg_pe_up_chsn']
                                   + log_dt['neg_pe_dw_chsn']))

        f2, ax2 = plt.subplots()

        limy2 = max(
            abs(scanobj.up_perc[0]), abs(scanobj.up_perc[-1]),
            abs(scanobj.down_perc[0]), abs(scanobj.down_perc[-1])) * 1.3

        ax2.plot(scanobj.flds, scanobj.up_perc, label='UP', color='crimson')
        ax2.plot(scanobj.flds, scanobj.down_perc, label='DOWN', color='green')
        ax2.set_title('E = {} eV, Epe = {}, T = {} +/- {} K,'.format(
            log_dt['exper_edge'], log_dt['setted_pedg'], log_dt['temp'],
            log_dt['temp_er']) + r'$\theta$ = {}°'.format(log_dt['angle']))
        ax2.legend()
        ax2.set_ylim(limy2, -limy2)
        ax2.axhline(y=0, color='darkgray')
        ax2.axvline(x=0, color='darkgray')
        ax2.set_xlabel('H (T)')
        ax2.set_ylabel('XMCD (%)')

        f2.suptitle('Normalized hysteresis - Scan ' + ' '.join(edge_used)
                    + ' '.join(pre_edge_used))
        plt.show()
        return f1, f2
    else:
        plt.show()
        return f1


def output_plot_t_aver(scanobj, log_dt):
    '''
    Create a final figure with plot reporting XMCD hysteresis point by
    point analized with time averaging.

    Parameters
    ----------
    scanobj : FieldScan object
        Contains results of X-Ray dichroism computations.

    log_dt : dict
        Collect data for logfile.

    Return
    ------
    Two pyplot figure objects: CR and CL data in the first, XMCD edge
    only in the sedond one.

    If pre-edge data are present also a third pyplot figure object with
    XMCD normalized with pre-edge data is returned.
    '''
    data = scanobj.xmcd
    log_tbl = log_dt['log_tbl']

    int_scn_num = [int(x) for x in log_tbl['scn_num']]

    low_scn = np.amin(int_scn_num)
    hgh_scn = np.amax(int_scn_num)
    time = '{:.1f}--{:.1f} s'.format(log_dt['t_scale'][0],
                                     log_dt['t_scale'][1])
    s_nums = 'scans {} - {} '.format(low_scn, hgh_scn)
    scanobj.subtitle = s_nums + time

    f1, ax1 = plt.subplots()

    ax1.plot(data['H'], data['CR'], label='CR')
    ax1.plot(data['H'], data['CL'], label='CL')
    ax1.set_xlabel('H (T)')
    ax1.set_ylabel('XAS (a.u.)')
    ax1.set_title(r'E = {} eV, T = {} K, $\theta$ = {}°'.format(
        log_dt['exper_edge'], log_dt['temp'], log_dt['angle']))
    ax1.legend()
    ax1.axhline(y=0, color='darkgray')
    ax1.axvline(x=0, color='darkgray')
    f1.suptitle('CR and CL\n{}'.format(scanobj.subtitle))

    f2, ax2 = plt.subplots()
    ax2.errorbar(data['H'], data['edge'], yerr=data['Dedge'],
                 label='Edge only')
    ax2.set_xlabel('H (T)')
    ax2.set_ylabel('XMCD (a.u.)')
    ax2.set_title(r'E = {} eV, T = {} K, $\theta$ = {}°'.format(
        log_dt['exper_edge'], log_dt['temp'], log_dt['angle']))
    ax2.legend()
    ax2.axhline(y=0, color='darkgray')
    ax2.axvline(x=0, color='darkgray')
    f2.suptitle('XMCD edge only\n{}'.format(scanobj.subtitle))

    if scanobj.pre_edge:
        f3, ax3 = plt.subplots()
        ax3.errorbar(data['H'], data['norm'], yerr=data['Dnorm'],
                     label='Normalized')
        ax3.set_xlabel('H (T)')
        ax3.set_ylabel('XMCD (%)')
        ax3.set_title(r'E = {} eV, Epe = {}, T = {} K, $\theta$ = {}°'.format(
            log_dt['exper_edge'], log_dt['setted_pedg'], log_dt['temp'],
            log_dt['angle']))
        ax3.legend()
        ax3.axhline(y=0, color='darkgray')
        ax3.axvline(x=0, color='darkgray')
        f3.suptitle('XMCD normalized\n{}'.format(scanobj.subtitle))

        plt.show()
        return f1, f2, f3
    else:
        plt.show()
        return f1, f2


def output_plot_t_split(guiobj, scanobj, log_dt):
    '''
    Create a final figure with plot reporting XMCD hysteresis point by
    point analized with time splitting.

    Parameters
    ----------
    guiobj : GUI object
        Provide GUI dialogues.

    scanobj : FieldScan object
        Contains results of X-Ray dichroism computations.

    log_dt : dict
        Collect data for logfile.

    Return
    ------
    Two pyplot figure objects: CR and CL data in the first, XMCD edge
    only in the sedond one.

    If pre-edge data are present also a third pyplot figure object with
    XMCD normalized with pre-edge data is returned.
    '''
    log_tbl = log_dt['log_tbl']

    int_scn_num = [int(x) for x in log_tbl['scn_num']]

    low_scn = np.amin(int_scn_num)
    hgh_scn = np.amax(int_scn_num)
    time = '{:.1f}--{:.1f} s'.format(log_dt['t_scale'][0],
                                     log_dt['t_scale'][1])
    s_nums = 'scans {} - {} '.format(low_scn, hgh_scn)
    scanobj.subtitle = s_nums + time

    num_times = guiobj.num_times(log_dt) - 1

    # Create colormap
    col = cm.jet(np.linspace(0, 1, num_times))
    c = 0

    # Group data by times
    t_gr = scanobj.xmcd.groupby('t')

    f1, (ax11, ax12) = plt.subplots(2, 1, sharex=True)
    f2, ax2 = plt.subplots()

    steps = scanobj.time_scale[::int(np.rint(len(scanobj.time_scale)
                                             / num_times)) + 1]
    for i in steps:
        data = t_gr.get_group(i)
        label = np.around(i, decimals=2)
        ax2.plot(data['H'], data['edge'], color=col[c], label=label)
        ax11.plot(data['H'], data['CR'], color=col[c], label=label)
        ax12.plot(data['H'], data['CL'], color=col[c], label=label)
        c += 1

    ax2.set_xlabel('H (T)')
    ax2.set_ylabel('XMCD (a.u.)')
    ax2.set_title(r'E = {} eV, T = {} K, $\theta$ = {}°'.format(
        log_dt['exper_edge'], log_dt['temp'], log_dt['angle']))
    ax2.legend()
    ax2.axhline(y=0, color='darkgray')
    ax2.axvline(x=0, color='darkgray')
    f2.suptitle('XMCD edge only\n{}'.format(scanobj.subtitle))

    ax12.set_xlabel('H (T)')
    ax11.set_ylabel('CR (a.u.)')
    ax11.set_title(r'E = {} eV, T = {} K, $\theta$ = {}°'.format(
        log_dt['exper_edge'], log_dt['temp'], log_dt['angle']))
    ax11.legend()
    ax11.axhline(y=0, color='darkgray')
    ax11.axvline(x=0, color='darkgray')
    ax12.set_ylabel('CL (a.u.)')
    ax12.legend()
    ax12.axhline(y=0, color='darkgray')
    ax12.axvline(x=0, color='darkgray')
    f1.suptitle('CR and CL\n{}'.format(scanobj.subtitle))

    if scanobj.pre_edge:
        c = 0
        f3, ax3 = plt.subplots()

        for i in steps:
            data = t_gr.get_group(i)
            label = np.around(i, decimals=2)
            ax3.plot(data['H'], data['norm'], color=col[c], label=label)
            c += 1

        ax3.set_xlabel('H (T)')
        ax3.set_ylabel('XMCD (%)')
        ax3.set_title(r'E = {} eV, Epe = {}, T = {} K, $\theta$ = {}°'.format(
            log_dt['exper_edge'], log_dt['setted_pedg'], log_dt['temp'],
            log_dt['angle']))
        ax3.legend()
        ax3.axhline(y=0, color='darkgray')
        ax3.axvline(x=0, color='darkgray')
        f3.suptitle('XMCD normalized\n{}'.format(scanobj.subtitle))

        plt.show()
        return f1, f2, f3
    else:
        plt.show()
        return f1, f2


def save_data_escan(confobj, guiobj, pos, neg, scanobj, log_dt, pos_ref,
                    neg_ref, scanobj_ref_norm, log_dt_ref):
    '''
    Save output data, logdata and final plots for energy scan analysis.

    Parameters
    ----------
    confobj : configuration obj.

    guiobj : GUI object
        Provide GUI dialogues.

    pos : ScanData object
        ScanData for positive scans (CR, LH, H-).

    neg : ScanData object
        ScanData for negative scans (CL, LV, H+).

    scanobj : EngyScan object
        Contains results of X-Ray dichroism computations.

    log_dt : dict
        Collect data for logfile.

    pos_ref : ScanData object
        ScanData for positive reference scans.

    neg_ref : ScanData object
        ScanData for negative reference scans.

    scanobj_ref_norm : EngyScan object
        Contains results of X-Ray dichroism computations normalized by
        reference sample.

    log_dt_ref : dict
        Collect data of reference and normalized dara for logfile.
    '''
    # Not normalized by reference data and plots
    out_data, col_nms, col_desc = output_fls_escan(guiobj, pos, neg, scanobj)
    f1, f2 = output_plot_escan(guiobj, pos, neg, scanobj, log_dt)
    logtxt = confobj.escan_logfl_creator(guiobj, log_dt)

    if confobj.ref_norm:
        guiobj.infile_ref = True  # For graphs and columns labelling
        out_norm_dt, col_norm_nms, col_norm_desc = output_fls_escan(
            guiobj, pos_ref, neg_ref, scanobj_ref_norm)
        f3, f4 = output_plot_escan(guiobj, pos_ref, neg_ref, scanobj_ref_norm,
                                   log_dt_ref)
        out_data = np.concatenate((out_data, out_norm_dt), axis=1)
        col_nms = col_nms.rstrip('\n') + ',' + col_norm_nms
        col_desc = col_desc + ',' + col_norm_desc

        logtxt += '\n\n'
        logtxt += 'Logs for ref scans and normalized data by refererence.'
        logtxt += '\n\n'
        logtxt += confobj.escan_logfl_creator(guiobj, log_dt_ref)

    default_nm = ('{}_{}_scan_{}-{}_{}K_{}T_{}_{}.dat'.format(
        log_dt['Edge_name'], guiobj.analysis,
        log_dt['log_tbl']['scn_num'].iloc[0].rstrip('(D)'),
        log_dt['log_tbl']['scn_num'].iloc[-1], log_dt['temp'],
        log_dt['field'], log_dt['angle'], guiobj.sense))

    out_nm = guiobj.outfile_name(default_nm)

    if out_nm is None:
        pass
    elif out_nm == '':
        pass
    else:
        with open(out_nm, 'w') as fl:
            fl.write(col_nms)
            fl.write(col_desc)
            np.savetxt(fl, out_data, fmt='%.7e', delimiter=',')

        f1_out_nm = out_nm.rstrip('dat') + 'png'
        f1.savefig(f1_out_nm)
        f2_out_nm = out_nm.rstrip('.dat') + '_int.png'
        f2.savefig(f2_out_nm)

        if confobj.ref_norm:
            f3_out_nm = out_nm.rstrip('dat') + '_norm_by_ref.png'
            f3.savefig(f3_out_nm)
            f4_out_nm = out_nm.rstrip('.dat') + '_norm_by_ref_int.png'
            f4.savefig(f4_out_nm)

        logfl_nm = out_nm.rstrip('dat') + 'log'

        with open(logfl_nm, 'w') as fl:
            fl.write(logfl_nm + '\n\n')
            fl.write(logtxt)


def save_data_hscan(confobj, guiobj, scanobj, log_dt):
    '''
    Save output data, logdata and final plots for scan field on the fly
    analysis.

    Parameters
    ----------
    confobj : configuration obj.

    guiobj : GUI object
        Provide GUI dialogues.

    scanobj : FieldScan object
        Contains results of X-Ray dichroism computations.

    log_dt : dict
        Collect data for logfile.
    '''
    out_data, col_nms, col_desc = output_fls_hscan(scanobj)

    if scanobj.pre_edge:
        f1, f2 = output_plot_hscan(scanobj, log_dt)
    else:
        f1 = output_plot_hscan(scanobj, log_dt)

    logtxt = confobj.hscan_logfl_creator(log_dt)

    default_nm = ('{}_Hyst_scan_{}-{}_{}K_{}T_{}_{}.dat'.format(
        log_dt['Edge_name'], log_dt['log_tbl']['scn_num'].iloc[0],
        log_dt['log_tbl']['scn_num'].iloc[-1], log_dt['temp'], log_dt['field'],
        log_dt['angle'], guiobj.sense))

    out_nm = guiobj.outfile_name(default_nm)

    if out_nm is None:
        pass
    elif out_nm == '':
        pass
    else:
        with open(out_nm, 'w') as fl:
            fl.write(col_nms)
            fl.write(col_desc)
            np.savetxt(fl, out_data, fmt='%.7e', delimiter=',')

        f1_out_nm = out_nm.removesuffix('.dat') + '_Only Edge.png'
        f1.savefig(f1_out_nm)

        if scanobj.pre_edge:
            f2_out_nm = out_nm.removesuffix('dat') + 'png'
            f2.savefig(f2_out_nm)

        logfl_nm = out_nm.removesuffix('dat') + 'log'

        with open(logfl_nm, 'w') as fl:
            fl.write(logfl_nm + '\n\n')
            fl.write(logtxt)


def save_data_tscan(confobj, guiobj, scanobj, log_dt):
    '''
    Save output data, logdata and final plots for scan field point by
    point analysis.

    Parameters
    ----------
    confobj : configuration obj.

    guiobj : GUI object
        Provide GUI dialogues.

    scanobj : FieldScan object
        Contains results of X-Ray dichroism computations.

    log_dt : dict
        Collect data for logfile.
    '''
    out_data, col_nms, col_desc = output_fls_ptbypt(guiobj, scanobj)

    if guiobj.analysis == 'hyst_t_aver':
        name = 'Hyst time aver'
        if scanobj.pre_edge:
            f1, f2, f3 = output_plot_t_aver(scanobj, log_dt)
        else:
            f1, f2 = output_plot_t_aver(scanobj, log_dt)
    if guiobj.analysis == 'hyst_t_split':
        name = 'Hyst time split'
        if scanobj.pre_edge:
            f1, f2, f3 = output_plot_t_split(guiobj, scanobj, log_dt)
        else:
            f1, f2 = output_plot_t_split(guiobj, scanobj, log_dt)

    default_nm = ('{}_{}_scan_{}-{}_{}-{}s_{}K_{}T_{}_{}.dat'.format(
        log_dt['Edge_name'], name, log_dt['log_tbl']['scn_num'].iloc[0],
        log_dt['log_tbl']['scn_num'].iloc[-1], log_dt['t_scale'][0],
        log_dt['t_scale'][1], log_dt['temp'], log_dt['field'], log_dt['angle'],
        guiobj.sense))

    logtxt = confobj.ptbypt_logfl_creator(log_dt)

    out_nm = guiobj.outfile_name(default_nm)

    if out_nm is None:
        pass
    elif out_nm == '':
        pass
    else:
        with open(out_nm, 'w') as fl:
            fl.write(col_nms)
            fl.write(col_desc)
            np.savetxt(fl, out_data, fmt='%.7e', delimiter=',')

            f1_out_nm = out_nm.removesuffix('.dat') + '_CR_CL.png'
            f1.savefig(f1_out_nm)
            f2_out_nm = out_nm.removesuffix('.dat') + '_Edge Only.png'
            f2.savefig(f2_out_nm)

        if scanobj.pre_edge:
            f3_out_nm = out_nm.removesuffix('dat') + 'png'
            f3.savefig(f3_out_nm)

        logfl_nm = out_nm.removesuffix('dat') + 'log'

        with open(logfl_nm, 'w') as fl:
            fl.write(logfl_nm + '\n\n')
            fl.write(logtxt)
