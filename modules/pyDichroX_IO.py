"""
pyDichroX_IO.py

Methods to manage input and output files for XMCD, XNCD, XNLD and hysteresys
data analysis.

Methods
--------
open_import_escan(guiobj, confobj)
    Open input files and import data for energy scan experiments.

scan_importer(guiobj, confobj, pos, neg)
    Manage scan importing.

dt_raw_import(confobj, data)
    Based on configuration file select data to import depending on sensing.
    It also normalize data by i0 if normalization is not provided in datafile.

escan_split(confobj, data):
    Split energy scan containing multiple scans into single scans.

set_scn_num(confobj, f_name, pos, neg)
    Associate an identifier for each scan.

separate_scans(guiobj, confobj, e_raw, dt_raw, scn_lbl, scn_num, ispol, lgrws, pos, neg)
    Separate scans and fill positive and negative ScanData objects with raw data
    and labels depending on analysis type

output_fls_escan(guiobj, pos, neg, scanobj)
    Organize output data and columns names.

output_plot_escan(guiobj, pos, neg, scanobj, log_dt)
    Create a plot final figure reporting averagaed positive and negative XAS and
    percentage X-Ray dichroism.

save_data_escan(confobj, guiobj, pos, neg, scanobj, log_dt, pos_ref, neg_ref,
    scanobj_ref_norm, log_dt_ref)
    Save output data, logdata and final plots.
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

import modules.pyDichroX_datatreat as dt


def open_import_escan(guiobj, confobj):
    '''
    Open input files for energy scan experiments (XMCD, XNCD, XNXD, XNLD),
    import data and groups them based on absorption coefficient for circularly
    polirized light or based on value of linear polarization for linearly
    polarized light.

    Collect experimental data for compiling output logfile.

    If confobj.ref_norm data from a reference sample are prompted 

    Parameters
    ----------
    guiobj : GUI object
        Provides GUI dialogs.

    confobj : Configurations object

    Returns
    -------
    ScanData object with positive scans (sigma+, CR, H- or LH scans, depending
        on analysis)

    ScanData object with negative scans (sigma-, CL, H+ or LV scans, depending
        on analysis)

    dict, contains log information
     . log_tbl : pd.Dataframes, columns depend on data provided by datalogfiles
     . Edge_name : edge name
     . Edge_en : edge energy
     . PreEdge_en : pre-edge energy
     . PostEdge_en : post-edge energy
     . angle : rotation angle of the sample
     . bm_angle : incidence beam angle
     . escale : list with energy scale info [e_min, e_max, num_points]
     . exper_edge : experimental edge energy
     . setted_pedg : setted pre-edge energy
     . setted_postedg : setted post-edge energy
     . pe_int : interpolated pre-edge
     . recal : bool, if True recalibrate energy scale
     . offset : energy scale calibration offset
     . pos_chs : list of positive scans chosed for analysis
     . neg_chs : list of negative scans chosed for analysis
     . pos_ej : edge jump for positive scans
     . pos_ej_int : edge jump for positive scans using interpolated pre-edge
     . neg_ej : edge jump for negative scans
     . neg_ej_int : edge jump for negative scans using interpolated pre-edge
     . xas_aver_ej : edge jump for avereaged XAS sptectrum
     . xas_aver_ej_int : edge jump for averaged XAS sptectrum using interpolated
            pre-edge

    If reference data are considered in confobj also a ScanData object with
    positive scans and a ScanData object with negative data from reference are
    returned. 
    '''
    # Set edge and pre-edge energies
    edge = guiobj.chs_edge()

    # dictionary collecting log data
    log_dt = {}

    # Insert experimental rotation angle of the sample
    angle, bm_angle = guiobj.ask_angle()

    # Ask for sample temperature if not provided by datalog file
    T = -100  # just initialize, if provided by log they will not be used
    H = -100
    if confobj.ask_for_T:
        T = guiobj.ask_T()
    # Ask for magnetic field if not provided by datalog file
    if confobj.ask_for_H:
        H = guiobj.ask_H()

    # Instantiation of objects for different polarization signs
    # pos will contain sigma+, CR, LH, H- polarization data
    # neg will contain sigma-, CL, LV, H+ polarization data
    # ref are referred to reference data. They will be filled only if
    # confobj.ref_norm is True
    pos = dt.ScanData()
    neg = dt.ScanData()

    pos_ref = dt.ScanData()
    neg_ref = dt.ScanData()

    # Attribute to select message for GUI
    guiobj.infile_ref = False

    log_tbl = scan_importer(guiobj, confobj, pos, neg, T, H)

    # dictionary collecting log data
    log_dt['log_tbl'] = log_tbl
    log_dt['temp'] = np.round(log_tbl['t'].mean(), 1)
    log_dt['field'] = np.round(log_tbl['field'].abs().mean(), 1)
    log_dt['Edge_name'] = edge[0]
    log_dt['Edge_en'] = float(edge[1])
    log_dt['PreEdge_en'] = float(edge[2])
    log_dt['PostEdge_en'] = float(edge[3])
    log_dt['angle'] = angle
    log_dt['bm_angle'] = bm_angle

    if confobj.ref_norm:
        # Attribute to select message for GUI
        guiobj.infile_ref = True

        log_ref_tbl = scan_importer(guiobj, confobj, pos_ref, neg_ref, T, H)

        log_dt['log_ref_tbl'] = log_ref_tbl

        guiobj.infile_ref = False

    return pos, neg, log_dt, pos_ref, neg_ref

def scan_importer(guiobj, confobj, pos, neg, T, H):
    '''
    Manage scan importing.

    Open input files for energy scan experiments (XMCD, XNCD, XNXD, XNLD),
    import data and groups them based on absorption coefficient for circularly
    polirized light or based on value of linear polarization for linearly
    polarized light.

    Collect experimental data for compiling output logfile.

    Parameters
    ----------
    guiobj : GUI object
        Provides GUI dialogs.

    confobj : Configurations object

    pos : ScanData object
        with positive scans (sigma+, CR, H- or LH scans, depending on analysis)

    neg : ScanData object 
        with negative scans (sigma-, CL, H+ or LV scans, depending on analysis)

    T : float
        temperature value manually inserted if not provided by data logfiles

    h : float
        field value manually inserted if not provided by data logfiles

    Returns
    -------
    Update and fill attributes of pos and neg ScanData objects

    Pandas DataFrame contains log information (depending on beamline)
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
                                   usecols=confobj.e_scn_cols(f_name))                
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

                # Mean and sign of magnetic field
                h_sgn = np.sign(lgrws['field'])
                # Light polarisation
                pol = lgrws['pol']

                e_raw, dt_raw = escan_split(confobj, data)

                for i in range(len(e_raw)):
                    # Set scan number from filename                
                    scn_num = set_scn_num(confobj, f_name, pos, neg)

                    # Check if magnetic field has changed: the first scan and
                    # the scans after the field has changed are usually rejected
                    # so a different label is provided for these scans.
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
                            continue  # continue if wrong file is found

                        separate_scans(guiobj, confobj, e_raw[i], dt_raw[i],
                            scn_lbl, scn_num, islv, lgrws, pos, neg)

                    else:  # if not XNLD check for circular polarisation
                        try:
                            iscr = confobj.cr_cond(pol)
                        except:
                            guiobj.wrongpol(scn_num, 'circular')
                            continue  # continue if wrong file is found

                        separate_scans(guiobj, confobj, e_raw[i], dt_raw[i],
                            scn_lbl, scn_num, iscr, lgrws, pos, neg)

                    log_tbl = log_tbl.append(lgrws, ignore_index=True)

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

def dt_raw_import(confobj, data):
    '''
    Based on configuration file select data to import depending on sensing.
    It also normalize data by i0 if normalization is not provided in datafile.

    Parameters
    ----------
    confobj : Configuration obj

    data : Pandas DataFrame
        DataFrame containing imported data from input file

    Return
    ------
    array with normalized data based on sensing reported in confobj
    '''
    if confobj.norm_curr:  # Normalize data
        if confobj.sense == 'TEY':
            dt_raw = data[confobj.it_escn]/data[confobj.i0_escn]
        else:
            dt_raw = data[confobj.if_escn]/data[confobj.if0_escn]
    else:  # Normalized data are provided
        if confobj.sense == 'TEY':
            dt_raw = data[confobj.iti0_escn]
        else:
            dt_raw = data[confobj.ifi0_escn]

    return dt_raw

def escan_split(confobj, data):
    '''
    Split energy scan containing multiple scans into single scans.

    Energy data trend in this case is supposed to have a sawtooth shape. To find
    when a scan ends and the next starts escan_split:
    . interpolate energy profile
    . computes the derivative of interpolation and its statistical mode
    . the points where the derivative is different from mode are the end-scan
      points

    If input data contain just one scan, return the single scan.

    Parameters
    ----------
    confobj : Configuration Object
    
    data : Pandas DataFrame
        DataFrame containing imported data from input file

    Returns
    -------
    Set of two list. The first list contains pd.Series with energy vlaues of
    splitted scans. The second lisc contains pd.Series with data values of
    splitted scans.
    If only one scan is present the two list will contain just one element each.
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
    der_mode = stats.mode(e_raw_der_round).mode

    st_idx = 0  # Energy scan start index

    # List with separated scans
    e_raw = []
    dt_raw = []

    for i in range(len(e_raw_der_round)):
        if e_raw_der_round[i] != der_mode:  # Scan ends here        
            stp_idx = i + 1
            e_raw.append(energy_raw.iloc[st_idx : stp_idx])
            dt_raw.append(data_raw.iloc[st_idx : stp_idx])
            
            st_idx = stp_idx
    # Append last scan, or the single scan in case no multiple scans are present
    e_raw.append(energy_raw.iloc[st_idx : ])
    dt_raw.append(data_raw.iloc[st_idx : ])

    return e_raw, dt_raw

def set_scn_num(confobj, f_name, pos, neg):
    '''
    Associate an identifier for each scan.
    This identifier is usually the scan number associated to datafile.
    If the input consists of more than one dataset and different scans coming
    from different set have the same number identifier is modified in order to
    avoid overwiritng data problems during data import.

    Parameters
    ----------
    confobj : configuration obj

    f_name : str
        scan name, used to exctract scan identifier

    pos : ScanData obj
        contains positive scan data

    neg : ScanData obj
        contains negative scan data

    Returns
    -------
    str, identifier for the scan
    '''
    scn_num = confobj.extract_num(f_name)
    scn_num_chk = scn_num

    suffx = 1  # suffix to be added if scan_id already used is found

    # if scan id is already used add suffix and check until no free id is found
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

def separate_scans(guiobj, confobj, e_raw, dt_raw, scn_lbl, scn_num, ispol,
    lgrws, pos, neg):
    '''
    Separate scans and fill positive and negative ScanData objects with raw data
    and labels depending on analysis type

    Parameters
    ----------
    confobj : configuration obj

    guiobj : GUI object
        Provides GUI dialogs

    e_raw : Pandas Series
        Series containing energy values of the scan

    dt_raw : Pandas Series
        Series containing imported data

    scn_lbl : str
        Scan label

    scn_num : str
        Scan number

    ispol : bool
        True for CR and LV polarisations
        False for CL and LH polarisations

    lgrws : dict
        Dictionary with log data collected from data logfile

    pos, neg : ScanData objects
        Collect raw data from positive and negative scans

    Return
    ------
    Fill pos and neg ScanData objects with raw data as well lgrws dictionary
    with scan labels

    '''
    
    if e_raw.iloc[0] > e_raw.iloc[-1]:
        e_raw_to_imp = np.flip(e_raw.to_numpy())
        dt_raw_to_imp = np.flip(dt_raw.to_numpy())
    else:
        e_raw_to_imp = e_raw.to_numpy()
        dt_raw_to_imp = dt_raw.to_numpy()

    rawen = pd.DataFrame()
    rawdt = pd.DataFrame() 

    # For XNLD, data are divided and different labels assigned based on
    # polarization sign.
    if guiobj.analysis in guiobj.type['xnld']:

        scn_lbl += ' H = {:.2f} T'.format(lgrws['field'])

        if ispol:  # LV data
            rawen['E' + scn_num] = e_raw_to_imp
            rawdt[scn_num] = dt_raw_to_imp
            neg.raw_imp = pd.concat([neg.raw_imp, rawen], axis=1)
            neg.raw_imp = pd.concat([neg.raw_imp, rawdt], axis=1)
            neg.label.append(scn_lbl)
            neg.idx.append(scn_num)
            neg.dtype = 'LV'
            lgrws['type'] = 'LV'
        else:  # LH data
            rawen['E' + scn_num] = e_raw_to_imp
            rawdt[scn_num] = dt_raw_to_imp
            pos.raw_imp = pd.concat([pos.raw_imp, rawen], axis=1)
            pos.raw_imp = pd.concat([pos.raw_imp, rawdt], axis=1)
            pos.label.append(scn_lbl)
            pos.idx.append(scn_num)
            pos.dtype = 'LH'
            lgrws['type'] = 'LH'
    # For XMCD, data are divided and different labels assigned based on sigma
    # sign.
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
            rawen['E' + scn_num] = e_raw_to_imp
            rawdt[scn_num] = dt_raw_to_imp
            pos.raw_imp = pd.concat([pos.raw_imp, rawen], axis=1)
            pos.raw_imp = pd.concat([pos.raw_imp, rawdt], axis=1)
            pos.label.append(scn_lbl)
            pos.idx.append(scn_num)
            pos.dtype = 'sigma^+'
        else:  # sigma - data
            rawen['E' + scn_num] = e_raw_to_imp
            rawdt[scn_num] = dt_raw_to_imp
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
            rawen['E' + scn_num] = e_raw_to_imp
            rawdt[scn_num] = dt_raw_to_imp
            pos.raw_imp = pd.concat([pos.raw_imp, rawen], axis=1)
            pos.raw_imp = pd.concat([pos.raw_imp, rawdt], axis=1)
            pos.label.append(scn_lbl)
            pos.idx.append(scn_num)
            pos.dtype = 'CR'
            lgrws['type'] = 'CR'
        else:  # CL data
            rawen['E' + scn_num] = e_raw_to_imp
            rawdt[scn_num] = dt_raw_to_imp
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
            rawen['E' + scn_num] = e_raw_to_imp
            rawdt[scn_num] = dt_raw_to_imp
            pos.raw_imp = pd.concat([pos.raw_imp, rawen], axis=1)
            pos.raw_imp = pd.concat([pos.raw_imp, rawdt], axis=1)
            pos.label.append(scn_lbl)
            pos.idx.append(scn_num)
            pos.dtype = 'H -'
        else:  # H+ data
            rawen['E' + scn_num] = e_raw_to_imp
            rawdt[scn_num] = dt_raw_to_imp
            neg.raw_imp = pd.concat([neg.raw_imp, rawen], axis=1)
            neg.raw_imp = pd.concat([neg.raw_imp, rawdt], axis=1)
            neg.label.append(scn_lbl)
            neg.idx.append(scn_num)
            neg.dtype = 'H +'

def output_fls_escan(guiobj, pos, neg, scanobj):
    '''
    Organize output data and columns names.

    Parameters
    ----------
    guiobj : GUI object
        Provide GUI dialogues

    pos : ScanData object
        ScanData for positive scans (CR, LH, H+ depending on experiment)

    neg : ScanData object
        ScanData for negative scans (CR, LH, H+ depending on experiment)

    scanobj : EngyScan object
        Contains results of X-Ray dichroism computations

    Return
    ------
    Numpy array with data to be saved
    col_nms str with column names
    col_desc str with column description
    '''
    # Collect data
    out_data = np.stack((scanobj.energy, scanobj.energycal, pos.aver,
                         neg.aver, scanobj.xd, scanobj.xd_aver, pos.norm,
                         neg.norm, scanobj.xd_pc, scanobj.pos_corr,
                         scanobj.neg_corr, scanobj.xd_pc_av_ej, pos.norm_int,
                         neg.norm_int, scanobj.xd_pc_int,
                         scanobj.xd_pc_av_ej_int), axis=1)
    if guiobj.infile_ref:
        ref = ' norm by ref'
    else:
        ref = ''
    # Output file column names
    col_nms = 'E_int{} (eV),'.format(ref)
    col_nms += 'E{} (eV),'.format(ref)
    col_nms += '{}_av{} (a.u.),'.format(pos.dtype, ref)
    col_nms += '{}_av{} (a.u.),'.format(neg.dtype, ref)
    col_nms += '{}{} (a.u.),'.format(guiobj.analysis, ref)
    col_nms += 'XAS avgd{} (a.u.),'.format(ref)
    col_nms += '{}oPE{} (a.u.),'.format(pos.dtype, ref)
    col_nms += '{}oPE{} (a.u.),'.format(neg.dtype, ref)
    col_nms += '{}{} (%),'.format(guiobj.analysis, ref)
    col_nms += '{} Corr{} (a.u.),'.format(pos.dtype, ref)
    col_nms += '{} Corr{} (a.u.),'.format(neg.dtype, ref)
    col_nms += '{} XAS avgd{} (%),'.format(guiobj.analysis, ref)
    col_nms += '{}oPE_int{} (a.u.),'.format(pos.dtype, ref)
    col_nms += '{}oPE_int{} (a.u.),'.format(neg.dtype, ref)
    col_nms += '{}_int{} (%),'.format(guiobj.analysis, ref)
    col_nms += '{} XAS avgd_int{} (%)\n'.format(guiobj.analysis, ref)

    # Output file column descriptions
    col_desc = 'Interpolated E,'
    col_desc += 'Calibrated E - offset: {} eV,'.format(scanobj.offset)
    col_desc += 'Averaged {},'.format(pos.dtype)
    col_desc += 'Averaged {},'.format(neg.dtype)
    col_desc += '{}_av - {}_av,'.format(neg.dtype, pos.dtype) 
    if guiobj.analysis in guiobj.type['xnld']:
        col_desc += '{}_av + ((2cos^2 - sin^2) {}_av)) / 3 cos^2,'.format(
            pos.dtype, neg.dtype)
    else:
        col_desc += '({}_av +  {}_av) / 2,'.format(neg.dtype, pos.dtype)
    col_desc += '{}_av/PE_av,'.format(pos.dtype)
    col_desc += '{}_av/PE_av,'.format(neg.dtype)
    if guiobj.analysis in guiobj.type['xnld']:
      
        col_desc += ('Normalized {} 300 * cos^2 *'.format(guiobj.analysis) +
            ' ({0}oPE - {1}oPE) / ({1}_EJNor'.format(neg.dtype, pos.dtype) +
            ' + (2cos^2 - sin^2) * {}_EJNor),'.format(neg.dtype))
    else:
        col_desc += ('Normalized {} 200 * '.format(guiobj.analysis) +
            '({0}oPE - {1}oPE) / ({1}_EJNorm'.format(neg.dtype, pos.dtype) +
            ' + {}_EJNorm),'.format(neg.dtype))
    col_desc += '{0}_av * XASavgdPE_av / {0}PE_av,'.format(pos.dtype)
    col_desc += '{0}_av * XASavgdPE_av / {0}PE_av,'.format(neg.dtype)
    col_desc += 'XAS avgd edge-jump norm,'
    col_desc += '{}_av/PE_int,'.format(pos.dtype)
    col_desc += '{}_av/PE_int,'.format(neg.dtype)
    if guiobj.analysis in guiobj.type['xnld']:
        col_desc += ('Normalized {} 300 * cos^2 *'.format(guiobj.analysis) +
            ' ({0}oPE_int - {1}oPE_int) / ({1}_EJNor_int'.format(neg.dtype,
            pos.dtype) + ' + (2cos^2 - sin^2) * {}_EJNor_int)'.format(
            neg.dtype))
    else:
        col_desc += ('Normalized {} 200 * '.format(guiobj.analysis) +
            '({0}oPE_int - {1}oPE_int) / ({1}_EJNor_int'.format(neg.dtype,
            pos.dtype) + ' + {}_EJNor_int)'.format(neg.dtype))

    return out_data, col_nms, col_desc

def output_plot_escan(guiobj, pos, neg, scanobj, log_dt):
    '''
    Create a plot final figure reporting averagaed positive and negative XAS and
    percentage X-Ray dichroism.

    Parameters
    ----------
    guiobj : GUI object
        Provide GUI dialogues

    pos : ScanData object
        ScanData for positive scans (CR, LH, H+ depending on experiment)

    neg : ScanData object
        ScanData for negative scans (CR, LH, H+ depending on experiment)

    scanobj : EngyScan object
        Contains results of X-Ray dichroism computations

    log_dt : dict
        Collect data for logfile

    Return
    ------
    Two pyplot figure objects one with data coming from computation of X-Ray
    dichroism using edge-jump and interpolated value of edge-jump respectively.
    '''
    if guiobj.infile_ref:
        ref = ' norm by ref'
    else:
        ref = ''

    f1, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
    if guiobj.analysis in guiobj.type['xnld']:
        ax1.plot(scanobj.energycal, neg.aver, color='darkorange',
                 label=neg.dtype)
        ax1.plot(scanobj.energycal, pos.aver, color='darkviolet',
                 label=pos.dtype)
    elif guiobj.analysis in guiobj.type['xmcd']:
        ax1.plot(scanobj.energycal, neg.aver,
                 color='blue', label=r'$\sigma^-$')
        ax1.plot(scanobj.energycal, pos.aver, color='red', label=r'$\sigma^+$')
    else:
        ax1.plot(scanobj.energycal, neg.aver, color='blue', label=neg.dtype)
        ax1.plot(scanobj.energycal, pos.aver, color='red', label=pos.dtype)
    ax1.set_ylabel('XAS{} (a.u.)'.format(ref))
    ax1.legend()
    if guiobj.analysis in guiobj.type['xnld']:
        ax2.plot(scanobj.energycal, scanobj.xd_pc, color='black',
                 label=guiobj.analysis)
        ax2.axhline(y=0, color='darkgray')
    else:
        ax2.plot(scanobj.energycal, scanobj.xd_pc, color='green',
                 label=guiobj.analysis)
        ax2.axhline(y=0, color='black')
    ax2.set_ylabel('{}{} (%)'.format(guiobj.analysis, ref))
    ax2.set_xlabel('E (eV)')
    ax2.legend()
    f1.suptitle('Edge : {:.2f}, PreEdge : {:.2f}, Edge-jump : {:.4f}\n'.format(
        log_dt['exper_edge'], log_dt['setted_pedg'], log_dt['xas_aver_ej']) +
        r'T = {} K, H = {} T, $\theta$ = {}°'.format(log_dt['temp'],
        log_dt['field'], log_dt['angle']))

    f2, (ax3, ax4) = plt.subplots(2, 1, sharex=True)
    if guiobj.analysis in guiobj.type['xnld']:
        ax3.plot(scanobj.energycal, neg.aver, color='darkorange',
                 label=neg.dtype)
        ax3.plot(scanobj.energycal, pos.aver, color='darkviolet',
                 label=pos.dtype)
    elif guiobj.analysis in guiobj.type['xmcd']:
        ax3.plot(scanobj.energycal, neg.aver,
                 color='blue', label=r'$\sigma^-$')
        ax3.plot(scanobj.energycal, pos.aver, color='red', label=r'$\sigma^+$')
    else:
        ax3.plot(scanobj.energycal, neg.aver, color='blue', label=neg.dtype)
        ax3.plot(scanobj.energycal, pos.aver, color='red', label=pos.dtype)
    ax3.set_ylabel('XAS{} (a.u.)'.format(ref))
    ax3.legend()
    if guiobj.analysis in guiobj.type['xnld']:
        ax4.plot(scanobj.energycal, scanobj.xd_pc_int, color='black',
                 label=guiobj.analysis)
        ax4.axhline(y=0, color='darkgray')
    else:
        ax4.plot(scanobj.energycal, scanobj.xd_pc_int, color='green',
                 label=guiobj.analysis)
        ax4.axhline(y=0, color='black')
    ax4.set_ylabel('{}{}_int (%)'.format(guiobj.analysis, ref))
    ax4.set_xlabel('E (eV)')
    ax4.legend()
    f2.suptitle('Edge : {:.2f}, PreEdge : {:.2f},'.format(log_dt['exper_edge'],
        log_dt['setted_pedg']) + ' PostEdge : {:.2f},'.format(
        log_dt['setted_postedg']) + ' Edge-jump : {:.4f}\n'.format(
        log_dt['xas_aver_ej_int']) +
        r'T = {} K, H = {} T, $\theta$ = {}°'.format(log_dt['temp'],
        log_dt['field'], log_dt['angle']))
    plt.show()

    return f1, f2

def save_data_escan(confobj, guiobj, pos, neg, scanobj, log_dt, pos_ref,
    neg_ref, scanobj_ref_norm, log_dt_ref):
    '''
    Save output data, logdata and final plots.

    Parameters
    ----------
    confobj : configuration obj

    guiobj : GUI object
        Provide GUI dialogues

    pos : ScanData object
        ScanData for positive scans (CR, LH, H+ depending on experiment)

    neg : ScanData object
        ScanData for negative scans (CR, LH, H+ depending on experiment)

    scanobj : EngyScan object
        Contains results of X-Ray dichroism computations

    log_dt : dict
        Collect data for logfile

    pos_ref : ScanData object
        ScanData for positive reference scans

    neg_ref : ScanData object
        ScanData for negative reference scans

    scanobj_ref_norm : EngyScan object
        Contains results of X-Ray dichroism computations normalized by reference
        sample

    log_dt_ref : dict
        Collect data of reference and normalized dara for logfile
    '''
    # Not normalized by reference data and plots
    out_data, col_nms, col_desc = output_fls_escan(guiobj, pos, neg, scanobj)
    f1, f2 = output_plot_escan(guiobj, pos, neg, scanobj, log_dt)
    logtxt = confobj.logfl_creator(log_dt)

    if confobj.ref_norm:
        guiobj.infile_ref = True  # For graphs and columns labelling
        out_norm_dt, col_norm_nms, col_norm_desc = output_fls_escan(guiobj,
            pos_ref, neg_ref, scanobj_ref_norm)
        f3, f4 = output_plot_escan(guiobj, pos_ref, neg_ref, scanobj_ref_norm,
            log_dt_ref)
        out_data = np.concatenate((out_data, out_norm_dt), axis=1)
        col_nms = col_nms.rstrip('\n') + ',' + col_norm_nms
        col_desc = col_desc + ',' + col_norm_desc

        logtxt += '\n\n'
        logtxt += 'Logs for reference scans and normalized data by refererence.'
        logtxt += '\n\n'
        logtxt += confobj.logfl_creator(log_dt_ref)

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