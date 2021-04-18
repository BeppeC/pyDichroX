"""
pyDichroX_IO.py

Methods to manage input and output files for XMCD, XNCD, XNLD and hysteresys
data analysis.

Methods
--------
log_scavenger(dataflnm, guiobj)
    Search in logfile for energies, field values, temperatures and sample
    position.

open_import_escan(guiobj, confobj)
    Open input files and import data for energy scan experiments.

logfl_creator(confobj, log_dt):
    Create string with log data to be saved in logfile

output_fls_escan(guiobj, pos, neg, scanobj)
    Save data and graphs for energy scan experiments.
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

import modules.pyDichroX_datatreat as dt


def log_scavenger(dataflnm, guiobj):
    '''
    Search for energies, field values, temperatures and sample position of a
    given datafile in related logfile.

    Parameters
    ----------
    dataflnm : datafile's name.
        The name of logfile is retrieved just changing in .log the extension of
        datafile, following SOLEIL convention.

    guiobj : GUI object
        Provides GUI dialogs.

    Returns
    -------
    dict:
     . mon_en : monocromator energy
     . field : magnetic field value
     . tb1 : sample temperature 1
     . tb2 : sample temperature 2
     . rz : sample rotation angle
     . tx : sample x position
     . tz : sample z position

    If a problem with opening logfile is encountered a message is printed on
    terminal and NaN values are returned.

    Notes
    -----
    Currently works only with Deimos @ Soleil logfiles
    '''
    # data log filename
    logfl = dataflnm.rstrip('txt') + 'log'

    try:
        with open(logfl, 'r', encoding='ISO-8859-1') as fl:
            logtx = fl.read()
            # separate paragraphs in logfile
            parlst = logtx.split('\n\n')
            # search in paragraphs the sections of interest
            for par in parlst:
                if 'Monochromator' in par:
                    # find line with energy and extract energy value
                    for ln in par.split('\n'):
                        if 'energy' in ln:
                            mon_en = float(ln.split(':')[1].strip(' eV'))
                if 'Sample magnetic field' in par:
                    # find line with field and extract field value
                    for ln in par.split('\n'):
                        if 'field ' in ln:
                            # field is logged only for energy scan, field
                            # absolute is recorded
                            field = abs(float(ln.split(':')[1].strip(
                                        ' TESLA')))
                if 'Sample temperature' in par:
                    # find line with sample temperature and extract TB values
                    for ln in par.split('\n'):
                        if '1_#1' in ln:
                            tb1 = float(ln.split('=')[1].strip(' K;'))
                        if '1_#2' in ln:
                            tb2 = float(ln.split('=')[1].strip(' K;'))
                if 'Position' in par:
                    # find lines with sample positions and extract
                    # positions values
                    for ln in par.split('\n'):
                        if 'exp1-mt_rz_#1' in ln:
                            rz = float(ln.split('=')[1].strip(' °;'))
                        if 'exp1-mt_tx_#1' in ln:
                            tx = float(ln.split('=')[1].strip(' mm;'))
                        if 'exp1-mt_tz.2_#1' in ln:
                            tz = float(ln.split('=')[1].strip(' mm;'))
        return {'mon_en': mon_en, 'field': field, 'tb1': tb1, 'tb2': tb2,
                'rz': rz, 'tx': tx, 'tz': tz}
    except:
        f_name = os.path.basename(logfl)
        guiobj.no_log(f_name)

        return {'mon_en': np.nan, 'field': np.nan, 'tb1': np.nan, 'tb2': np.nan,
                'rz': np.nan, 'tx': np.nan, 'tz': np.nan}


def open_import_escan(guiobj, confobj):
    '''
    Open input files for energy scan experiments (XMCD, XNCD, XNXD, XNLD),
    import data and groups them based on absorption coefficient for circularly
    polirized light or based on value of linear polarization for linearly
    polarized light.

    Collect experimental data for compiling output logfile (Currently only for
    Deimos @ Soleil beamline).

    Parameters
    ----------
    guiobj : GUI object
        Provides GUI dialogs.

    confobj : Configurations object

    Returns
    -------
    ScandData object with positive scans (sigma+, CR, H- or LH scans, depending
        on analysis)

    ScandData object with negative scans (sigma-, CL, H+ or LV scans, depending
        on analysis)

    dict, contains log information
     . log_tbl : pd.Dataframes, columns are:
         ...... Only if scan log files are present
         . mon_en : float, monocromator energy
         . field : float, magnetic field value
         . tb1 : float, sample temperature 1
         . tb2 : float, sample temperature 2
         . rz : float, sample rotation angle
         . tx : float, sample x position
         . tz : float, sample z position
         ...........................................
         . scan_num : str, scan number (dummies scan are highlighted)
         . type : str, scan polarization
     . Edge_name : edge name
     . Edge_en : edge energy
     . PreEdge_en : pre-edge energy
     . PostEdge_en : post-edge energy
     . angle : rotation angle of the sample
     . bm_angle : incidence beam angle
     ...... Only if the are no scan log files
     . temp : sample temperature
     . field : magnetic field

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
    '''
    # Set edge and pre-edge energies
    edge = guiobj.chs_edge()

    # dictionary collecting log data
    log_dt = {}

    # If no logfile search in configuaration ask for temperature and field
    if not confobj.log:
        temp, field = guiobj.ask_temp_field()
        log_dt['temp'] = temp
        log_dt['field'] = field

    # Insert experimental rotation angle of the sample
    angle, bm_angle = guiobj.ask_angle()

    # Instantiation of objects for different polarization signs
    # pos will contain sigma+, CR, LH, H- polarization data
    # neg will contain sigma-, CL, LV, H+ polarization data
    pos = dt.ScanData()
    neg = dt.ScanData()

    while True:
        # log table with scan log details
        lgrws = {}
        log_tbl = pd.DataFrame()

        in_sets = guiobj.in_dtst()  # Import data sets

        for dtset in in_sets:
            chk_sgn_h = 1.  # To check a field change

            for file in dtset:
                data = pd.read_csv(file, delim_whitespace=True)

                e_raw = data[confobj.energy]

                # Select data based on sensing in configuration file
                if confobj.sense == 'TEY':
                    dt_raw = data[confobj.iti0_escn]
                else:
                    dt_raw = data[confobj.ifi0_escn]

                # Mean and sign of magnetic field
                h_mean = data[confobj.field_escn].mean()
                h_sgn = np.sign(h_mean)
                # Mean of light phase
                phi_mean = data[confobj.phase_escn].mean()

                # Extract scan number from filename
                f_name = os.path.basename(file)
                scn_num = confobj.extract_num(f_name)

                # Create label to tag scans.

                # Check if magnetic field has changed: the first scan and
                # the scans after the field has changed are usually rejected
                # so a different label is provided for these scans.
                if (dtset.index(file) == 0) or (chk_sgn_h != h_sgn):
                    scn_lbl = scn_num + ' Dummy scan'
                    lgrws['scan_num'] = scn_num + '(D)'
                    chk_sgn_h = h_sgn
                else:
                    scn_lbl = scn_num
                    lgrws['scan_num'] = scn_num

                # For XMCD, data are divided and different labels assigned
                # based on sigma sign.
                if guiobj.case in guiobj.type['xmcd']:
                    if confobj.cr_cond(phi_mean):
                        scn_lbl += ' CR, H = {:.2f} T'.format(h_mean)
                        lgrws['type'] = 'CR'
                    else:
                        scn_lbl += ' CL, H = {:.2f} T'.format(h_mean)
                        lgrws['type'] = 'CL'

                    sigma_sgn = h_sgn * confobj.phi_sgn

                    if sigma_sgn > 0:  # sigma + data
                        pos.raw_imp['E' + scn_num] = e_raw
                        pos.raw_imp[scn_num] = dt_raw
                        pos.label.append(scn_lbl)
                        pos.idx.append(scn_num)
                        pos.dtype = 'sigma^+'
                    else:  # sigma - data
                        neg.raw_imp['E' + scn_num] = e_raw
                        neg.raw_imp[scn_num] = dt_raw
                        neg.label.append(scn_lbl)
                        neg.idx.append(scn_num)
                        neg.dtype = 'sigma^-'
                # For XNCD, data are divided and different labels assigned
                # based on polarizaion sign.
                elif guiobj.case in guiobj.type['xncd']:
                    scn_lbl += ' H = {:.2f} T'.format(h_mean)

                    if confobj.cr_cond(phi_mean):  # CR data
                        pos.raw_imp['E' + scn_num] = e_raw
                        pos.raw_imp[scn_num] = dt_raw
                        pos.label.append(scn_lbl)
                        pos.idx.append(scn_num)
                        pos.dtype = 'CR'
                        lgrws['type'] = 'CR'
                    else:  # CL data
                        neg.raw_imp['E' + scn_num] = e_raw
                        neg.raw_imp[scn_num] = dt_raw
                        neg.label.append(scn_lbl)
                        neg.idx.append(scn_num)
                        neg.dtype = 'CL'
                        lgrws['type'] = 'CL'
                # For XNXD, data are divided and different labels assigned
                # based on magnetic field sign.
                elif guiobj.case in guiobj.type['xnxd']:
                    if confobj.cr_cond(phi_mean):
                        scn_lbl += ' CR, H = {:.2f} T'.format(h_mean)
                        lgrws['type'] = 'CR'
                    else:
                        scn_lbl += ' CL, H = {:.2f} T'.format(h_mean)
                        lgrws['type'] = 'CL'

                    if h_sgn < 0:  # H- data
                        pos.raw_imp['E' + scn_num] = e_raw
                        pos.raw_imp[scn_num] = dt_raw
                        pos.label.append(scn_lbl)
                        pos.idx.append(scn_num)
                        pos.dtype = 'H -'
                    else:  # H+ data
                        neg.raw_imp['E' + scn_num] = e_raw
                        neg.raw_imp[scn_num] = dt_raw
                        neg.label.append(scn_lbl)
                        neg.idx.append(scn_num)
                        neg.dtype = 'H +'
                # For XNLD, data are divided and different labels assigned
                # based on polarization sign.
                elif guiobj.case in guiobj.type['xnld']:
                    scn_lbl += ' H = {:.2f} T'.format(h_mean)

                    if confobj.lv_cond(np.around(phi_mean, 1)):  # LV data
                        neg.raw_imp['E' + scn_num] = e_raw
                        neg.raw_imp[scn_num] = dt_raw
                        neg.label.append(scn_lbl)
                        neg.idx.append(scn_num)
                        neg.dtype = 'LV'
                        lgrws['type'] = 'LV'
                    else:  # LH data
                        pos.raw_imp['E' + scn_num] = e_raw
                        pos.raw_imp[scn_num] = dt_raw
                        pos.label.append(scn_lbl)
                        pos.idx.append(scn_num)
                        pos.dtype = 'LH'
                        lgrws['type'] = 'LH'

                if confobj.log:
                    lgrws.update(log_scavenger(file, guiobj))
                log_tbl = log_tbl.append(lgrws, ignore_index=True)

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

    # dictionary collecting log data
    log_dt['log_tbl'] = log_tbl
    if confobj.log:
        log_dt['temp'] = np.round(
            ((log_tbl['tb1'].mean() + log_tbl['tb2'].mean()) / 2), 1)
        log_dt['field'] = np.round(log_tbl['field'].mean(), 1)
    log_dt['Edge_name'] = edge[0]
    log_dt['Edge_en'] = float(edge[1])
    log_dt['PreEdge_en'] = float(edge[2])
    log_dt['PostEdge_en'] = float(edge[3])
    log_dt['angle'] = angle
    log_dt['bm_angle'] = bm_angle

    return pos, neg, log_dt

def logfl_creator(confobj, log_dt):
    '''
    Create string with log data to be saved in logfile

    Parameters
    ----------
    confobj : Configurations object

    log_dt : dictionary with log data

    Returns
    -------
    str, data formatted to be saved in logfile
    '''
    logtxt = ''
    log_tbl = log_dt['log_tbl']

    if confobj.log:
        logtxt += 'Sample temperature\n'
        logtxt += 'TB1 : {} +/- {} K\n'.format(log_tbl['tb1'].mean(),
                                               log_tbl['tb1'].std())
        logtxt += 'TB2 : {} +/- {} K\n\n'.format(log_tbl['tb2'].mean(),
                                                 log_tbl['tb2'].std())
        logtxt += 'Magnetic field {} +/- {} T\n\n'.format(
                                log_tbl['field'].mean(), log_tbl['field'].std())
        logtxt += 'Sample position\n'
        logtxt += 'Rz : {} +/- {} °\n'.format(log_tbl['rz'].mean(),
                                              log_tbl['rz'].std())
        logtxt += 'Tx : {} +/- {} mm\n'.format(log_tbl['tx'].mean(),
                                               log_tbl['tx'].std())
        logtxt += 'Tz : {} +/- {} mm\n\n'.format(log_tbl['tz'].mean(),
                                                 log_tbl['tz'].std())
    else:
        logtxt += 'Temperature : {} K\n'.format(log_dt['temp'])
        logtxt += 'Magnetic field {} T\n\n'.format(log_dt['field'])

    logtxt += 'Setted angle : {}°\n\n'.format(log_dt['angle'])

    logtxt += 'Input scans\n'
    for i in range(len(log_tbl)):
        logtxt += '{} ({}), '.format(log_tbl['scan_num'].iloc[i],
                                     log_tbl['type'].iloc[i])
    logtxt += '\n\n'
    logtxt += 'Positive selected scans\n'
    for i in log_dt['pos_chs']:
        logtxt += '{}, '.format(i)
    logtxt += '\n\n'
    logtxt += 'Negative selected scans\n'
    for i in log_dt['neg_chs']:
        logtxt += '{}, '.format(i)
    logtxt += '\n\n'
    logtxt += 'Edge used : {}\n'.format(log_dt['Edge_name'])
    logtxt += 'Edge energy used : {} eV - tabulated {} eV\n'.format(
        log_dt['exper_edge'], log_dt['Edge_en'])
    logtxt += 'Pre-edge energy : {} eV\n'.format(log_dt['setted_pedg'])
    logtxt += 'Post-edge energy : {} eV\n'.format(log_dt['setted_postedg'])
    logtxt += 'Recalibration : {}\n'.format(log_dt['recal'])
    logtxt += 'Offset : {} eV'.format(log_dt['offset'])

    logtxt += '\n\n'
    logtxt += 'Edge jump for positive scans : {}\n'.format(log_dt['pos_ej'])
    logtxt += 'Edge jump for positive scans - int.d pre-edge: {}\n'.format(
        log_dt['pos_ej_int'])
    logtxt += 'Edge jump for negative scans : {}\n'.format(log_dt['neg_ej'])
    logtxt += 'Edge jump for negative scans - int.d pre-edge: {}\n'.format(
        log_dt['neg_ej_int'])
    logtxt += 'Edge jump for Avgd XAS spectrum : {}\n'.format(
        log_dt['xas_aver_ej'])
    logtxt += 'Edge jump for Avgd XAS sptectrum - int.d pre-edge: {}\n'.format(
        log_dt['xas_aver_ej_int'])

    return logtxt

def output_fls_escan(confobj, guiobj, pos, neg, scanobj, log_dt):
    '''
    Organize data and save them in output file. Also plot and save graphs.

    Parameters
    ----------
    guiobj : GUI object
        Provide GUI dialogues

    pos : ScanData object
        ScanData for positive scans (CR, LH, H+ depending on experiment)

    neg : ScanData object
        ScanData for negative scans (CR, LH, H+ depending on experiment)

    scanobj : EngyScan object
    '''
    # Collect data
    out_data = np.stack((scanobj.energy, scanobj.energycal, pos.aver,
                        neg.aver, scanobj.xd, scanobj.xd_aver, pos.norm,
                        neg.norm, scanobj.xd_pc, scanobj.pos_corr,
                        scanobj.neg_corr, scanobj.xd_pc_av_ej, pos.norm_int,
                        neg.norm_int, scanobj.xd_pc_int,
                        scanobj.xd_pc_av_ej_int), axis=1)
    # Output file column names
    col_nms = 'E_int (eV),'
    col_nms += 'E (eV),'
    col_nms += '{}_av (a.u.),'.format(pos.dtype)
    col_nms += '{}_av (a.u.),'.format(neg.dtype)
    col_nms += '{} (a.u.),'.format(guiobj.analysis)
    col_nms += 'ISO (a.u.),'
    col_nms += '{}oPE (a.u.),'.format(pos.dtype)
    col_nms += '{}oPE (a.u.),'.format(neg.dtype)
    col_nms += '{} (%),'.format(guiobj.analysis)
    col_nms += '{} Corr (a.u.),'.format(pos.dtype)
    col_nms += '{} Corr (a.u.),'.format(neg.dtype)
    col_nms += '{} Iso (%),'.format(guiobj.analysis)
    col_nms += '{}oPE_int (a.u.),'.format(pos.dtype)
    col_nms += '{}oPE_int (a.u.),'.format(neg.dtype)
    col_nms += '{}_int (%),'.format(guiobj.analysis)
    col_nms += '{} Iso_int (%)\n'.format(guiobj.analysis)

    # Output file column descriptions
    col_desc = 'Interpolated E,'
    col_desc += 'Calibrated E - offset: {} eV,'.format(scanobj.offset)
    col_desc += 'Averaged {},'.format(pos.dtype)
    col_desc += 'Averaged {},'.format(neg.dtype)
    col_desc += '{}_av - {}_av,'.format(neg.dtype, pos.dtype)
    if guiobj.case in guiobj.type['xnld']:
        col_desc += '{}_av + ((2cos^2 - sin^2) {}_av)) / 3 cos^2,'.format(
            pos.dtype, neg.dtype)
    else:
        col_desc += '({}_av +  {}_av) / 2,'.format(neg.dtype, pos.dtype)
    col_desc += '{}_av/PE_av,'.format(pos.dtype)
    col_desc += '{}_av/PE_av,'.format(neg.dtype)
    if guiobj.case in guiobj.type['xnld']:
        col_desc += ('Normalized {} 300 * cos^2 *'.format(guiobj.analysis) +
            ' ({0}oPE - {1}oPE) / ({1}_EJNor'.format(neg.dtype, pos.dtype) +
            ' + (2cos^2 - sin^2) * {}_EJNor),'.format(neg.dtype))
    else:
        col_desc += ('Normalized {} 200 * '.format(guiobj.analysis) +
            '({0}oPE - {1}oPE) / ({1}_EJNorm'.format(neg.dtype, pos.dtype) +
            ' + {}_EJNorm),'.format(neg.dtype))
    col_desc += '{0}_av * IsoPE_av / {0}PE_av,'.format(pos.dtype)
    col_desc += '{0}_av * IsoPE_av / {0}PE_av,'.format(neg.dtype)
    col_desc += 'Iso edge-jump norm,'
    col_desc += '{}_av/PE_int,'.format(pos.dtype)
    col_desc += '{}_av/PE_int,'.format(neg.dtype)
    if guiobj.case in guiobj.type['xnld']:
        col_desc += ('Normalized {} 300 * cos^2 *'.format(guiobj.analysis) +
            ' ({0}oPE_int - {1}oPE_int) / ({1}_EJNor_int'.format(neg.dtype,
            pos.dtype) + ' + (2cos^2 - sin^2) * {}_EJNor_int),'.format(
            neg.dtype))
    else:
        col_desc += ('Normalized {} 200 * '.format(guiobj.analysis) +
            '({0}oPE_int - {1}oPE_int) / ({1}_EJNor_int'.format(neg.dtype,
            pos.dtype) + ' + {}_EJNor_int),'.format(neg.dtype))

    # Final plot
    f1, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
    if guiobj.case in guiobj.type['xnld']:
        ax1.plot(scanobj.energycal, neg.aver, color='darkorange',
            label=neg.dtype)
        ax1.plot(scanobj.energycal, pos.aver, color='darkviolet',
            label=pos.dtype)
    elif guiobj.case in guiobj.type['xmcd']:
        ax1.plot(scanobj.energycal, neg.aver,
                 color='blue', label=r'$\sigma^-$')
        ax1.plot(scanobj.energycal, pos.aver, color='red', label=r'$\sigma^+$')
    else:
        ax1.plot(scanobj.energycal, neg.aver, color='blue', label=neg.dtype)
        ax1.plot(scanobj.energycal, pos.aver, color='red', label=pos.dtype)
    ax1.set_ylabel('XAS (a.u.)')
    ax1.legend()
    if guiobj.case in guiobj.type['xnld']:
        ax2.plot(scanobj.energycal, scanobj.xd_pc, color='black',
            label=guiobj.analysis)
        ax2.axhline(y=0, color='darkgray')
    else:
        ax2.plot(scanobj.energycal, scanobj.xd_pc, color='green',
            label=guiobj.analysis)
        ax2.axhline(y=0, color='black')
    ax2.set_ylabel('{} (%)'.format(guiobj.analysis))
    ax2.set_xlabel('E (eV)')
    ax2.legend()
    f1.suptitle('Edge : {:.2f}, PreEdge : {:.2f}, Edge-jump : {:.4f}\n'.format(
        log_dt['exper_edge'], log_dt['setted_pedg'], log_dt['xas_aver_ej']) +
        r'T = {} K, H = {} T, $\theta$ = {}°'.format(log_dt['temp'],
        log_dt['field'], log_dt['angle']))

    f2, (ax3, ax4)=plt.subplots(2, 1, sharex=True)
    if guiobj.case in guiobj.type['xnld']:
        ax3.plot(scanobj.energycal, neg.aver, color='darkorange',
            label=neg.dtype)
        ax3.plot(scanobj.energycal, pos.aver, color='darkviolet',
            label=pos.dtype)
    elif guiobj.case in guiobj.type['xmcd']:
        ax3.plot(scanobj.energycal, neg.aver,
                 color='blue', label=r'$\sigma^-$')
        ax3.plot(scanobj.energycal, pos.aver, color='red', label=r'$\sigma^+$')
    else:
        ax3.plot(scanobj.energycal, neg.aver, color='blue', label=neg.dtype)
        ax3.plot(scanobj.energycal, pos.aver, color='red', label=pos.dtype)
    ax3.set_ylabel('XAS (a.u.)')
    ax3.legend()
    if guiobj.case in guiobj.type['xnld']:
        ax4.plot(scanobj.energycal, scanobj.xd_pc_int, color='black',
            label=guiobj.analysis)
        ax4.axhline(y=0, color='darkgray')
    else:
        ax4.plot(scanobj.energycal, scanobj.xd_pc_int, color='green',
            label=guiobj.analysis)
        ax4.axhline(y=0, color='black')
    ax4.set_ylabel('{}_int (%)'.format(guiobj.analysis))
    ax4.set_xlabel('E (eV)')
    ax4.legend()
    f2.suptitle('Edge : {:.2f}, PreEdge : {:.2f},'.format(log_dt['exper_edge'],
        log_dt['setted_pedg']) + ' PostEdge : {:.2f},'.format(
        log_dt['setted_postedg']) + ' Edge-jump : {:.4f}\n'.format(
        log_dt['xas_aver_ej_int']) + 
        r'T = {} K, H = {} T, $\theta$ = {}°'.format(log_dt['temp'],
            log_dt['field'], log_dt['angle']))
    plt.show()

    default_nm=('{}_{}_scan_{}-{}_{}K_{}T_{}_{}.dat'.format(log_dt['Edge_name'],
        guiobj.analysis, log_dt['log_tbl']['scan_num'].iloc[0].rstrip('(D)'),
        log_dt['log_tbl']['scan_num'].iloc[-1], log_dt['temp'],
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

        logfl_nm = out_nm.rstrip('dat') + 'log'

        logtxt = logfl_creator(confobj, log_dt)

        with open(logfl_nm, 'w') as fl:
            fl.write(logfl_nm + '\n\n')
            fl.write(logtxt)
