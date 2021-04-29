"""
pyDichroX_IO.py

Methods to manage input and output files for XMCD, XNCD, XNLD and hysteresys
data analysis.

Methods
--------
open_import_escan(guiobj, confobj)
    Open input files and import data for energy scan experiments.

output_fls_escan(confobj, guiobj, pos, neg, scanobj)
    Save data and graphs for energy scan experiments.

set_scn_num(confobj, f_name, pos, neg)
    Associate an identifier for each scan.
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
         . mon_en : float, monocromator energy
         . pol : int, polarisation identifier
         . field : float, magnetic field value
         . tb1 : float, sample temperature 1
         . tb2 : float, sample temperature 2
         . rz : float, sample rotation angle
         . tx : float, sample x position
         . tz : float, sample z position
         . scan_num : str, scan number (dummies scan are highlighted)
         . type : str, scan polarization
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
    '''
    # Set edge and pre-edge energies
    edge = guiobj.chs_edge()

    # dictionary collecting log data
    log_dt = {}

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

        in_sets = guiobj.in_dtst(confobj)  # Import data sets

        for dtset in in_sets:
            chk_sgn_h = 1.  # To check a field change

            for file in dtset:
                data = pd.read_csv(file, sep=confobj.sep,
                                   usecols=confobj.e_scn_cols)
                f_name = os.path.basename(file)
                try:
                    lgrws.update(confobj.log_scavenger(file))
                except:
                    logfn = confobj.single_lognm(f_name)
                    guiobj.no_log(logfn, confobj.nologmess)
                    continue

                e_raw = data[confobj.energy]

                # Select data based on sensing in configuration file
                # Normalize data if it/i0 is not directly provided
                if confobj.norm_curr:
                    if confobj.sense == 'TEY':
                        dt_raw = data[confobj.it_escn]/data[confobj.i0_escn]
                    else:
                        dt_raw = data[confobj.if_escn]/data[confobj.if0_escn]
                else:
                    if confobj.sense == 'TEY':
                        dt_raw = data[confobj.iti0_escn]
                    else:
                        dt_raw = data[confobj.ifi0_escn]

                # Mean and sign of magnetic field
                field = lgrws['field']
                h_sgn = np.sign(field)
                # Light polarisation
                pol = lgrws['pol']

                # Set scan number from filename                
                scn_num = set_scn_num(confobj, f_name, pos, neg)

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
                    try:
                        iscr = confobj.cr_cond(pol)
                    except:
                        guiobj.wrongpol(scn_num, 'circular')
                        continue  # continue if wrong file is found

                    if iscr:
                        scn_lbl += ' CR, H = {:.2f} T'.format(field)
                        lgrws['type'] = 'CR'
                    else:
                        scn_lbl += ' CL, H = {:.2f} T'.format(field)
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

                    log_tbl = log_tbl.append(lgrws, ignore_index=True)

                # For XNCD, data are divided and different labels assigned
                # based on polarizaion sign.
                elif guiobj.case in guiobj.type['xncd']:
                    try:
                        iscr = confobj.cr_cond(pol)
                    except:
                        guiobj.wrongpol(scn_num, 'circular')
                        continue  # continue if wrong file is found

                    scn_lbl += ' H = {:.2f} T'.format(field)

                    if iscr:  # CR data
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

                    log_tbl = log_tbl.append(lgrws, ignore_index=True)

                # For XNXD, data are divided and different labels assigned
                # based on magnetic field sign.
                elif guiobj.case in guiobj.type['xnxd']:
                    try:
                        iscr = confobj.cr_cond(pol)
                    except:
                        guiobj.wrongpol(scn_num, 'circular')
                        continue  # continue if wrong file is found

                    if iscr:
                        scn_lbl += ' CR, H = {:.2f} T'.format(field)
                        lgrws['type'] = 'CR'
                    else:
                        scn_lbl += ' CL, H = {:.2f} T'.format(field)
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

                    log_tbl = log_tbl.append(lgrws, ignore_index=True)

                # For XNLD, data are divided and different labels assigned
                # based on polarization sign.
                elif guiobj.case in guiobj.type['xnld']:
                    try:
                        islv = confobj.lv_cond(pol)
                    except:
                        guiobj.wrongpol(scn_num, 'linear')
                        continue  # continue if wrong file is found

                    scn_lbl += ' H = {:.2f} T'.format(field)

                    if islv:  # LV data
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

    return pos, neg, log_dt


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

    f2, (ax3, ax4) = plt.subplots(2, 1, sharex=True)
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

    default_nm = ('{}_{}_scan_{}-{}_{}K_{}T_{}_{}.dat'.format(
        log_dt['Edge_name'], guiobj.analysis,
        log_dt['log_tbl']['scan_num'].iloc[0].rstrip('(D)'),
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

        logtxt = confobj.logfl_creator(log_dt)

        with open(logfl_nm, 'w') as fl:
            fl.write(logfl_nm + '\n\n')
            fl.write(logtxt)


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
            scn_num_chk = scan_num + '_{}'.format(suffx)
            suffx += 1
        else:
            break

    while True:
        if scn_num_chk in neg.raw_imp.columns:
            scn_num_chk = scan_num + '_{}'.format(suffx)
            suffx += 1
        else:
            break

    return scn_num_chk