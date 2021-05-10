"""
pyDichroX.py

A program to analyse XAS data.

"""

# Copyright (C) Giuseppe Cucinotta.
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

import modules.pyDichroX_IO as io
import modules.pyDichroX_gui as pdxgui
import modules.pyDichroX_datatreat as dt
import modules.pyDichroX_cfgmanager as cfgman


if __name__ == '__main__':

    # Set configurations
    conf = cfgman.open_config()
    
    while True:
        # Create GUI object for data analysis
        gui = pdxgui.GUI(conf)        

        # Analysis type
        gui.chs_analysis()    
        
        # Import data
        pos, neg, log_dt, pos_ref, neg_ref = io.open_import_escan(gui, conf)

        # Create common energy scale
        e_scale = dt.e_scale(gui, pos, neg, log_dt, pos_ref, neg_ref)

        # Analize data
        if conf.ref_norm:
            # First compute analysis of data not normalized by reference
            # Empty ScanData object, data must not be normalized at first
            log_dt_ref_norm = log_dt         
            e_scan = dt.EngyScan(gui, conf, e_scale, pos, neg, log_dt)            
            # Compute analysis normalizing pos and neg ScanData
            gui.infile_ref = True  # For gui messages
            e_scan_norm = dt.EngyScan(gui, conf, e_scale, pos_ref, neg_ref,
                log_dt_ref_norm, pos, neg)
            gui.infile_ref = False 
        else:
            e_scan = dt.EngyScan(gui, conf, e_scale, pos, neg, log_dt)
            e_scan_norm = []
            log_dt_ref_norm = ''

        # Show and save results
        io.save_data_escan(conf, gui, pos, neg, e_scan, log_dt, pos_ref,
            neg_ref, e_scan_norm, log_dt_ref_norm)

        cont = pdxgui.ask_continue()

        if cont == None:
            break
        if not cont:
            break
