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
        pos, neg, log_dt = io.open_import_escan(gui, conf)

        # Analize data
        e_scan = dt.EngyScan(pos, neg, gui, log_dt)

        # Show and save results
        io.output_fls_escan(conf, gui, pos, neg, e_scan, log_dt)

        cont = pdxgui.ask_continue()

        if cont == None:
            break
        if not cont:
            break
