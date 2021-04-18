'''
pyDichroX_cfgmanager

Import beamline file configurations.
'''

# Copyright (C) Giuseppe Cucinotta.
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

import os
import importlib.util
import sys

import modules.pyDichroX_gui as pdxgui


def open_config():
    '''
    Open and read configuration file. If no configuration file is found
    it is created.

    Returns
    -------
    configurations object
    '''
    # look for configuration files
    cfg_files = []
    for file in os.listdir('conf/'):
        if (file.endswith(".py")) and (file != '__init__.py'):
            cfg_files.append(file)

    if cfg_files == []:
        pdxgui.no_config()  # if no config file exit program
    elif len(cfg_files) == 1:  # if one config file import it
        name = cfg_files[0]
        path = 'conf/' + name
        spec = importlib.util.spec_from_file_location(name, path)
        conf = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(conf)

    else:  # if more than one config file ask for choose
        sel_cfg = pdxgui.set_config(cfg_files)
        path = 'conf/' + sel_cfg
        spec = importlib.util.spec_from_file_location(sel_cfg, path)
        conf = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(conf)

    # Create configuration object
    cfg = conf.Configuration()

    return cfg
