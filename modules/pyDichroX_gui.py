"""
pyDichroX_GUI.py

Class and methods to provide user interfaces for XMCD, XNCD, XNLD and
hysteresys data analysis.

Classes
--------
GUI : Provides methods to manage GUI.

Methods
-------
ask_continue()
    GUI dialogue to ask for continuing with other analysis.

ask_quit(title, mes):
    GUI dialogue to ask if quit or not the program.

no_config():
    GUI dialogue for no presence of configuration file.

set_config(cfg_list):
    GUI dialogue to select configuration file.
"""

# Copyright (C) Giuseppe Cucinotta.
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

import sys
import easygui as eg
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import modules.pyDichroX_datatreat as dt


class GUI:
    '''
    Provides methods to manage GUI.

    Attributes
    ----------
    type : dict
        Links each type of measurements (keys) to corresponding a_case (values)

    sense : dict
        Links the type of detection (TEY or Fluorescence) to corresponding 
        a_cases

    title : str
        Title of windows GUI

    analysis : str
        Name of analysis type.

    case : int
        Identifier for data analysis.
        . 0 : XMCD
        . 1 : XNCD
        . 2 : XNLD
        . 3 : XNXD
        . 4 : Hysteresis

    interactive : bool
        True for interactive mode false otherwise

    Methods
    -------    
    chs_analysis()
        GUI dialogue to choose the data analysis to perform.

    ipt_fls()
        GUI dialogue to select data input files.

    add_set()
        GUI dialogue to add further sets of XMCD/XNCD/XNLD scans.

    not_enough_fls(a_case, pol)
        GUI warning that not enough files for a given polarization have been
        supplied.    

    def chs_scns(choices)
        GUI dialogue to select the scans to be averaged from a list of labels.

    confirm_choice()
        Provides a GUI to let user confirms his choice and continue or make
        a different choice.

    int_pnts(e_num)
        Provides a GUI to let the user set the number of points used for
        energy scan interpolations.

    ask_temp_field()
        GUI dialogue to ask for experimental temperature and magnetic field.

    ask_angle()
        Provides a GUI to ask for experimental angle of the sample respect the 
        beam.

    chs_edge()
        Provides a GUI to set the edge and pre-edge energies for energy scan
        analysis.

    edge(sel_edg)
        Allows user set the edge and pre-edge energies from experimental data
        in order to perform energy calibration and normalizations.

    outfile_name(self, default_nm):
        GUI dialogue to choose output file name
    '''

    def __init__(self, confobj):
        '''
        At instantiation the title of windows, and the attributes type,
        interactive, log and sense are created.

        Parameters
        ----------
        confobj : configuration object
        '''
        self.type = {'hyst': [4],  # case values for hysteresis
                     'xmcd': [0],  # case values for XMCD
                     'xncd': [1],  # case values for XNCD
                     'xnld': [2],  # case values for XNLD
                     'xnxd': [3]}  # case values for XNXD

        self.title = 'pyDichroX'

        self.interactive = confobj.interactive
        self.log = confobj.log
        self.sense = confobj.sense

    def chs_analysis(self):
        '''
        GUI dialogue to choose the data analysis to perform. It also set window
        title and case attributes.

        Returns
        -------
        set title attribure (str) 
        and case attribute (int):

        . 0 : XMCD
        . 1 : XNCD
        . 2 : XNLD
        . 3 : XNXD
        . 4 : Hysteresis
        '''
        self.title = 'pyDichroX'

        msg = 'Select the analysis you want to perform.'
        #choices = ['XMCD', 'XNCD', 'XNLD', 'XNXD', 'Hysteresis']
        choices = ['XMCD', 'XNCD', 'XNLD', 'XNXD']

        self.title += ' {} analysis'.format(self.sense)

        while True:
            a_choice = eg.choicebox(msg=msg, title=self.title, choices=choices)

            if a_choice is None:
                ask_quit(self.title)
            else:
                break

        self.analysis = a_choice
        a_case = choices.index(a_choice)
        self.case = a_case

        self.title = '{} {} analysis'.format(self.analysis, self.sense)

    def sel_edg_fls(self):
        '''
        GUI dialogue to select or create edge-list file.

        Returns
        -------
        str, the name of the edge-list file.
        '''
        # Current directory and txt files are setted as default.
        msg = ("Choose to open an existen edge list file or create a new one.")
        choices = ("Open a file", "Create a new file")

        default = 'edge files/*.txt'

        while True:
            chs = eg.boolbox(msg, self.title, choices)
            if chs:
                msg2 = "Choose the edge list file"
                f_nm = eg.fileopenbox(msg2, self.title, default=default)
            elif chs is False:
                msg2 = ("A new edge list file will be created.\n" +
                        "Choose the directory and the filename.")
                f_nm = eg.filesavebox(msg2, self.title, default=default)
            elif chs is None:
                ask_quit(self.title, 1)
                continue

            if f_nm:
                break
            else:
                ask_quit(self.title, 1)
                continue

        return f_nm

    def chs_edge(self):
        '''
        GUI dialogue to set the edge and pre-edge energies for energy scan
        analysis.
        The edge and pre-edge energies can be choosen from a file, if present.
        Otherwise the file is created.
        If the needed edge is missing it can be added and saved to the file.

        Returns
        -------
        list (str) containing for the selcted edge:
        - [0] : name of the edge
        - [1] : edge energy
        - [2] : pre-edge energy
        - [3] : post-edge energy

        Notes
        -----
        The file edge_list.txt containing the edges' data must be a
        four-columns csv, with comma ',' separated values. First line
        contains the headers of the columns, namely 'Name', 'Edge Energy',
        'Pre-edge Energy' and 'Post-edge Energy'.
        '''
        if self.case in self.type['hyst']:
            # hysteresis analysis does not require edge file
            return ['not provided', 'not provided', 'not provided',
                    'not provided']
        else:
            # File lisitng edge and pre-edge energies
            edg_filename = self.sel_edg_fls()
            try:
                edg_lst = pd.read_csv(edg_filename, sep=',')
            except:  # if the file is empty
                with open(edg_filename, 'w') as f:
                    # Write file headers with column names
                    f.write('Name,Edge Energy,Pre-edge Energy,Post-edge Energy')
                edg_lst = pd.read_csv(edg_filename, sep=',')

            msg = 'Choose the edge you want to use.'
            # Create list of edges shown to user
            edges = []
            for i in edg_lst['Name']:
                edges.append(i)
            edges.append('Add')

            chsd_edge = eg.choicebox(msg, self.title, edges)

            # Check that a choice has been made
            while True:
                errmsg = ''
                if chsd_edge is None:
                    errmsg = '\nPlease choose an edge or add a new one.'
                if not errmsg:
                    break
                chsd_edge = eg.choicebox(msg + errmsg, self.title, edges)

            # If user select Add, a new edge is added to the edge list
            if chsd_edge == 'Add':
                msg = 'Add a new edge.'
                field_nms = ['Name', 'Edge Energy', 'Pre-edge Energy',
                             'Post-edge Energy']
                field_vals = eg.multenterbox(msg, self.title, field_nms)

                # Check that enetered values are valid
                while True:
                    errmsg = ''
                    for val in field_vals:
                        if not val.strip():
                            errmsg += '\n"{}" is a required field.\n\n'.format(
                                val)
                    try:
                        float(field_vals[1])
                    except:
                        errmsg += ('\nOnly numerical values are accepted for' +
                                   ' Edge energy.\n\n')
                    try:
                        float(field_vals[2])
                    except:
                        errmsg += ('\nOnly numerical values are accepted for' +
                                   ' Pre-Edge energy.\n\n')
                    try:
                        float(field_vals[3])
                    except:
                        errmsg += ('\nOnly numerical values are accepted for' +
                                   ' Post-Edge energy.\n\n')

                    if field_vals[0] in edg_lst['Name']:
                        errmsg += ('\nThere\'s already an Edge named ' +
                                   '{}.\n'.format(field_vals[0]) +
                                   'Please choose another name.')

                    if not errmsg:
                        break

                    field_vals = eg.multenterbox(msg + errmsg, self.title,
                                                 field_nms)

                # Add the new edge to the list
                edg_lst.loc[len(edg_lst)] = field_vals
                edg_lst.to_csv(edg_filename, index=False, mode='w+')

                return field_vals
            else:
                # Select row in DataFrame
                # values.tolist() returns a list of selected rows, each of them
                # is on turn a list of the values in the rows.
                sel_edg = edg_lst[edg_lst['Name'] == chsd_edge]
                sel_edg = sel_edg.values.tolist()[0]
                return sel_edg

    def set_edges(self, sel_edg, exper_edge, x, y1, y2, y2int):
        '''
        Allow setting the edge energy from experimental data, pre-edge energy,
        post-edge energy and energy range for pre-edge average. It also compute
        a linear interpolation considering pre-edge and post-edge values to
        take into account of baseline effects.

        Parameters
        ----------
        sel_edg : list
            Contains data of the edges considered
            sel_edg[0] expected edge energy
            sel_edg[1] pre-edge energy
            sel_edg[2] post-edge energy

        exper_edge : float
            experimental edge energy

        x : array
            energy scale of spectrum

        y1 : array
            xd spectrum data

        y2 : array
            xd_average spectrum data

        y2int : UnivariateSpline
            interpolation of spectrum y2

        Returns
        -------
        list :
        [new_vals[1], new_vals[2], new_vals[4], pe_int, new_vals[5]]
        - [0] exper_edge : experimental edge energy
        - [1] e_pe : pre-edge energy
        - [2] e_poste : post-edge energy
        - [3] pe_wdt : number of points representing the half-width energy
                range used for pre-edge average
        #- [2] pe_wdt : number of points representing the half-width energy
        #        range used for pre-edge average
        #- [3] pe_int : interpolated value of pre-edge
        - [4] cal : bool

        Notes
        -----

        '''
        # Initializes number of points for half-width averaging interval
        # of pre-edge energy to 4.
        pe_wdt = 4

        msg = 'Set values for edge and pre-edge energies.'
        # Field names for enterbox
        fn_exptd = 'Expected edge energy'
        fn_exper = 'Experimental edge energy'
        fn_pe = 'Pre-edge energy'
        fn_pse = 'Post-edge energy'
        fn_pewdt = ('Number of points of half-width interval for ' +
                    'pre-edge average')
        fn_cal = ('Recalbrate energy scale (Y(y) / N(n)) \n' +
                  'Choose \'Y\' if you want recalibrate the energy scale ' +
                  'with the expected edge energy.')
        accepted = ['y', 'yes', 'n', 'no']
        field_nms = [fn_exptd, fn_exper, fn_pe, fn_pse, fn_pewdt, fn_cal]

        init_vals = [sel_edg[0], exper_edge,
                     sel_edg[1], sel_edg[2], pe_wdt, 'N']

        chk_ok = False  # While loop control

        while not chk_ok:
            new_vals = eg.multenterbox(msg, self.title, field_nms, init_vals)

            # Check input fields are ok
            while True:
                errmsg = ''
                if new_vals is None:
                    ask_quit(self.title)
                else:                    
                    for i in range(len(field_nms)):
                        if not new_vals[i].strip():
                            errmsg += ('{} is a required field\n\n'.format(
                                       field_nms[i]))
                            break
                        if i < (len(field_nms) - 1):
                            try:
                                float(new_vals[i])
                            except:
                                errmsg += ('Please insert only numbers ' +
                                           'in {}'.format(field_nms[i]))
                                break
                        # Check that enetered energy values are included in the
                        # energy range
                        if (i < (len(field_nms) - 2)) and (float(new_vals[i])
                                                           <= x[0]):
                            errmsg = + ('Please choose for {}'.format(
                                field_nms[i]) + ' a value included' +
                                ' in the considered energy range')
                        if (i < (len(field_nms) - 2)) and (float(new_vals[i])
                                                           >= x[-1]):
                            errmsg = + ('Please choose for {}'.format(
                                field_nms[i]) + ' a value included' +
                                ' in the considered energy range')
                        if (i == len(field_nms) - 1) and (new_vals[i].lower()
                                                          not in accepted):
                            errmsg += ('Please insert only \'Y\', \'y\', ' +
                                        '\'N\' or \'n\' in Recalbrate field.')
                    if not errmsg:
                        break
                new_vals = eg.multenterbox(msg + errmsg, self.title,
                                           field_nms, init_vals)

            for i in range(len(field_nms) - 2):
                new_vals[i] = float(new_vals[i])
            # Range for pre-edge averaging
            pe_idx = np.argmin(np.abs(x - new_vals[2]))

            lpe_idx = pe_idx - int(new_vals[4])  # Left endpoint range index
            # If average interval extends over energy range shrink it
            if lpe_idx < 0:
                lpe_e = x[0]
                new_vals[4] = pe_idx
            else:
                lpe_e = x[lpe_idx]

            rpe_idx = pe_idx + int(new_vals[4])  # Right endpoint range index
            if rpe_idx >= len(x):
                rpe_e = x[-1]
                new_vals[4] = len(x) - 1 - pe_idx
            else:
                rpe_e = x[rpe_idx]

            # Updates init_vals with the setted values - except expected edge
            # value - for the next loop
            init_vals[1] = new_vals[1]
            init_vals[2] = new_vals[2]
            init_vals[3] = new_vals[3]
            init_vals[4] = new_vals[4]
            init_vals[5] = new_vals[5]

            # Compute linear interpolation of baseline considering pre-edge and
            # post-edge energies

            # Pre-edge and post-edge points on interpolated curve.
            x_int = [init_vals[2], init_vals[3]]
            y_int = [y2int(init_vals[2]), y2int(init_vals[3])]

            pe_int = dt.lin_interpolate(x_int, y_int, new_vals[1])

            # Plot spectrum with position of edge expected and measured
            # together with pre-edges range selected for average
            fig, ax1 = plt.subplots()
            ax1.set_xlabel('E (eV)')
            ax1.set_ylabel(self.analysis + ' (a.u.)',
                           color='black')
            ax1.tick_params(axis='y', labelcolor='black')
            ax1.plot(x, y1, color='black')
            ax1.axvline(x=new_vals[1], color='red', label='Experimental edge')
            ax1.axvline(x=new_vals[0], color='blue', label='Expected edge')
            ax1.axvspan(lpe_e, rpe_e, color='mistyrose')
            ax1.axvline(x=new_vals[2], color='coral', label='Pre-edge energy')
            ax1.axvline(x=new_vals[3], color='plum', label='Post-edge energy')
            
            ax1.axhline(y=0, color='black')
            ax1.legend()
            # Plot averaged xas spectrum in a second y-axis
            ax2 = ax1.twinx()  # second axes that shares the same x-axis
            ax2.set_ylabel('Averaged XAS (a.u.)', color='pink')
            ax2.plot(x, y2, color='pink')
            ax2.plot(x_int, y_int, color='deepskyblue')
            ax2.plot(new_vals[1], pe_int, marker='x', color='deepskyblue',
                     label='Linear baseline and pre-edge interpolated')
            ax2.tick_params(axis='y', labelcolor='pink')
            ax2.legend()

            fig.suptitle(self.title)

            fig.tight_layout()
            plt.show()

            # Loop control
            chk_ok = self.confirm_choice()

        if new_vals[5] in accepted[:2]:
            new_vals[5] = True
        else:
            new_vals[5] = False

        #return [new_vals[1], new_vals[2], new_vals[4], pe_int, new_vals[5]]
        return [new_vals[1], new_vals[2], new_vals[3], new_vals[4], new_vals[5]]

    def ask_temp_field(self):
        '''
        GUI dialogue to ask for experimental temperature and magnetic field.

        Returns
        -------
        str, sample temperature, for log purpose
        str, magnetic field, for log purpose.
        '''
        msg = 'Enter the sample temperature (K).'
        temp = eg.enterbox(msg, self.title)
        if temp is None:
            temp = ' '

        msg = 'Enter the magnetic field value (T).'
        field = temp = eg.enterbox(msg, self.title)
        if field is None:
            field = ' '

        return temp, field

    def ask_angle(self):
        '''
        GUI dialogue to ask for experimental angle of the sample respect to the
        X-Ray beam.

        Returns
        -------
        float, the inserted angle value, for log purpose
        float, the complementary angle of the entered one in radiants.
        '''
        msg = 'Enter the rotation angle of the sample.'
        in_ang = 0

        new_angle = eg.enterbox(msg, self.title, in_ang)

        # Error messages for wrong inputs
        err_empty = '\nPlease provide an angle.\n'
        err_num = '\nOnly numbers are accepted.\n'

        while True:
            errmsg = ''
            if new_angle is None:  # If Cancel is pressed exit program
                ask_quit(self.title)
                errmsg = ' '
            if not new_angle:
                errmsg += err_empty
            else:
                try:
                    float(new_angle)
                except:
                    errmsg += err_num
            if not errmsg:
                break
            new_angle = eg.enterbox(msg + errmsg, self.title, in_ang)

        # Return input angle value (for log purpose) and complementary angle
        # (the angle between X-Ray beam and the sample sufrace normal) in
        # radians for XNLD computations
        return new_angle, np.deg2rad(90 - float(new_angle))

    def sel_ifls(self):
        '''
        GUI dialogue to select input files.

        Returns
        -------
        List (str) with the filenames to be opened.
        For hysteresis analysis returns a two elements nested list, the first
        element being a list with edge scans filenames and the second one a
        list with pre-edge scans filenames.
        '''
        # Current directory and txt files are setted as default.

        # Non hysteresis scans.
        if self.case not in self.type['hyst']:
            msg = ("Choose data files")
            f_nms = eg.fileopenbox(msg=msg, title=self.title, multiple=True,
                                   default='*.txt')
            return f_nms

        # Hysteresis scans.
        else:
            msg = ("Choose edge scans")
            f_edg_nms = eg.fileopenbox(msg=msg, title=self.title,
                                       multiple=True, default='*.txt')
            msg = ("Choose pre-edge scans")
            f_pedg_nms = eg.fileopenbox(msg=msg, title=self.title,
                                        multiple=True, default='*.txt')

            # Pre-edge scans are not mandatory, if Cancel is pressed of windows
            # is closed an empty list is passed
            if f_pedg_nms is None:
                f_pedg_nms = []

            return f_edg_nms, f_pedg_nms

    def add_set(self):
        '''
        GUI dialogue to ask for adding further data set.

        Returns
        -------
        bool, True to add another set of data, False otherwise.
        '''
        msg = ("Dou you want to add another set of data?\n" +
               "(N.B. These data will be averaged with the set of data " +
               "already inserted)")

        add = eg.ynbox(msg, self.title)

        if add is None:
            return False  # If window is closed don't add further scans
        else:
            return add

    def in_dtst(self):
        '''
        GUI dialogue to select input files data sets.

        Returns
        -------
        nested list (str), each element of the list being a list with the
        filenames of a data set:
        [[dtset11, dtset12, ...], [dtset21, dtset22, ...], ...].

        In case of hysteresis each added dataset is represented by a tuple
        whose first element is a list including the filenames of edge scans and
        the second being a list including the filenames of pre-edge scans:
        [([dtset1edg1, dtset1edg2, ...], [dtset1pedg1, dtset1pedg2, ...]),
        ([dtset2edg1, dtset2edg2, ...], [dtset2pedg1, dtset2pedg2, ...]), ...].
        '''
        dataset = []
        # File selection
        while True:
            dataset.append(self.sel_ifls())
            if self.case not in self.type['hyst']:
                if None in dataset:
                    dataset.pop()
                    ask_quit(self.title, 1)
                else:
                    break
            else:
                if None in dataset[0][0]:
                    dataset[0][0].pop()
                    ask_quit(self.title, 2)
                else:
                    break
        # Additional data selection
        while True:
            if self.add_set():
                add = self.sel_ifls()

                if self.case not in self.type['hyst']:
                    if add is None:
                        break
                    else:
                        dataset.append(add)
                else:
                    if add[0] is None:
                        add = ([], add[1])
                    if not(add[0] or add[1]):
                        break
                    dataset.append(add)
            else:
                break

        return dataset

    def not_enough_fls(self, pol):
        '''
        GUI warning that not enough files for a given polarization have been
        supplied.

        Parameters
        ----------
        pol : bool
            True for sigma+, CR and LH
            False for sigma-, CL and LV.

        Returns
        -------
        bool, True to choose again input files, False to quit program.
        '''
        if self.case in self.type['xmcd']:
            if pol:
                data_type = 'sigma +'
            else:
                data_type = 'sigma -'
        elif a_case in self.type['xnld']:
            if pol:
                data_type = 'LH'
            else:
                data_type = 'LV'
        else:
            if pol:
                data_type = 'CR'
            else:
                data_type = 'CL'

        msg = ("You did not provide enough {} ".format(data_type) +
               "files to continue data analysis.\n\n" +
               "Do you want to choose files again or do you want to quit?")

        return eg.ccbox(msg=msg, title=self.title, choices=('Choose again',
                                                'Quit'), cancel_choice='Quit')

    def no_log(self, fname):
        '''
        Message box indicating a missing logfile.

        Parameters
        ----------
        fname : filename related to the missing logfile
        '''
        msg = "Logfile {} not found".format(fname)

        eg.msgbox(msg, self.title)

    def e_num_pnts(self, e_num):
        '''
        GUI dialogue to set the number of points used for energy scan
        interpolations.

        Parameters
        ----------
        e_num : int
            Default number

        Returns
        -------
        int, the number of points for scan interpolation.
        '''
        msg = ("Please enter the number of points you want to consider " +
               "for data interpolations.\n\nDefault value is given by the " +
               "average number of data points.")

        while True:
            num_points = eg.integerbox(msg, self.title, default=e_num,
                                       upperbound=100000)
            if num_points is None:
                ask_quit(self.title)
            else:
                break

        return num_points

    def chs_scns(self, choices):
        '''
        GUI dialogue to select the scans to be averaged from a list of labels.

        Parameters
        ----------
        choices : list (str)
            Labels of scans.

        Returns
        -------
        list (str), the scan numbers of choiced scans.
        '''
        msg = "Select the scans you want to average."
        poss_chs = choices
        choiced = eg.multchoicebox(msg, self.title, poss_chs)

        # Check that user has made a choice
        while True:
            errmsg = ""
            if choiced is None:  # If Cancel is pressed exit program
                ask_quit(self.title)
            if not choiced:
                errmsg += ("\nPlease, choose at least one scan")
            if not errmsg:
                break
            choiced = eg.multchoicebox(msg + errmsg, self.title, poss_chs)

        # Retrieves scan numbers from labels
        for i in range(len(choiced)):
            choiced[i] = choiced[i].split(maxsplit=1)[0]

        return choiced

    def confirm_choice(self):
        '''
        GUI dialogue to confirm a choice and continue or make a different
        choice.

        Returns
        -------
        eg.boolbox obj
        '''
        msg = ("Do you confirm your choice?")
        buttons = ["Yes, continue.", "No, make another choice."]

        # It returns in turn bool True for 'Yes' and False for 'No'
        return eg.boolbox(msg, self.title, buttons)

    def outfile_name(self, default_nm):
        '''
        GUI dialogue to choose output file name

        Parameters
        ----------
        default_nm : str
            default value

        Returns
        -------
        str, file name
        '''
        msg = ''
        fname = eg.filesavebox(msg=None, title=self.title, default=default_nm,
                               filetypes='*.dat')

        return fname
        #if fname is None:
        #    return default_nm
        #else:
        #    return fname

def ask_continue():
    '''
    GUI dialogue to ask for continuing with other analysis.
    '''
    title = 'pyDichroX'
    msg = 'Do you want continue doing other analysis?'

    return eg.ynbox(msg=msg, title=title)

def ask_quit(title, mes=0):
    '''
    GUI dialogue to ask if quit or not the program.
    Quit program if Yes or windows is closed, otherwhise pass.

    Parameters
    ----------
    title : str
        pop-up window's title
    mes : int
        Identifiers for the message
    '''
    id_msg = {0: '',
              1: 'You didn\'t choose any file\n\n',
              2: 'Edge scans are mandatory\n\n'}

    msg = id_msg[mes] + 'Do you really want to quit?'
    quit = eg.ynbox(msg=msg, title=title)

    if quit or quit is None:
        sys.exit(0)
    else:
        pass

def no_config():
    '''
    GUI dialogue for no presence of configuration file.
    Quit program.
    '''
    title = 'pyDichroX'
    msg = 'No configuration file found, pyDichroX will quit.'
    eg.msgbox(msg=msg, title=title)

    sys.exit(0)

def set_config(cfg_list):
    '''
    GUI dialogue to select configuration file.

    Parameters
    ----------
    cfg_list : list
        list of configuration files

    Returns
    -------
    Configuration file choosed
    '''
    title = 'pyDichroX'
    msg = 'Select the configuration file.'

    while True:
        cfg = eg.choicebox(msg=msg, title=title, choices=cfg_list)

        if cfg is None:
            ask_quit(title)
        else:
            break

    return cfg
