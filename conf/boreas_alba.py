# Configurations for Deimos beamline at Soleil Synchrotron, Paris (France)
import os
import numpy as np
import pandas as pd
import modules.pyDichroX_gui as pdxgui

class Configuration():
    '''
    Set program configurations based on data provenance

    Attributes
    ----------
    default_ext : str
        default file extension

    interactive : bool
        if True interactive mode is setted, if False not interactive mode is
        setted
    
    sep : str
        column separator in datafiles

    sense : str
        type of sensing used

    list_analysis : list
        list of supported analysis for the beamline

    scanlog_nms : list
        list of scanlog filenames associated at different dataset

    scanlog_cnt : int
        counter to run through scanlog_nms

    singlescan : bool
        True if each datafile contains a single scan data
        False if more than one scan is prensent in a single datafile

    norm_curr : bool
        True if it/i0 is not provided by data
        False if it/i0 is provided by data

    energy : str
        datafile column name for energy data

    it_escn : str
        datafile column name for it data in energy scan experiments

    i0_escn : str
        datafile column name for i0 data in energy scan experiments

    phi_sgn : int
        sign assigned to CR (+1) and CL (-1) for the discrimination of sigma+
        from sigma-

    Methods
    -------
    cr_cond(x)
        Set condition to discriminate for right and left circular polarizations.

    lv_cond(x)
        Set condition to discriminate for vertical and horizontal linear
        polarizations.

    exctract_num(f_name)
        Exctract scan number from file name.

    scanlog_fname():
        Collect logfile associated to scandata.

    def single_lognm(self, dataflnm):
        Reconstruct name of datalog file.
        Not used for Boreas @ Alba.

    log_scavenger(dataflnm, guiobj)
        Search in logfile for energies, field values, temperatures and sample
        position.

    logfl_creator(log_dt):
        Create string with log data to be saved in logfile
    '''

    def __init__(self):
        '''
        Instatiate object setting all the attributes
        '''
        # Default file extension
        self.default_ext = '*.dat'

        # True for interactive mode program execution
        self.interactive = True

        self.sep = '\t'

        # Set 'TEY' for total electron yeld experiments

        # Currently only TEY is computed for Boreas @ Alba
        #-------------------------------------------------
        self.sense = 'TEY'

        # List of of performed analysis
        self.list_analysis = ['XMCD', 'XNCD', 'XNLD', 'XNXD']

        # Attributes for logfiles
        self.scanlog_nms = []
        self.scanlog_cnt = 0

        # One scan per datafile
        self.singlescan = True

        # no it/i0 is provided
        self.norm_curr = True

        # Columns assignemt
        # Columns for energy scans
        self.energy = 'energy_mono_corrected'  # column with energy data
        self.it_escn = 'adc2_i3'  # it data - TEY
        self.i0_escn = 'adc2_i2'  # i0 data - TEY

        # Energy scan colums list to be imported
        self.e_scn_cols = [self.energy, self.it_escn, self.i0_escn]

    def cr_cond(self, x):
        '''
        Set condition to discriminate for right and left circular polarizations.

        Parameters
        ----------
        x : float polarisation identifier
            CR = 0.78...
            CL = -0.78...

        Returns
        -------
        bool, True if CR, False if CL
        '''
        x_trunc = np.trunc(x * 100)
        if x_trunc == 78:
            self.phi_sgn = 1
            return True
        elif x_trunc == -78:
            self.phi_sgn = -1
            return False
        else:
            raise Exception()

    def lv_cond(self, x):
        '''
        Set condition to discriminate for vertical and horizontal linear
        polarizations.

        Parameters
        ----------
        x : float polarisation identifier
            LV = 1.57...
            LH = 0

        Returns
        -------
        bool, True if LV, False if LH
        '''
        x_trunc = np.trunc(x)
        if np.trunc(x) == 0 :
            return False
        elif np.trunc(x * 100) == 157:
            return True
        else:
            raise Exception()

    def extract_num(self, f_name):
        '''
        Exctract scan number from file name.

        Parameters
        ----------
        f_name : str
            filename

        Returns
        -------
        str, scan-number
        '''
        scn_num = f_name.rstrip('.dat')

        return scn_num

    def scanlog_fname(self, guiobj):
        '''
        Collect logfile associated to scandata.
        Boreas provide a unique file with datalogs for all the scans in a
        dataset.

        Parameters
        ----------
        guiobj : GUI object
            Provides GUI dialogs.

        Retrun
        ------
        Fill the class attribute scanlog_nms with the file names of datalogs.
        '''
        while True:
            self.scanlog_nms.append(guiobj.ask_logfn())

            if None in self.scanlog_nms:
                self.scanlog_nms.pop()
                pdxgui.ask_quit(guiobj.title, 1)
            else:
                break

    def single_lognm(self, dataflnm):
        '''
        Reconstruct name of datalog file.
        Used in case of single logfile associated to single datafile, so logfile
        name is associated to scanfile name.

        Not used for Boreas @ Alba.

        Parameters
        ----------
        dataflnm : str
            name of datafile associated to logfile

        Return
        ------
        str, name of logfile associated to dataflnm
        '''
        return 'Not used for Boreas.'

    def log_scavenger(self, dataflnm):
        '''
        Search for energies, field values, temperatures and sample position of a
        given datafile in related logfile.

        Parameters
        ----------
        dataflnm : datafile's name.
            The name of logfile is retrieved just changing in .log the extension
            of datafile, following SOLEIL convention.

        Returns
        -------
        dict:
         . mon_en : monocromator energy
         . pol : polarisation identifier
         . field : magnetic field value
         . polarization : id code for different polarisations ()
         . tb1 : sample temperature 1
         . tb2 : sample temperature 2
         . rz : sample rotation angle
         . tx : sample x position
         . tz : sample z position
        '''
        # retrive data log filename
        logfl = self.scanlog_nms[self.scanlog_cnt]

        # Message for no log presence - Not used for Boreas @ Alba
        self.nologmess = ''

        # scan number id from data filename
        scannum = self.extract_num(os.path.basename(dataflnm))

        try:
            with open(logfl, 'r', encoding='ISO-8859-1') as fl:
                logtx = fl.read()
                # separate paragraphs in logfile
                parlst = logtx.split('\n\n')
                # search in paragraphs the section related to current scan
                for par in parlst:
                    if ('#S ' + scannum) in par:
                        # Collects fieldnames preceeded by #O and values
                        # preceeded by #P
                        headnms = []
                        vals = []          
                        for ln in par.split('\n'):
                            if '#O' in ln:
                                line_nm = ln.split('  ')
                                line_nm[0] = line_nm[0].split(' ')[1]
                                headnms.extend(line_nm)
                            if '#P' in ln:
                                vals.extend(ln.split(' ')[1:])          

                        mon_en = float(vals[headnms.index(
                            'ideu71_motor_energy')])
                        pol = float(vals[headnms.index(
                            'ideu71_motor_polarization')])
                        t = float(vals[headnms.index('xmcd_temp_1K')])
                        field = float(vals[headnms.index('magnet_y')])

            return {'mon_en': mon_en, 'pol': pol, 'field': field, 't':t}
        except:            
            raise Exception()

    def logfl_creator(self, log_dt):
        '''
        Create string with log data to be saved in logfile

        Parameters
        ----------
        log_dt : dictionary with log data

        Returns
        -------
        str, data formatted to be saved in logfile
        '''
        logtxt = ''
        log_tbl = log_dt['log_tbl']

        logtxt += 'Sample temperature\n'
        logtxt += 'T : {} +/- {} K\n'.format(log_tbl['t'].mean(),
                                               log_tbl['t'].std())
        logtxt += 'Magnetic field {} +/- {} T\n\n'.format(
            log_tbl['field'].abs().mean(), log_tbl['field'].abs().std())
        
        logtxt += 'Sample position\n'
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
        logtxt += 'Edge jump for positive scans : {}\n'.format(
            log_dt['pos_ej'])
        logtxt += 'Edge jump for positive scans - int.d pre-edge: {}\n'.format(
            log_dt['pos_ej_int'])
        logtxt += 'Edge jump for negative scans : {}\n'.format(
            log_dt['neg_ej'])
        logtxt += 'Edge jump for negative scans - int.d pre-edge: {}\n'.format(
            log_dt['neg_ej_int'])
        logtxt += 'Edge jump for Avgd XAS spectrum : {}\n'.format(
            log_dt['xas_aver_ej'])
        logtxt += ('Edge jump for Avgd XAS sptectrum' +
                   ' - int.d pre-edge: {}\n'.format(log_dt['xas_aver_ej_int']))

        return logtxt