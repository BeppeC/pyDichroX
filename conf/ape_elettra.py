# Configurations for APE beamline at Elettra Synchrotron, Trieste (Italy)
import os
import numpy as np
import modules.pyDichroX_gui as pdxgui

class Configuration():
    '''
    Set program configurations based on data provenance

    Attributes
    ----------
    default_only_ext : str
        default datafile exension

    default_ext : str
        default datafile extension (mask)

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

    norm_curr : bool
        True if it/i0 is not provided by data
        False if it/i0 is provided by data

    ask_for_T : bool
        True if no information on sample temperature is provided by datalog file
        False if sample temperature is provided in datalog file

    ask_for_H : bool
        True if no information on magnetic field is provided by datalog file
        False if mangetic field is provided in datalog file
    
    ref_norm : bool
        True if spectra must be normalized by a reference spectrum
        False otherwise
        
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
    e_scn_cols(self, f_name):
        Assing column names for energy scans based on beamline settings.

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
        self.default_only_ext = '.txt'
        self.default_ext = '*.txt'  # mask

        # True for interactive mode program execution
        self.interactive = True

        self.sep = '\t'

        # Set 'TEY' for total electron yeld experiments

        # Currently only TEY is computed for APE @ Elettra
        #-------------------------------------------------
        self.sense = 'TEY'

        # List of of performed analysis
        self.list_analysis = ['XNCD', 'XNLD']

        # Attributes for logfiles
        self.scanlog_nms = []
        self.scanlog_cnt = 0

        # it/i0 is provided
        self.norm_curr = False

        # APE does NOT provide log information for sample T
        self.ask_for_T = True

        # APE does NOT provide log information for magnetic field
        self.ask_for_H = True

        # Normalizaion by reference scans
        self.ref_norm = True

    def e_scn_cols(self, f_name):
        '''
        Assing column names for energy scans based on beamline settings.

        Parameters
        ----------
        f_name : str
            data filename, some beamline have filename in column names,
            APE case.

        Return
        ------
        list of column names to be imprted
        '''
        f_nm = f_name.rstrip(self.default_only_ext)
        self.energy = 'Energy_{}'.format(f_nm)  # column with energy data
        self.iti0_escn = 'Is/Im (#)_1_{}'.format(f_nm)  # it/i0 data - TEY

        # Energy scan colums list to be imported
        return [self.energy, self.iti0_escn]

    def cr_cond(self, x):
        '''
        Set condition to discriminate for right and left circular polarizations.

        Parameters
        ----------
        x : float polarisation identifier
            CR = -21...
            CL = 21...

        Returns
        -------
        bool, True if CR, False if CL
        '''
        x_trunc = np.trunc(x)
        if x_trunc == -21:
            self.phi_sgn = 1
            return True
        elif x_trunc == 21:
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
            LV = -30...
            LH = 0

        Returns
        -------
        bool, True if LV, False if LH
        '''
        x_trunc = np.trunc(x)
        if x_trunc == 0 :
            return False
        elif x_trunc == -30:
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
        scn_num = f_name.rstrip(self.default_only_ext)

        return scn_num

    def scanlog_fname(self, guiobj):
        '''
        Collect logfile associated to scandata.
        Not needed at APE, simply pass

        Parameters
        ----------
        guiobj : GUI object
            Provides GUI dialogs.
        '''
        pass

    def single_lognm(self, dataflnm):
        '''
        Reconstruct name of datalog file.
        Used in case of single logfile associated to single datafile, so logfile
        name is associated to scanfile name.

        Parameters
        ----------
        dataflnm : str
            name of datafile associated to logfile

        Return
        ------
        str, name of logfile associated to dataflnm
        '''
        return dataflnm.rstrip(self.default_only_ext) + '.meta'

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
         . pol : polarisation identifier
         . x : sample x position
         . y : sample y position
         . z : sample z position
        '''
        # data log filename
        logfl = self.single_lognm(dataflnm)

        # Message for no log presence
        self.nologmess = 'Related datascan will be ignored.'

        try:
            with open(logfl, 'r') as fl:
                for ln in fl:
                    if 'FRONTEND/MACHINE/MACHINE_1/PHA092' in ln:
                        pol = float(ln.split('=')[1].rstrip(' mm\n'))
                    if 'x' in ln:
                        x = float(ln.split('=')[1])
                    if 'y' in ln:
                        y = float(ln.split('=')[1])
                    if 'z' in ln:
                        z = float(ln.split('=')[1])

            return {'pol': pol, 'x': x, 'y': y, 'z':z}
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
        logtxt += 'T : {} K\n'.format(log_tbl['t'].mean())
        logtxt += 'Magnetic field {} T\n\n'.format(
            log_tbl['field'].abs().mean())
        
        logtxt += 'Sample position\n'
        logtxt += 'x : {} +/- {}\n'.format(log_tbl['x'].mean(),
            log_tbl['x'].std())
        logtxt += 'y : {} +/- {}\n'.format(log_tbl['y'].mean(),
            log_tbl['y'].std())
        logtxt += 'z : {} +/- {}\n\n'.format(log_tbl['z'].mean(),
            log_tbl['z'].std())
        logtxt += 'Setted angle : {}°\n\n'.format(log_dt['angle'])

        logtxt += 'Input scans\n'
        for i in range(len(log_tbl)):
            logtxt += '{} ({}), '.format(log_tbl['scn_num'].iloc[i],
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