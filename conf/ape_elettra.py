"""
Configurations for APE Beamline at Elettra Synchrotron, Trieste (Italy)
"""

import numpy as np

class Configuration():
    '''
    Set program configurations based on data provenance.

    Attributes
    ----------
    default_only_ext : str
        default datafile exension.

    default_ext : str
        default datafile extension (mask).

    filetypes : list of str
        list of masks for other input data files (mainly designed to
        include the single cumulative file from which data must be
        extracted).

    interactive : bool
        True interactive mode is setted
        False not interactive mode is setted.
    
    sep : str
        column separator in datafiles.

    sense : str
        type of sensing used.

    list_analysis : list
        list of supported analysis for the beamline.

    scanlog_nms : list
        list of scanlog filenames associated at different dataset.

    scanlog_cnt : int
        counter to run through scanlog_nms.

    norm_curr : bool
        True if it/i0 is not provided by data
        False if it/i0 is provided by data.

    ask_for_T : bool
        True if no information on sample temperature is provided by
            datalog file
        False if sample temperature is provided in datalog file.

    ask_for_H : bool
        True if no information on magnetic field is provided by datalog
            file
        False if mangetic field is provided in datalog file.
    
    ref_norm : bool
        True if spectra must be normalized by a reference spectrum
        False otherwise.
        
    energy : str
        datafile's column name for energy data.

    iti0 : str
        datafile's column name for normalized it/i0 - TEY data.

    phi_sgn : int
        sign assigned to CR (+1) and CL (-1) for the discrimination of
        sigma+ from sigma-.

    Methods
    -------
    scn_cols(guiobj, f_name):
        Assign column names for columns to be imported based on beamline
        settings.

    cr_cond(x)
        Set condition to discriminate for right and left circular
        polarizations.

    lv_cond(x)
        Set condition to discriminate for vertical and horizontal linear
        polarizations.

    exctract_num(f_name)
        Exctract scan number from file name.

    scanlog_fname():
        Collect logfile associated to scandata.

    def single_lognm(self, dataflnm):
        Reconstruct name of datalog file.

    log_scavenger(dataflnm, guiobj)
        Search in logfile for energies, field values, temperatures and
        sample position.

    escan_logfl_creator(log_dt)
        Create string with log data to be saved in logfile for energy
        scans analysis.

    hscan_logfl_creator(log_dt)
        Create string with log data to be saved in logfile for field
        scans analysis. - NOT YET PROVIDED FOR APE.

    ptbypt_logfl_creator(log_dt)
        Create string with log data to be saved in logfile for
        hysteresis point by point analysis.
        NOT YET PROVIDED FOR APE.
    '''

    def __init__(self):
        '''
        Instatiate object setting all the attributes.
        '''
        # Default file extension
        self.default_only_ext = '.txt'
        self.default_ext = '*.txt'  # mask
        self.filetypes = []  # leave empty - not needed

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

    def scn_cols(self, guiobj, f_name):
        '''
        Assing column names for energy scans based on beamline settings.

        Parameters
        ----------
        guiobj : GUI object
            Provides GUI dialogs.

        f_name : str
            data filename, some beamline have filename in column names,
            APE case.

        Return
        ------
        list of column names to be imprted.
        '''
        if guiobj.analysis in guiobj.type['hyst']:
            # Yet not considered hysteresis analysis for this beamline
            pass
        else:
            f_nm = f_name.rstrip(self.default_only_ext)
            # column with energy data
            self.energy = 'Energy_{}'.format(f_nm)
            # it/i0 data - TEY
            self.iti0 = 'Is/Im (#)_1_{}'.format(f_nm)

            # Energy scan colums list to be imported
            return [self.energy, self.iti0]

    def cr_cond(self, x):
        '''
        Set condition to discriminate for right and left circular
        polarizations.

        Parameters
        ----------
        x : float polarisation identifier
            CR = -21...
            CL = 21...

        Returns
        -------
        bool, True if CR, False if CL.
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
        bool, True if LV, False if LH.
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
            filename.

        Returns
        -------
        str, scan-number.
        '''
        scn_num = f_name.rstrip(self.default_only_ext)

        return scn_num

    def scanlog_fname(self, guiobj):
        '''
        Collect logfile associated to scandata.
        Not needed at APE, simply pass.

        Parameters
        ----------
        guiobj : GUI object
            Provides GUI dialogs.
        '''
        pass

    def single_lognm(self, dataflnm):
        '''
        Reconstruct name of datalog file.
        Used in case of single logfile associated to single datafile, so
        logfile name is associated to scanfile name.

        Parameters
        ----------
        dataflnm : str
            name of datafile associated to logfile.

        Return
        ------
        str, name of logfile associated to dataflnm.
        '''
        return dataflnm.rstrip(self.default_only_ext) + '.meta'

    def log_scavenger(self, dataflnm):
        '''
        Search for energies, field values, temperatures and sample
        position of a given datafile in related logfile.

        Parameters
        ----------
        dataflnm : datafile's name.
            The name of logfile is retrieved just changing in .log the
            extension of datafile, following SOLEIL convention.

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

    def escan_logfl_creator(self, log_dt):
        '''
        Create string with log data to be saved in logfile.

        Parameters
        ----------
        guiobj : GUI object
            Provide GUI dialogues.

        log_dt : dictionary with log data.

        Returns
        -------
        str, data formatted to be saved in logfile.
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

        logtxt += 'Weight energy for positive scans: {} eV\n'.format(
                                                                log_dt['PWEn'])
        logtxt += 'Weight energy for negative scans: {} eV\n\n'.format(
                                                                log_dt['NWEn'])

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
        if guiobj.bsl_int:
            logtxt += 'Baseline interpolated with ArpLS method\n\n'
        else:
            logtxt += ('Baseline linearly interpolated considering values at' +
                    ' pre-edge and post-edge energies\n\n')
        logtxt += 'Edge used : {}\n'.format(log_dt['Edge_name'])
        logtxt += 'Edge energy used : {} eV - tabulated {} eV\n'.format(
            log_dt['exper_edge'], log_dt['Edge_en'])
        logtxt += 'Pre-edge energy : {} eV\n'.format(log_dt['setted_pedg'])
        if guiobj.bsl_int:
            logtxt += 'Smooth parameter : {:f}\n'.format(log_dt['lambda']/1E7)
        else:
            logtxt += 'Post-edge energy : {} eV\n'.format(
                                                    log_dt['setted_postedg'])
        logtxt += 'Recalibration : {}\n'.format(log_dt['recal'])
        logtxt += 'Offset : {} eV'.format(log_dt['offset'])

        logtxt += '\n\n'
        logtxt += 'Edge jump for positive scans : {} +/- {}\n'.format(
            log_dt['pos_ej'], log_dt['pos_ej_er'])
        logtxt += 'Edge jump for positive scans - int.d pre-edge: {} '.format(
            log_dt['pos_ej_int']) + '+/- {}\n'.format(log_dt['pos_ej_int_er'])
        logtxt += 'Edge jump for negative scans : {} +/-{}\n'.format(
            log_dt['neg_ej'], log_dt['neg_ej_er'])
        logtxt += 'Edge jump for negative scans - int.d pre-edge: {} '.format(
            log_dt['neg_ej_int']) + '+/- {}\n'.format(log_dt['neg_ej_int_er'])
        logtxt += 'Edge jump for Avgd XAS spectrum : {} +/- {}\n'.format(
            log_dt['xas_aver_ej'], log_dt['xas_aver_ej_er'])
        logtxt += ('Edge jump for Avgd XAS sptectrum' +
                   ' - int.d pre-edge: {} +/- {}\n'.format(
                    log_dt['xas_aver_ej_int'], log_dt['xas_aver_ej_int_er']))

        return logtxt

    def hscan_logfl_creator(self, log_dt):
        '''
        Create string with log data to be saved in logfile for field
        scans analysis.

        Parameters
        ----------
        log_dt : dictionary with log data.

        Returns
        -------
        Not yet provided for Ape, currently does nothing.
        '''
        pass

    def ptbypt_logfl_creator(self, log_dt):
        '''
        Create string with log data to be saved in logfile for
        hysteresis point by point analysis.

        Parameters
        ----------
        log_dt : dictionary with log data.

        Returns
        -------
        Not yet provided for Ape, currently does nothing.
        '''
        pass