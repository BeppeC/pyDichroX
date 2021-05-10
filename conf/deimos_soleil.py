# Configurations for Deimos beamline at Soleil Synchrotron, Paris (France)

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

    scanlog_cnt : int
        counter to run through scanlog_nms
        Not used for Deimos

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

    iti0_escn : str
        datafile column name for it/i0 data in energy scan experiments

    ifi0_escn : str
        datafile column name for if/f0 data in energy scan experiments

    field_hyst : str
        datafile column name for magnetic field data in field scan experiments

    iti0_hyst : str
        datafile column name for it/i0 data in field scan experiments

    ifi0_hyst : str
        datafile column name for if/f0 data in field scan experiments

    phi_sgn : int
        sign assigned to CR (+1) and CL (-1) for the discrimination of sigma+
        from sigma-

    Methods
    -------
    e_scn_cols(f_name='')
        Assing column names for energy scans based on beamline settings.

    hyst_scn_cols(f_name='')
        Assing column names for hysteresis scans based on beamline settings.

    cr_cond(x)
        Set condition to discriminate for right and left circular polarizations.

    lv_cond(x)
        Set condition to discriminate for vertical and horizontal linear
        polarizations.

    exctract_num(f_name)
            Exctract scan number from file name.

    scanlog_fname():
        Collect logfile associated to scandata. Not needed for Deimos.

    single_lognm(dataflnm):
        Reconstruct name of datalog file.

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
        self.default_only_ext = '.txt'
        self.default_ext = '*.txt'  # mask

        # True for interactive mode program execution
        self.interactive = True

        self.sep = '\s+'

        # Set 'TEY' for total electron yeld experiments
        # Set 'Fluo' for total fluorescence experiments
        self.sense = 'TEY'

        # List of of performed analysis
        self.list_analysis = ['XMCD', 'XNCD', 'XNLD', 'XNXD']

        # Attributes for logfiles - present but not used
        self.scanlog_cnt = 0

        # it/i0 is provided
        self.norm_curr = False

        # Deimos provide log information for sample T
        self.ask_for_T = False

        # Deimos provide log information for magnetic field
        self.ask_for_H = False

        # Normalizaion by reference scans
        self.ref_norm = False

    def e_scn_cols(self, f_name):
        '''
        Assing column names for energy scans based on beamline settings.

        Parameters
        ----------
        f_name : str
            data filename, some beamline have filename in column names,
            NOT Deimos case.

        Return
        ------
        list of column names to be imprted
        '''
        self.energy = 'data_01'  # column with energy data
        self.iti0_escn = 'data_07'  # it/i0 data - TEY

        self.ifi0_escn = 'data_12'  # if/if0 data - Fluorescence

        # Energy scan colums list to be imported
        return [self.energy, self.iti0_escn, self.ifi0_escn]

    def hyst_scn_cols(self, f_name):
        '''
        Assing column names for hysteresis scans based on beamline settings.

        Parameters
        ----------
        f_name : str
            data filename, some beamline have filename in column names,
            NOT Deimos case.

        Return
        ------
        list of column names to be imprted
        '''
        self.field_hyst = 'data_02'  # magnetic field data
        self.iti0_hyst = 'data_08'  # it/i0 data - TEY

        # self.phase_hyst = 'data_09'  # phase data

        self.ifi0_hyst = 'data_12'  # if/if0 data - Fluorescence

        # Hysteresis scan colums list to be imported
        return [self.field_hyst, self.iti0_hyst, self.ifi0_hyst]

    def cr_cond(self, x):
        '''
        Set condition to discriminate for right and left circular polarizations.

        Parameters
        ----------
        x : int polarisation identifier
            CR id = 4
            CL id = 3

        Returns
        -------
        bool, True if CR, False if CL
        '''
        if x == 4:
            self.phi_sgn = 1
            return True
        elif x == 3:
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
        x : int polarisation identifier
            LV id = 2
            LH id = 1

        Returns
        -------
        bool, True if LV, False if LH
        '''
        if x == 2:
            return True
        elif x == 1:
            return False
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
        scn_num = f_name.lstrip('scan_').rstrip(self.default_only_ext)

        return scn_num

    def scanlog_fname(self, guiobj):
        '''
        Collect logfile associated to scandata.
        Not needed at Deimos, simply pass

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
        return dataflnm.rstrip(self.default_only_ext) + '.log'

    def log_scavenger(self, dataflnm):
        '''
        Search for energies, field values, temperatures and sample position of a
        given datafile in related logfile.

        Parameters
        ----------
        dataflnm : datafile's name.
            The name of logfile is retrieved just changing in .log the extension
            of datafile, following SOLEIL convention.

        guiobj : GUI object
            Provides GUI dialogs.

        Returns
        -------
        dict:
         . mon_en : monocromator energy
         . pol : polarisation identifier
         . field : magnetic field value
         . tb1 : sample temperature 1
         . tb2 : sample temperature 2
         . rz : sample rotation angle
         . tx : sample x position
         . tz : sample z position
        '''
        # data log filename
        logfl = self.single_lognm(dataflnm)

        # Message for no log presence
        self.nologmess = 'Related datascan will be ignored.'

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
                    if 'Ondulator HU52' in par:
                        # find polarisation id
                        for ln in par.split('\n'):
                            if 'polarisation' in ln:
                                pol = int(ln.split(':')[1])
                    if 'Sample magnetic field' in par:
                        # find line with field and extract field value
                        for ln in par.split('\n'):
                            if 'field ' in ln:
                                field = float(ln.split(':')[1].strip(' TESLA'))
                    if 'Sample temperature' in par:
                        # find line with sample temperature and take TB values
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
            t = (tb1 + tb2) / 2
            return {'mon_en': mon_en, 'pol': pol, 'field': field, 'tb1': tb1,
                    'tb2': tb2, 't': t, 'rz': rz, 'tx': tx, 'tz': tz}
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
        logtxt += 'TB1 : {} +/- {} K\n'.format(log_tbl['tb1'].mean(),
                                               log_tbl['tb1'].std())
        logtxt += 'TB2 : {} +/- {} K\n\n'.format(log_tbl['tb2'].mean(),
                                                 log_tbl['tb2'].std())
        logtxt += 'Magnetic field {} +/- {} T\n\n'.format(
            log_tbl['field'].abs().mean(), log_tbl['field'].abs().std())
        logtxt += 'Sample position\n'
        logtxt += 'Rz : {} +/- {} °\n'.format(log_tbl['rz'].mean(),
                                              log_tbl['rz'].std())
        logtxt += 'Tx : {} +/- {} mm\n'.format(log_tbl['tx'].mean(),
                                               log_tbl['tx'].std())
        logtxt += 'Tz : {} +/- {} mm\n\n'.format(log_tbl['tz'].mean(),
                                                 log_tbl['tz'].std())

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