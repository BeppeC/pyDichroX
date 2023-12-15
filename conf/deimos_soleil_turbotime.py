"""
Configurations for DEIMOS Beamline at Soleil Synchrotron, Paris (France)
"""

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

    spl_brnch : bool
        True for hysteresis branches splitted
        False otherwise

    scanlog_cnt : int
        counter to run through scanlog_nms.
        Not used for Deimos.

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

    ifi0 : str
        datafile's column name for normalized if/f0 - fluorescence data.

    field : str
        datafile's column name for magnetic field data.

    time : str
        datafile's column name for acquisition timestamps.

    i0 : str
        datafile's column name for i0 data.

    it : str
        datafile's column name for not normalized it - TEY data.
    
    if1 : str
        datafile's column name for not normalized if - fluo data.
    
    if0 : str
        datafile's column name for if0 - fluorescence data.

    phi_sgn : int
        sign assigned to CR (+1) and CL (-1) for the discrimination of
        sigma+ from sigma-.

    Methods
    -------
    scn_cols(f_name='')
        Assign column names for columns to be imported based on beamline
        settings.

    hyst_scn_cols(f_name='')
        Assing column names for hysteresis scans based on beamline
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
        Collect logfile associated to scandata. Not needed for Deimos.

    single_lognm(dataflnm):
        Reconstruct name of datalog file.

    log_scavenger(dataflnm, guiobj)
        Search in logfile for energies, field values, temperatures and
        sample position.

    escan_logfl_creator(log_dt)
        Create string with log data to be saved in logfile for energy
        scans analysis.

    hscan_logfl_creator(log_dt)
        Create string with log data to be saved in logfile for field
        scans analysis.

    ptbypt_logfl_creator(log_dt)
        Create string with log data to be saved in logfile for
        hysteresis point by point analysis.
    '''

    def __init__(self):
        '''
        Instatiate object setting all the attributes.
        '''
        self.default_only_ext = '.txt'
        self.default_ext = '*.txt'  # mask
        self.filetypes = []  #  leave empty - not needed

        # True for interactive mode program execution
        self.interactive = True

        self.sep = '\s+'

        # Set 'TEY' for total electron yeld experiments
        # Set 'Fluo' for total fluorescence experiments
        self.sense = 'TEY'

        # List of of performed analysis
        self.list_analysis = ['XMCD', 'XNCD', 'XNLD', 'XNXD',
                            'XMCD Hysteresis on the fly',
                            'XMCD Hysteresis point by point - time average',
                            'XMCD Hysteresis point by point - time split']

        # at Deimos hysteresis are usually collected with scans covering
        # the whole branch
        self.spl_brnch = False
        
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

    def scn_cols(self, guiobj, f_name):
        '''
        Assign column names for columns to be imported based on beamline 
        settings.

        Parameters
        ----------
        guiobj : GUI object
            Provides GUI dialogs.

        f_name : str
            data filename, some beamline have filename in column names,
            NOT Deimos case.

        Return
        ------
        list of column names to be imprted.
        '''
        if guiobj.analysis in guiobj.type['hyst']:
            # Columns for hysteresis on fly collected with TurboHyst
            if guiobj.analysis == 'hyst_fly':
                self.field = 'data_02'  # magnetic field data
                self.iti0 = 'data_08'  # it/i0 data - TEY

                self.ifi0 = 'data_12'  # if/if0 data - Fluorescence

                self.time = 'abs_time'  # timestamps

                # They should not be used
                self.i0 = 'data_06'  # i0 data - TEY
                self.it = 'data_07'  # it data - TEY
                self.if1 = 'data_11'  # if data - Fluorescence
                self.if0 = 'data_06'  # if0 data - Fluorescence

                # Hysteresis scan colums list to be imported
                return [self.field, self.iti0, self.ifi0, self.time]

            else:
                # Columns for hysteresis point by point collected with
                # TurboTimeScan
                self.field = 'data_09'  # magnetic field data
                self.iti0 = 'data_07'  # it/i0 data - TEY

                # self.phase_hyst = 'data_09'  # phase data

                self.ifi0 = 'data_12'  # if/if0 data - Fluorescence

                self.time = 'abs_time'  # timestamps

                # They should not be used
                self.i0 = 'data_05'  # i0 data - TEY
                self.it = 'data_06'  # it data - TEY
                self.if1 = 'data_11'  # if data - Fluorescence
                self.if0 = 'data_05'  # if0 data - Fluorescence

                # Hysteresis scan colums list to be imported
                return [self.field, self.iti0, self.ifi0, self.time]
        else:
            # columns for energy scan experiments
            self.energy = 'data_01'  # column with energy data
            self.iti0 = 'data_07'  # it/i0 data - TEY

            self.ifi0 = 'data_12'  # if/if0 data - Fluorescence

            # They should not be used
            self.i0 = 'data_05'  # i0 data - TEY
            self.it = 'data_06'  # it data - TEY
            self.if1 = 'data_11'  # if data - Fluorescence
            self.if0 = 'data_05'  # if0 data - Fluorescence

            # Energy scan colums list to be imported
            return [self.energy, self.iti0, self.ifi0]

    def cr_cond(self, x):
        '''
        Set condition to discriminate for right and left circular
        polarizations.

        Parameters
        ----------
        x : int polarisation identifier
            CR id = 4
            CL id = 3

        Returns
        -------
        bool, True if CR, False if CL.
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
        bool, True if LV, False if LH.
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
            filename.

        Returns
        -------
        str, scan-number.
        '''
        scn_num = f_name.lstrip('scan_').rstrip(self.default_only_ext)

        return scn_num

    def scanlog_fname(self, guiobj):
        '''
        Collect logfile associated to scandata.
        Not needed at Deimos, just pass.

        Parameters
        ----------
        guiobj : GUI object
            Provides GUI dialogs.
        '''
        pass

    def single_lognm(self, dataflnm):
        '''
        Reconstruct name of datalog file.
        Used in case of single logfile associated to single datafile,
        so logfile name is associated to scanfile name.

        Parameters
        ----------
        dataflnm : str
            name of datafile associated to logfile.

        Return
        ------
        str, name of logfile associated to dataflnm.
        '''
        return dataflnm.rstrip(self.default_only_ext) + '.log'

    def log_scavenger(self, dataflnm):
        '''
        Search for energies, field values, temperatures and sample
        position of a given datafile in related logfile.

        Parameters
        ----------
        dataflnm : datafile's name.
            The name of logfile is retrieved just changing in .log the
            extension of datafile, following SOLEIL convention.

        guiobj : GUI object
            Provides GUI dialogs.

        Returns
        -------
        dict:
         . mon_en : monocromator energy
         . pol : polarisation identifier
         . field : magnetic field value
         . t : sample temperature
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
                            if ('field ' in ln) and ('TESLA' in ln):
                                field = float(ln.split(':')[1].strip(' TESLA'))
                    if 'Sample temperature' in par:
                        # find line with sample temperature and take TB
                        # values
                        for ln in par.split('\n'):
                            # nov 2023: 
                            # 1_#1 is no more present in Deimos logfiles
                            # if '1_#1' in ln:
                            #    tb1 = float(ln.split('=')[1].strip(' K;'))
                            if '1_#2' in ln:
                            # t once was t2 and the average was returned
                                t = float(ln.split('=')[1].strip(' K;'))
                    if 'Position' in par:
                        # find lines with sample positions and extract
                        # positions values
                        for ln in par.split('\n'):
                            if 'exp1-mt_rz_#2' in ln:
                                rz = float(ln.split('=')[1].strip(' °;'))
                            if 'exp1-mt_tx_#2' in ln:
                                tx = float(ln.split('=')[1].strip(' mm;'))
                            if 'exp1-mt_tz.2_#2' in ln:
                                tz = float(ln.split('=')[1].strip(' mm;'))
            return {'mon_en': mon_en, 'pol': pol, 'field': field, 't': t,
                    'rz': rz, 'tx': tx, 'tz': tz}
        except:            
            raise Exception()

    def escan_logfl_creator(self, guiobj, log_dt):
        '''
        Create string with log data to be saved in logfile for energy
        scans analysis.

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
        logtxt += 'T : {} +/- {} K\n\n'.format(log_tbl['t'].mean(),
                                               log_tbl['t'].std())
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
        logtxt += '\n'
        logtxt += ('Dichroism @ edge energy : {} +/- {} %\n'.format(
                    log_dt['edge_pc'], log_dt['edge_pc_er']))
        logtxt += ('Dichroism @ edge energy int. : {} +/- {} %\n'.format(
                    log_dt['edge_pc_int'], log_dt['edge_pc_int_er']))

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
        str, data formatted to be saved in logfile.
        '''
        logtxt = ''
        log_tbl = log_dt['log_tbl']

        logtxt += 'Sample temperature\n'
        logtxt += 'T : {} +/- {} K\n\n'.format(log_tbl['t'].mean(),
                                               log_tbl['t'].std())
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
        logtxt += 'CR Up selected scans\n'
        for i in log_dt['pos_up_chsn']:
            logtxt += '{}, '.format(i)
        logtxt += '\n\n'
        logtxt += 'CR Down selected scans\n'
        for i in log_dt['pos_dw_chsn']:
            logtxt += '{}, '.format(i)
        logtxt += '\n\n'
        logtxt += 'CL Up selected scans\n'
        for i in log_dt['neg_up_chsn']:
            logtxt += '{}, '.format(i)
        logtxt += '\n\n'
        logtxt += 'CL Down selected scans\n'
        for i in log_dt['neg_dw_chsn']:
            logtxt += '{}, '.format(i)
        logtxt += '\n\n'
        if log_dt['pos_pe_up_chsn']:
            logtxt += 'CR Up pre-edge selected scans\n'
            for i in log_dt['pos_pe_up_chsn']:
                logtxt += '{}, '.format(i)
            logtxt += '\n\n'
        if log_dt['pos_pe_dw_chsn']:
            logtxt += 'CR Down pre-edge selected scans\n'
            for i in log_dt['pos_pe_dw_chsn']:
                logtxt += '{}, '.format(i)
            logtxt += '\n\n'
        if log_dt['neg_pe_up_chsn']:
            logtxt += 'CL Up pre-edge selected scans\n'
            for i in log_dt['neg_pe_up_chsn']:
                logtxt += '{}, '.format(i)
            logtxt += '\n\n'
        if log_dt['neg_pe_dw_chsn']:
            logtxt += 'CL Down pre-edge selected scans\n'
            for i in log_dt['neg_pe_dw_chsn']:
                logtxt += '{}, '.format(i)
            logtxt += '\n\n'
        logtxt += 'Edge used : {}\n'.format(log_dt['Edge_name'])
        logtxt += 'Edge energy used : {} eV - tabulated {} eV\n'.format(
            log_dt['exper_edge'], log_dt['Edge_en'])
        logtxt += 'Pre-edge energy : {} eV\n'.format(log_dt['setted_pedg'])
        
        return logtxt

    def ptbypt_logfl_creator(self, log_dt):
        '''
        Create string with log data to be saved in logfile for
        hysteresis point by point analysis.

        Parameters
        ----------
        log_dt : dictionary with log data.

        Returns
        -------
        str, data formatted to be saved in logfile.
        '''
        logtxt = ''
        log_tbl = log_dt['log_tbl']

        logtxt += 'Sample temperature\n'
        logtxt += 'T : {} +/- {} K\n\n'.format(log_tbl['t'].mean(),
                                               log_tbl['t'].std())
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
        
        logtxt += 'Time window : {} -- {} s\n'.format(log_dt['t_scale'][0],
                                                        log_dt['t_scale'][1])
        logtxt += 'Number of time steps : {}\n'.format(log_dt['t_scale'][2])
        logtxt += 'Time step length : {} s\n\n'.format(log_dt['t_scale'][3])

        logtxt += 'Edge used : {}\n'.format(log_dt['Edge_name'])
        logtxt += 'Edge energy used : {} eV - tabulated {} eV\n'.format(
            log_dt['exper_edge'], log_dt['Edge_en'])
        logtxt += 'Pre-edge energy : {} eV\n'.format(log_dt['setted_pedg'])
        
        return logtxt