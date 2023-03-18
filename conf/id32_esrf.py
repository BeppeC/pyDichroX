"""
Configurations for ID-32 Beamline at ESRF Synchrotron, Grenoble (France)
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
        Not used for id32.

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
        Collect logfile associated to scandata. Not needed for id32.

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

    data_ex(in_file, out_dir)
        Extract single data scans and log from cumulative data file and
        save them in given directory.
    '''

    def __init__(self):
        '''
        Instatiate object setting all the attributes.
        '''
        self.default_only_ext = '.dat'
        self.default_ext = '*.dat'  # mask
        # Include .spec cumulative files in open file dialogue window
        self.filetypes = ['*.spec']

        # True for interactive mode program execution
        self.interactive = True

        self.sep = '\s+'

        # Set 'TEY' for total electron yeld experiments
        # Set 'Fluo' for total fluorescence experiments
        self.sense = 'TEY'

        # List of of performed analysis
        self.list_analysis = ['XMCD', 'XNCD', 'XNLD', 'XNXD',
                              'XMCD Hysteresis on the fly']

        # ID32 hysteresis scan are collected splitting branches with
        # positive and negative fields
        self.spl_brnch = True

        # Attributes for logfiles - present but not used
        self.scanlog_cnt = 0

        # it/i0 is provided
        self.norm_curr = False

        # id32 provide log information for sample T
        self.ask_for_T = False

        # id32 provide log information for magnetic field
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
            NOT id32 case.

        Return
        ------
        list of column names to be imprted.
        '''
        if guiobj.analysis in guiobj.type['hyst']:
            # Columns for hysteresis on fly
            if guiobj.analysis == 'hyst_fly':
                self.field = 'field'  # magnetic field data
                self.iti0 = 'itratio'  # it/i0 data - TEY

                self.ifi0 = 'ifratio'  # if/if0 data - Fluorescence

                self.time = 'Epoch'  # timestamps

                # Hysteresis scan colums list to be imported
                return [self.field, self.iti0, self.ifi0, self.time]

            else:
                # Columns for hysteresis point by point - yet not
                # considered for this bemaline
                pass
        else:
            # columns for energy scan experiments
            self.energy = 'zap_energyM'  # column with energy data
            self.iti0 = 'zap_itratio'  # it/i0 data - TEY

            self.ifi0 = 'zap_ifratio'  # if/if0 data - Fluorescence

            # They should not be used
            self.i0 = 'zap_i0'  # i0 data - TEY
            self.it = 'zap_it'  # it data - TEY
            self.if1 = 'zap_ifluo'  # if data - Fluorescence
            self.if0 = 'zap_i0'  # if0 data - Fluorescence

            # Energy scan colums list to be imported
            return [self.energy, self.iti0, self.ifi0]

    def cr_cond(self, x):
        '''
        Set condition to discriminate for right and left circular
        polarizations.

        Parameters
        ----------
        x : int polarisation identifier
            CR id = 1
            CL id = -1

        Returns
        -------
        bool, True if CR, False if CL.
        '''
        if x == 1:
            self.phi_sgn = 1
            return True
        elif x == -1:
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
            LV id = 1
            LH id = -1

        Returns
        -------
        bool, True if LV, False if LH.
        '''
        if x == 1:
            return True
        elif x == -1:
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
        scn_num = f_name.lstrip('S_').rstrip(self.default_only_ext)

        return scn_num

    def scanlog_fname(self, guiobj):
        '''
        Collect logfile associated to scandata.
        Not needed at ESRF, just pass.

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
         . tb1 : sample temperature 1
         . tb2 : sample temperature 2
         . t : average of tb1 and tb2
         . rz : sample rotation angle
         . tx : sample x position
         . tz : sample z position
        '''
        # data log filename
        logfl = self.single_lognm(dataflnm)

        # Message for no log presence
        self.nologmess = 'Related datascan will be ignored.'

        try:
            with open(logfl, 'r') as fl:
                for line in fl:
                    sp_line = line.split(':')
                    name = sp_line[0]
                    value = sp_line[1].rstrip('\n')
                    # find polarisation id
                    if name == 'pol':
                        pol = int(value)
                    # find energy value
                    if name == 'energy':
                        mon_en = float(value)
                    # find magnetic field value
                    if name == 'magnet':
                        field = float(value)
                    # find temperature value
                    # id-32 provides two temperatures readings:
                    # tset which is the setted temperature and
                    # corresponds to the heater temperature
                    # tsam which is the temperature measured near the
                    # sample but it is not measured during the scans
                    # acquisition, so it is not considered here
                    if name == 'tset':
                        t = float(value)
                    # find sample positions
                    if name == 'sz':
                        tz = float(value)
                    if name == 'sy':
                        tx = float(value)
                    if name == 'srot':
                        rz = float(value)
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
        logtxt += 'SRot : {} +/- {} °\n'.format(log_tbl['rz'].mean(),
                                                log_tbl['rz'].std())
        logtxt += 'Sy : {} +/- {} mm\n'.format(log_tbl['tx'].mean(),
                                               log_tbl['tx'].std())
        logtxt += 'Sz : {} +/- {} mm\n\n'.format(log_tbl['tz'].mean(),
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
            logtxt += ('Baseline linearly interpolated considering values at'
                       + ' pre-edge and post-edge energies\n\n')
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
        logtxt += 'SRot : {} +/- {} °\n'.format(log_tbl['rz'].mean(),
                                                log_tbl['rz'].std())
        logtxt += 'Sy : {} +/- {} mm\n'.format(log_tbl['tx'].mean(),
                                               log_tbl['tx'].std())
        logtxt += 'Sz : {} +/- {} mm\n\n'.format(log_tbl['tz'].mean(),
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
        logtxt += 'SRot : {} +/- {} °\n'.format(log_tbl['rz'].mean(),
                                                log_tbl['rz'].std())
        logtxt += 'Sy : {} +/- {} mm\n'.format(log_tbl['tx'].mean(),
                                               log_tbl['tx'].std())
        logtxt += 'Sz : {} +/- {} mm\n\n'.format(log_tbl['tz'].mean(),
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

    def data_ex(self, in_file, out_dir):
        '''
        Extract single data scans and log from cumulative data file and
        save them in given directory

        Parameters
        ----------
        if_file : str cumulative data filename
        out_dir : str output directory name
        '''
        # Open input file
        with open(in_file) as fin:
            # loop through the file
            for line in fin:
                if line.split():
                    sline = line.split()
                    if sline[0] == '#S':
                        scn_num = sline[1]
                        # if first scan just open the files
                        if scn_num == '1':
                            # set log and data file names
                            log_fn = out_dir + 'S_' + scn_num + '.log'
                            scn_fn = out_dir + 'S_' + scn_num + '.dat'

                            f_log = open(log_fn, 'w+')
                            f_dat = open(scn_fn, 'w+')
                        # if not first scan first close previous files
                        else:
                            f_log.close()
                            f_dat.close()
                            # set log and data file names
                            log_fn = out_dir + 'S_' + scn_num + '.log'
                            scn_fn = out_dir + 'S_' + scn_num + '.dat'

                            f_log = open(log_fn, 'w+')
                            f_dat = open(scn_fn, 'w+')
                    # P1 line contains the phase
                    elif sline[0] == '#P1':
                        # At ESRF
                        # LH = -0.15
                        # LV = 43
                        # CL ca. -26 (-25 -- -27)
                        # CR ca. 26 (25 -- 27)
                        # P1 value will be L or C depending on analysis
                        # which is chosen by user
                        # LH id = -1
                        # LV id = 1
                        # CL id = -1
                        # CR id = 1
                        fl_id = float(sline[1])
                        if fl_id < 10:
                            f_log.write('pol:-1\n')
                        else:
                            f_log.write('pol:1\n')
                    # P6 line contains energy
                    elif sline[0] == '#P6':
                        f_log.write('energy:{}\n'.format(sline[1]))
                    elif sline[0] == '#P13':
                        f_log.write('magnet:{}\n'.format(sline[4]))
                    # P14 line contains motors positions
                    elif sline[0] == '#P14':
                        f_log.write('sz:{}\n'.format(sline[1]))
                        f_log.write('srot:{}\n'.format(sline[2]))
                        f_log.write('sy:{}\n'.format(sline[3]))
                    # C A: is heater temperature -> setted
                    elif sline[0] == '#C' and sline[1] == 'A:':
                        f_log.write('tset:{}\n'.format(sline[2]))
                    # C D: is sample temperature but it is always 0
                    # during measurements
                    elif sline[0] == '#C' and sline[1] == 'D:':
                        f_log.write('tsam:{}\n'.format(sline[2]))
                    # L line contains data headers
                    elif sline[0] == '#L':
                        # write header of datafile
                        f_dat.write(line.lstrip('#L '))
                    # line not starting with # contains data
                    elif line[0] != '#':
                        f_dat.write(line)
                    else:
                        continue
                else:
                    continue
            f_log.close()
            f_dat.close()
