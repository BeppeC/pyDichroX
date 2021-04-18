# Configurations for Deimos beamline at Soleil Synchrotron, Paris (France)

class Configuration():
    '''
    Set program configurations based on data provenance

    Attributes
    ----------
    interactive : bool
        if True interactive mode is setted, if False not interactive mode is
        setted

    log : bool
        Ture for exctract data from scan logfiles

    sense : str
        type of sensing used

    energy : str
        datafile column name for energy data

    iti0_escn : str
        datafile column name for it/i0 data in energy scan experiments

    phase_escn : str
        datafile column name for phase data in energy scan experiments

    field_escn : str
        datafile column name for magnetic field data in energy scan experiments

    ifi0_escn : str
        datafile column name for if/f0 data in energy scan experiments

    field_hyst : str
        datafile column name for magnetic field data in field scan experiments

    iti0_hyst : str
        datafile column name for it/i0 data in field scan experiments

    phase_hyst : str
        datafile column name for phase data in field scan experiments

    ifi0_hyst : str
        datafile column name for if/f0 data in field scan experiments

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

    exctract_num(f_name):
            Exctract scan number from file name.
    '''
    def __init__(self):
        '''
        Instatiate object setting all the attributes
        '''

        # True for interactive mode program execution
        self.interactive = True

        # True for exctract data from scan logfiles
        self.log = True

        # set 'TEY' for total electron yeld experiments
        # set 'Fluo' for total fluorescence experiments
        self.sense = 'TEY'

        # columns assignemt
        self.energy = 'data_01'  # column with energy data
        # columns for energy scans
        self.iti0_escn = 'data_07'  # it/i0 data - TEY
        self.phase_escn = 'data_08'  # phase data
        self.field_escn = 'data_09'  # magnetic field data
        self.ifi0_escn = 'data_12'  # if/if0 data - Fluorescence
        # columns for hysteresis scans
        self.field_hyst = 'data_02'  # magnetic field data
        self.iti0_hyst = 'data_08'  # it/i0 data - TEY
        self.phase_hyst = 'data_09'  # phase data
        self.ifi0_hyst = 'data_12'  # if/if0 data - Fluorescence

    def cr_cond(self, x):
        '''
        Set condition to discriminate for right and left circular polarizations.

        Returns
        -------
        bool, True if CR, False if CL
        '''
        if x > 0:
            self.phi_sgn = 1
            return True
        else:
            self.phi_sgn = -1
            return False


    def lv_cond(self, x):
        '''
        Set condition to discriminate for vertical and horizontal linear
        polarizations.

        Returns
        -------
        bool, True if LV, False if LH
        '''
        if x:
            return True
        else:
            return False

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
        scn_num = f_name.lstrip('scan_').rstrip('.txt')

        return scn_num