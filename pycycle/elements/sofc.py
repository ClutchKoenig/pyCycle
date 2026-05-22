

class SOFC(pyc.Element):
    """
    
    """
    def initalize(self):
        self.options.declare('segments', default = ('Segment1'), types=(list,tuple), 
                             desc= 'List of MEA Segments in the Model')
        self.options.declare('N_seg')
        self.options.declare('statics', default=True, 
                             desc='If True, calculate static properties')
        self.options.declare('fuel_type', default='H2', 
                             desc='Type of fuel')
        self.default_des_od_conns=
        [
            #tbd
        ]
        super().initialize()

    def pyc_setup_output_ports(self):
        thermo_method = self.options['thermo_method']
        thermo_data = self.options['thermo_data']
        fuel_type = self.options['fuel_type']

        self.thermo_add_comp = ThermoAdd(method = thermo_method, mix_mode = 'reactant',
                                         thermo_kwargs = {'spec': thermo_data,
                                                          'inflow_composition': self.Fl_I_data['Fl_I'],
                                                          'mix_composition': fuel_type})
        self.copy_flow(self.thermo_add_comp, 'Fl_O')
