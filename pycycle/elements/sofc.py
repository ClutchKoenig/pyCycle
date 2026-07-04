import openmdao.api as om
import pycycle.api as pyc
from pycycle.api import ThermoAdd
from sofc_api  import SOFCThermoAdd, SegmentSOFC 
from pycycle.flow_in import FlowIn
import numpy as np
from pycycle.constants import CEA_AIR_COMPOSITION

class SOFC(pyc.Element):
    """
    
    """
    def initialize(self):
        self.options.declare('segments', default = ('Segment1'), types=(list,tuple), 
                             desc= 'List of MEA Segments in the Model')
        self.options.declare('N_active',    default = 10,   desc= 'Number of Segments with current Flow and chemical reactions')
        self.options.declare('N_passive',   default = 7,    desc= 'Number of segments without current and chemical reactions')
        self.options.declare('in_composition_anode',    default = )
        self.options.declare('in_composition_cathode',  default = CEA_AIR_COMPOSITION)
        # self.options.declare('statics', default=True, 
        #                      desc='If True, calculate static properties')

        self.options.declare('fuel_type', default='H2', 
                             desc='Type of fuel')
        self.default_des_od_conns = 
        [
            #tbd
        ]
        self.segments: list[str] = []
        self.N_seg_tot = self.options['N_passive'] + self.options['N_active'] + self.options['N_passive']
        self.segments =  [f"segment{i}" for i in range(self.options['N_passive'] + self.options['N_active'] + self.options['N_passive'])]
        self.segments_active =#TODO
        super().initialize()

    def pyc_setup_output_ports(self):
        spec = self.options['thermo_data']

        thermo_data = self.options['thermo_data']

        anode_rxn = SOFCThermoAdd(reaction_type='anode',
                                  spec =  spec,
                                  inflow_composition = self.Fl_I_data['Fl_I_an'])
        
        cathode_rxn = SOFCThermoAdd(reaction_type = 'cathode',
                                    spec = spec,
                                    inflow_composition = self.Fl_I_data['Fl_I_cat'])
        
        self.init_output_flow('Fl_O_an', anode_rxn)
        self.init_output_flow('Fl_O_cat', cathode_rxn)
        
        
    
    def setup(self):
        N_active = self.options['N_active']
        N_passive = self.options['N_passive']

        flow_in_an  = FlowIn('Fl_I_an')
        flow_in_cat = FlowIn('Fl_I_cat') 
        self.add_subsystem('inflow_an', flow_in_an, 
                           promotes = ['Fl_I_an:tot:*', 'Fl_I_an:stat:*'])
        self.add_subsystem('inflow_cat', flow_in_cat,
                           promotes = ['Fl_I_cat:tot:*', 'Fl_I_cat:stat:*'])
        # NOTE: I made a big error in the way the compositions are currently passed into the segements,
        #       or, to put it more precisly, they are not getting passed at all. Every segment gets the 
        #       same composition. 
        # 
        for i, name in enumerate(self.segments):
            if i < (N_passive - 1):
                self.add_subsystem(name, SegmentSOFC(type='passive',
                                                     N_segments             = self.N_seg_tot,
                                                     anode_composition      = self.Fl_I_data['Fl_I_an'],
                                                     cathode_composition    = self.Fl_I_data['Fl_I_cat'],
                                                     spec                   = self.options['thermo_data']))
                
            elif (N_passive - 1) < i <= (N_passive + N_active - 1): 
                self.add_subsystem(name, SegmentSOFC(type='active',
                                                     N_segments             = self.N_seg_tot,
                                                     anode_composition      = self.Fl_I_data['Fl_I_an'],
                                                     cathode_composition    = self.Fl_I_data['Fl_I_cat'],
                                                     spec                   = self.options['thermo_data']),
                                                     promotes_input=['n_cell'])        
            elif i > (N_passive + N_active - 1):
                self.add_subsystem(name, SegmentSOFC(type='passive',
                                                     N_segments             = self.N_seg_tot,
                                                     anode_composition      = self.Fl_I_data['Fl_I_an'],
                                                     cathode_composition    = self.Fl_I_data['Fl_I_cat'],
                                                     spec                   = self.options['thermo_data']))
        for i, name in enumerate(self.segments[:-1]):
            if i == len(self.segments):
                break
            
            self.connect(f'{name}.W_out_A',         f'{self.segments[i + 1]}.W_in_A')
            self.connect(f'{name}.W_out_C',         f'{self.segments[i + 1]}.W_in_C')
            self.connect(f'{name}.h_out_A',         f'{self.segments[i + 1]}.h_in_A')
            self.connect(f'{name}.T_A_out',         f'{self.segments[i + 1]}.T_A_in')
            self.connect(f'{name}.T_C_out',         f'{self.segments[i + 1]}.T_C_in')
            self.connect(f'{name}.x_out_A',         f'{self.segments[i + 1]}.x_in_A')
            self.connect(f'{name}.x_out_C',         f'{self.segments[i + 1]}.x_in_C')

            self.connect(f'{name}.composition_out_A',         f'{self.segments[i + 1]}.composition_in_A')
            self.connect(f'{name}.composition_out_C',         f'{self.segments[i + 1]}.composition_in_C')

            self.connect(f'{name}.n_out_A',          f'{self.segments[i + 1]}.n_in_A')
            self.connect(f'{name}.n_out_C',          f'{self.segments[i + 1]}.n_in_C')

            self.connect(f'{name}.p_A_out',         f'{self.segments[i + 1]}.P_A_in')
            self.connect(f'{name}.p_C_out',         f'{self.segments[i + 1]}.P_C_in')

            self.connect(f'{name}.T_PEN',           f'{self.segments[i + 1]}.T_PEN_left')
            self.connect(f'{name}.T_IC',           f'{self.segments[i + 1]}.T_IC_left')


            self.connect(f'{self.segments[i + 1]}.T_PEN',       f'{name}.T_PEN_right')
            self.connect(f'{self.segments[i + 1]}.T_IC',        f'{name}.T_IC_right')
            


        # Flowstation connections
        # Anode:
        self.connect('Fl_I_an:stat:W',              'segment0.W_in_A')
        self.connect('Fl_I_an:tot:h',               'segment0.h_in_A')
        self.connect('Fl_I_an:tot:T',               'segment0.T_A_in')
        self.connect('Fl_I_an:tot:composition',     'segment0.composition_in_A')
        self.connect('Fl_I_an:tot:P',               'segment0.')
        # Temperature Connections
        self.connect('T_PEN_left_0',                'segment0.T_PEN_left')
        self.connect('T_PEN_right_0',               'segment0.T_PEN_right')# Should be the last segment if I am thinking correctly
        self.connect('T_IC_right_0',                'segment0.T_IC_right')# Should be the last segment if I am thinking correctly
        self.connect('T_IC_left_0',                 'segment0.T_IC_left')
        
        # Cathode:
        self.connect('Fl_I_cat:stat:W',              'segment0.W_in_C')
        self.connect('Fl_I_cat:tot:h',               'segment0.h_in_C')
        self.connect('Fl_I_cat:tot:T',               'segment0.T_C_in')
        self.connect('Fl_I_cat:tot:composition',     'segment0.composition_in_C')
        self.connect('Fl_I_cat:tot:P',               'segment0.')

        self.add_subsystem('i_sum',         om.ExecComp(f'i_sum = {i}' for ))

        self.add_subsystem('Ucell',         om.ExecComp('Ucell=Vcell'))

        balances = self.add_subsystem('balances', om.BalanceComp())
        for i, segment in enumerate(self.segments):
            balances.add_balance(f'{segment}_I', val = 1, units='A', eq_units='V' )
            self.connect(f'balances.{segment}_I',        f'{segment}.I')
            self.connect(f'{segment}.V_cell',           f'balances.lhs:{segment}_I')
            self.connect('U_cell',                      f'balances.rhs:{segment}_I')


        balances.add_balance('U_cell', val= 1, units='V', 
                             eq_units='A',
                             lower=0.3, upper=1.2,
                             use_mult=False, #mult_val=-1
                             )
        self.connect('balances.U_cell',         'U_cell')
        self.connect('i_sum.I_total',           'balances.lhs:U_cell')
        self.connect('I',                       'balances.rhs:U_cell')
        
    

class i_sum(om.ExplicitComponent):
    def initialize(self):
        self.options.declare('N_active', default=10, desc='Number of active Segments')
    
    def setup(self):
        #TODO: finish
