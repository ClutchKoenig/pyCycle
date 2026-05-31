import openmdao.api as om
import pycycle.api as pyc
from pycycle.api import ThermoAdd
from sofc_api  import SOFCThermoAdd, SegmentSOFC 
from pycycle.flow_in import FlowIn
import numpy as np


class SOFC(pyc.Element):
    """
    
    """
    def initalize(self):
        self.options.declare('segments', default = ('Segment1'), types=(list,tuple), 
                             desc= 'List of MEA Segments in the Model')
        self.options.declare('N_seg_active', default= 10, desc= 'Number of Segments with current Flow and chemical reactions')
        self.options.declare('N_seg_passive', default = 7, desc= 'Number of segments without current and chemical reactions')

        # self.options.declare('statics', default=True, 
        #                      desc='If True, calculate static properties')

        self.options.declare('fuel_type', default='H2', 
                             desc='Type of fuel')
        self.default_des_od_conns = 
        [
            #tbd
        ]
        self.segments: list[str] = []
        self.N_seg_tot = self.options['N_passive'] + self.options['N_active'] + self.options['N_passive'])
        self.segments =  [f"segment_{i}" for i in range(self.options['N_passive']) + self.options['N_active'] + self.options['N_passive'])]
        super().initialize()

    def pyc_setup_output_ports(self):
        spec = self.options['thermo_data']

        thermo_data = self.options['thermo_data']

        anode_rxn = SOFCThermoAdd(reaction_type='anode',
                                                   spec =  spec,
                                                   inflow_composition = self.Fl_I_dat['Fl_I_an']
                                                   )
        cathode_rxn = SOFCThermoAdd(reaction_type = 'cathode',
                                    spec = spec,
                                    inflow_composition = self.Fl_I_data['Fl_I_cat'])
        self.init_output_flow('Fl_O_an', anode_rxn)
        self.init_output_flow('Fl_O_cat', cathode_rxn)
        
        
    
    def setup(self):
        N_active = self.options['N_seg_active']
        N_passive = self.options['N_seg_passive']

        flow_in_an = FlowIn('Fl_I_an')
        flow_in_cat = FlowIn('Fl_I_cat') 
        self.add_subsystem('inflow_an', flow_in_an, 
                           promotes=['Fl_I_an:tot:*', 'Fl_I_an:stat:*'])
        self.add_subsystem('inflow_cat', flow_in_cat,
                           promotes=['Fl_I_cat:tot:*', 'Fl_I_cat:stat:*'])
        
        for i, name in enumerate(N_passive + N_active + N_passive):
            if i < (N_passive - 1):
                self.add_subsystem(name, SegmentSOFC(type='passive',
                                                     N_seg = self.N_seg_tot),
                                                     anode_composition= self.Fl_I_data['Fl_I_an'],
                                                     cathode_composition= self.Fl_I_data['Fl_I_cat'],
                                                     spec = self.options['thermo_data'])
                
            elif (N_passive - 1) < i < (N_passive + N_active - 1): 
                self.add_subsystem(name, SegmentSOFC(type='active',
                                                     N_seg = self.N_seg_tot),
                                                     anode_composition= self.Fl_I_data['Fl_I_an'],
                                                     cathode_composition= self.Fl_I_data['Fl_I_cat'],
                                                     spec = self.options['thermo_data'])        
            elif i > (N_passive + N_active - 1):
                self.add_subsystem(name, SegmentSOFC(type='passive',
                                                     N_seg = self.N_seg_tot),
                                                     anode_composition= self.Fl_I_data['Fl_I_an'],
                                                     cathode_composition= self.Fl_I_data['Fl_I_cat'],
                                                     spec = self.options['thermo_data'],
                                                     promotes_input=['n_cell'])
        # Flowstation connections
        self.connect('Fl_I_an:stat:W',              'Segment0.W_in_A')
        self.connect('Fl_I_an:tot:h',               'Segment0.h_in_A')
        self.connect('Fl_I_an:tot:T',               'Segment0.T_A_in')
        self.connect('Fl_I_an:tot:composition',     'Segment0.composition_in_A')
        self.connect('Fl_I_an:tot:P',               'Segment0.')
        
        self.connect('T_PEN_left_0',                'Segment0.T_PEN_left')
        self.connect('T_PEN_right_0',               'Segment0.T_PEN_right')
        self.connect('T_IC_right_0',                'Segment0.T_IC_right')
        self.connect('T_IC_left_0',                 'Segment0.T_IC_left')

        self.add_subsystem('i_sum',         om.ExecComp(f'i_sum = {i}' for ))

        balances = self.add_subsystem('balance', om.BalanceComp())
        balances.add_balance('I', units='A', eq_units='V')
        self.connect('balance.I', )

        balances.add_balance('I')