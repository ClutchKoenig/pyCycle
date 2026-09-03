import openmdao.api as om
from pycycle.element_base import Element
from pycycle.sofc_api import SOFCThermoAdd, SegmentSOFC
from pycycle.flow_in import FlowIn
import numpy as np
from pycycle.thermo.thermo import Thermo
from pycycle.thermo.cea.species_data import janaf, Properties
from pycycle.constants import CEA_AIR_COMPOSITION
import pycycle.sofc_api as sofc



class SOFC(Element):
    """
    Full SOFC-Model
    ------------------------------------------------
    Inputs:
        I                   Total current

    """
    def initialize(self):
        # Must run first: Element.initialize() also declares 'thermo_data'
        # (default=False). A later declare() of the same key overwrites the
        # earlier one, so calling super() first lets SOFC's own default
        # (janaf, below) win instead of being clobbered back to False.
        super().initialize()

        self.options.declare('N_active',    default = 10,   desc= 'Number of Segments with current Flow and chemical reactions')
        self.options.declare('N_passive',   default = 7,    desc= 'Number of segments without current and chemical reactions')
        self.options.declare('in_composition_anode',    default = {'H': 1})
        self.options.declare('in_composition_cathode',  default = CEA_AIR_COMPOSITION)
        self.options.declare('thermo_data', default=janaf)
        self.options.declare('T_PEN_BC_left', default=973.0,
                             desc='Hardwired PEN temperature on the open (left) face of segment0 [K]')
        self.options.declare('T_IC_BC_left', default=973.0,
                             desc='Hardwired IC temperature on the open (left) face of segment0 [K]')

        self.options.declare('T_PEN_BC_right', default=973.0,
                             desc='Hardwired PEN temperature on the open (right) face of segment[-1] [K]')
        self.options.declare('T_IC_BC_right', default=973.0,
                             desc='Hardwired IC temperature on the open (right) face of segment[-1] [K]')
        # self.options.declare('statics', default=True, 
        #                      desc='If True, calculate static properties')

        self.options.declare('fuel_type', default='H2', 
                             desc='Type of fuel')
        self.default_des_od_conns = [
            #tbd
        ]

    def _seed_default_inflow_composition(self):
        """
        Fall back to the in_composition_* options when nothing upstream has
        populated Fl_I_data (e.g. SOFC used as a standalone/root element with
        no element connected to its Fl_I_an/Fl_I_cat ports). No-op once
        Fl_I_data has already been set, e.g. by the pyCycle graph wiring.
        """
        self.Fl_I_data.setdefault('Fl_I_an', self.options['in_composition_anode'])
        self.Fl_I_data.setdefault('Fl_I_cat', self.options['in_composition_cathode'])

    def pyc_setup_output_ports(self):
        self._seed_default_inflow_composition()
        spec = self.options['thermo_data']

        thermo_data = self.options['thermo_data']

        anode_rxn = SOFCThermoAdd(reaction_type='anode',
                                  spec =  spec,
                                  inflow_composition = self.Fl_I_data['Fl_I_an'])
        
        cathode_rxn = SOFCThermoAdd(reaction_type = 'cathode',
                                    spec = spec,
                                    inflow_composition = self.Fl_I_data['Fl_I_cat'])
        
        # init_output_flow() only calls .output_port_data() for pyCycle's own
        # ThermoAdd; SOFCThermoAdd doesn't subclass it, so it must be called
        # explicitly here or Fl_O_data ends up holding the component instance.
        self.init_output_flow('Fl_O_an', anode_rxn.output_port_data())
        self.init_output_flow('Fl_O_cat', cathode_rxn.output_port_data())
        
        
    
    def setup(self):
        self.segments: list[str] = []
        self.N_seg_tot = self.options['N_passive'] + self.options['N_active'] + self.options['N_passive']
        self.segments =  [f"segment{i}" for i in range(self.options['N_passive'] + self.options['N_active'] + self.options['N_passive'])]
        self.segments_active =[f"segment{i+self.options['N_passive']}" for i in range(self.options['N_active'])]

        self._seed_default_inflow_composition()
        spec = self.options['thermo_data']
        N_active = self.options['N_active']
        N_passive = self.options['N_passive']

        flow_in_an  = FlowIn(fl_name='Fl_I_an')
        flow_in_cat = FlowIn(fl_name='Fl_I_cat')
        self.add_subsystem('inflow_an', flow_in_an, 
                           promotes = ['Fl_I_an:tot:*', 'Fl_I_an:stat:*'])
        self.add_subsystem('inflow_cat', flow_in_cat,
                           promotes = ['Fl_I_cat:tot:*', 'Fl_I_cat:stat:*'])
        anode_comp = dict(self.Fl_I_data['Fl_I_an'])
        cathode_comp = dict(self.Fl_I_data['Fl_I_cat'])
        # NOTE: I made a big error in the way the compositions are currently passed into the segements,
        #       or, to put it more precisly, they are not getting passed at all. Every segment gets the 
        #       same composition. 
        # 
        # Input Definition for Segment0
        anode_in_props = Properties(spec, init_elements=anode_comp)
        cathode_in_props = Properties(spec, init_elements=cathode_comp)

        anode_thermo_in = Thermo(mode='total_TP', fl_name='Fl_I_an:tot',
                                 method='CEA',
                                 thermo_kwargs={'composition':anode_comp, 'spec':spec})
        cathode_thermo_in = Thermo(mode='total_TP', fl_name='Fl_I_cat:tot',
                                 method='CEA',
                                 thermo_kwargs={'composition':cathode_comp, 'spec':spec})
        
        self.add_subsystem('anode_thermo_in', anode_thermo_in,
                           promotes_inputs=[('T',   'T_A_in'),
                                            ('P',   'p_A_in')],
                           promotes_outputs=[('Fl_I_an:tot:h',  'h_A_in'),
                                             ('Fl_I_an:tot:S',  'S_A_in'),
                                             ('Fl_I_an:tot:Cp', 'Cp_A_in'),
                                             ('Fl_I_an:tot:Cv', 'Cv_A_in')])
        self.add_subsystem('cathode_thermo_in', cathode_thermo_in,
                   promotes_inputs=[('T',   'T_C_in'),
                                    ('P',   'p_C_in')],
                   promotes_outputs=[('Fl_I_cat:tot:h',  'h_C_in'),
                                     ('Fl_I_cat:tot:S',  'S_C_in'),
                                     ('Fl_I_cat:tot:Cp', 'Cp_C_in'),
                                     ('Fl_I_cat:tot:Cv', 'Cv_C_in')])
        
        self.add_subsystem('anode_species_in',
                            sofc.SpeciesFlowCalc(thermo=anode_in_props,
                                                 electrode='anode',
                                                 composition=anode_comp),
                            promotes_inputs=[('W', 'W_in_A')],
                            promotes_outputs=[('x_i',    'x_in_A'),
                                              ('n_flow', 'n_in_A')])
        
        self.add_subsystem('cathode_species_in',
                           sofc.SpeciesFlowCalc(thermo=cathode_in_props,
                                                electrode='cathode',
                                                composition=cathode_comp),
                           promotes_inputs=[('W', 'W_in_C')],
                           promotes_outputs=[('x_i',    'x_in_C'),
                                             ('n_flow', 'n_in_C')])

        self.connect('Fl_I_an:stat:W',   'W_in_A')
        self.connect('Fl_I_cat:stat:W',  'W_in_C')
        self.connect('Fl_I_an:tot:T',    'T_A_in')
        self.connect('Fl_I_an:tot:P',    'p_A_in')
        self.connect('Fl_I_cat:tot:T',   'T_C_in')
        self.connect('Fl_I_cat:tot:P',   'p_C_in')

        for i, name in enumerate(self.segments):
            if i <= (N_passive - 1):
                seg_type = 'passive'
            elif (N_passive - 1) < i <= (N_passive + N_active - 1): 
                seg_type = 'active'
            elif i > (N_passive + N_active - 1):
                seg_type = 'passive'

            seg = SegmentSOFC(type=seg_type,
                            N_segments             = self.N_seg_tot,
                            anode_composition      = anode_comp,
                            cathode_composition    = cathode_comp,
                            spec                   = self.options['thermo_data'])

            if seg_type == 'active': 
                self.add_subsystem(name, seg, promotes_inputs=['n_cell', ('p_A_out', 'p_A'), ('p_C_out', 'p_C')])
            else:
                self.add_subsystem(name, seg, promotes_inputs=[('p_A_out', 'p_A'), ('p_C_out', 'p_C')])

            anode_comp = SOFCThermoAdd(reaction_type='anode', inflow_composition=anode_comp).output_port_data()
            cathode_comp = SOFCThermoAdd(reaction_type='cathode', inflow_composition=cathode_comp).output_port_data()
            
        # Segement0.W_in_A/C 
        for i, name in enumerate(self.segments[:-1]):            
            self.connect(f'{name}.W_out_A',         f'{self.segments[i + 1]}.W_in_A')
            self.connect(f'{name}.W_out_C',         f'{self.segments[i + 1]}.W_in_C')
            self.connect(f'{name}.h_A_out',         f'{self.segments[i + 1]}.h_A_in')
            self.connect(f'{name}.h_C_out',         f'{self.segments[i + 1]}.h_C_in')
            self.connect(f'{name}.T_A_out',         f'{self.segments[i + 1]}.T_A_in')
            self.connect(f'{name}.T_C_out',         f'{self.segments[i + 1]}.T_C_in')
            self.connect(f'{name}.x_out_A',         f'{self.segments[i + 1]}.x_in_A')
            self.connect(f'{name}.x_out_C',         f'{self.segments[i + 1]}.x_in_C')

            self.connect(f'{name}.composition_out_A',         f'{self.segments[i + 1]}.composition_in_A')
            self.connect(f'{name}.composition_out_C',         f'{self.segments[i + 1]}.composition_in_C')

            self.connect(f'{name}.n_out_A',          f'{self.segments[i + 1]}.n_in_A')
            self.connect(f'{name}.n_out_C',          f'{self.segments[i + 1]}.n_in_C')

            # Now connected using promotion across same
            # self.connect(f'{name}.p_A_out',         f'{self.segments[i + 1]}.p_A_out')
            # self.connect(f'{name}.p_C_out',         f'{self.segments[i + 1]}.p_C_out')

            self.connect(f'{name}.T_PEN',           f'{self.segments[i + 1]}.T_PEN_left')
            self.connect(f'{name}.T_IC',           f'{self.segments[i + 1]}.T_IC_left')


            self.connect(f'{self.segments[i + 1]}.T_PEN',       f'{name}.T_PEN_right')
            self.connect(f'{self.segments[i + 1]}.T_IC',        f'{name}.T_IC_right')
            


        # Initial Flowstation connections
        # Anode:
        self.connect('Fl_I_an:stat:W',              'segment0.W_in_A')
        self.connect('h_A_in',                      'segment0.h_A_in')
        self.connect('x_in_A',                      'segment0.x_in_A')
        self.connect('n_in_A',                      'segment0.n_in_A')
        self.connect('Fl_I_an:tot:T',               'segment0.T_A_in')
        self.connect('Fl_I_an:tot:composition',     'segment0.composition_in_A')
        self.connect('Fl_I_an:tot:P',               'segment0.p_A_out')
        # Temperature Connections
        # segment0.T_PEN_right / T_IC_right are already wired from segment1
        # by the inter-segment loop above. Only the open left face (no
        # left-neighbour segment) needs a boundary value.
        self.set_input_defaults('segment0.T_PEN_left', val=self.options['T_PEN_BC_left'], units='K')
        self.set_input_defaults('segment0.T_IC_left',  val=self.options['T_IC_BC_left'],  units='K')
        
        self.set_input_defaults(f'{self.segments[-1]}.T_PEN_right', val=self.options['T_PEN_BC_right'], units='K')
        self.set_input_defaults(f'{self.segments[-1]}.T_IC_right',  val=self.options['T_IC_BC_right'],  units='K')

        # Cathode:
        self.connect('Fl_I_cat:stat:W',              'segment0.W_in_C')
        self.connect('h_C_in',                       'segment0.h_C_in')
        self.connect('x_in_C',                       'segment0.x_in_C')
        self.connect('n_in_C',                       'segment0.n_in_C')
        self.connect('Fl_I_cat:tot:T',               'segment0.T_C_in')
        self.connect('Fl_I_cat:tot:composition',     'segment0.composition_in_C')
        self.connect('Fl_I_cat:tot:P',               'segment0.p_C_out') # only for consistency no channel pressure drop is assumed

        # Output Flowstation Connections
        # Last segment's outlet state -> pyCycle Fl_O_an/Fl_O_cat ports.
        # (mirrors the single-segment pattern documented in SOFCThermoAdd's
        # docstring: Thermo(mode='total_TP') for tot:*, direct connect for
        # stat:W since total_TP has no W input/output.)
        anode_out = Thermo(mode='total_TP', fl_name='Fl_O_an:tot',
                           method='CEA',
                           thermo_kwargs={'composition': self.Fl_O_data['Fl_O_an'],
                                          'spec': spec})
        cathode_out = Thermo(mode='total_TP', fl_name='Fl_O_cat:tot',
                             method='CEA',
                             thermo_kwargs={'composition': self.Fl_O_data['Fl_O_cat'],
                                            'spec': spec})

        self.add_subsystem('anode_out',   anode_out,   promotes_outputs=['Fl_O_an:*'])
        self.add_subsystem('cathode_out', cathode_out, promotes_outputs=['Fl_O_cat:*'])

        self.connect(f'{self.segments[-1]}.composition_out_A', 'anode_out.composition')
        self.connect(f'{self.segments[-1]}.T_A_out',            'anode_out.T')
        self.connect('p_A',                                     'anode_out.P')
        self.connect(f'{self.segments[-1]}.W_out_A',            'Fl_O_an:stat:W')

        self.connect(f'{self.segments[-1]}.composition_out_C', 'cathode_out.composition')
        self.connect(f'{self.segments[-1]}.T_C_out',            'cathode_out.T')
        self.connect('p_C',                                     'cathode_out.P')
        self.connect(f'{self.segments[-1]}.W_out_C',            'Fl_O_cat:stat:W')

        # Equipotential condition: every active segment sits on the same
        # electrical interconnect, so they all share one terminal voltage
        # (balances.U_cell). Each segment's current is unknown, solved so
        # that segment.V_cell matches that shared voltage.
        segment_current_inputs = {f'I_{segment}': {'units': 'A'} for segment in self.segments_active}
        i_sum_expr = 'I_total = ' + ' + '.join(segment_current_inputs)
        self.add_subsystem('i_sum', om.ExecComp(i_sum_expr, I_total={'units': 'A'}, **segment_current_inputs))

        balances = self.add_subsystem('balances', om.BalanceComp(), promotes_inputs=[('rhs:U_cell', 'I')])

        for segment in self.segments_active:
            balances.add_balance(f'{segment}_I', val = 1, units='A', eq_units='V' )
            self.connect(f'balances.{segment}_I',        f'{segment}.I')
            self.connect(f'{segment}.V_cell',            f'balances.lhs:{segment}_I')
            self.connect('balances.U_cell',              f'balances.rhs:{segment}_I')
            self.connect(f'balances.{segment}_I',        f'i_sum.I_{segment}')

        # Sum-of-currents condition: U_cell is unknown, solved so that the
        # sum of all active-segment currents matches the total stack
        # current I (promoted from balances' rhs:U_cell above).
        balances.add_balance('U_cell', val= 1, units='V',
                             eq_units='A',
                             lower=0.3, upper=1.2,
                             use_mult=False, #mult_val=-1
                             )
        self.connect('i_sum.I_total',   'balances.lhs:U_cell')

        super().setup()