import openmdao.api as om
from pycycle.element_base import Element
from .api import SOFCThermoAdd, SegmentSOFC
from pycycle.flow_in import FlowIn
import numpy as np
from pycycle.thermo.thermo import Thermo
from pycycle.thermo.cea.species_data import janaf, Properties
from pycycle.constants import CEA_AIR_COMPOSITION
from . import api as sofc



class SOFC(Element):
    """
    Full SOFC-Model: N_passive + N_active + N_passive SegmentSOFC
    segments (anode/cathode flow, PEN/IC energy balances, electrochemistry)
    enforces an equipotential condition
    (every active segment shares one terminal voltage, solved so that the
    sum of segment currents matches the total stack current).

    Boundary variables promoted to element level
    ----------------------------------------------
    Inputs:
        Fl_I_an:tot:T, Fl_I_an:tot:P           anode inlet temperature/pressure [K]/[Pa]
        Fl_I_an:tot:composition                anode inlet element composition [mol/g]
        Fl_I_an:stat:W                         anode inlet mass flow [kg/s]
        Fl_I_cat:tot:T, Fl_I_cat:tot:P         cathode inlet temperature/pressure [K]/[Pa]
        Fl_I_cat:tot:composition               cathode inlet element composition [mol/g]
        Fl_I_cat:stat:W                        cathode inlet mass flow [kg/s]
        I                                       total stack current [A] (promoted from
                                                 U_balance.rhs:U_cell)
        n_cell                                  number of cells (shared across active segments) [-]
        length_seg_active, length_seg_passive   length of one active / passive segment along the
                                                 flow (conduction distance dx) [m]
        A_seg_active                            active area of one active segment (current
                                                 density i = I_seg / A_seg_active) [m**2]
        height_cell, n_cell,
        thickness_PEN, lambda_PEN,
        thickness_IC, lambda_IC,
        width_IC, height_ch, N_IC_walls          PEN/IC conduction geometry, shared by every
                                                 segment's Conduction group

    Outputs:
        Fl_O_an:tot:*, Fl_O_an:stat:W           anode outlet flow station (h/S/Cp/Cv/T/P, composition, W)
        Fl_O_cat:tot:*, Fl_O_cat:stat:W         cathode outlet flow station (h/S/Cp/Cv/T/P, composition, W)
        U_cell                                  shared terminal voltage across all active
                                                 segments, solved by the equipotential balance [V]
        P_el                                    stack electrical power U_cell * I * n_cell [W]
        H2_utilization_stack,
        O2_utilization_stack                    fuel/oxygen utilization relative to the stack
                                                 inlet flows [-]

    Options
    -------
    N_active, N_passive                     number of active/passive segments (N_passive segments
                                             are placed on both the inlet and outlet side)
    in_composition_anode/cathode            inlet element-basis composition dicts, used as a
                                             fallback when no upstream element sets Fl_I_data
    thermo_data                             CEA thermo spec (default janaf)
    BC_type_left, BC_type_right              'T' or 'Q' boundary type for the stack's two open
                                             faces (segment0's left face, segments[-1]'s right face)
    T_PEN_BC_left, T_IC_BC_left               hardwired PEN/IC temperatures on the open left face,
                                             used when BC_type_left=='T' [K]
    T_PEN_BC_right, T_IC_BC_right             same, for the open right face (BC_type_right=='T') [K]
    Q_cond_PEN_left, Q_cond_IC_left           hardwired heat flux into the open left face,
                                             used when BC_type_left=='Q' [W]
    Q_cond_PEN_right, Q_cond_IC_right         hardwired heat flux into the open right face,
                                             used when BC_type_right=='Q' [W]
    fuel_type                               type of fuel (default 'H2')
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
        self.options.declare('thermo_data',             default=janaf)
        
        self.options.declare('BC_type_left',  default='T', values=['T', 'Q'],
                             desc="Boundary type shared by PEN and IC on the open left face (segment0): "
                                  "'T' hardwires temperatures T_PEN_BC_left/T_IC_BC_left. " \
                                  "'Q' hardwires the heat fluxes Q_cond_PEN_left/Q_cond_IC_left "
                                  "(default 0 W / adiabatic).")
        self.options.declare('BC_type_right', default='Q', values=['T', 'Q'],
                             desc='Same as BC_type_left, for the open right face (segments[-1]), using '
                                  'T_PEN_BC_right/T_IC_BC_right or Q_cond_PEN_right/Q_cond_IC_right.')
        self.options.declare('T_PEN_BC_left', default=973.0,
                             desc='Hardwired PEN temperature on the open (left) face of segment0 [K]. '
                                  'Used when BC_type_left=="T".')
        self.options.declare('T_IC_BC_left', default=973.0,
                             desc='Hardwired IC temperature on the open (left) face of segment0 [K]. '
                                  'Used when BC_type_left=="T".')
        
        self.options.declare('T_PEN_BC_right', default=973.0,
                             desc='Hardwired PEN temperature on the open (right) face of segment[-1] [K]. '
                                  'Used when BC_type_right=="T".')
        self.options.declare('T_IC_BC_right', default=973.0,
                             desc='Hardwired IC temperature on the open (right) face of segment[-1] [K]. '
                                  'Used when BC_type_right=="T".')
        self.options.declare('Q_cond_PEN_left', default=0,
                             desc='Hardwired PEN heat flux on the left face of segment0 [W]')
        self.options.declare('Q_cond_IC_left', default=0,
                             desc='Hardwired IC heat flux on the left face of segment0 [W]')
        self.options.declare('Q_cond_PEN_right', default=0,
                             desc='Hardwired PEN heat flux on the right face of segment[-1] [W]')
        self.options.declare('Q_cond_IC_right', default=0,
                             desc='Hardwired IC heat flux on the right face of segment[-1] [W]')
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
        #
        # Input Definition for Segment0
        #
        # segment0 has no upstream segment to inherit x_in_A/n_in_A from, so it
        # needs its own scratch Properties/Thermo/SpeciesFlowCalc block here.
        # That block must be sized on the same OUTLET-basis element set that
        # segment0 itself uses internally for anode_bulk.x_in etc. (all species
        # calcs inside a segment use the outlet basis -- see sofc_segment.py),
        # not the raw inlet dict: SOFCThermoAdd always adds 'O' to the anode
        # outlet element set (sofc_thermo_add.py, oxygen ions arrive from the
        # cathode even before any current flows), defaulted to exactly 0.0, so
        # this correctly represents "no reaction has happened yet" rather than
        # some nonzero placeholder. Cathode's output_port_data() doesn't add
        # anything, so this is a no-op for cathode -- done anyway for symmetry.
        anode_comp0   = SOFCThermoAdd(reaction_type='anode',   inflow_composition=anode_comp).output_port_data()
        cathode_comp0 = SOFCThermoAdd(reaction_type='cathode', inflow_composition=cathode_comp).output_port_data()

        anode_in_props = Properties(spec, init_elements=anode_comp0)
        cathode_in_props = Properties(spec, init_elements=cathode_comp0)

        anode_thermo_in = Thermo(mode='total_TP', fl_name='Fl_I_an:tot',
                                 method='CEA',
                                 thermo_kwargs={'composition':anode_comp0, 'spec':spec})
        cathode_thermo_in = Thermo(mode='total_TP', fl_name='Fl_I_cat:tot',
                                 method='CEA',
                                 thermo_kwargs={'composition':cathode_comp0, 'spec':spec})
        
        # NOTE: Fl_I_an:tot:T/P, Fl_I_an:stat:W, Fl_I_an:tot:composition (and the
        # cathode equivalents) are themselves promoted INPUTS (from `inflow_an`/
        # `inflow_cat` above) meant to be fed externally (via pyc_connect_flow, or
        # a plain prob.set_val when SOFC is used standalone). OpenMDAO's connect()
        # only goes output->input, so these aliases are promoted directly to the
        # same name as the real port instead of being wired with self.connect()
        # (which is invalid input->input and was never actually validated before,
        # since Problem.setup() alone doesn't check connections -- that only
        # happens in final_setup()/run_model()).
        self.add_subsystem('anode_thermo_in', anode_thermo_in,
                           promotes_inputs=[('T',   'Fl_I_an:tot:T'),
                                            ('P',   'Fl_I_an:tot:P')],
                           promotes_outputs=[('Fl_I_an:tot:h',  'h_A_in'),
                                             ('Fl_I_an:tot:S',  'S_A_in'),
                                             ('Fl_I_an:tot:Cp', 'Cp_A_in'),
                                             ('Fl_I_an:tot:Cv', 'Cv_A_in')])
        self.add_subsystem('cathode_thermo_in', cathode_thermo_in,
                   promotes_inputs=[('T',   'Fl_I_cat:tot:T'),
                                    ('P',   'Fl_I_cat:tot:P')],
                   promotes_outputs=[('Fl_I_cat:tot:h',  'h_C_in'),
                                     ('Fl_I_cat:tot:S',  'S_C_in'),
                                     ('Fl_I_cat:tot:Cp', 'Cp_C_in'),
                                     ('Fl_I_cat:tot:Cv', 'Cv_C_in')])

        self.add_subsystem('anode_species_in',
                            sofc.SpeciesFlowCalc(thermo=anode_in_props,
                                                 electrode='anode',
                                                 composition=anode_comp0),
                            promotes_inputs=[('W', 'Fl_I_an:stat:W')],
                            promotes_outputs=[('x_i',    'x_in_A'),
                                              ('n_flow', 'n_in_A')])

        self.add_subsystem('cathode_species_in',
                           sofc.SpeciesFlowCalc(thermo=cathode_in_props,
                                                electrode='cathode',
                                                composition=cathode_comp0),
                           promotes_inputs=[('W', 'Fl_I_cat:stat:W')],
                           promotes_outputs=[('x_i',    'x_in_C'),
                                             ('n_flow', 'n_in_C')])

        # n and n_moles are not promoted through the Thermo group boundary, so
        # they need a direct connect (same pattern as sofc_segment.py's
        # anode_thermo_out.base_thermo.n -> anode_species_out.n_i) -- without
        # this, *_species_in.n_i/n_moles silently stay at their component
        # defaults (uniform 1.0 / 0.01) instead of real CEA-equilibrium values.
        self.connect('anode_thermo_in.base_thermo.n',        'anode_species_in.n_i')
        self.connect('anode_thermo_in.base_thermo.n_moles',  'anode_species_in.n_moles')
        self.connect('cathode_thermo_in.base_thermo.n',       'cathode_species_in.n_i')
        self.connect('cathode_thermo_in.base_thermo.n_moles', 'cathode_species_in.n_moles')

        self.segment_types = {}
        for i, name in enumerate(self.segments):
            if i <= (N_passive - 1):
                seg_type = 'passive'
            elif (N_passive - 1) < i <= (N_passive + N_active - 1):
                seg_type = 'active'
            elif i > (N_passive + N_active - 1):
                seg_type = 'passive'
            self.segment_types[name] = seg_type

            seg_bc_kwargs = {}
            if i == 0:
                seg_bc_kwargs['left_bc'] = self.options['BC_type_left']
            if i == self.N_seg_tot - 1:
                seg_bc_kwargs['right_bc'] = self.options['BC_type_right']

            seg = SegmentSOFC(type=seg_type,
                            N_segments             = self.N_seg_tot,
                            anode_composition      = anode_comp,
                            cathode_composition    = cathode_comp,
                            spec                   = self.options['thermo_data'],
                            **seg_bc_kwargs)

            seg_promotes = [('p_A_out', 'Fl_I_an:tot:P'), ('p_C_out', 'Fl_I_cat:tot:P')]
            if seg_type == 'active':
                seg_promotes += [('A_seg', 'A_seg_active')]
            # conduction geometry only exists if the segment builds a link
            # (not on a segment0 with a 'Q' left face and no right link)
            has_link = 'T' in (seg.options['left_bc'], seg.options['right_bc'])
            if has_link:
                seg_promotes += [('length_seg', f'length_seg_{seg_type}'), 'height_cell',
                                 'thickness_PEN', 'lambda_PEN',
                                 'thickness_IC', 'lambda_IC',
                                 'width_IC', 'height_ch', 'N_IC_walls']
            if has_link or seg_type == 'active':
                seg_promotes.append('n_cell')
            if i == 0:
                # segment0 has no upstream segment to inherit its inlet state
                # from -- promote its inlet ports directly to the real Fl_I_an/
                # Fl_I_cat flow ports (same reasoning as the pressure above).
                seg_promotes += [('W_in_A', 'Fl_I_an:stat:W'),
                                  ('T_A_in', 'Fl_I_an:tot:T'),
                                  ('composition_in_A', 'Fl_I_an:tot:composition'),
                                  ('W_in_C', 'Fl_I_cat:stat:W'),
                                  ('T_C_in', 'Fl_I_cat:tot:T'),
                                  ('composition_in_C', 'Fl_I_cat:tot:composition')]

            self.add_subsystem(name, seg, promotes_inputs=seg_promotes)

            anode_comp = SOFCThermoAdd(reaction_type='anode', inflow_composition=anode_comp).output_port_data()
            cathode_comp = SOFCThermoAdd(reaction_type='cathode', inflow_composition=cathode_comp).output_port_data()
            
         
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


            if self.segment_types[self.segments[i + 1]] == 'active':
                self.connect(f'{name}.n_out_A',          f'{self.segments[i + 1]}.n_in_A')
                self.connect(f'{name}.n_out_C',          f'{self.segments[i + 1]}.n_in_C')

            # Now connected using promotion across same
            # self.connect(f'{name}.p_A_out',         f'{self.segments[i + 1]}.p_A_out')
            # self.connect(f'{name}.p_C_out',         f'{self.segments[i + 1]}.p_C_out')

            # Conduction across this interface is computed once, by the segment on
            # the right, as its own left link (T_PEN_left = this segment's T_PEN).
            # The flux entering *this* segment through that same link is the
            # link's Q_cond_PEN_left_out -- cross-wired here instead of being
            # recomputed on this side.
            self.connect(f'{name}.T_PEN', f'{self.segments[i + 1]}.T_PEN_left')
            self.connect(f'{name}.T_IC',  f'{self.segments[i + 1]}.T_IC_left')

            self.connect(f'{self.segments[i + 1]}.Q_cond_PEN_left_out', f'{name}.Q_cond_PEN_right')
            self.connect(f'{self.segments[i + 1]}.Q_cond_IC_left_out',  f'{name}.Q_cond_IC_right')


        self.connect('h_A_in',                      'segment0.h_A_in')
        self.connect('x_in_A',                      'segment0.x_in_A')
        # n_in_A only exists on segment0 if it's active (see the n_in_A/n_in_C
        # note in the inter-segment loop above).
        if self.segment_types['segment0'] == 'active':
            self.connect('n_in_A',                  'segment0.n_in_A')
        # Boundary conditions on the two open faces (segment0's left face,
        # segments[-1]'s right face). Every interior face is already wired by the
        # inter-segment loop above (interior right faces are 'neighbour'). 'T':
        # the boundary conduction link gets a hardwired neighbour temperature.
        # 'Q': the segment's flux boundary input Q_*_BC_left/right gets the
        # hardwired flux (0 W = adiabatic).
        last = self.segments[-1]
        if self.options['BC_type_left'] == 'T':
            self.set_input_defaults('segment0.T_PEN_left', val=self.options['T_PEN_BC_left'], units='K')
            self.set_input_defaults('segment0.T_IC_left',  val=self.options['T_IC_BC_left'],  units='K')
        else:
            self.set_input_defaults('segment0.Q_PEN_BC_left', val=self.options['Q_cond_PEN_left'], units='W')
            self.set_input_defaults('segment0.Q_IC_BC_left',  val=self.options['Q_cond_IC_left'],  units='W')

        if self.options['BC_type_right'] == 'T':
            self.set_input_defaults(f'{last}.T_PEN_right', val=self.options['T_PEN_BC_right'], units='K')
            self.set_input_defaults(f'{last}.T_IC_right',  val=self.options['T_IC_BC_right'],  units='K')
        else:
            self.set_input_defaults(f'{last}.Q_PEN_BC_right', val=self.options['Q_cond_PEN_right'], units='W')
            self.set_input_defaults(f'{last}.Q_IC_BC_right',  val=self.options['Q_cond_IC_right'],  units='W')

        # Cathode: same reasoning as the anode block above.
        self.connect('h_C_in',                       'segment0.h_C_in')
        self.connect('x_in_C',                       'segment0.x_in_C')
        if self.segment_types['segment0'] == 'active':
            self.connect('n_in_C',                  'segment0.n_in_C')

        anode_out = Thermo(mode='total_TP', fl_name='Fl_O_an:tot',
                           method='CEA',
                           thermo_kwargs={'composition': self.Fl_O_data['Fl_O_an'],
                                          'spec': spec})
        cathode_out = Thermo(mode='total_TP', fl_name='Fl_O_cat:tot',
                             method='CEA',
                             thermo_kwargs={'composition': self.Fl_O_data['Fl_O_cat'],
                                            'spec': spec})

        self.add_subsystem('anode_out',   anode_out,   promotes_inputs=[('P', 'Fl_I_an:tot:P')],
                                                        promotes_outputs=['Fl_O_an:*'])
        self.add_subsystem('cathode_out', cathode_out, promotes_inputs=[('P', 'Fl_I_cat:tot:P')],
                                                        promotes_outputs=['Fl_O_cat:*'])

        self.connect(f'{self.segments[-1]}.composition_out_A', 'anode_out.composition')
        self.connect(f'{self.segments[-1]}.T_A_out',            'anode_out.T')

        self.connect(f'{self.segments[-1]}.composition_out_C', 'cathode_out.composition')
        self.connect(f'{self.segments[-1]}.T_C_out',            'cathode_out.T')


        self.add_subsystem('anode_out_W',
                           om.ExecComp('W_out = W_in', W_out={'units': 'kg/s'}, W_in={'units': 'kg/s'}),
                           promotes_outputs=[('W_out', 'Fl_O_an:stat:W')])
        self.connect(f'{self.segments[-1]}.W_out_A', 'anode_out_W.W_in')

        self.add_subsystem('cathode_out_W',
                           om.ExecComp('W_out = W_in', W_out={'units': 'kg/s'}, W_in={'units': 'kg/s'}),
                           promotes_outputs=[('W_out', 'Fl_O_cat:stat:W')])
        self.connect(f'{self.segments[-1]}.W_out_C', 'cathode_out_W.W_in')

        # Equipotential condition: every active segment sits on the same
        # electrical interconnect, so they all share one terminal voltage
        # (U_cell, from U_balance below). Each segment's current is unknown,
        # solved so that segment.V_cell matches that shared voltage.
        segment_current_inputs = {f'I_{segment}': {'units': 'A'} for segment in self.segments_active}
        i_sum_expr = 'I_total = ' + ' + '.join(segment_current_inputs)
        self.add_subsystem('i_sum', om.ExecComp(i_sum_expr, I_total={'units': 'A'}, **segment_current_inputs))

        balances = self.add_subsystem('balances', om.BalanceComp(),
                                      promotes_inputs=[(f'rhs:{segment}_I', 'U_cell')
                                                       for segment in self.segments_active])

        for segment in self.segments_active:
            balances.add_balance(f'{segment}_I', val = 1, units='A', eq_units='V' )
            self.connect(f'balances.{segment}_I',        f'{segment}.I')
            self.connect(f'{segment}.V_cell',            f'balances.lhs:{segment}_I')
            self.connect(f'balances.{segment}_I',        f'i_sum.I_{segment}')


        U_balance = self.add_subsystem('U_balance', om.BalanceComp(),
                                       promotes_inputs=[('rhs:U_cell', 'I')],
                                       promotes_outputs=['U_cell'])
        U_balance.add_balance('U_cell', val= 1, units='V',
                              eq_units='A',
                              lower=0.3, upper=1.5,
                              use_mult=False)
        self.connect('i_sum.I_total',   'U_balance.lhs:U_cell')

        self.add_subsystem('power', om.ExecComp('P_el = U_cell * I * n_cell',
                                                P_el={'units': 'W'},
                                                U_cell={'units': 'V'},
                                                I={'units': 'A'}),
                           promotes=['*'])

        self.add_subsystem('H2_util_stack',
                           sofc.SpeciesUtilization(thermo=anode_in_props, species='H2'),
                           promotes_inputs=[('n_in', 'n_in_A'), ('x_in', 'x_in_A'), 'I', 'n_cell'],
                           promotes_outputs=[('H2_utilization',  'H2_utilization_stack'),
                                             ('H2_consumed_mol', 'H2_consumed_stack')])
        self.add_subsystem('O2_util_stack',
                           sofc.SpeciesUtilization(thermo=cathode_in_props, species='O2'),
                           promotes_inputs=[('n_in', 'n_in_C'), ('x_in', 'x_in_C'), 'I', 'n_cell'],
                           promotes_outputs=[('O2_utilization',  'O2_utilization_stack'),
                                             ('O2_consumed_mol', 'O2_consumed_stack')])
        # I feeds U_balance, power and both utilizations; one shared default
        self.set_input_defaults('I', val=1.0, units='A')

        newton = self.nonlinear_solver = om.NewtonSolver(solve_subsystems=True)
        newton.options['maxiter'] = 25
        newton.options['atol'] = 1e-6
        newton.options['rtol'] = 1e-6
        newton.options['iprint'] = 2
        newton.options['reraise_child_analysiserror'] = False
        newton.linesearch = om.ArmijoGoldsteinLS(iprint=0)
        self.linear_solver = om.DirectSolver()

        super().setup()