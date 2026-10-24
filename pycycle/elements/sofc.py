import openmdao.api as om
from pycycle.element_base import Element
from pycycle.sofc_api import SOFCThermoAdd, SegmentSOFC, HeatConductionStructure
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
        self.options.declare('thermo_data',             default=janaf)
        
        self.options.declare('BC_type_left',  default='T', values=['T', 'Q'],
                             desc="Boundary type shared by PEN and IC on the open left face (segment0): "
                                  "'T' hardwires temperatures T_PEN_BC_left/T_IC_BC_left. " \
                                  "'Q' leaves the conductive flux at its component default (0 W / adiabatic)"
                                  "set a nonzero flux via prob.set_val('<path>.segment0.Q_cond_PEN_left', ...) "
                                  "(and/or Q_cond_IC_left) after setup() if needed.")
        self.options.declare('BC_type_right', default='Q', values=['T', 'Q'],
                             desc='Same as BC_type_left, for the open right face (segments[-1]), using '
                                  'T_PEN_BC_right/T_IC_BC_right or Q_cond_PEN_right/Q_cond_IC_right.')
        self.options.declare('T_PEN_BC_left', default=973.0,
                             desc='Hardwired PEN temperature on the open (left) face of segment0 [K]. '
                                  'Used when BC_type_left=="T".')
        self.options.declare('T_IC_BC_left', default=973.0,
                             desc='Hardwired IC temperature on the open (left) face of segment0 [K]. '
                                  'Used when BC_type_left=="T".')
        
        self.options.declare('Q_cond_PEN_right', default=0,
                             desc='Hardwired PEN heat flux on the right face of segment[-1] [W]')
        self.options.declare('Q_cond_PEN_right', default=0,
                                     desc='Hardwired PEN heat flux on the right face of segment[-1] [W]')
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

            ic_cond = HeatConductionStructure(struture='IC', N_segments=self.N_seg_tot)
            pen_cond = HeatConductionStructure(struture='PEN', N_segments=self.N_seg_tot)

            cond_promotes   = ['length_cell', 'height_cell','thickness_struct', 'lambda_struct', 'n_cell']
            #ic_promotes     = [('T_cell_left',  'T_IC_left'),
                             #('T_cell_right', 'T_IC')]
            #pen_promotes    = [('T_cell_left',  'T_IC_left'), ('T_cell_right', 'T_PEN')]
            #ic_promotes     = ic_promotes.append(cond_promotes)
            #pen_promotes    = pen_promotes.append(cond_promotes)

            if i==0:
                 cond_promotes += [('T_cell_left',  'T_boundary_left')]
            # elif i == self.N_seg_tot - 1:
            #     cond_promotes += [('T_cell_right')]
            self.add_subsystem(f'{name}_IC_Conduction', ic_cond, promotes_inputs=cond_promotes)
            self.add_subsystem(f'{name}_PEN_Conduction', pen_cond, promotes_inputs=cond_promotes) 

            seg_promotes = [('p_A_out', 'Fl_I_an:tot:P'), ('p_C_out', 'Fl_I_cat:tot:P')]
             
            if seg_type == 'active':
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

            self.connect(f'{name}.T_PEN',           f'{self.segments[i + 1]}_PEN_Conduction.T_PEN_right')
            self.connect(f'{name}.T_IC',            f'{self.segments[i + 1]}_IC_Conduction.T_IC_right')

            self.connect(f'{name}.T_PEN',           f'{self.segments[i + 1]}.T_PEN_right')
            self.connect(f'{name}.T_IC',            f'{self.segments[i + 1]}.T_IC_right')

            
            # Conduction across this interface is computed once, by the segment on
            # the right, as its own left link (T_PEN_left=this segment's T_PEN,
            # T_PEN=its own state). The flux entering *this* segment from that same
            # link is that link's Q_cond_PEN_left_out -- cross-wired here instead of
            # being independently recomputed on this side.
            if i == self.N_seg_tot - 1:
                self.set_input_defaults(f'{self.segments[i]}.Q_cond_PEN_right', val= 0, units='W')
                self.set_input_defaults(f'{self.segments[i]}.Q_cond_IC_right', val= 0, units='W')
                
            self.connect(f'{self.segments[i + 1]}_PEN_Conduction.Q_left', f'{name}.Q_cond_PEN_right')
            self.connect(f'{self.segments[i + 1]}_IC_Conduction.Q_left',  f'{name}.Q_cond_IC_right')

            self.connect(f'{name}_PEN_Conduction.Q_right',      f'{name}.Q_cond_PEN_right')
            self.connect(f'{name}_IC_Conduction.Q_right',      f'{name}.Q_cond_IC_right')

            self.connect(f'{name}_PEN_Conduction.Q_left')

        # Initial Flowstation connections
        # Anode: W_in_A/T_A_in/composition_in_A/p_A_out are promoted directly to
        # Fl_I_an:stat:W/Fl_I_an:tot:T/Fl_I_an:tot:composition/Fl_I_an:tot:P above
        # (segment0 has no upstream segment to connect() those from -- and those
        # Fl_I_an:* names are themselves inputs, so connect() couldn't be used
        # even if it had one). h_A_in/x_in_A/n_in_A ARE real outputs (from
        # anode_thermo_in/anode_species_in above), so those connects are valid.
        self.connect('h_A_in',                      'segment0.h_A_in')
        self.connect('x_in_A',                      'segment0.x_in_A')
        # n_in_A only exists on segment0 if it's active (see the n_in_A/n_in_C
        # note in the inter-segment loop above).
        if self.segment_types['segment0'] == 'active':
            self.connect('n_in_A',                  'segment0.n_in_A')
        # Boundary conditions on the two open faces (segment0's left face,
        # segments[-1]'s right face). segment0.Q_cond_PEN_right / Q_cond_IC_right
        # (and, symmetrically, segments[-1].Q_cond_PEN_left / Q_cond_IC_left) are
        # already wired by the inter-segment loop above -- only the truly open
        # faces need a boundary value here, and its type ('T' vs 'Q') is what
        # SegmentSOFC was built with above (left_bc for segment0, right_bc for
        # segments[-1]). For 'Q', nothing is set here: Q_cond_PEN/IC_left(right)
        # already default to 0 W (adiabatic) on PENEnergyBalance/ICEnergyBalance;
        # set a nonzero flux with prob.set_val(...) after setup() if needed.
        if self.options['BC_type_left'] == 'T':
            self.set_input_defaults('segment0.T_PEN_left', val=self.options['T_PEN_BC_left'], units='K')
            self.set_input_defaults('segment0.T_IC_left',  val=self.options['T_IC_BC_left'],  units='K')

        if self.options['BC_type_right'] == 'Q':
            self.set_input_defaults(f'{self.segments[-1]}.T_PEN_right', val=self.options['T_PEN_BC_right'], units='K')
            self.set_input_defaults(f'{self.segments[-1]}.T_IC_right',  val=self.options['T_IC_BC_right'],  units='K')

        # Cathode: same reasoning as the anode block above.
        self.connect('h_C_in',                       'segment0.h_C_in')
        self.connect('x_in_C',                       'segment0.x_in_C')
        if self.segment_types['segment0'] == 'active':
            self.connect('n_in_C',                  'segment0.n_in_C')

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

        self.add_subsystem('anode_out',   anode_out,   promotes_inputs=[('P', 'Fl_I_an:tot:P')],
                                                        promotes_outputs=['Fl_O_an:*'])
        self.add_subsystem('cathode_out', cathode_out, promotes_inputs=[('P', 'Fl_I_cat:tot:P')],
                                                        promotes_outputs=['Fl_O_cat:*'])

        self.connect(f'{self.segments[-1]}.composition_out_A', 'anode_out.composition')
        self.connect(f'{self.segments[-1]}.T_A_out',            'anode_out.T')

        self.connect(f'{self.segments[-1]}.composition_out_C', 'cathode_out.composition')
        self.connect(f'{self.segments[-1]}.T_C_out',            'cathode_out.T')

        # Fl_O_*:stat:W has no natural source: total_TP mode has no W input/output,
        # so there's nothing to promote it from -- a tiny pass-through gives it a
        # real output to connect() the last segment's mass flow into.
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
        # (balances.U_cell). Each segment's current is unknown, solved so
        # that segment.V_cell matches that shared voltage.
        segment_current_inputs = {f'I_{segment}': {'units': 'A'} for segment in self.segments_active}
        i_sum_expr = 'I_total = ' + ' + '.join(segment_current_inputs)
        self.add_subsystem('i_sum', om.ExecComp(i_sum_expr, I_total={'units': 'A'}, **segment_current_inputs))

        # rhs:{segment}_I is promoted to the same name as U_cell's own output
        # (below) rather than connected to it -- connect() requires two
        # different subsystems, and both live on this same `balances`
        # BalanceComp instance. Promoting an input and the (single) output that
        # feeds it to one shared name, at this parent level, is the standard
        # way to broadcast one component's own solved value back into another
        # input on itself.
        balance_input_promotes = [('rhs:U_cell', 'I')]
        balance_input_promotes += [(f'rhs:{segment}_I', 'U_cell') for segment in self.segments_active]
        balances = self.add_subsystem('balances', om.BalanceComp(),
                                      promotes_inputs=balance_input_promotes,
                                      promotes_outputs=['U_cell'])

        for segment in self.segments_active:
            balances.add_balance(f'{segment}_I', val = 1, units='A', eq_units='V' )
            self.connect(f'balances.{segment}_I',        f'{segment}.I')
            self.connect(f'{segment}.V_cell',            f'balances.lhs:{segment}_I')
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

        # Segment-to-segment coupling (T_PEN/T_IC and cross-wired Q_cond chains)
        # plus the equipotential/sum-of-currents `balances` above form real
        # algebraic loops across the whole stack -- each SegmentSOFC converges
        # itself internally, but nothing above drove that outer loop to
        # convergence.
        #
        # NonlinearBlockGS was tried first (matching the originally documented
        # intent of an outer GS sweep across segments) but doesn't work here:
        # it has no mechanism to invert an ImplicitComponent's residual (the
        # `balances` BalanceComp above), so it just leaves I_{segment}/U_cell
        # sitting at their initial guesses forever while the segments'
        # thermal states genuinely converge -- and since NLBGS's reported
        # residual is dominated by those converging quantities, it reports
        # "Converged" while silently never having solved the electrical
        # balance at all (confirmed empirically: i_sum.I_total stayed at its
        # initial-guess sum instead of matching the requested `I`). A flattened
        # Newton solve, which linearizes and correctly inverts implicit states
        # wherever they live in the tree, is needed instead.
        newton = self.nonlinear_solver = om.NewtonSolver(solve_subsystems=True)
        newton.options['maxiter'] = 25
        newton.options['atol'] = 1e-6
        newton.options['rtol'] = 1e-6
        newton.options['iprint'] = 2
        newton.options['reraise_child_analysiserror'] = False
        newton.linesearch = om.ArmijoGoldsteinLS(iprint=0)
        self.linear_solver = om.DirectSolver()

        super().setup()