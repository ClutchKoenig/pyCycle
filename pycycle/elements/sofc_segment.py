import numpy as np
import openmdao.api as om

from pycycle.thermo.thermo import Thermo
from pycycle.thermo.cea.species_data import Properties, janaf
import pycycle.sofc_api as sofc


class SegmentSOFC(om.Group):
    """
    One spatial segment of the discretized SOFC model.

    Boundary variables promoted to segment level
    ---------------------------------------------
    Inputs:
        W_in_A, W_in_C                 inlet mass flows [kg/s]
        composition_in_A/C             inlet element composition [mol/g]
        h_A_in, h_C_in                 inlet specific enthalpies [J/kg]
        T_A_in, T_C_in                 inlet temperatures [K]
        p_A_out, p_C_out               channel pressure (no drop assumed) [Pa]
        x_in_A, x_in_C                 inlet mole fraction vectors [-]
        n_in_A, n_in_C                 inlet total molar flows [mol/s]
        I                              segment current [A]  (active only)
        T_PEN_left, T_PEN_right        neighbor PEN temperatures [K]
        T_IC_left,  T_IC_right         neighbor IC  temperatures [K]
        n_cell                         number of cells [-]
        (geometric parameters promoted via Conduction)

    Outputs:
        W_out_A, W_out_C               outlet mass flows [kg/s]
        composition_out_A/C            outlet element composition [mol/g]
        h_A_out, h_C_out               outlet specific enthalpies [J/kg]
        T_A_out, T_C_out               outlet temperatures (implicit states) [K]
        p_A_out, p_C_out               outlet pressures [Pa]
        x_out_A, x_out_C              outlet mole fraction vectors [-]
        n_out_A, n_out_C               outlet total molar flows [mol/s]
        T_PEN, T_IC                    solid temperatures (implicit states) [K]
        V_cell                         cell voltage [V]  (active only)
        H2_utilization, O2_utilization fuel/air utilization [-]  (active only)
    """

    def initialize(self):
        self.options.declare('type', default='passive',
                             values=['passive', 'active'],
                             desc='switches mass-flow calculation and energy balances')
        self.options.declare('N_segments', default=10,
                             desc='Number of segments for current type')
        self.options.declare('anode_composition',
                             desc='Element dict for anode inlet, e.g. {"H": 1.0}')
        self.options.declare('cathode_composition',
                             desc='Element dict for cathode inlet, e.g. {"N": 0.76, "O": 0.24}')
        self.options.declare('spec', default=janaf, recordable=False)

    def setup(self):
        spec         = self.options['spec']
        anode_comp   = self.options['anode_composition']
        cathode_comp = self.options['cathode_composition']
        segment_type = self.options['type']
        N_seg        = self.options['N_segments']

        # Instantiate reaction components before add_subsystem so
        # output_port_data() is available to build downstream Thermo objects.
        anode_rxn   = sofc.SOFCThermoAdd(reaction_type='anode',   spec=spec,
                                          inflow_composition=anode_comp)
        cathode_rxn = sofc.SOFCThermoAdd(reaction_type='cathode', spec=spec,
                                          inflow_composition=cathode_comp)

        # Outlet element composition dicts — used to build Properties and Thermo.
        anode_out_comp   = anode_rxn.output_port_data()    # {'H':..., 'O':...}
        cathode_out_comp = cathode_rxn.output_port_data()  # {'N':..., 'O':...}

        # Properties objects for SpeciesFlowCalc / BulkComposition / SpeciesUtilization.
        # All species-level calculations inside the segment use the outlet basis so
        # that x_in_A / x_out_A always have the same length.
        # Stored as instance attributes so callers can query species indices/sizes
        # after prob.setup() without duplicating SOFCThermoAdd instantiation.
        anode_out_props   = self.anode_out_props   = Properties(spec, init_elements=anode_out_comp)
        cathode_out_props = self.cathode_out_props = Properties(spec, init_elements=cathode_out_comp)

        # Species indices needed for src_indices connections to Electrochemistry.
        H2_idx  = (anode_out_props.products.index('H2')
                   if 'H2'  in anode_out_props.products  else None)
        H2O_idx = (anode_out_props.products.index('H2O')
                   if 'H2O' in anode_out_props.products  else None)
        O2_idx  = (cathode_out_props.products.index('O2')
                   if 'O2'  in cathode_out_props.products else None)

        # Outlet Thermo components (total_TP: T and P are inputs, h/cp/S/Cv outputs).
        anode_thermo_out = Thermo(mode='total_TP', fl_name='Fl_O_an:tot',
                                  method='CEA',
                                  thermo_kwargs={'composition': anode_out_comp, 'spec': spec})
        cathode_thermo_out = Thermo(mode='total_TP', fl_name='Fl_O_cat:tot',
                                    method='CEA',
                                    thermo_kwargs={'composition': cathode_out_comp, 'spec': spec})

        # ==============================================================
        # Add subsystems
        # ==============================================================

        # --- Electrochemical reactions ---
        self.add_subsystem('anode_rxn', anode_rxn,
                           promotes_inputs=[('Fl_I:stat:W',          'W_in_A'),
                                            ('Fl_I:tot:composition',  'composition_in_A'),
                                            'I'],
                           promotes_outputs=[('composition_out', 'composition_out_A'),
                                             ('Wout',            'W_out_A'),
                                             'H2_consumed', 'H2O_produced'])

        self.add_subsystem('cathode_rxn', cathode_rxn,
                           promotes_inputs=[('Fl_I:stat:W',          'W_in_C'),
                                            ('Fl_I:tot:composition',  'composition_in_C'),
                                            'I'],
                           promotes_outputs=[('composition_out', 'composition_out_C'),
                                             ('Wout',            'W_out_C'),
                                             'O2_consumed'])

        # --- Outlet thermodynamic state ---
        # T input comes from ChannelEnergyBalance (implicit T_A_out / T_C_out).
        # P is passed straight through (no pressure drop per segment).
        self.add_subsystem('anode_thermo_out', anode_thermo_out,
                           promotes_inputs=[('T', 'T_A_out'),
                                            ('P', 'p_A_out')],
                           promotes_outputs=[('Fl_O_an:tot:h',  'h_A_out'),
                                             ('Fl_O_an:tot:S',  'S_A_out'),
                                             ('Fl_O_an:tot:Cp', 'Cp_A_out'),
                                             ('Fl_O_an:tot:Cv', 'Cv_A_out')])

        self.add_subsystem('cathode_thermo_out', cathode_thermo_out,
                           promotes_inputs=[('T', 'T_C_out'),
                                            ('P', 'p_C_out')],
                           promotes_outputs=[('Fl_O_cat:tot:h',  'h_C_out'),
                                             ('Fl_O_cat:tot:S',  'S_C_out'),
                                             ('Fl_O_cat:tot:Cp', 'Cp_C_out'),
                                             ('Fl_O_cat:tot:Cv', 'Cv_C_out')])

        # --- Species molar flows at outlet (outlet composition basis) ---
        self.add_subsystem('anode_species_out',
                           sofc.SpeciesFlowCalc(thermo=anode_out_props,
                                                electrode='anode',
                                                composition=anode_out_comp),
                           promotes_inputs=[('W', 'W_out_A')],
                           promotes_outputs=[('x_i',    'x_out_A'),
                                             ('n_flow', 'n_out_A')])

        self.add_subsystem('cathode_species_out',
                           sofc.SpeciesFlowCalc(thermo=cathode_out_props,
                                                electrode='cathode',
                                                composition=cathode_out_comp),
                           promotes_inputs=[('W', 'W_out_C')],
                           promotes_outputs=[('x_i',    'x_out_C'),
                                             ('n_flow', 'n_out_C')])

        # --- Bulk mole fractions (average of segment inlet and outlet) ---
        # x_in_A / x_in_C are external inputs supplied by the parent SOFC element
        # from the previous segment's x_out (same outlet-basis size).
        self.add_subsystem('anode_bulk',
                           sofc.BulkComposition(thermo=anode_out_props),
                           promotes_inputs=[('x_in',  'x_in_A'),
                                            ('x_out', 'x_out_A')],
                           promotes_outputs=[('x_bulk', 'x_bulk_A')])

        self.add_subsystem('cathode_bulk',
                           sofc.BulkComposition(thermo=cathode_out_props),
                           promotes_inputs=[('x_in',  'x_in_C'),
                                            ('x_out', 'x_out_C')],
                           promotes_outputs=[('x_bulk', 'x_bulk_C')])

        # --- Heat transfer ---
        self.add_subsystem('Convection', sofc.HeatConvection(),
                           promotes_inputs=['W_in_A',  'W_out_A',
                                            'W_in_C',  'W_out_C',
                                            'h_A_in',  'h_A_out',
                                            'h_C_in',  'h_C_out',
                                            'T_A_in',  'T_A_out',
                                            'T_C_in',  'T_C_out',
                                            'T_PEN',   'T_IC'],
                           promotes_outputs=['Q_conv_A',     'Q_conv_C',
                                             'Q_conv_PEN_A', 'Q_conv_PEN_C',
                                             'Q_conv_IC_A',  'Q_conv_IC_C'])

        self.add_subsystem('Conduction', sofc.HeatConduction(N_segments=N_seg),
                           promotes=['*'])

        # --- Energy balances ---
        self.add_subsystem('PEN', sofc.PENEnergyBalance(seg_type=segment_type),
                           promotes=['*'])

        self.add_subsystem('IC', sofc.ICEnergyBalance(),
                           promotes=['*'])

        self.add_subsystem('Cathode', sofc.ChannelEnergyBalance(electrode='cathode'),
                           promotes_inputs=[('Q_conv_channel', 'Q_conv_C'),
                                            ('Q_conv_PEN',     'Q_conv_PEN_C'),
                                            ('Q_conv_IC',      'Q_conv_IC_C')],
                           promotes_outputs=[('T_channel', 'T_C_out')])

        self.add_subsystem('Anode', sofc.ChannelEnergyBalance(electrode='anode'),
                           promotes_inputs=[('Q_conv_channel', 'Q_conv_A'),
                                            ('Q_conv_PEN',     'Q_conv_PEN_A'),
                                            ('Q_conv_IC',      'Q_conv_IC_A')],
                           promotes_outputs=[('T_channel', 'T_A_out')])

        # --- Active-only: electrochemistry and utilization ---
        if segment_type == 'active':
            self.add_subsystem('Electrochemistry',
                               sofc.ElectroChemistry(N_segments=N_seg),
                               promotes_inputs=['I', 'A', 'T_PEN',
                                                ('P_cat', 'p_C_out')],
                               promotes_outputs=['V_cell', 'Qdot_chem',
                                                 'U_Nernst', 'U_OCV'])

            self.add_subsystem('H2_Utilization',
                               sofc.SpeciesUtilization(thermo=anode_out_props,
                                                       species='H2'),
                               promotes_inputs=[('n_in', 'n_in_A'),
                                                ('x_in', 'x_in_A'),
                                                'I', 'n_cell'],
                               promotes_outputs=['H2_utilization', 'H2_consumed_mol'])

            self.add_subsystem('O2_Utilization',
                               sofc.SpeciesUtilization(thermo=cathode_out_props,
                                                       species='O2'),
                               promotes_inputs=[('n_in', 'n_in_C'),
                                                ('x_in', 'x_in_C'),
                                                'I', 'n_cell'],
                               promotes_outputs=['O2_utilization', 'O2_consumed_mol'])

        # ==============================================================
        # Internal connections
        # ==============================================================

        # Reaction composition → outlet Thermo
        self.connect('composition_out_A',   'anode_thermo_out.composition')
        self.connect('composition_out_C', 'cathode_thermo_out.composition')

        # Outlet Thermo base_thermo → SpeciesFlowCalc
        # (n and n_moles are not promoted through the Thermo group boundary)
        self.connect('anode_thermo_out.base_thermo.n',        'anode_species_out.n_i')
        self.connect('anode_thermo_out.base_thermo.n_moles',  'anode_species_out.n_moles')
        self.connect('cathode_thermo_out.base_thermo.n',       'cathode_species_out.n_i')
        self.connect('cathode_thermo_out.base_thermo.n_moles', 'cathode_species_out.n_moles')

        # Bulk mole fractions → Electrochemistry (scalar extraction via src_indices)
        if segment_type == 'active':
            if H2_idx is not None:
                self.connect('x_bulk_A', 'Electrochemistry.x_H2',
                             src_indices=[H2_idx])
            if H2O_idx is not None:
                self.connect('x_bulk_A', 'Electrochemistry.x_H2O',
                             src_indices=[H2O_idx])
            if O2_idx is not None:
                self.connect('x_bulk_C', 'Electrochemistry.x_O2',
                             src_indices=[O2_idx])
        
        # Solver Setup:
        newton = self.nonlinear_solver = om.NewtonSolver()
        newton.options['rtol'] = 1e-8
        newton.options['atol'] = 1e-8
        newton.options['iprint'] = 2
        newton.options['maxiter'] = 150
        newton.options['solve_subsystems'] = True
        newton.options['max_sub_solves'] = 50
        newton.options['reraise_child_analysiserror'] = True
             
        ls = newton.linesearch = om.ArmijoGoldsteinLS()
        ls.options['maxiter'] = 20
        ls.options['rho'] = 0.75
        ls.options['print_bound_enforce']=True

        self.linear_solver = om.DirectSolver()

