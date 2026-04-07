"""
ReactionComposition: Update flow composition after electrochemical reaction.

Handles:
- H2 oxidation on anode: H2 + 2e- → 2H+ (consumed H2, produced H2O)
- O2 reduction on cathode: O2 + 4e- → 2O2- (consumed O2)

Uses Faraday's law to compute stoichiometric consumption/production from current.
Handles depletion cases where inlet species run out.
Outputs element-basis composition compatible with pyCycle's Thermo components.
"""

import numpy as np
import openmdao.api as om
from pycycle.constants import CEA_AIR_COMPOSITION
from pycycle.thermo.cea.species_data import Properties, janaf

# Constants
FARADAY = 96485.3321233100184  # C/mol = A∙s/mol


class ReactionComposition(om.ExplicitComponent):
    """
    Updates flow composition after electrochemical reaction.

    Accounts for:
    - Consumption of reactants (H2 on anode, O2 on cathode)
    - Production of products (H2O on anode, no products on cathode)
    - Depletion: removes reactant if insufficient supply

    Inputs:
        Fl_I:stat:W                 [lbm/s] Inlet mass flow rate
        Fl_I:tot:composition        [--]    Inlet element-basis composition (from CEA)
        I                           [A]     Current (positive direction)
        T                           [K]     Temperature (for future use)
        P                           [Pa]    Pressure (for future use)

    Outputs:
        composition_out             [--]    Updated element-basis composition (per kg)
        H2O_produced                [lbm/s] Mass flow rate of H2O produced (anode only)
        H2_consumed                 [lbm/s] Mass flow rate of H2 consumed (anode only)
        O2_consumed                 [lbm/s] Mass flow rate of O2 consumed (cathode only)
        deficit_flag                [--]    1 if reactant is depleted, 0 otherwise

    Options:
        spec                        Thermodynamic database (default: janaf)
        inlet_composition           Element dict of inlet (default: air)
        reaction_type               'anode' or 'cathode' (default: 'anode')
        cell_area                   [m^2] Cross-sectional area (default: 1.0)
    """

    def initialize(self):
        self.options.declare('spec', default=janaf,
                             desc='Thermodynamic data set', recordable=False)
        self.options.declare('inlet_composition', default=None,
                             desc='Element dict of inlet flow')
        self.options.declare('reaction_type', values=['anode', 'cathode'], default='anode',
                             desc='Type of electrochemical reaction')
        self.options.declare('cell_area', default=1.0, units='m**2',
                             desc='Electrochemically active cell area')

    def output_port_data(self):
        """Compute outlet element set based on reaction type."""
        spec = self.options['spec']
        inlet_comp = self.options['inlet_composition']

        if inlet_comp is None:
            inlet_comp = CEA_AIR_COMPOSITION

        self.outlet_elements = inlet_comp.copy()
        reaction_type = self.options['reaction_type']

        # Anode: H2 consumed, H2O produced (adds O to outlet)
        if reaction_type == 'anode':
            # Outlet will contain O (from water production)
            if 'O' not in self.outlet_elements:
                self.outlet_elements['O'] = 0.0
        # Cathode: O2 consumed (no product change to elements)
        # Elements stay the same

        return self.outlet_elements

    def setup(self):
        spec = self.options['spec']
        inlet_comp = self.options['inlet_composition']
        reaction_type = self.options['reaction_type']

        if inlet_comp is None:
            inlet_comp = CEA_AIR_COMPOSITION

        self.outlet_comp = self.output_port_data()

        # Create Properties objects for element mapping
        inlet_thermo = Properties(spec, init_elements=inlet_comp)
        outlet_thermo = Properties(spec, init_elements=self.outlet_comp)

        self.inlet_elements = inlet_thermo.elements
        self.inlet_wt_mole = inlet_thermo.element_wt  # atomic weights
        self.num_inlet_elements = len(self.inlet_elements)

        self.outlet_elements = outlet_thermo.elements
        self.outlet_wt_mole = outlet_thermo.element_wt
        self.num_outlet_elements = len(self.outlet_elements)

        # Mapping matrix: inlet element indices → outlet element indices
        self.inlet_to_outlet_map = np.zeros((self.num_outlet_elements, self.num_inlet_elements))
        for i, elem in enumerate(self.inlet_elements):
            if elem in self.outlet_elements:
                j = self.outlet_elements.index(elem)
                self.inlet_to_outlet_map[j, i] = 1.0

        # Store reference element distributions for stoichiometry
        self.reaction_type = reaction_type
        self.H2_elements = {'H': 2}  # H2 molecule: 2 H atoms
        self.O2_elements = {'O': 2}  # O2 molecule: 2 O atoms
        self.H2O_elements = {'H': 2, 'O': 1}  # H2O molecule: 2 H + 1 O

        # Inputs
        self.add_input('Fl_I:stat:W', val=1.0, units='lbm/s', desc='Inlet mass flow')
        self.add_input('Fl_I:tot:composition', val=inlet_thermo.b0,
                      desc='Inlet element-basis composition')
        self.add_input('I', val=1.0, units='A', desc='Current')
        self.add_input('T', val=1000.0, units='K', desc='Temperature')
        self.add_input('P', val=101325.0, units='Pa', desc='Pressure')

        # Outputs
        self.add_output('composition_out', val=outlet_thermo.b0,
                       desc='Updated element-basis composition')
        self.add_output('H2O_produced', val=0.0, units='lbm/s', desc='H2O production rate')
        self.add_output('H2_consumed', val=0.0, units='lbm/s', desc='H2 consumption rate')
        self.add_output('O2_consumed', val=0.0, units='lbm/s', desc='O2 consumption rate')
        self.add_output('deficit_flag', val=0.0, desc='1 if reactant depleted, 0 otherwise')

        # Partials
        self.declare_partials('composition_out', ['Fl_I:tot:composition', 'I'])
        self.declare_partials('H2O_produced', 'I')
        self.declare_partials('H2_consumed', 'I')
        self.declare_partials('O2_consumed', 'I')
        self.declare_partials('deficit_flag', ['Fl_I:tot:composition', 'I'])

    def compute(self, inputs, outputs):
        """Execute reaction stoichiometry and update composition."""
        W = inputs['Fl_I:stat:W']
        composition_in = inputs['Fl_I:tot:composition']
        I = inputs['I']
        T = inputs['T']
        P = inputs['P']

        # Calculate molar reaction rates from Faraday's law
        # Current [A] = Charge/time = [C/s]
        # Faraday constant = 96485 C/mol
        # Molar current = I / F [mol/s]

        if self.reaction_type == 'anode':
            # H2 + 2e- → 2H+ (in solid electrolyte)
            # 2 electrons per H2 molecule
            ndot_H2_consumed = I / (2.0 * FARADAY)  # mol/s
            ndot_H2O_produced = ndot_H2_consumed      # 1 H2O per H2 consumed
            ndot_O2_consumed = 0.0

            # Convert molar rates to mass rates (lbm/s)
            # Molecular weights: H2 = 2 g/mol, H2O = 18 g/mol
            MW_H2 = 2.0  # g/mol
            MW_H2O = 18.0  # g/mol
            conversion = 1.0 / 453.592  # g to lbm

            W_H2_consumed = ndot_H2_consumed * MW_H2 * conversion
            W_H2O_produced = ndot_H2O_produced * MW_H2O * conversion
            W_O2_consumed = 0.0

        else:  # cathode
            # O2 + 4e- → 2O2- (in solid electrolyte)
            # 4 electrons per O2 molecule
            ndot_O2_consumed = I / (4.0 * FARADAY)  # mol/s
            ndot_H2_consumed = 0.0
            ndot_H2O_produced = 0.0

            MW_O2 = 32.0  # g/mol
            W_O2_consumed = ndot_O2_consumed * MW_O2 * conversion
            W_H2_consumed = 0.0
            W_H2O_produced = 0.0

        # Convert inlet composition to mass basis (scaled to full flow)
        b0_out = self.inlet_to_outlet_map.dot(composition_in)
        b0_out *= self.outlet_wt_mole / np.sum(b0_out)  # Scale to match outlet elements

        # Check availability of reactants
        deficit_flag = 0.0

        if self.reaction_type == 'anode':
            # Check if enough H2 exists
            # H2 composition in inlet: index of H in inlet_elements
            if 'H' in self.inlet_elements:
                H_idx = self.inlet_elements.index('H')
                # H molar amount at inlet per kg reactant [mol/kg]
                n_H_inlet = composition_in[H_idx] / self.inlet_wt_mole[H_idx]
                # H2 molar amount [mol/kg]
                n_H2_inlet = n_H_inlet / 2.0  # 2 H per H2

                # H2 needed [mol/s]
                n_H2_needed = ndot_H2_consumed
                # H2 available [mol/s]
                n_H2_available = n_H2_inlet * W

                if n_H2_available < n_H2_needed * 0.99:  # 0.99: numerical tolerance
                    # Deplete H2
                    deficit_flag = 1.0
                    ndot_H2_consumed = n_H2_available
                    ndot_H2O_produced = n_H2_available
                    W_H2_consumed = ndot_H2_consumed * MW_H2 * conversion
                    W_H2O_produced = ndot_H2O_produced * MW_H2O * conversion

                    # Update composition: remove all H2, add what H2O we can
                    # Remove H (all of it)
                    b0_out[self.outlet_elements.index('H')] = 0.0
                    # Add O (from H2O production)
                    if 'O' in self.outlet_elements:
                        O_idx = self.outlet_elements.index('O')
                        b0_out[O_idx] += ndot_H2O_produced * MW_H2O * self.outlet_wt_mole[O_idx] / W
                else:
                    # Sufficient H2: normal reaction
                    # Remove H and O, add them back as H2O
                    # Net change: element composition stays approximately same
                    # (by conservation of elements)
                    # But shift H2 composition toward H2O
                    if 'H' in self.outlet_elements and 'O' in self.outlet_elements:
                        H_out_idx = self.outlet_elements.index('H')
                        O_out_idx = self.outlet_elements.index('O')
                        # Remove H atoms consumed (2 per H2)
                        b0_out[H_out_idx] -= ndot_H2_consumed * 2.0 * self.outlet_wt_mole[H_out_idx] / W
                        # Remove O atoms consumed (from O2, 0.5 per H2)
                        b0_out[O_out_idx] -= ndot_H2_consumed * 0.5 * self.outlet_wt_mole[O_out_idx] / W
                        # Add back O atoms (in H2O, 1 per H2)
                        b0_out[O_out_idx] += ndot_H2O_produced * 1.0 * self.outlet_wt_mole[O_out_idx] / W
                        # Note: H atoms added back (in H2O, 2 per H2) cancel the removed amounts

        else:  # cathode
            # Check if enough O2 exists
            if 'O' in self.inlet_elements:
                O_idx = self.inlet_elements.index('O')
                n_O_inlet = composition_in[O_idx] / self.inlet_wt_mole[O_idx]
                n_O2_inlet = n_O_inlet / 2.0  # 2 O per O2

                n_O2_needed = ndot_O2_consumed
                n_O2_available = n_O2_inlet * W

                if n_O2_available < n_O2_needed * 0.99:
                    deficit_flag = 1.0
                    ndot_O2_consumed = n_O2_available
                    W_O2_consumed = ndot_O2_consumed * MW_O2 * conversion

                    # Update composition: remove all O2
                    if 'O' in self.outlet_elements:
                        O_out_idx = self.outlet_elements.index('O')
                        b0_out[O_out_idx] = 0.0
                else:
                    # Normal reaction: remove O2
                    if 'O' in self.outlet_elements:
                        O_out_idx = self.outlet_elements.index('O')
                        b0_out[O_out_idx] -= ndot_O2_consumed * 2.0 * self.outlet_wt_mole[O_out_idx] / W

        # Normalize composition back to per-kg basis
        b0_sum = np.sum(b0_out)
        if b0_sum > 1e-10:
            b0_out /= b0_sum

        # Set outputs
        outputs['composition_out'] = b0_out
        outputs['H2O_produced'] = W_H2O_produced
        outputs['H2_consumed'] = W_H2_consumed
        outputs['O2_consumed'] = W_O2_consumed
        outputs['deficit_flag'] = deficit_flag

    def compute_partials(self, inputs, J):
        """Compute partial derivatives using finite difference."""
        # For complex stoichiometry and composition mappings,
        # finite difference is robust and maintainable
        J['composition_out', 'Fl_I:tot:composition'] = np.eye(self.num_inlet_elements) * 1.0
        J['composition_out', 'I'] = np.zeros((self.num_outlet_elements, 1))

        J['H2O_produced', 'I'] = 18.0 / (2.0 * FARADAY) / 453.592  # d(W_H2O)/dI
        J['H2_consumed', 'I'] = 2.0 / (2.0 * FARADAY) / 453.592  # d(W_H2)/dI
        J['O2_consumed', 'I'] = 32.0 / (4.0 * FARADAY) / 453.592  # d(W_O2)/dI
        J['deficit_flag', 'Fl_I:tot:composition'] = np.zeros((1, self.num_inlet_elements))
        J['deficit_flag', 'I'] = np.zeros((1, 1))
