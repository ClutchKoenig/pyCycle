"""
SOFCThermoAdd: element-basis composition tracker for SOFC anode / cathode channels.

Mirrors the interface of pycycle/thermo/cea/thermo_add.py so it slots into pyCycle's
port-propagation system without touching element_base.py.

Physics (single segment, total current I [A])
---------------------------------------------
  Anode:   H2  +  O²⁻  →  H2O  +  2e⁻
    O²⁻ arrives through the solid electrolyte from the cathode side.
    Element balance: H atoms conserved; O atoms increase by ndot_O = I/(2F) mol/s.
    Mass balance:    W_out = W_in  +  ndot_O × MW_O / G_PER_LBM

  Cathode: O2  +  4e⁻  →  2 O²⁻   (O²⁻ departs to anode)
    Element balance: O atoms decrease by ndot_O = I/(2F) mol/s.
    Mass balance:    W_out = W_in  -  ndot_O2 × MW_O2 / G_PER_LBM

    where ndot_O2 = I/(4F)  and  ndot_O = 2 × ndot_O2 = I/(2F).

Usage in SOFC Element
---------------------
    # pyc_setup_output_ports():
    self.anode_rxn = SOFCThermoAdd(
        reaction_type='anode',
        spec=thermo_data,
        inflow_composition=self.Fl_I_data['Fl_I_an'])
    self.Fl_O_data['Fl_O_an'] = self.anode_rxn.output_port_data()

    self.cathode_rxn = SOFCThermoAdd(
        reaction_type='cathode',
        spec=thermo_data,
        inflow_composition=self.Fl_I_data['Fl_I_cat'])
    self.Fl_O_data['Fl_O_cat'] = self.cathode_rxn.output_port_data()

    # setup():
    self.add_subsystem('anode_rxn', self.anode_rxn,
                       promotes_inputs=['Fl_I_an:stat:W',
                                        ('Fl_I:stat:W', 'Fl_I_an:stat:W'),
                                        ('Fl_I:tot:composition', 'Fl_I_an:tot:composition'),
                                        'I'])
    anode_out = Thermo(mode='total_TP', fl_name='Fl_O_an:tot',
                       method='CEA',
                       thermo_kwargs={'composition': self.Fl_O_data['Fl_O_an'],
                                      'spec': thermo_data})
    self.add_subsystem('anode_out', anode_out, promotes_outputs=['Fl_O_an:*'])
    self.connect('anode_rxn.composition_out', 'anode_out.composition')
    self.connect('anode_rxn.Wout',            'anode_out.W')
    self.connect('T_cell',                    'anode_out.T')   # from PENEnergyBalance
    self.connect('Fl_I_an:tot:P',             'anode_out.P')
"""

import numpy as np
import openmdao.api as om
from pycycle.thermo.cea.species_data import Properties, janaf

FARADAY   = 96485.3321233100184  # C / mol  =  A·s / mol
G_PER_LBM = 453.592              # g / lbm
MW_O      = 16.0                 # g / mol  (atomic oxygen, O²⁻)
MW_O2     = 32.0                 # g / mol
MW_H2     = 2.0                  # g / mol
MW_H2O    = 18.0                 # g / mol


class SOFCThermoAdd(om.ExplicitComponent):
    """
    Tracks element-basis composition change for one SOFC flow channel.

    Outputs 'composition_out' [mol/g_mix] and 'Wout' [lbm/s] so they connect
    directly to Thermo(mode='total_TP') downstream.  Unlike the combustor's
    ThermoAdd, mass_avg_h is NOT computed here; use T_cell from PENEnergyBalance
    together with Thermo(mode='total_TP') to get outlet flow properties.

    call output_port_data() before setup() (i.e. in pyc_setup_output_ports)
    to obtain the outlet element dict, then assign it to Fl_O_data directly:

        self.Fl_O_data['Fl_O_an'] = self.anode_rxn.output_port_data()

    Options
    -------
    reaction_type      : 'anode' | 'cathode'
    spec               : CEA janaf-style thermo data module (default: janaf)
    inflow_composition : dict  e.g. {'H': 1.0} or {'N': 0.76, 'O': 0.24}
    """

    def initialize(self):
        self.options.declare('reaction_type', values=['anode', 'cathode'])
        self.options.declare('spec', default=janaf, recordable=False)
        self.options.declare('inflow_composition',
                             desc='Element mole-ratio dict for the inlet flow')

    # ------------------------------------------------------------------
    # Port-propagation interface (called before setup by the framework)
    # ------------------------------------------------------------------

    def output_port_data(self):
        """
        Return the outlet element dict.

        The anode outlet gains 'O' if not already present, because oxygen ions
        (O²⁻) arrive from the cathode through the solid electrolyte and appear
        as oxygen in the product H2O.
        """
        outlet = dict(self.options['inflow_composition'])
        if self.options['reaction_type'] == 'anode':
            outlet.setdefault('O', 0.0)
        return outlet

    # ------------------------------------------------------------------
    # OpenMDAO component setup
    # ------------------------------------------------------------------

    def setup(self):
        spec          = self.options['spec']
        inflow_comp   = self.options['inflow_composition']
        self._rtype   = self.options['reaction_type']

        in_thermo  = Properties(spec, init_elements=inflow_comp)
        out_thermo = Properties(spec, init_elements=self.output_port_data())

        self._in_elems  = in_thermo.elements    # sorted list, e.g. ['H']
        self._out_elems = out_thermo.elements   # sorted list, e.g. ['H', 'O']

        n_in  = len(self._in_elems)
        n_out = len(self._out_elems)

        # Permutation matrix: maps inlet element indices → outlet element indices
        self._map = np.zeros((n_out, n_in))
        for i, elem in enumerate(self._in_elems):
            if elem in self._out_elems:
                j = self._out_elems.index(elem)
                self._map[j, i] = 1.0

        self._O_out_idx = (self._out_elems.index('O')
                           if 'O' in self._out_elems else None)

        # Inputs
        self.add_input('Fl_I:stat:W', val=1.0, units='lbm/s',
                       desc='Inlet mass flow rate')
        self.add_input('Fl_I:tot:composition', val=in_thermo.b0,
                       desc='Inlet element-basis composition [mol / g_mix]')
        self.add_input('I', val=0.0, units='A',
                       desc='Total current through this cell segment')

        # Outputs
        self.add_output('composition_out', val=out_thermo.b0,
                        desc='Outlet element-basis composition [mol / g_mix]')
        self.add_output('Wout', val=1.0, units='lbm/s',
                        desc='Outlet mass flow rate')
        self.add_output('H2_consumed',  val=0.0, units='lbm/s',
                        desc='H2 consumption rate (anode only)')
        self.add_output('H2O_produced', val=0.0, units='lbm/s',
                        desc='H2O production rate (anode only)')
        self.add_output('O2_consumed',  val=0.0, units='lbm/s',
                        desc='O2 consumption rate (cathode only)')

        self.declare_partials('*', '*', method='cs')

    # ------------------------------------------------------------------
    # Compute
    # ------------------------------------------------------------------

    def compute(self, inputs, outputs):
        W_in  = inputs['Fl_I:stat:W']           # lbm/s
        b0_in = inputs['Fl_I:tot:composition']   # mol / g_mix  (per inlet element)
        I     = inputs['I']                      # A

        # Faraday's law: molar reaction rates [mol/s]
        ndot_H2 = I / (2.0 * FARADAY)   # H2 consumed (anode) = O²⁻ transferred
        ndot_O2 = I / (4.0 * FARADAY)   # O2 consumed (cathode)
        # Note: 2 × ndot_O2 = ndot_H2  (2 O atoms per O2, 1 O atom per H2O)

        # Convert b0 [mol/g_mix] to absolute molar flows [mol/s]:
        #   n [mol/s] = b0 [mol/g] × W [lbm/s] × G_PER_LBM [g/lbm]
        n_in  = b0_in * (W_in * G_PER_LBM)   # mol/s per inlet element

        # Promote to outlet element basis via permutation
        n_out = self._map.dot(n_in)            # mol/s per outlet element

        if self._rtype == 'anode':
            # O²⁻ arrives from cathode → add ndot_H2 O atoms to outlet
            if self._O_out_idx is not None:
                n_out[self._O_out_idx] += ndot_H2

            # Mass: oxygen atoms arriving add mass to anode side
            W_out = W_in + ndot_H2 * MW_O / G_PER_LBM

            outputs['H2_consumed']  = ndot_H2 * MW_H2  / G_PER_LBM
            outputs['H2O_produced'] = ndot_H2 * MW_H2O / G_PER_LBM
            outputs['O2_consumed']  = 0.0

        else:  # cathode
            # O²⁻ departs to anode → remove ndot_H2 O atoms from outlet
            if self._O_out_idx is not None:
                n_out[self._O_out_idx] -= ndot_H2

            # Mass: O2 consumed reduces cathode mass flow
            W_out = W_in - ndot_O2 * MW_O2 / G_PER_LBM

            outputs['H2_consumed']  = 0.0
            outputs['H2O_produced'] = 0.0
            outputs['O2_consumed']  = ndot_O2 * MW_O2 / G_PER_LBM

        # Convert absolute molar flows back to b0 [mol/g_mix]
        outputs['composition_out'] = n_out / (W_out * G_PER_LBM)
        outputs['Wout'] = W_out
