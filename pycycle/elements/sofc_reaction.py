import numpy as np
import openmdao.api as om
from pycycle.element_base import Element
from pycycle.thermo.thermo import Thermo, ThermoAdd
from pycycle.thermo.cea.species_data import Properties, janaf

from pycycle.flow_in import FlowIn

RM = 8.3145                    # J/(mol·K)
FARADAY = 96485.3321233100184  # C/mol = As/mol
P_REF   = 101325.0             # Pa  (1 atm reference pressure for Nernst)


class NernstThermo(om.ExplicitComponent):
    """
    Computes Nernst voltage, thermoneutral voltage, and PEN reaction heat
    using NASA polynomial data from pyCycle's CEA thermo — no polynomial fit.

    Reaction:  H2 + 0.5 O2 → H2O,  n_e = 2

    E_OCV   = -ΔG°(T) / (2F)                              [V]
    E_Nernst= E_OCV + (R·T/2F)·ln(x_H2·√(x_O2·P/P_ref) / x_H2O)  [V]
    V_tn    = -ΔH°(T) / (2F)                              [V]
    Qdot_chem = V_tn · I                                   [W]  (total reaction enthalpy rate)

    The net heat deposited in the PEN is: Qdot_chem - V_cell·I = (V_tn - V_cell)·I

    Options
    -------
    spec : janaf-style thermo data module (default: janaf)
    """

    def initialize(self):
        self.options.declare('spec', default=janaf, recordable=False)

    def setup(self):
        # Build a Properties object that covers H and O elements so that
        # H2, O2, and H2O are all present in the products list.
        props = Properties(self.options['spec'], init_elements={'H': 2, 'O': 1})
        self._props    = props
        self._idx_H2   = props.products.index('H2')
        self._idx_O2   = props.products.index('O2')
        self._idx_H2O  = props.products.index('H2O')

        self.add_input('T',      val=1000., units='K',   desc='Reaction temperature (T_cell)')
        self.add_input('I',      val=0.0,   units='A',   desc='Total current')
        self.add_input('x_H2',   val=0.5,   units=None,  desc='H2 mole fraction at anode TPB')
        self.add_input('x_H2O',  val=0.5,   units=None,  desc='H2O mole fraction at anode TPB')
        self.add_input('x_O2',   val=0.21,  units=None,  desc='O2 mole fraction at cathode TPB')
        self.add_input('P',      val=101325., units='Pa', desc='Operating pressure')

        self.add_output('E_OCV',     val=1.0,  units='V', desc='Standard Nernst voltage -ΔG°/(2F)')
        self.add_output('E_Nernst',  val=1.0,  units='V', desc='Nernst voltage with concentration correction')
        self.add_output('V_tn',      val=1.25, units='V', desc='Thermoneutral voltage -ΔH°/(2F)')
        self.add_output('Qdot_chem', val=0.0,  units='W', desc='Total reaction enthalpy rate V_tn·I')

        self.declare_partials('*', '*', method='cs')

    def compute(self, inputs, outputs):
        T      = inputs['T']
        I      = inputs['I']
        x_H2   = inputs['x_H2']
        x_H2O  = inputs['x_H2O']
        x_O2   = inputs['x_O2']
        P      = inputs['P']

        H0 = self._props.H0(T)   # dimensionless H_j / (R·T), shape (num_prod,)
        S0 = self._props.S0(T)   # dimensionless S_j / R,     shape (num_prod,)

        # ΔH°(T) and ΔS°(T) for H2 + 0.5 O2 → H2O  [J/mol]
        dH0 = H0[self._idx_H2O] - H0[self._idx_H2] - 0.5 * H0[self._idx_O2]
        dS0 = S0[self._idx_H2O] - S0[self._idx_H2] - 0.5 * S0[self._idx_O2]

        dH = dH0 * RM * T          # J/mol  (H0 = H/(R·T)  →  H = H0·R·T)
        dG = dH - T * dS0 * RM     # J/mol  (G = H - T·S,  S = S0·R)

        outputs['V_tn']      = -dH / (2.0 * FARADAY)
        outputs['E_OCV']     = -dG / (2.0 * FARADAY)
        outputs['E_Nernst']  = outputs['E_OCV'] + (RM * T / (2.0 * FARADAY)) * np.log(
                                    x_H2 * np.sqrt(x_O2 * P / P_REF) / x_H2O)
        outputs['Qdot_chem'] = outputs['V_tn'] * I

class VoltageCalc(om.ExplicitComponent):
    """
    Calculates the voltage of a SOFC MEA segment using 
    the Area specific ohmic resistance overpotential 'eta_ASR and Nernst Potential 'E'
    V_cell = E_nernst - eta_ASR 
    """
    def setup(self):
        self.add_input('U_Nernst', val=1.23, units='V', desc='Voltage considering concentration of gases')
        self.add_input('eta_asr', val=0.3, units='V', desc='Overpotential Voltage - Losses due to combined "area specific resistance (ASR)"')
        #self.add_input('i', val=1, units='A/m**2', desc=)
        
        self.add_output('V_cell', val=0.93, units='V', desc='Voltage')
        self.declare_partials('V_cell', ['U_Nernst', 'eta_asr'])

    def compute(self, inputs, outputs):
        U_Nernst = inputs['U_Nernst']
        eta = inputs['eta_asr']
        outputs['V_cell'] = U_Nernst - eta

    def compute_partials(self, inputs, J):
        U_Nernst = inputs['U_Nernst']
        eta = inputs['eta_asr']

        J['V_cell', 'U_Nernst'] = 1
        J['V_cell', 'eta_asr'] = -1


class AreaSpecificResistanceOverpotential(om.ExplicitComponent):
    """
    ASR = (exp((1/T - 1/T_ref_asr) * E_a_asr * F / RM) * R_ref_asr + R_c_asr) * ASR_nom_750
    eta = ASR * i
    """
    def setup(self):
        self.add_input('T_PEN',     val=1001.9, units='K', desc='Temperature at and in the electrochemical active area')
        self.add_input('T_ref_asr', val=800,    units='K', desc='Reference Temperature at which Correlation of ASR is calibrated with')
        self.add_input('E_a_asr',   val=0.2,       units='V', desc='Reference activation energy')
        self.add_input('R_ref_asr', val=,       units='ohm', desc='Reference asr resistance')
        self.add_input('R_c_asr',   val=,       units='ohm', desc='Reference asr resistance')
        self.add_input('ASR_nom_750',val=,      units='ohm*m**2', desc=)
        self.add_input('i',         val=,       units='A/m/m', desc='Current density')
        # Einheiten inkonsistent? Was hat ASR_nom für ne Einheit?
        # eta_asr = A/ m**2 * ohm * ASR (ohm*m**2?)
        self.add_output('ASR', val=, units='ohm*m**2', desc=)
        self.add_output('eta_asr', val=, units='V', desc='')

        self.declare_partials('ASR', ['T_PEN', 'T_ref_asr', 'E_a_asr', 'R_ref_asr', 'R_c_asr', 'ASR_nom_750'])
        self.declare_partials('eta_asr', ['T_PEN', 'i', 'T_ref_asr', 'E_a_asr', 'R_ref_asr', 'R_c_asr', 'ASR_nom_750'])

    def compute(self, inputs, outputs):
        T = inputs['T_PEN']
        T_ref = inputs['T_ref_asr']
        E_a = inputs['E_a_asr']
        R_ref = inputs['R_ref_asr']
        R_c = inputs['R_c_asr']
        ASR_nom = inputs['ASR_nom_750']
        i = inputs['i']

        ASR = (np.exp((1/T - 1/T_ref) * E_a * FARADAY / RM) * R_ref + R_c) * ASR_nom

        outputs['ASR'] = ASR
        outputs['eta_asr'] = ASR * i

    def compute_partials(self, inputs, J):
        T = inputs['T_PEN']
        T_ref = inputs['T_ref_asr']
        E_a = inputs['E_a_asr']
        R_ref = inputs['R_ref_asr']
        R_c = inputs['R_c_asr']
        ASR_nom = inputs['ASR_nom_750']
        i = inputs['i']
        exp_term = np.exp((1/T - 1/ T_ref) * E_a * FARADAY / RM)
                          
        J['ASR', 'T_PEN'] = (- (1/T**2) * E_a * FARADAY / RM * exp_term * R_ref) * ASR_nom
        J['ASR', 'T_ref_asr'] = ( (1/T_ref**2) * E_a * FARADAY / RM  * exp_term * R_ref)* ASR_nom
        J['ASR', 'E_a_asr'] = (1/T - 1/ T_ref) * FARADAY / RM * exp_term * R_ref * ASR_nom
        J['ASR', 'R_ref_asr'] = exp_term * ASR_nom
        J['ASR', 'R_c_asr'] = ASR_nom
        J['ASR', 'ASR_nom_750'] = exp_term * R_ref + R_c

        J['eta_asr', 'T_PEN'] =             J['ASR', 'T'] * i
        J['eta_asr', 'i'] =             (exp_term * R_ref + R_c) * ASR_nom
        J['eta_asr', 'T_ref_asr'] =     J['ASR', 'T_ref_asr'] * i
        J['eta_asr', 'E_a_asr'] =       J['ASR', 'E_a_asr'] * i
        J['eta_asr', 'R_ref_asr'] =     J['ASR', 'R_ref_asr'] * i
        J['eta_asr', 'R_c_asr'] =       J['ASR', 'R_c_asr'] * i
        J['eta_asr', 'ASR_nom_750'] =   J['ASR', 'ASR_nom_750'] * i


class NernstPotential(om.ExplicitComponent):
    """
    Calculates the NernstPotential
    In the loss term there seems to be 2 terms that cancel in the ssimullink mdoel pfuel/p0, is that maybe an error??
    """
    def setup(self):
        self.add_input('GR_ecr_a', val=1, units='J/mol/K/K', desc='1st Fitting Parameter for gibbs reaction enthalpy')
        self.add_input('GR_ecr_b', val=1, units='J/mol/K', desc='2nd Fitting Parameter for gibbs reaction enthalpy')
        self.add_input('GR_ecr_c', val=1, units='J/mol', desc='3rd Fitting Parameter for gibbs reaction enthalpy')

        self.add_input('T', val=1000., units='K', desc= 'Temperature')

        self.add_input('p0', val=10e5, units='Pa', desc='static pressure of flow')
        self.add_input('p_fuel', val=10e5, units='Pa', desc='static pressure of fuel flow')
        self.add_input('p_air', val=10e5, units='Pa', desc='static pressure of air flow')

        self.add_input('a_H2', val=0.5, units=None, desc='Activity of H2')
        self.add_input('a_O2', val=0.1, units=None, desc='Activity of O2')
        self.add_input('a_H2O', val=0.1, units=None, desc='Activity of H2O')


        self.add_output('E_Nernst', val=0.2, units='V', desc='Nernst Potential')
        self.add_output('E_OCV', val = 1.23, units='V',desc='Open circuit voltage Nernst Potential' )

        self.declare_partials('E_Nernst', ['T', 'p0','p_air', 'a_H2O', 'a_O2', 'a_H2', 'GR_ecr_a', 'GR_ecr_b', 'GR_ecr_c'])
        self.declare_partials('E_OCV', ['T', 'GR_ecr_a', 'GR_ecr_b', 'GR_ecr_c'])
        

    def compute(self, inputs, outputs):
        # E_nernst = E0 - E_loss
        # with E0 = -(GR_ecr,a * T² + GR_ecr,b * T + GR_ecr,c) / (2F)
        # with E_loss = Rm * T / (2F) * ln( a_H2O / (a_H2 * sqrt(a_O2 * p_air / p_theta)))
        GR_ecr_a = inputs['GR_ecr_a']
        GR_ecr_b = inputs['GR_ecr_b']
        GR_ecr_c = inputs['GR_ecr_c']
        a_H2 = inputs['a_H2']
        a_H2O = inputs['a_H2O']
        a_O2 = inputs['a_O2']
        p_air = inputs['p_air']
        p_fuel = inputs['p_fuel']
        p0 = inputs['p0']
        T = inputs['T']

        outputs['E_OCV'] = E0 = -(GR_ecr_a * T**2 + GR_ecr_b * T + GR_ecr_c) / (2*FARADAY)
        E_losses = RM * T / (2*FARADAY) * np.log(a_H2O / (a_H2 * np.sqrt(a_O2 * p_air / p0)))
        outputs['E_Nernst'] = E0 - E_losses

    def compute_partials(self, inputs , J):
        GR_ecr_a = inputs['GR_ecr_a']
        GR_ecr_b = inputs['GR_ecr_b']
        GR_ecr_c = inputs['GR_ecr_c']
        a_H2 = inputs['a_H2']
        a_H2O = inputs['a_H2O']
        a_O2 = inputs['a_O2']
        p_air = inputs['p_air']
        p_fuel = inputs['p_fuel']
        p0 = inputs['p0']
        T = inputs['T']

        ln_term = np.log(a_H2O / (a_H2 * np.sqrt(a_O2 * p_air / p0)))

        J['E_Nernst', 'T'] =       - RM / (2*FARADAY) * ln_term - (2*GR_ecr_a * T + GR_ecr_b) / (2 * FARADAY)
        J['E_Nernst', 'p0'] =      - T * RM / (2*FARADAY) * 1 / (2*p0)
        J['E_Nernst', 'a_H2O'] =   - T * RM / (2*FARADAY) * 1 / a_H2O
        J['E_Nernst', 'a_H2'] =    - T * RM / (2*FARADAY) * (-1 / a_H2)
        J['E_Nernst', 'a_O2'] =    - T * RM / (2*FARADAY) * (- 1 / (2*a_O2))
        J['E_Nernst', 'p_air'] =   - T * RM / (2*FARADAY) * (- 1 / (2*p_air))
        J['E_Nernst', 'GR_ecr_a'] = - T**2 / (2*FARADAY)
        J['E_Nernst', 'GR_ecr_b'] = - T / (2*FARADAY)
        J['E_Nernst', 'GR_ecr_c'] = - 1 / (2*FARADAY)

        J['E_OCV', 'T'] =           - (2*GR_ecr_a * T + GR_ecr_b) / (2 * FARADAY)
        J['E_OCV', 'GR_ecr_a'] =    - T**2 / (2*FARADAY)
        J['E_OCV', 'GR_ecr_b'] =    - T / (2*FARADAY)
        J['E_OCV', 'GR_ecr_c'] =    - 1 / (2*FARADAY)


class SpeciesPartialPressure(om.ExplicitComponent):
    """
    Extract species composition from CEA and calculate partial pressures/activities.
    Deprecated and not used anymore
    """
    
    def initialize(self):
        self.options.declare('thermo_obj', desc='CEA thermo object', recordable=False)
        self.options.declare('target_species', default=['H2', 'O2', 'N2', 'H2O', 'CO2'],
                            types=list, desc='Species to track')
    
    def setup(self):
        species_list = self.options['target_species']
        thermo_obj = self.options['thermo_obj']
        
        # The CEA internal products list
        all_species = thermo_obj.products  # All 15+ species CEA knows about
        
        # Create a mapping from species name to index
        self.species_idx = {}
        for species in species_list:
            if species in all_species:
                self.species_idx[species] = all_species.index(species)
        
        # Inputs: Get the mole fractions from CEA's ChemEq component
        self.add_input('n_species', val=np.ones(len(all_species)), 
                      desc='Mole fractions of ALL species from CEA')
        self.add_input('n_total', val=1.0, 
                      desc='Total moles')
        self.add_input('P_total', val=1.0, units='bar',
                      desc='Total pressure')
        
        # Outputs: Partial pressures and activities
        for species in species_list:
            if species in self.species_idx:
                self.add_output(f'P_{species}', val=0.1, units='bar',
                               desc=f'Partial pressure of {species}')
                self.add_output(f'a_{species}', val=0.1, 
                               desc=f'Activity of {species}')
        
        self.declare_partials('*', '*', method='fd')
    
    def compute(self, inputs, outputs):
        n_species = inputs['n_species']
        n_total = inputs['n_total']
        P_total = inputs['P_total']
        P_ref = 1.0  # bar
        
        # Calculate mole fractions
        x_i = n_species / n_total
        
        # Calculate partial pressures
        for species, idx in self.species_idx.items():
            P_i = x_i[idx] * P_total
            a_i = P_i / P_ref  # Activity for ideal gas
            
            outputs[f'P_{species}'] = P_i
            outputs[f'a_{species}'] = a_i
            

class CurrentDensityCalc(Element):
    """
    Computes current density i = I / (A / N_seg) = I * N_seg / A

    Inputs:
        I     : Total current [A]
        A     : Cross-sectional area [m**2]
        N_seg : Number of segments [-] (integer constant, default=8)

    Output:
        i     : Current density [A/m**2]
    """
    def initialize(self):
        self.options.declare('N_seg', default=10)

    def setup(self):
        self.add_input('I', val=100,    units='A',      desc='Current in the segment. Only > 0 for active segments')
        self.add_input('A', val=,       units='m**2',   desc='Cross-sectional area of current flux - perpendicular to FLow direction') # promote input

        
        self.add_output('i', val=, units='A/m/m',desc='Current density in segment')
        self.declare_partials('i', ['I','A'])

    def compute(self, inputs, outputs):
        I = inputs['I']
        A = inputs['A']
        N_seg = self.options['N_seg']
        outputs['i'] = I / (A / N_seg)

    def compute_partials(self, inputs, J):
        I = inputs['I']
        A = inputs['A']
        N_seg = self.options['N_seg']
        J['i', 'I'] = N_seg / A
        J['i', 'A'] = - N_seg * I / (A**2 )
        

class segment_calc

class ElectroChemistry(om.Group):
    """
    Assembly that models the electrochemical calculations
    -------------
    Design
    -------------
        inputs
        --------
        T
        p

        x_h2_in
        x_h2_out
        
        x_o2_in
        x_o2_out


        outputs
        --------
        V_cell
        ASR
        E_nernst
    """
    def initialize(self):
        self.options.declare('N_segments', default=10)
    def setup(self):
        N_seg = self.options['N_segments']
        self.add_subsystem('current_density', CurrentDensityCalc(N_seg=N_seg),
                           promotes_inputs=['I', 'A'],
                           promotes_outputs=['i'])

        self.add_subsystem('ASR', AreaSpecificResistanceOverpotential(),
                           promotes_inputs=['i', 'T_PEN'],
                           promotes_outputs=['ASR', 'eta_asr'])

        self.add_subsystem('Nernst_Calc', NernstThermo(),
                           promotes_inputs=['I', ('T', 'T_PEN'), ('P', 'P_cat')],
                           promotes_outputs=[('E_Nernst', 'U_Nernst'), ('E_OCV', 'U_OCV')])
        
        self.add_subsystem('voltage', VoltageCalc(),
                           promotes=['*'])
        
