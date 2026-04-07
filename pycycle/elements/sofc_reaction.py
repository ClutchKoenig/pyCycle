import numpy as np
import openmdao.api as om
from pycycle.element_base import Element
from pycycle.thermo.thermo import Thermo, ThermoAdd

from pycycle.flow_in import FlowIn

RM = 8.3145                    # J/ (mol K)
FARADAY =  96485.3321233100184 # C/mol = As/mol

class VoltageCalc(om.ExplicitComponent):
    """
    Calculates the voltage of a SOFC MEA segment using 
    the Area specific ohmic resistance overpotential 'eta_ASR and Nernst Potential 'E'
    V_cell = E_nernst - eta_ASR 
    """
    def setup(self):
        self.add_input('U_Nernst', val=, units='V', desc=)
        self.add_input('eta_asr', val=, units='V', desc=)
        #self.add_input('i', val=1, units='A/m**2', desc=)
        
        self.add_output('V_cell', val=1.23, units='V', desc='Voltage')
        self.declare_partials('V_cell', ['U_Nernst', 'eta_asr'])

    def compute(self, inputs, outputs):
        U_Nernst = inputs['U_Nernst']
        eta = inputs['eta_asr']
        outputs['V_cell'] = U_Nernst - eta

    def compute_partials(self, inputs, J):
        U_Nernst = inputs['U_Nernst']
        eta = inputs['eta_asr']

        J['V_cell', 'U_Nernst'] = 1
        J['V_cell', 'eta_asr'] = 1


class AreaSpecificResistanceOverpotential(om.ExplicitComponent):
    """
    ASR = (exp((1/T - 1/T_ref_asr) * E_a_asr * F / RM) * R_ref_asr + R_c_asr) * ASR_nom_750
    eta = ASR * i
    """
    def setup(self):
        self.add_input('T', val=, units='K', desc=)
        self.add_input('T_ref_asr', val=, units='K', desc=)
        self.add_input('E_a_asr', val=, units='V', desc=)
        self.add_input('R_ref_asr', val=, units='ohm', desc=)
        self.add_input('R_c_asr', val=, units='ohm', desc=)
        self.add_input('ASR_nom_750', val=, units='ohm*m**2', desc=)
        self.add_input('i', val=, units='A/m/m', desc='Current density')
        # Einheiten inkonsistent? Was hat ASR_nom für ne Einheit?
        # eta_asr = A/ m**2 * ohm * ASR (ohm*m**2?)
        self.add_output('ASR', val=, units='ohm*m**2', desc=)
        self.add_output('eta_asr', val=, units='V', desc='')

        self.declare_partials('ASR', ['T', 'T_ref_asr','E_a_asr' , 'R_ref_asr', 'R_c_asr', 'ASR_nom_750'])
        self.declare_partials('eta_asr', ['T', 'i', 'T_ref_asr','E_a_asr' , 'R_ref_asr', 'R_c_asr', 'ASR_nom_750'])

    def compute(self, inputs, outputs):
        T = inputs['T']
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
        T = inputs['T']
        T_ref = inputs['T_ref_asr']
        E_a = inputs['E_a_asr']
        R_ref = inputs['R_ref_asr']
        R_c = inputs['R_c_asr']
        ASR_nom = inputs['ASR_nom_750']
        i = inputs['i']
        exp_term = np.exp((1/T - 1/ T_ref) * E_a * FARADAY / RM)
                          
        J['ASR', 'T'] = (- (1/T**2) * E_a * FARADAY / RM * exp_term * R_ref) * ASR_nom
        J['ASR', 'T_ref_asr'] = ( (1/T_ref**2) * E_a * FARADAY / RM  * exp_term * R_ref)* ASR_nom
        J['ASR', 'E_a_asr'] = (1/T - 1/ T_ref) * FARADAY / RM * exp_term * R_ref * ASR_nom
        J['ASR', 'R_ref_asr'] = exp_term * ASR_nom
        J['ASR', 'R_c_asr'] = ASR_nom
        J['ASR', 'ASR_nom_750'] = exp_term * R_ref + R_c

        J['eta_asr', 'T'] =             J['ASR', 'T'] * i
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
    def setup(self):
        self.add_input('I', val=, units='', desc=)
        self.add_input('A', val=, units='', desc=) # promote input
        self.add_input('N_seg', val=8,units=None)
        
        self.add_output('i', val=, units='',desc=)
        self.declare_partials('i', ['I','A','N_seg'])

    def compute(self, inputs, outputs):
        I = inputs['I']
        A = inputs['A']
        N_seg = inputs['N_seg']
        outputs['i'] = I / (A / N_seg)

    def compute_partials(self, inputs, J):
        I = inputs['I']
        A = inputs['A']
        N_seg = inputs['N_seg']
        J['i', 'I'] = N_seg / A
        J['i', 'A'] = - N_seg * I / (A**2 )
        J['i', 'N_seg'] = I / A

class segment_calc

class ElectroChemistry(Element):
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

class SOFC(Element):
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



class PENEnergyBalance(om.ImplicitComponent):
    """
    Steady-state energy balance for the PEN layer.
    Residual: Qdot_conv_an + Qdot_conv_cat - Qdot_chem - P_elec + Qdot_cond_left + Qdot_cond_right = 0
    T_cell is the implicit state solved by the Newton solver.
    """

    def setup(self):
        self.add_input('Qdot_conv_PEN_an',    val=0.0, units='W', desc='Convective heat from anode to PEN')
        self.add_input('Qdot_conv_PEN_cat',   val=0.0, units='W', desc='Convective heat from cathode to PEN')
        self.add_input('Qdot_chem',           val=0.0, units='W', desc='Chemical heat release')
        self.add_input('Qdot_conduct_PEN_left',  val=0.0, units='W', desc='Conductive heat flux, left')
        self.add_input('Qdot_conduct_PEN_right', val=0.0, units='W', desc='Conductive heat flux, right')
        self.add_input('I',      val=1.0,  units='A',   desc='Current')
        self.add_input('V_cell', val=0.7,  units='V',   desc='Cell voltage')
        self.add_input('N_cell', val=1.0,  units=None,  desc='Number of cells')

        self.add_output('T_cell',        val=1001.9, units='K', desc='PEN temperature')
        self.add_output('PEN_residuum',  val=0.0,   units='W', desc='Energy balance residual')

        self.declare_partials('T_cell', ['Qdot_conv_PEN_an', 'Qdot_conv_PEN_cat', 'Qdot_chem',
                                         'Qdot_conduct_PEN_left', 'Qdot_conduct_PEN_right',
                                         'I', 'V_cell', 'N_cell'])
        self.declare_partials('T_cell', 'T_cell', val=0.0)
        self.declare_partials('PEN_residuum', ['Qdot_conv_PEN_an', 'Qdot_conv_PEN_cat', 'Qdot_chem',
                                               'Qdot_conduct_PEN_left', 'Qdot_conduct_PEN_right',
                                               'I', 'V_cell', 'N_cell'])

    def apply_nonlinear(self, inputs, outputs, residuals):
        P_elec = inputs['I'] * inputs['V_cell'] * inputs['N_cell']
        energy_balance = (inputs['Qdot_conv_PEN_an'] +
                          inputs['Qdot_conv_PEN_cat'] -
                          inputs['Qdot_chem'] -
                          P_elec +
                          inputs['Qdot_conduct_PEN_left'] +
                          inputs['Qdot_conduct_PEN_right'])

        residuals['T_cell'] = energy_balance
        outputs['PEN_residuum'] = energy_balance

    def linearize(self, inputs, outputs, J):
        J['T_cell', 'Qdot_conv_PEN_an']      = 1.0
        J['T_cell', 'Qdot_conv_PEN_cat']     = 1.0
        J['T_cell', 'Qdot_chem']             = -1.0
        J['T_cell', 'Qdot_conduct_PEN_left'] = 1.0
        J['T_cell', 'Qdot_conduct_PEN_right']= 1.0
        J['T_cell', 'I']                     = -inputs['V_cell'] * inputs['N_cell']
        J['T_cell', 'V_cell']                = -inputs['I'] * inputs['N_cell']
        J['T_cell', 'N_cell']                = -inputs['I'] * inputs['V_cell']

        # PEN_residuum shares the same derivatives
        for key in ['Qdot_conv_PEN_an', 'Qdot_conv_PEN_cat', 'Qdot_chem',
                    'Qdot_conduct_PEN_left', 'Qdot_conduct_PEN_right', 'I', 'V_cell', 'N_cell']:
            J['PEN_residuum', key] = J['T_cell', key]

class ICEnergyBalance(om.ImplicitComponent):
    """
    Steady-state energy balance for the interconnect
    R = Qdot_conv_IC_an + Qdot_conv_IC_cat -  Q_dot_loss + Qdot_conduct_IC_left + Qdot_conduct_IC_right 
    """
    def setup(self):
        self.add_input('Qdot_conv_IC_an', val=0, units='W', desc= 'Convective heat from anode to IC')
        self.add_input('Qdot_conv_IC_cat', val=0, units='W', desc= 'Convective heat from cathode to IC')
        self.add_input('Q_dot_loss', val=0, units='W') # need a seperate component to calculate q_loss using Power and percent heat loss
        self.add_input('Qdot_conduct_IC_left', val=0, units='W')
        self.add_input('Qdot_conduct_IC_right', val=0, units='W')

        self.add_output('T_cell', val=1000, units='K')
        self.add_output('IC_residuum', val=0, units='W')
        self.declare_partials('T_cell', ['Qdot_conv_IC_an', 'Qdot_conv_IC_cat', 
                                         'Q_dot_loss', 'Qdot_conduct_IC_left', 'Qdot_conduct_IC_right'])
        self.declare_partials('T_cell', 'T_cell', val=0)
        self.declare_partials('IC_residuum', ['Qdot_conv_IC_an', 'Qdot_conv_IC_cat', 
                                              'Q_dot_loss', 'Qdot_conduct_IC_left', 'Qdot_conduct_IC_right'])
        
    def apply_nonlinear(self, inputs, outputs, residuals):
        energy_balance = (inputs['Qdot_conv_IC_an'] + inputs['Qdot_conv_IC_cat'] - 
                          inputs['Q_dot_loss'] + 
                          inputs['Qdot_conduct_IC_left'] + inputs['Qdot_conduct_IC_right'] )
        residuals['T_cell'] = energy_balance
        outputs['IC_residuum'] = energy_balance

    def linearize(self, inputs, outputs, J):
        J['T_cell', 'Qdot_conv_IC_an'] = 1
        J['T_cell', 'Qdot_conv_IC_cat'] = 1
        J['T_cell', 'Q_dot_loss'] = - 1
        J['T_cell', 'Qdot_conduct_IC_left'] = 1
        J['T_cell', 'Qdot_conduct_IC_left'] = 1
        
        for key in ['Qdot_conv_IC_an', 'Qdot_conv_IC_cat', 'Q_dot_loss', 
                    'Qdot_conduct_IC_left', 'Qdot_conduct_IC_right']:
            J['IC_residuum', key] = J['T_cell', key]