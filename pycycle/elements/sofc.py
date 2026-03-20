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
        self.add_input('T', val=, units='', desc=)
        self.add_input('I', val=, units='', desc=)
        self.add_input('', val=, units='', desc=)
        self.add_input('', val=, units='', desc=)
        self.add_input('', val=, units='', desc=)
        self.add_input('', val=, units='', desc=)
        self.add_input('', val=, units='', desc=)
        self.add_input('', val=, units='', desc=)


        self.add_output('V_cell', val=, units='V', desc='Voltage')

        self.declare_partials()

    def compute(self, inputs, outputs):
        




    def compute_partials(self, inputs, outputs):




class NernstPotential(om.ExplicitComponent):
    """
    Calculates the NernstPotential
    """
    def setup(self):
        self.add_input('GR_ecr_a', val=, units='J/mol/K/K', desc='1st Fitting Parameter for gibbs reaction enthalpy')
        self.add_input('GR_ecr_b', val=, units='J/mol/K', desc='2nd Fitting Parameter for gibbs reaction enthalpy')
        self.add_input('GR_ecr_c', val=, units='J/mol', desc='3rd Fitting Parameter for gibbs reaction enthalpy')

        self.add_input('T', val=1000., units='K', desc= 'Temperature')
        self.add_input('p', val=10e5, units='Pa', desc='static pressure of flow')
        self.add_input('a_H2', val=0.5, units=None, desc='Activity of H2')
        self.add_input('a_O2', val=0.1, units=None, desc='Activity of O2')
        self.add_input('a_H2O', val=0.1, units=None, desc='Activity of H2O')

        self.add_input('i', val=1, units='A/cm/cm',desc= 'current density' )

        self.add_output('E_Nerst', val=1.23, units='V', desc='Nernst Potential')

        self.declare_partials('E_Nernst', ['T', 'p', 'a_H2O', 'a_O2', 'a_H2O'])


    def compute(self, inputs, outputs):
        # E_nernst = E0 - E_loss
        # with E0 = -(GR_ecr,a * T² + GR_ecr,b * T + GR_ecr,c) / (2F)
        # with E_loss = Rm * T / (2F) * ln( a_H2O / (a_H2 * sqrt(a_O2 * p_air / p_theta)))

        
        i = inputs['i']
        T = inputs['T']
        p = inputs['p']
        outputs['E_Nernst'] = 


    def compute_partials(self, inputs , J):
        J['E_Nernst', 'T'] = 
        J['E_Nernst', 'p'] =
        J['E_Nernst', 'a_H2O'] =
        J['E_Nernst', 'a_H2'] =
        J['E_Nernst', 'a_O2'] =
        # Insert partial derivss

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
            

class CurrentDensity(Element):
    """


    """
    def setup(self):
        self.add_input('I', val=, units='', desc=)
    

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