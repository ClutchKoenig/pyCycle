import numpy as np
import openmdao.api as om

from pycycle.element_base import Element
from pycycle.thermo.thermo import Thermo, ThermoAdd

from pycycle.flow_in import FlowIn

RM = 8.3145                    # J/ (mol K)
FARADAY =  96485.3321233100184 # C/mol = As/mol

class SpeciesUtilization(om.ExplicitComponent):
    """
    Calculates the Oxygen/Fuel utilization
    """
    def initialize(self):
        self.options.declare('thermo', desc='thermodynamic data object', recordable=False)
    def setup(self):
        thermo = self.options['thermo']
        self._H2idx = thermo.products.index('H2') if 'H2' in thermo.products else None
        self._O2idx = thermo.products.index('O2') if 'O2' in thermo.products else None
        idx = np.arange(thermo.num_prod)
        self.add_input('n_in', units='mol/s', desc='Total molar flow')
        self.add_input('x_in', units= None, val=np.ones(thermo.num_prod), desc='Vector of molar fractions for each species')
        self.add_input('I',     units='A')
        self.add_input('n_cell', units=None)

        self.add_output('O2_utilization', units=None)
        self.add_output('H2_utilization', units=None)

        self.add_output('O2_consumed_mol', units='mol/s')
        self.add_output('H2_consumed_mol', units='mol/s')

        self.declare_partials('O2_utilization', ['n_in', 'I', 'n_cell'])
        self.declare_partials('O2_utilization', 'x_in', rows=[0], cols=[self._O2idx])
        self.declare_partials('H2_utilization', ['n_in', 'I', 'n_cell'])
        self.declare_partials('H2_utilization', 'x_in', rows=[0], cols=[self._H2idx])

        self.declare_partials('O2_consumed_mol', 'n_in', val=0)
        self.declare_partials('O2_consumed_mol', 'x_in', val=0)
        self.declare_partials('H2_consumed_mol', 'n_in', val=0)
        self.declare_partials('H2_consumed_mol', 'x_in', val=0)

        self.declare_partials('O2_consumed_mol', ['I', 'n_cell'])
        self.declare_partials('H2_consumed_mol', ['I', 'n_cell'])

    def compute(self, inputs, outputs):
        n_in = inputs['n_in']
        x_in = inputs['x_in']
        n_cell= inputs['n_cell']
        I = inputs['I']
        
        O2_used = I * n_cell / (4*FARADAY)
        H2_used = I * n_cell / (2*FARADAY)

        outputs['O2_consumed_mol'] = O2_used
        outputs['H2_consumed_mol'] = H2_used
        outputs['O2_utilization'] = (O2_used) / (n_in * x_in[self._O2idx])
        outputs['H2_utilization'] = (H2_used) / (n_in * x_in[self._H2idx])

    def compute_partials(self, inputs, J):
        n_in = inputs['n_in']
        x_in = inputs['x_in']
        n_cell= inputs['n_cell']
        I = inputs['I']
        
        O2_used = I * n_cell / (4*FARADAY)
        H2_used = I * n_cell / (2*FARADAY)

        J['O2_consumed_mol', 'I']        = n_cell / (4*FARADAY)
        J['O2_consumed_mol', 'n_cell']   = I / (4*FARADAY)
        J['O2_utilization', 'I']     = n_cell / (4*FARADAY) / (n_in * x_in[self._O2idx])
        J['O2_utilization', 'n_cell']= I / (4*FARADAY) / (n_in * x_in[self._O2idx])
        J['O2_utilization', 'x_in']  = - O2_used / (n_in * x_in[self._O2idx]**2)
        J['O2_utilization', 'n_in']  = - O2_used / (n_in**2 * x_in[self._O2idx])
        J['H2_consumed_mol', 'I']        = n_cell / (2*FARADAY)
        J['H2_consumed_mol', 'n_cell']   = I / (2*FARADAY)
        J['H2_utilization','I']      = n_cell / (2*FARADAY) / (n_in * x_in[self._H2idx])
        J['H2_utilization','n_cell'] = I / (2*FARADAY) / (n_in * x_in[self._H2idx])
        J['H2_utilization','x_in']   = -H2_used / (n_in * x_in[self._H2idx]**2)
        J['H2_utilization','n_in']   = -H2_used / (n_in**2 * x_in[self._H2idx])

class MassSanityCheck(om.ExplicitComponent):
    def setup(self):
        self.add_input('W_in_A' , val=1.0, units='kg/s', desc='Inlet massflow W [kg/s]')
        self.add_input('W_in_C' , val=1.0, units='kg/s', desc='Inlet massflow W [kg/s]')
        self.add_input('W_out_A', val=1.0, units='kg/s', desc='Outlet massflow W [kg/s]')
        self.add_input('W_out_C', val=1.0, units='kg/s', desc='Outlet massflow W [kg/s]')

        self.add_output('W_res', val=0., units='kg/s')



class BulkComposition(om.ExplicitComponent):
    def initialize(self):
        self.options.declare('thermo', desc='thermodynamic data object', recordable=False)
    def setup(self):
        thermo = self.options['thermo']
        self._n = thermo.num_prod
        self._H2idx = thermo.products.index('H2') if 'H2' in thermo.products else None
        self._O2idx = thermo.products.index('O2') if 'O2' in thermo.products else None
        self._H2Oidx = thermo.products.index('H2O') if 'H2O' in thermo.products else None
        idx = np.arange(self._n)

        self.add_input('x_in', val=np.ones(self._n), units=None)
        self.add_input('x_out',val=np.ones(self._n), units=None)

        self.add_output('x_bulk', val=np.ones(self._n), units=None)
        self.declare_partials('x_bulk', 'x_in',     rows=idx, cols=idx)
        self.declare_partials('x_bulk',  'x_out',    rows=idx, cols=idx)

    def compute(self, inputs, outputs):
        outputs['x_bulk'] = (inputs['x_in'] + inputs['x_out']) / 2

    def compute_partials(self, inputs, J):
        J['x_bulk', 'x_in']     = 0.5
        J['x_bulk', 'x_out']    = 0.5