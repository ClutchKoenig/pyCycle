import numpy as np
import openmdao.api as om

from pycycle.element_base import Element
from pycycle.thermo.thermo import Thermo, ThermoAdd

from pycycle.flow_in import FlowIn

RM = 8.3145                    # J/ (mol K)
FARADAY =  96485.3321233100184 # C/mol = As/mol

class SpeciesUtilization(om.ExplicitComponent):
    """
    Calculates utilization for a single species (H2 or O2).

    Instantiate twice in the parent group — once with species='H2' connected
    to the anode stream, once with species='O2' connected to the cathode stream.

    Options
    -------
    thermo  : Properties   species data for the relevant electrode
    species : 'H2' | 'O2'
    """
    def initialize(self):
        self.options.declare('thermo', desc='thermodynamic data object', recordable=False)
        self.options.declare('species', values=['H2', 'O2'])

    def setup(self):
        thermo  = self.options['thermo']
        species = self.options['species']

        self._sp_idx = (thermo.products.index(species)
                        if species in thermo.products else None)
        n = thermo.num_prod

        self.add_input('n_in',   units='mol/s', desc='Total molar flow of this electrode stream')
        self.add_input('x_in',   units=None, val=np.ones(n), desc='Mole fraction vector for this stream')
        self.add_input('I',      units='A')
        self.add_input('n_cell', units=None)

        self.add_output(f'{species}_utilization',  units=None)
        self.add_output(f'{species}_consumed_mol', units='mol/s')

        util = f'{species}_utilization'
        cons = f'{species}_consumed_mol'

        self.declare_partials(util, ['n_in', 'I', 'n_cell'])
        if self._sp_idx is not None:
            self.declare_partials(util, 'x_in', rows=[0], cols=[self._sp_idx])

        self.declare_partials(cons, 'n_in', val=0)
        self.declare_partials(cons, 'x_in', val=0)
        self.declare_partials(cons, ['I', 'n_cell'])

    def compute(self, inputs, outputs):
        species = self.options['species']
        n_in   = inputs['n_in']
        x_in   = inputs['x_in']
        n_cell = inputs['n_cell']
        I      = inputs['I']

        F = 4 * FARADAY if species == 'O2' else 2 * FARADAY
        n_consumed = I * n_cell / F

        outputs[f'{species}_consumed_mol'] = n_consumed
        if self._sp_idx is not None:
            outputs[f'{species}_utilization'] = n_consumed / (n_in * x_in[self._sp_idx])

    def compute_partials(self, inputs, J):
        species = self.options['species']
        n_in   = inputs['n_in']
        x_in   = inputs['x_in']
        n_cell = inputs['n_cell']
        I      = inputs['I']

        F = 4 * FARADAY if species == 'O2' else 2 * FARADAY
        n_consumed = I * n_cell / F

        cons = f'{species}_consumed_mol'
        util = f'{species}_utilization'

        J[cons, 'I']      = n_cell / F
        J[cons, 'n_cell'] = I / F

        if self._sp_idx is not None:
            x_sp = x_in[self._sp_idx]
            J[util, 'I']      = n_cell / F / (n_in * x_sp)
            J[util, 'n_cell'] = I / F / (n_in * x_sp)
            J[util, 'x_in']   = -n_consumed / (n_in * x_sp**2)
            J[util, 'n_in']   = -n_consumed / (n_in**2 * x_sp)

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