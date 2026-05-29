import numpy as np

import openmdao.api as om
from pycycle.element_base import Element
from pycycle.thermo.thermo import ThermoAdd
from pycycle.thermo.cea.species_data import Properties, janaf
from pycycle.flow_in import FlowIn



class PENEnergyBalance(om.ImplicitComponent):
    """
    Steady-state energy balance for the PEN layer.
    Residual: Qdot_conv_an + Qdot_conv_cat - Qdot_chem - P_elec + Qdot_cond_left + Qdot_cond_right = 0
    T_cell is the implicit state solved by the Newton solver.
    """
    def initialize(self):
        self.options.declare('seg_type', values=['passive', 'active'])

    def setup(self):

        self.add_input('Qdot_conv_PEN_an',    val=0.0, units='W', desc='Convective heat from anode to PEN')
        self.add_input('Qdot_conv_PEN_cat',   val=0.0, units='W', desc='Convective heat from cathode to PEN')
        
        self.add_input('Qdot_conduct_PEN_left',  val=0.0, units='W', desc='Conductive heat flux, left')
        self.add_input('Qdot_conduct_PEN_right', val=0.0, units='W', desc='Conductive heat flux, right')
        
        self.add_input('N_cell', val=1.0,  units=None,  desc='Number of cells')

        self.add_output('T_PEN', val=1001.9, units='K', desc='PEN temperature')
        #self.add_output('PEN_residuum',  val=0.0,   units='W', desc='Energy balance residual')
        passive_inputs = ['Qdot_conv_PEN_an', 'Qdot_conv_PEN_cat', 
                        'Qdot_conduct_PEN_left', 'Qdot_conduct_PEN_right']
        
        seg_type = self.options['seg_type']

        if seg_type == 'active':
            self.add_input('Qdot_chem',val=0.0, units='W', desc='Chemical heat release')
            self.add_input('I',        val=1.0, units='A', desc='Current')
            self.add_input('V_cell',   val=0.7, units='V', desc='Cell voltage')
            active_inputs = ['Qdot_chem', 'I', 'V_cell', 'N_cell']
            self.declare_partials('T_PEN', active_inputs)

        self.declare_partials('T_PEN', passive_inputs)
        self.declare_partials('T_PEN', 'T_PEN', val=0.0)
        # self.declare_partials('PEN_residuum', ['Qdot_conv_PEN_an', 'Qdot_conv_PEN_cat', 'Qdot_chem',
        #                                        'Qdot_conduct_PEN_left', 'Qdot_conduct_PEN_right',
        #                                        'I', 'V_cell', 'N_cell'])

    def apply_nonlinear(self, inputs, outputs, residuals):
        seg_type = self.options['seg_type']

        Q_conv_an = inputs['Qdot_conv_PEN_an']
        Q_conv_cat = inputs['Qdot_conv_PEN_cat']
        Q_cond_left = inputs['Qdot_conduct_PEN_left']
        Q_cond_right = inputs['Qdot_conduct_PEN_right']

        if seg_type == 'active':
            P_elec = inputs['I'] * inputs['V_cell'] * inputs['N_cell']
            Q_chem = inputs['Qdot_chem']

        elif seg_type == 'passive':
            P_elec = 0
            Q_chem = 0

        energy_balance = ( Q_conv_cat + Q_conv_an - Q_chem - P_elec + Q_cond_left +Q_cond_right)
        residuals['T_PEN'] = energy_balance

        #outputs['PEN_residuum'] = energy_balance

    def linearize(self, inputs, outputs, J):
        seg_type = self.options['seg_type']

        J['T_PEN', 'Qdot_conv_PEN_an']       = 1.0
        J['T_PEN', 'Qdot_conv_PEN_cat']      = 1.0
        J['T_PEN', 'Qdot_conduct_PEN_left']  = 1.0
        J['T_PEN', 'Qdot_conduct_PEN_right'] = 1.0

        if seg_type == 'active':
            J['T_PEN', 'N_cell']    = -inputs['I'] * inputs['V_cell']
            J['T_PEN', 'V_cell']    = -inputs['I'] * inputs['N_cell']
            J['T_PEN', 'I']         = -inputs['V_cell'] * inputs['N_cell']
            J['T_PEN', 'Qdot_chem'] = -1.0

        # # PEN_residuum shares the same derivatives
        # for key in ['Qdot_conv_PEN_an', 'Qdot_conv_PEN_cat', 'Qdot_chem',
        #             'Qdot_conduct_PEN_left', 'Qdot_conduct_PEN_right', 'I', 'V_cell', 'N_cell']:
        #     J['PEN_residuum', key] = J['T_cell', key]



class ICEnergyBalance(om.ImplicitComponent):
    """
    Steady-state energy balance for the interconnect
    R = Qdot_conv_IC_an + Qdot_conv_IC_cat -  Q_dot_loss + Qdot_conduct_IC_left + Qdot_conduct_IC_right
    
    For the IC the segment type doesnt matter. No chemical reaction or Current Flows.
    I think a source term for joule heating since current is flowing in the IC. 
    Q_joule = I**2 * A_cell * ASR_IC ? 
    """
    def initalize(self):
        self.options.declare('seg_type', default='passive', values=['active', 'passive'])
    def setup(self):
        seg_type = self.options['seg_type']

        self.add_input('Qdot_conv_IC_an', val=0, units='W', desc= 'Convective heat from anode to IC')
        self.add_input('Qdot_conv_IC_cat', val=0, units='W', desc= 'Convective heat from cathode to IC')
        self.add_input('Q_dot_loss', val=0, units='W') # need a seperate component to calculate q_loss using Power and percent heat loss
        self.add_input('Qdot_conduct_IC_left', val=0, units='W')
        self.add_input('Qdot_conduct_IC_right', val=0, units='W')

        self.add_output('T_IC', val=1000, units='K')
        self.add_output('IC_residuum', val=0, units='W')
        self.declare_partials('T_IC', ['Qdot_conv_IC_an', 'Qdot_conv_IC_cat', 
                                         'Q_dot_loss', 'Qdot_conduct_IC_left', 'Qdot_conduct_IC_right'])
        self.declare_partials('T_IC', 'T_cell', val=0)
        self.declare_partials('IC_residuum', ['Qdot_conv_IC_an', 'Qdot_conv_IC_cat', 
                                              'Q_dot_loss', 'Qdot_conduct_IC_left', 'Qdot_conduct_IC_right'])
        
    def apply_nonlinear(self, inputs, outputs, residuals):
        energy_balance = (inputs['Qdot_conv_IC_an'] + inputs['Qdot_conv_IC_cat'] - 
                          inputs['Q_dot_loss'] + 
                          inputs['Qdot_conduct_IC_left'] + inputs['Qdot_conduct_IC_right'] )
        residuals['T_IC'] = energy_balance
        outputs['IC_residuum'] = energy_balance

    def linearize(self, inputs, outputs, J):
        J['T_IC', 'Qdot_conv_IC_an'] = 1
        J['T_IC', 'Qdot_conv_IC_cat'] = 1
        J['T_IC', 'Q_dot_loss'] = - 1
        J['T_IC', 'Qdot_conduct_IC_left'] = 1
        J['T_IC', 'Qdot_conduct_IC_left'] = 1
        
        for key in ['Qdot_conv_IC_an', 'Qdot_conv_IC_cat', 'Q_dot_loss', 
                    'Qdot_conduct_IC_left', 'Qdot_conduct_IC_right']:
            J['IC_residuum', key] = J['T_IC', key]

class ChannelEnergyBalance(om.ImplicitComponent):
    """
    Steady-state energy balance for the anode or cathode flow channel.

    Residual:
        R = Q_conv_channel - Q_conv_PEN - Q_conv_IC [- Q_loss] = 0

    where Q_conv_channel = ṁ·(h_out - h_in)  is supplied by heat_convection_channel,
    and Q_conv_PEN / Q_conv_IC are the wall-to-gas convective heat fluxes.

    T_channel is the implicit state.  It does not appear directly in R — the
    sensitivity dR/dT_channel = ṁ·cp is seen by the group-level Newton solver
    through the connected Thermo and heat_convection_channel components.

    Options
    -------
    electrode : 'anode' | 'cathode'
    loss      : bool    if True, a heat-loss term Q_loss must be connected
    """

    def initialize(self):
        self.options.declare('electrode', values=['anode', 'cathode'])
        self.options.declare('loss', default=False, desc='Include external heat loss term')

    def setup(self):
        self.add_input('Q_conv_channel', val=0.0, units='W',
                       desc='Enthalpy flux absorbed by channel gas: ṁ·(h_out - h_in)')
        self.add_input('Q_conv_PEN',     val=0.0, units='W',
                       desc='Convective heat from PEN wall to channel gas')
        self.add_input('Q_conv_IC',      val=0.0, units='W',
                       desc='Convective heat from IC wall to channel gas')

        if self.options['loss']:
            self.add_input('Q_loss',     val=0.0, units='W',
                           desc='External heat loss from channel')

        self.add_output('T_channel', val=1000.0, units='K',
                        desc='Channel outlet temperature (implicit state)')

        active_inputs = ['Q_conv_channel', 'Q_conv_PEN', 'Q_conv_IC']
        if self.options['loss']:
            active_inputs.append('Q_loss')

        self.declare_partials('T_channel', active_inputs)
        self.declare_partials('T_channel', 'T_channel', val=0.0)

    def apply_nonlinear(self, inputs, outputs, residuals):
        R = (inputs['Q_conv_channel']
             - inputs['Q_conv_PEN']
             - inputs['Q_conv_IC'])

        if self.options['loss']:
            R -= inputs['Q_loss']

        residuals['T_channel'] = R

    def linearize(self, inputs, outputs, J):
        J['T_channel', 'Q_conv_channel'] =  1.0
        J['T_channel', 'Q_conv_PEN']     = -1.0
        J['T_channel', 'Q_conv_IC']      = -1.0
        J['T_channel', 'T_channel']      =  0.0  # sensitivity via group Newton only

        if self.options['loss']:
            J['T_channel', 'Q_loss']     = -1.0
    def linearize(self, inputs, outputs, J):        
        J['T_channel', 'Q_conv_channel'] = 1
        J['T_channel', 'Q_conv_PEN'] = -1
        J['T_channel', 'Q_conv_IC'] = -1
        J['T_channel', 'T_channel'] = 0

class SpeciesFlowCalc(om.ExplicitComponent):
    """
    This class is not really a mass balance. The consistency of molar flows is guaranteed by construction.
    A sanity check will be implemented:

        sanity: total system mass in = total mass out
        assert np.isclose(W_an_in + W_cat_in, W_an_out + W_cat_out, rtol=1e-8)

    This class rather takes inputs from anode_out/cathode_out.base_thermo.n and computes the molar fractions xi for
    each species of the gas channel.

    Outputs should connect to the 
    """
    def initialize(self):
        self.options.declare('thermo', desc='thermodynamic data object', recordable=False)
        self.options.declare('spec', default=janaf, recordable=False)
        self.options.declare('electrode', values=['anode', 'cathode'])

        self.options.declare('composition', desc='Element mole-ratio dict for the inlet flow')
    def setup(self):
        thermo = self.options['thermo']
        spec = self.options['spec']

        self._n = thermo.num_prod
        num_elem = thermo.num_element
        products = thermo.products

        idx = np.arange(self._n)

        self.add_input('n_i', val= np.ones(self._n), units=None, 
                       desc= 'array containing all n [mol/g] for all possible species ' \
                       'given the channel elemental composition.') # ChemEq doesnt have units declared for the moles :/
        self.add_input('n_moles', val= 1e-2, units=None, desc='Total molar flow')
        self.add_input('W', val=1.0, units='g/s')

        self.add_output('n_flow', val=1e-2, units= 'mol/s')
        self.add_output('x_i', val=np.ones(self._n), units=None,
                        desc='array containing all mole fractions')
        
        self.declare_partials('x_i', 'n_i',     rows=idx, cols=idx)
        self.declare_partials('x_i', 'n_moles', rows=idx, cols=np.zeros(self._n, dtype=int))
        self.declare_partials('x_i', 'W',       val=0)
        
        self.declare_partials('n_flow', 'n_i', val= 0)
        self.declare_partials('n_flow', 'n_moles')
        self.declare_partials('n_flow', 'W')

    def compute(self, inputs, outputs):
        n_i = inputs['n_i']
        n = inputs['n_moles']
        W = inputs['W']
        x_i = n_i / n

        outputs['n_flow'] = n * W
        outputs['x_i'] = x_i

    
    def compute_partials(self, inputs, J):
        J['x_i', 'n_i'] = 1.0 / inputs['n_moles']# Since the shape is already defined in declare_partials this part is unneccessary_ * np.ones(self._n)
        J['x_i', 'n_moles'] = - inputs['n_i'] / inputs['n_moles']**2

        J['n_flow', 'n_moles']  = inputs['W']
        J['n_flow', 'W']        = inputs['n_moles']