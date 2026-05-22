import numpy as np

import openmdao as om
from pycycle.element_base import Element
from pycycle.thermo.thermo import ThermoAdd

from pycycle.flow_in import FlowIn



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