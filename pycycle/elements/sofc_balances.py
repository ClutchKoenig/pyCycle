import numpy as np

import openmdao.api as om
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

class heat_convection_channel(om.ExplicitComponent):
    """
    Class responsible for calculating the convective heat transfer 
    inside the Anode or Cathode flow channels. I.e. handles:
        - Anode or Cathode heat convection.

    Outputs Q_conv_anode [W] or Q_conv_cathode [W]

    Options
    -------
    electrode       : 'anode' | 'cathode'
    """
    def initialize(self):
        self.options.declare('electrode', values=['anode', 'cathode'])
    def setup(self):
        self.add_input('W_channel', units='kg/s')
        self.add_input('')


class heat_convection_electrode(om.ExplicitComponent):
    """
    Class capable of handling:
        
        - Heat convection to Anode or Cathode.

    Outputs Q_conv_PEN [W], Q_conv_IC [W], directly inputs to the Energy balances
    PENEnergyBalance(), ICEnergyBalance()

    Options 
    -------
    electrode       : 'anode'   | 'cathode'
    structure       : 'PEN'     | 'IC'  
    """
    def initialize(self):
        self.options.declare('electrode', values=['anode', 'cathode'])
        self.options.declare('structure', values=['PEN', 'IC'])

    def setup(self):
        self.add_input('T_channel_in',  units='K')
        self.add_input('T_channel_out', units='K',      desc='Outlet temperature of channel mixture - result of Implicit calculation in EnergyBalance()')
        self.add_input('T_struc_out',   units='K',      desc='Outlet temperature of either PEN or IC')

        self.add_input('Nu',            units=None,     desc = 'Nusselt Number for selected electrode')
        self.add_input('d_hyd',         units='m',      desc='characteristic length Channel')
        self.add_input('A',             units='m**2',   desc='Area between PEN and channel')
        self.add_input('lambda',        units='W/m/K',  desc='Thermal conductivity of channel mixture - obtained by thermal_conductivity()')
        self.add_input('A_ic_side',     units='m**2',   desc='Area between IC and channel')
        # TODO: look into the MR_an_aus switch

        self.add_output('Q_conv_',      units='W',      desc='Convective Heat transfer from selected structure to electrode flow channel')

        self.add_output('T_bulk',   units='K',          desc='Bulk mixture temperature. Needed for thermal_conductivity()')
        self.declare_partials('Q_conv_', ['T_channel_out', 'T_struc_out' ,'Nu', 'lambda', 'd_hyd', 'A', 'A_ic_side' ]) # T_fuel partial is =0 : left undeclared
        self.declare_partials('T_bulk', 'T_channel_in', val = 0.5)
        self.declare_partials('T_bulk', 'T_channel_out', val=0.5)

    def compute(self, inputs, outputs):
        electrode = self.options['electrode']
        structure = self.options['structure']
               
        T_bulk = (inputs['T_channel_in'] + inputs['T_channel_out']) / 2
        outputs['T_bulk'] = T_bulk
        alpha =  inputs['Nu'] * inputs['lambda'] / inputs['d_hyd']
        dT = inputs['T_channel_out'] - inputs['T_struc_out']

        if structure == 'PEN':
            Q_conv_ = alpha * inputs['A'] * dT

        elif structure == 'IC':
            Q_conv_ = alpha * (inputs['A'] + inputs['A_ic_side']) * dT

        outputs['Q_conv_'] = Q_conv_

    def compute_partials(self, inputs, J):
        structure = self.options['structure']

        alpha =  inputs['Nu'] * inputs['lambda'] / inputs['d_hyd']
        dT = inputs['T_channel_out'] - inputs['T_struc_out']

        if structure == 'PEN':
            J['Q_conv_', 'T_channel_out']   = alpha * inputs['A'] 
            J['Q_conv_', 'T_struc_out']     = -alpha * inputs['A'] 
            J['Q_conv_', 'Nu']              = inputs['lambda'] / inputs['d_hyd'] * inputs['A'] * dT
            J['Q_conv_', 'd_hyd']           = - inputs['Nu'] * inputs['lambda'] / inputs['d_hyd']**2 * inputs['A'] * dT
            J['Q_conv_', 'A']               = alpha * dT
            J['Q_conv_', 'A_ic_side']       = 0
            J['Q_conv_',  'lambda']          = inputs['Nu'] / inputs['d_hyd'] * inputs['A'] * dT 

        elif structure == 'IC':
            A_combined = inputs['A'] + inputs['A_ic_side']
            J['Q_conv_', 'T_channel_out']   = alpha * (A_combined)
            J['Q_conv_', 'T_struc_out']     = -alpha * (A_combined)
            J['Q_conv_', 'Nu']              = inputs['lambda'] / inputs['d_hyd'] * (A_combined) * dT
            J['Q_conv_', 'd_hyd']           = - inputs['Nu'] * inputs['lambda'] / inputs['d_hyd']**2 * (A_combined) * dT
            J['Q_conv_', 'A']               = alpha * dT
            J['Q_conv_', 'A_ic_side']       = alpha * dT
            J['Q_conv_', 'lambda']          = inputs['Nu'] / inputs['d_hyd'] * A_combined * dT



class thermal_conductivity(om.ExplicitComponent):
    """
    Calculates the thermal conductivity lambda for a mixture with the 
    Mason-Saxena coefficents F12/F21 using the Wassilijewa rule
    
    lambda_mix = (x_H2 * lambda_H2) / (x_H2 + x_H2O * F12)
                +(x_H2O * lambda_H2O) / (x_H2O + x_H2 * F21) 
    Pure components are fitted linearly by lambda_H2 = a_H2 * T + b_H2

    F12 and F21 require dynamic viscosities which are fitted using a 
    4th-order polymomial
    """
    def initialize(self):
        self.options.declare('mixture', values=['air', 'H2/H2O'])

    def setup(self):
        self.add_input('T', units='K', desc='bulk mixture temperature - obtained from heat')
        
        self.add_input('x1', units=None,
                       desc='component 1 of mixture (O2/H2)')
        self.add_input('x2', units=None,
                       desc='component 2 of mixture (N2/H2O)')
        self.add_output('lambda_mix', units='W/m/K', 
                        desc='Thermal conductivity of declared mixture')
        self.declare_partials('lambda_mix', ['T', 'x1', 'x2'])


    def compute(self, inputs, outputs):
        mixture = self.options['mixture']
        if mixture == 'air':
            lambda_coeff1 = [a, b]
            eta_coeff1    = [A,B,C,D,E]
            lambda_coeff2 = [a, b]
            eta_coeff2    = [A,B,C,D,E]
            M_1 = 31.99880 # g/mol O2
            M_2 = 28.014    # g/mol N2

        elif mixture =='H2/H2O':
            lambda_coeff1 = [a, b]
            eta_coeff1   = [A,B,C,D,E]
            lambda_coeff2 = [a, b]
            eta_coeff2   = [A,B,C,D,E]
            M_1 = 2.016     # g/mol H2
            M_2 = 18.015    # g/mol H2O
        else: 
            raise ValueError("mixture must be 'air' or 'H2/H2O', but '{}' was given.".format(mixture))

        T = inputs['T']
        x1 = inputs['x1']
        x2 = inputs['x2']

        lambda_1 =  a1 * T + b1
        lambda_2 = a2 * T + b2

        eta_1 = A1 + B1*T + C1* T**2 + D1 * T**3 + E1 * T**4
        eta_2 = A2 + B2*T + C2* T**2 + D2 * T**3 + E2 * T**4

        F12 = (1 + np.sqrt(eta_1/eta_2) * (M_2/M_1)**(1/4) )**2 / np.sqrt(8 * (1 + M_1 / M_2))
        F21 = (1 + np.sqrt(eta_2/eta_1) * (M_1/M_2)**(1/4) )**2 / np.sqrt(8 * (1 + M_2 / M_1))
        
        outputs['lambda_mix'] = (x1 * lambda_1) / (x1 + x2 * F12) + (x2 * lambda_2) / (x2 + x1 * F21)
    
    def compute_partials(self, inputs, J):
        mixture = self.options['mixture']
        if mixture == 'air':
            # TODO: find coefficents in literature
            lambda_coeff1 = [a, b]
            eta_coeff1    = [A,B,C,D,E]
            lambda_coeff2 = [a, b]
            eta_coeff2    = [A,B,C,D,E]
            M_1 = 31.99880 # g/mol O2
            M_2 = 28.014    # g/mol N2

        elif mixture =='H2/H2O':
            # TODO: find coefficents in literature
            lambda_coeff1 = [a, b]
            eta_coeff1   = [A,B,C,D,E]
            lambda_coeff2 = [a, b]
            eta_coeff2   = [A,B,C,D,E]
            M_1 = 2.016     # g/mol H2
            M_2 = 18.015    # g/mol H2O
        else: 
            raise ValueError("mixture must be 'air' or 'H2/H2O', but '{}' was given.".format(mixture))

        T = inputs['T']
        x1 = inputs['x1']
        x2 = inputs['x2']

        lambda_1 =  a1 * T + b1 
        lambda_2 = a2*T + b2
        eta_1 = A1 + B1*T + C1* T**2 + D1 * T**3 + E1 * T**4
        eta_2 = A2 + B2*T + C2* T**2 + D2 * T**3 + E2 * T**4
        F12 = (1 + np.sqrt(eta_1/eta_2) * (M_2/M_1)**(1/4) )**2 / np.sqrt(8 * (1 + M_1 / M_2))
        F21 = (1 + np.sqrt(eta_2/eta_1) * (M_1/M_2)**(1/4) )**2 / np.sqrt(8 * (1 + M_2 / M_1))
        F12_denom = np.sqrt(8 * (1 + M_1 / M_2))
        F21_denom =np.sqrt(8 * (1 + M_2 / M_1))

        Deriv_denom1 = (x1 + x2 * F12)
        Deriv_denom2 = (x2 + x1 * F21)

        dL_dLam1 = x1 / Deriv_denom1
        dL_dLam2 = x2 / Deriv_denom2
        dLam1_dT = a1
        dLam2_dT = a2

        deta1_dT = B1 + 2*C1 * T + 3*D1 * T**2 + 4*E1* T**3
        deta2_dT = B2 + 2*C2 * T + 3*D2 * T**2 + 4*E2* T**3
        
        dF12_deta1 = 2 * (1 + np.sqrt(eta_1/eta_2) * (M_2/M_1)**(1/4) ) / F12_denom  * (1/ (2 * np.sqrt(eta_1 * eta_2)) * (M_2/M_1)**(1/4))
        dF12_deta2 = 2 * (1 + np.sqrt(eta_1/eta_2) * (M_2/M_1)**(1/4) ) / F12_denom * (- 1/2 * np.sqrt(eta_1) / np.sqrt(eta_2**3) * (M_2/M_1)**(1/4))

        dF21_deta1 = 2 * (1 + np.sqrt(eta_2/eta_1) * (M_1/M_2)**(1/4) ) / F21_denom * (- 1/2 * np.sqrt(eta_2) / np.sqrt(eta_1**3) * (M_1/M_2)**(1/4))
        dF21_deta2 = 2 * (1 + np.sqrt(eta_2/eta_1) * (M_1/M_2)**(1/4) ) / F21_denom * (1/ (2 * np.sqrt(eta_2 * eta_1)) * (M_1/M_2)**(1/4))

        dF12_dT = (dF12_deta1 * deta1_dT
                   + dF12_deta2 * deta2_dT)
        
        dF21_dT = (dF21_deta1 * deta1_dT
                   + dF21_deta2 * deta2_dT)
        dL_dF12 = -x1*x2*lambda_1 / Deriv_denom1**2
        dL_dF21 = - x1*x2*lambda_2 / Deriv_denom2**2

        J['lambda_mix', 'T'] = (dL_dLam1 * dLam1_dT
                                + dL_dLam2 * dLam2_dT
                                + dL_dF12 * dF12_dT
                                + dL_dF21 * dF21_dT) 
        J['lambda_mix', 'x1']= lambda_1 / (Deriv_denom1) - x1*lambda_1*1 /(Deriv_denom1)**2 - x2*lambda_2*F21 / Deriv_denom2**2
        J['lambda_mix', 'x2']= lambda_2 / (Deriv_denom2) - x2*lambda_2*1 /(Deriv_denom2)**2 - x1*lambda_1*F12 / Deriv_denom1**2