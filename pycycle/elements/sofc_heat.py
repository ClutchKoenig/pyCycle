import openmdao.api as om
import numpy as np
from pycycle.thermo.thermo import ThermoAdd
from pycycle.thermo.cea.species_data import Properties, janaf

class HeatConvection(om.Group):
    """
    Group that handles convective heat transfer of all kinds.
    """
    #def initialize(self):
        #self.options.declare('N_seg', default=10)
    def setup(self):
        #N_seg = self.options['N_seg']

        self.add_subsystem('T_bulk', om.ExecComp(['T_bulk_A = (T_A_in + T_A_out) / 2',
                                                  'T_bulk_C = (T_C_in + T_C_out) / 2'],
                                                  T_A_in={'units':'K'}, T_C_in={'units':'K'}, T_bulk_A={'units':'K'},
                                                  T_C_out={'units':'K'}, T_A_out = {'units':'K'}, T_bulk_C={'units':'K'}),
                                                  promotes=['*'])
        
        self.add_subsystem('anode', HeatConvectionChannel(electrode='anode'),
                           promotes_inputs=[('W_in', 'W_in_A'), ('W_out', 'W_out_A'),
                                            ('h_in', 'h_A_in'), ('h_out', 'h_A_out')],
                           promotes_outputs=[('Q_conv_ch', 'Q_conv_A')])
        self.add_subsystem('cathode', HeatConvectionChannel(electrode='cathode'),
                           promotes_inputs=[('W_in', 'W_in_C'), ('W_out', 'W_out_C'),
                                            ('h_in', 'h_C_in'), ('h_out', 'h_C_out')],
                           promotes_outputs=[('Q_conv_ch', 'Q_conv_C')])
        
        self.add_subsystem('lambda_anode', ThermalConductivityMixture(mixture='H2/H2O'),
                           promotes_outputs=[('lambda_mix', 'lambda_A')])
        self.add_subsystem('lambda_cathode', ThermalConductivityMixture(mixture='air'),
                           promotes_outputs=[('lambda_mix', 'lambda_C')])
        
        
        self.add_subsystem('PEN_A', HeatConvectionElectrode(electrode='anode', structure='PEN'),
                           promotes_inputs = [('lambda', 'lambda_A'),
                                             ('T_channel_in', 'T_A_in'),
                                             ('T_channel_out', 'T_A_out')],
                           promotes_outputs =[('Q_conv_', 'Q_conv_PEN_A')]) 
                                             
        
        self.add_subsystem('PEN_C', HeatConvectionElectrode(electrode='cathode', structure='PEN'),
                           promotes_inputs = [('lambda', 'lambda_C'),
                                             ('T_channel_in', 'T_C_in'),
                                             ('T_channel_out', 'T_C_out')],
                           promotes_outputs=[('Q_conv_', 'Q_conv_PEN_C')])  
                                             
        self.add_subsystem('IC_A', HeatConvectionElectrode(electrode='anode', structure='IC'),
                           promotes_inputs = [('lambda', 'lambda_A'),
                                             ('T_channel_in', 'T_A_in'),
                                             ('T_channel_out', 'T_A_out')],
                           promotes_outputs=[('Q_conv_', 'Q_conv_IC_A')])
        
        self.add_subsystem('IC_C', HeatConvectionElectrode(electrode='cathode', structure='IC'),
                           promotes_inputs = [('lambda', 'lambda_C'),
                                             ('T_channel_in', 'T_C_in'),
                                             ('T_channel_out', 'T_C_out')],
                           promotes_outputs=[('Q_conv_', 'Q_conv_IC_C')])

        # Connect thermodynamic variables: 
        self.connect('T_bulk_A','lambda_anode.T')
        self.connect('T_bulk_C','lambda_cathode.T')
        # Not neccessary since they are already connected by promoted namespace
        # self.connect('lambda_A', 'PEN_A.lambda')
        # self.connect('lambda_C', 'PEN_C.lambda')
        # self.connect('lambda_A', 'IC_A.lambda')
        # self.connect('lambda_C', 'IC_C.lambda')


class HeatConvectionChannel(om.ExplicitComponent):
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
        self.add_input('W_in', units='kg/s',   desc= 'Mass flow of channel')
        self.add_input('h_in',      units='J/kg', )
        self.add_input('h_out',     units='J/kg')
        self.add_input('W_out', units='kg/s')

        self.add_output('Q_conv_ch', units='W')
        self.declare_partials('Q_conv_ch', ['W_in', 'W_out', 'h_in', 'h_out'])

    def compute(self, inputs, outputs):
        outputs['Q_conv_ch'] = inputs['W_out'] * inputs['h_out'] - inputs['W_in'] * inputs['h_in']

    def compute_partials(self, inputs, J):
        J['Q_conv_ch', 'W_out'] =  inputs['h_out']
        J['Q_conv_ch', 'h_out'] =  inputs['W_out']
        J['Q_conv_ch', 'W_in']  = -inputs['h_in']
        J['Q_conv_ch', 'h_in']  = -inputs['W_in']
        

class HeatConvectionElectrode(om.ExplicitComponent):
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



class ThermalConductivityMixture(om.ExplicitComponent):
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
        self.add_input('T', units='K', val=1000, 
                       desc='bulk mixture temperature - obtained from heat')
        
        self.add_input('x1', units=None, val = 0.5,
                       desc='component 1 of mixture (O2/H2)')
        self.add_input('x2', units=None, val = 0.5,
                       desc='component 2 of mixture (N2/H2O)')
        self.add_output('lambda_mix', units='W/m/K', val = 0.026,#for air
                        desc='Thermal conductivity of declared mixture')
        self.declare_partials('lambda_mix', ['T', 'x1', 'x2'])


    def compute(self, inputs, outputs):
        mixture = self.options['mixture']
        if mixture == 'air':
            lambda_coeff1 = {'a': 0, 'b': 6.116e-05, 'c': 1.128e-02}
            eta_coeff1    = {'A': -0.10257e-5,
                             'B':    0.92625e-7,
                             'C':   -0.80657e-10,
                             'D':    0.05113e-12,
                             'E':   -0.01295e-15}
            
            lambda_coeff2 = {'a': -0.0000000137, 'b': 0.0000740423, 'c':0.0052068912}
            eta_coeff2    = {'A':     -0.01020e-5,
                             'B':        0.74785e-7,
                             'C':       -0.59037e-10,
                             'D':        0.03230e-12,
                             'E':       -0.00673e-15}
            M_1 = 31.99880 # g/mol O2
            M_2 = 28.014    # g/mol N2

        elif mixture =='H2/H2O':
            lambda_coeff1 = {'a1': 3.799454545e-4, 'b1': 0.0804581818} #[a1, b1]

            eta_coeff1   = {'A':    0.18024e-5, 
                            'B':    0.27174e-7,
                            'C':   -0.13395e-10,
                            'D':    0.00585e-12,
                            'E':   -0.00104e-15}
            
            #[A,B,C,D,E]
            lambda_coeff2 = {'a2': -0.0086593279, 'b2': 0.0000747671, 'c2': 0.0000000294}

            eta_coeff2 = {'A':   -3.01414856711936e-06,
                          'B':    4.076228959276e-08}
            
            M_1 = 2.016     # g/mol H2
            M_2 = 18.015    # g/mol H2O
        else: 
            raise ValueError("mixture must be 'air' or 'H2/H2O', but '{}' was given.".format(mixture))

        T = inputs['T']
        x1 = inputs['x1']
        x2 = inputs['x2']

        thermal_cond = lambda x, coeffs: np.polyval(list(coeffs.values()), x)
        viscosity    = lambda x, coeffs: np.polyval(list(reversed(coeffs.values())), x)

        lambda_1 = thermal_cond(T, lambda_coeff1)
        lambda_2 = thermal_cond(T, lambda_coeff2)
        eta_1 = viscosity(T, eta_coeff1)
        eta_2 = viscosity(T, eta_coeff2)
        

        F12 = (1 + np.sqrt(eta_1/eta_2) * (M_2/M_1)**(1/4) )**2 / np.sqrt(8 * (1 + M_1 / M_2))
        F21 = (1 + np.sqrt(eta_2/eta_1) * (M_1/M_2)**(1/4) )**2 / np.sqrt(8 * (1 + M_2 / M_1))
        
        outputs['lambda_mix'] = (x1 * lambda_1) / (x1 + x2 * F12) + (x2 * lambda_2) / (x2 + x1 * F21)
    
    def compute_partials(self, inputs, J):
        mixture = self.options['mixture']
        if mixture == 'air':
            lambda_coeff1 = {'a': 0, 'b': 6.116e-05, 'c': 1.128e-02}
            eta_coeff1    = {'A': -0.10257e-5,
                             'B':    0.92625e-7,
                             'C':   -0.80657e-10,
                             'D':    0.05113e-12,
                             'E':   -0.01295e-15}
            
            lambda_coeff2 = {'a': -0.0000000137, 'b': 0.0000740423, 'c':0.0052068912}
            eta_coeff2    = {'A':     -0.01020e-5,
                             'B':        0.74785e-7,
                             'C':       -0.59037e-10,
                             'D':        0.03230e-12,
                             'E':       -0.00673e-15}
            M_1 = 31.99880 # g/mol O2
            M_2 = 28.014    # g/mol N2

        elif mixture =='H2/H2O':
            lambda_coeff1 = {'a1': 3.799454545e-4, 'b1': 0.0804581818} #[a1, b1]

            eta_coeff1   = {'A':    0.18024e-5, 
                            'B':    0.27174e-7,
                            'C':   -0.13395e-10,
                            'D':    0.00585e-12,
                            'E':   -0.00104e-15}
            
            #[A,B,C,D,E]
            lambda_coeff2 = {'a2': -0.0086593279, 'b2': 0.0000747671, 'c2': 0.0000000294}

            eta_coeff2 = {'A':   -3.01414856711936e-06,
                          'B':    4.076228959276e-08}
            
            M_1 = 2.016     # g/mol H2
            M_2 = 18.015    # g/mol H2O
        else: 
            raise ValueError("mixture must be 'air' or 'H2/H2O', but '{}' was given.".format(mixture))

        T = inputs['T']
        x1 = inputs['x1']
        x2 = inputs['x2']

        thermal_cond      = lambda x, coeffs: np.polyval(list(coeffs.values()), x)
        viscosity         = lambda x, coeffs: np.polyval(list(reversed(coeffs.values())), x)
        d_thermal_cond_dT = lambda x, coeffs: np.polyval(np.polyder(list(coeffs.values())), x)
        d_viscosity_dT    = lambda x, coeffs: np.polyval(np.polyder(list(reversed(coeffs.values()))), x)

        lambda_1 = thermal_cond(T, lambda_coeff1)
        lambda_2 = thermal_cond(T, lambda_coeff2)
        eta_1    = viscosity(T, eta_coeff1)
        eta_2    = viscosity(T, eta_coeff2)
        F12 = (1 + np.sqrt(eta_1/eta_2) * (M_2/M_1)**(1/4) )**2 / np.sqrt(8 * (1 + M_1 / M_2))
        F21 = (1 + np.sqrt(eta_2/eta_1) * (M_1/M_2)**(1/4) )**2 / np.sqrt(8 * (1 + M_2 / M_1))
        F12_denom = np.sqrt(8 * (1 + M_1 / M_2))
        F21_denom =np.sqrt(8 * (1 + M_2 / M_1))

        Deriv_denom1 = (x1 + x2 * F12)
        Deriv_denom2 = (x2 + x1 * F21)

        dL_dLam1 = x1 / Deriv_denom1
        dL_dLam2 = x2 / Deriv_denom2
        dLam1_dT = d_thermal_cond_dT(T, lambda_coeff1)
        dLam2_dT = d_thermal_cond_dT(T, lambda_coeff2)
        deta1_dT = d_viscosity_dT(T, eta_coeff1)
        deta2_dT = d_viscosity_dT(T, eta_coeff2)
        
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



class HeatConductionStructure(om.ExplicitComponent):
    """
    Lateral heat conduction (along flow direction) through the PEN or IC
    structure between two adjacent segments.

    Fourier (PEN, single cross-section):
        Q = lambda * (h_cell * delta) * n_cell * (N_seg / l_cell) * dT

    Fourier (IC, two cross-sections in series per Simulink):
        R1 = (l_cell/N_seg) / (lambda * h_cell * delta * n_cell)
        R2 = (l_cell/N_seg) / (lambda * width_IC * delta_ch * n_cell * N_IC_walls)
        Q  = dT / (R1 + R2)

    Sign convention (consistent with energy balance):
        Q_cond_struc_right = +Q   (heat entering the right-neighbour segment)
        Q_cond_struc_left  = -Q   (heat entering the left-neighbour segment)
        with dT = T_cell_left - T_cell_right

    Options
    -------
    structure  : 'PEN' | 'IC'
    N_segments : int   number of discretisation segments (default 10)
    """
    # TODO: check if partials are setup correctly!!!!
    def initialize(self):
        self.options.declare('structure', values=['PEN', 'IC'])
        self.options.declare('N_segments', default=10, desc='Number of discretisation segments')

    def setup(self):
        struct = self.options['structure']

        self.add_input('T_cell_left',      val=1000.0, units='K',     desc='Temperature of left-neighbour segment')
        self.add_input('T_cell_right',     val=1000.0, units='K',     desc='Temperature of right-neighbour segment')
        self.add_input('length_cell',      val=0.1,    units='m',     desc='Total active cell length')
        self.add_input('height_cell',      val=0.1,    units='m',     desc='Cell height (conduction cross-section)')
        self.add_input('thickness_struct', val=2e-4,   units='m',     desc='Thickness of PEN or IC layer')
        self.add_input('lambda_struct',    val=2.0,    units='W/m/K', desc='Thermal conductivity of PEN or IC')
        self.add_input('n_cell',           val=1.0,    units=None,    desc='Number of cells in stack')

        if struct == 'IC':
            self.add_input('width_IC',    val=1e-3, units='m',  desc='IC channel width')
            self.add_input('thickness_ch',val=1e-3, units='m',  desc='Channel height (IC side-wall area)')
            self.add_input('N_IC_walls',  val=2.0,  units=None, desc='Number of IC side-wall faces per channel')

        self.add_output('Q_cond_struc_right', val=0.0, units='W', desc='Heat flux into right-neighbour segment')
        self.add_output('Q_cond_struc_left',  val=0.0, units='W', desc='Heat flux into left-neighbour segment')

        common = ['T_cell_left', 'T_cell_right', 'length_cell', 'height_cell',
                  'thickness_struct', 'lambda_struct', 'n_cell']
        active = common + (['width_IC', 'thickness_ch', 'N_IC_walls'] if struct == 'IC' else [])

        self.declare_partials('Q_cond_struc_right', active)
        self.declare_partials('Q_cond_struc_left',  active)

    def compute(self, inputs, outputs):
        struct     = self.options['structure']
        N          = self.options['N_segments']
        dT         = inputs['T_cell_left'] - inputs['T_cell_right']
        l_cell     = inputs['length_cell']
        h_cell     = inputs['height_cell']
        delta      = inputs['thickness_struct']
        lambda_str = inputs['lambda_struct']
        n_cell     = inputs['n_cell']
        dx         = l_cell / N  # segment length

        if struct == 'PEN':
            A = h_cell * delta
            Q = lambda_str * A * n_cell / dx * dT
        else:
            width_IC  = inputs['width_IC']
            delta_ch  = inputs['thickness_ch']
            n_walls   = inputs['N_IC_walls']
            A1 = h_cell * delta
            A2 = width_IC * delta_ch
            R1 = dx / (lambda_str * A1 * n_cell)
            R2 = dx / (lambda_str * A2 * n_cell * n_walls)
            Q  = dT / (R1 + R2)

        outputs['Q_cond_struc_right'] =  Q
        outputs['Q_cond_struc_left']  = -Q

    def compute_partials(self, inputs, J):
        struct     = self.options['structure']
        N          = self.options['N_segments']
        dT         = inputs['T_cell_left'] - inputs['T_cell_right']
        l_cell     = inputs['length_cell']
        h_cell     = inputs['height_cell']
        delta      = inputs['thickness_struct']
        lambda_str = inputs['lambda_struct']
        n_cell     = inputs['n_cell']
        dx         = l_cell / N

        if struct == 'PEN':
            A = h_cell * delta
            G = lambda_str * A * n_cell / dx   # Q = G * dT

            dQ_dTleft  =  G
            dQ_dTright = -G
            dQ_dh      = lambda_str * delta  * n_cell / dx * dT
            dQ_ddelta  = lambda_str * h_cell * n_cell / dx * dT
            dQ_dlambda = A * n_cell / dx * dT
            dQ_dn      = lambda_str * A / dx * dT
            dQ_dl      = -lambda_str * A * n_cell * N / l_cell**2 * dT  # dx=l/N → G~N/l

        else:
            width_IC  = inputs['width_IC']
            delta_ch  = inputs['thickness_ch']
            n_walls   = inputs['N_IC_walls']

            A1      = h_cell * delta
            A2      = width_IC * delta_ch
            R1      = dx / (lambda_str * A1 * n_cell)
            R2      = dx / (lambda_str * A2 * n_cell * n_walls)
            R_total = R1 + R2
            Q       = dT / R_total

            dQ_dTleft  =  1.0 / R_total
            dQ_dTright = -1.0 / R_total

            # ∂Q/∂R_total = -Q / R_total
            dQ_dR = -Q / R_total

            # ∂R/∂param using R1 = dx/(lambda*A1*n), so ∂R1/∂x = -R1/x for any multiplicative x
            dQ_dh      = dQ_dR * (-R1 / h_cell)
            dQ_ddelta  = dQ_dR * (-R1 / delta)
            dQ_dlambda = dQ_dR * (-(R1 + R2) / lambda_str)
            dQ_dn      = dQ_dR * (-(R1 + R2) / n_cell)
            dQ_dl      = dQ_dR * ( (R1 + R2) / l_cell)  # dx=l/N → R~l → ∂R/∂l = R/l
            dQ_dwidth  = dQ_dR * (-R2 / width_IC)
            dQ_dch     = dQ_dR * (-R2 / delta_ch)
            dQ_dnwalls = dQ_dR * (-R2 / n_walls)

        # Q_right = +Q → same sign; Q_left = -Q → flipped sign
        for out, sign in [('Q_cond_struc_right', 1.0), ('Q_cond_struc_left', -1.0)]:
            J[out, 'T_cell_left']      = sign * dQ_dTleft
            J[out, 'T_cell_right']     = sign * dQ_dTright
            J[out, 'height_cell']      = sign * dQ_dh
            J[out, 'thickness_struct'] = sign * dQ_ddelta
            J[out, 'lambda_struct']    = sign * dQ_dlambda
            J[out, 'n_cell']           = sign * dQ_dn
            J[out, 'length_cell']      = sign * dQ_dl
            if struct == 'IC':
                J[out, 'width_IC']     = sign * dQ_dwidth
                J[out, 'thickness_ch'] = sign * dQ_dch
                J[out, 'N_IC_walls']   = sign * dQ_dnwalls


class HeatConduction(om.Group):
    """
    Organizes the heat conduction of the IC and PEN in Flow direction of the Channels. 
    """
    def initialize(self):
        self.options.declare('N_segments', default=10)

    def setup(self):
        N_segments = self.options['N_segments']
        self.add_subsystem('IC_Conduction', HeatConductionStructure(structure='IC', N_segments=N_segments),
                           promotes_inputs=[('T_cell_left',  'T_IC_left'),
                                            ('T_cell_right', 'T_IC_right')],  #TODO: promote geometric parameters
                           promotes_outputs=[('Q_cond_struc_right', 'Q_cond_IC_right'), 
                                             ('Q_cond_struc_left',  'Q_cond_IC_left')])
        self.add_subsystem('PEN_Conduction', HeatConductionStructure(structure='PEN', N_segments=N_segments),
                           promotes_inputs=[('T_cell_left',  'T_PEN_left'),
                                            ('T_cell_right', 'T_PEN_right')],  #TODO: promote geometric parameters
                           promotes_outputs=[('Q_cond_struc_right', 'Q_cond_PEN_right'),
                                             ('Q_cond_struc_left',  'Q_cond_PEN_left')])
                            #TODO: get all geometric parameters promoted to this level and higher!!
        self.add_subsystem('Q_conduc_left', om.ExecComp('Qc_left = Q_cond_IC_left + Q_cond_PEN_left',
                                                        Q_cond_IC_left= {'val':-1000, 'units':'W'},
                                                        Q_cond_PEN_left= {'val':-1000, 'units': 'W'},
                                                        Qc_left = {'val': -2000, 'units': 'W'}),
                                                        promotes_outputs=['Qc_left'])
        self.add_subsystem('Q_conduc_right', om.ExecComp('Qc_right = Q_cond_IC_right + Q_cond_PEN_right',
                                                        Q_cond_IC_right= {'val':1000, 'units':'W'},
                                                        Q_cond_PEN_right= {'val':1000, 'units': 'W'},
                                                        Qc_right = {'val': 2000, 'units': 'W'}),
                                                        promotes_outputs=['Qc_right'])
        
        self.connect('Q_cond_IC_left',  'Q_conduc_left.Q_cond_IC_left')
        self.connect('Q_cond_IC_right', 'Q_conduc_right.Q_cond_IC_right')

        self.connect('Q_cond_PEN_left',  'Q_conduc_left.Q_cond_PEN_left')
        self.connect('Q_cond_PEN_right', 'Q_conduc_right.Q_cond_PEN_right')

        self.linear_solver = om.LinearRunOnce()

