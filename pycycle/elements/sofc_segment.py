import numpy as np
import openmdao.api as om

from pycycle.element_base import Element 
from pycycle.thermo.thermo import Thermo, ThermoAdd

from pycycle.flow_in import FlowIn
import pycycle.api as pyc
import pycycle.sofc_api as sofc


class sofc_segment(om.Group):
    """
    Represents 1 segment of the discretized SOFC Model. 
    Can be either part of a passive assembly or an active part with reactions involved.
    """

    def initialize(self):
        self.options.declare('type', default='passive', 
                             values=['passive', 'active'], 
                             desc ='switches the massflow calculation and energy balances')
        self.options.declare('N_segments', default=10, desc='Number of segments for current type')
        self.options

    def setup(self):
        """ """ 
        segment_type = self.options['type']
        N_seg = self.options['N_segments']
        
        if segment_type == 'active':
            self.add_subsystem('Electrochemistry', sofc.ElectroChemistry(N_segments= N_seg))
            self.connect('Electrochemistry.V_cell', 'PEN.V_cell')
            self.connect('Electrochemistry.Qdot_chem', 'PEN.Qdot_chem')

        self.add_subsystem('Convection', sofc.HeatConvection(),
                           promotes_inputs=['W_A', 'W_C', 'h_A_in' 'h_C_in',
                                            'T_A_in', 'T_A_out',
                                            'T_C_in', 'T_C_out'],
                           promotes_outputs=['Q_conv_A', 'Q_conv_C',
                                             'Q_conv_PEN_C', 'Q_conv_PEN_A',
                                             'Q_conv_IC_C', 'Q_conv_IC_A'])

        self.add_subsystem('Conduction', sofc.HeatConduction(N_segments=N_seg),
                           promotes_inputs=[], #TODO: get all geometric parameters promoted to this level and higher!!
                           promotes_outputs=['Qc_left', 'Qc_right'])

        self.add_subsystem('PEN',   sofc.PENEnergyBalance(), 
                           promotes_inputs=['N_cell'],
                           promotes_outputs=['T_PEN'])
        self.add_subsystem('IC',            sofc.ICEnergyBalance(), 
                           promotes_inputs=[''],
                           promotes_outputs=['T_IC'])
        
        self.add_subsystem('Cathode',       sofc.ChannelEnergyBalance(), 
                           promotes_inputs=[''])
        self.add_subsystem('Anode',         sofc.ChannelEnergyBalance(), promotes_inputs=[''])
        self.add_subsystem('molar_fractions',  sofc.ChannelMassBalance(), promotes_inputs=[''])

