import numpy as np
import openmdao.api as om

from pycycle.element_base import Element 
from pycycle.thermo.thermo import Thermo, THermoAdd

from pycycle.flow_in import FlowIn


class sofc_segment(om.ExplicitComponent):
    """
    Represents 1 segment of the discretized SOFC Model. 
    Can be either part of a passive assembly or an active part with reactions involved.

    """
    def initialize(self):
        self.options.declare('type', default='passive', desc ='switches the massflow calculation and energy balances')
    
    def output_port_data(self):
        """Compute the output (element set?) anode and cathode"""
        inlet_comp = self.options['inlet_composition']

        return inlet_comp

    def setup(self):
        """ """ 
        segment_type = self.options['type']




        self.add_subsystem('PEN_balance', pen_balance(), promotes_inputs=['V'])
        self.add_subsystem('IC_balance', ic_balance(), promotes_inputs=[''])
        self.add_subsystem('cathode', cathode_balance(), promotes_inputs=[''])
        self.add_subsystem('anode', anode_balance(), promotes_inputs=[''])
        self.add_subsystem('mass_balance', mass_balance(), promotes_inputs=[''])
        if segment_type == 'active':
            self.add_subsystem('reaction', ElectroChemistry(), promotes_inputs['V'])
        return