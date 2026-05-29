import numpy as np
import openmdao as om

from pycycle.element_base import Element
from pycycle.thermo.thermo import Thermo, ThermoAdd

from pycycle.flow_in import FlowIn

RM = 8.3145                    # J/ (mol K)
FARADAY =  96485.3321233100184 # C/mol = As/mol

class utilization(om.ExplicitComponent):
    """
    Calculates the Oxygen/Fuel utilization
    """

    def setup(self):
        self.add_input('n', units=None, desc='The base_thermo.thermo.n mole vector n [mol/g] ')
        