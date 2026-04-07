import numpy as np
import openmdao as om

from pycycle.element_base import Element
from pycycle.thermo.thermo import Thermo, ThermoAdd

from pycycle.flow_in import FlowIn

RM = 8.3145                    # J/ (mol K)
FARADAY =  96485.3321233100184 # C/mol = As/mol
