import openmdao.api as om
import numpy as np
from pycycle.elements.sofc.api import HeatConvection, HeatConvectionChannel, HeatConvectionElectrode, HeatConductionStructure, HeatConduction

# Input vars
Nu = 4.9
h_in_A = 100
h_out_A = 200
W_in_A = 1
W_out_A = 1


model = om.group()
model.add_subsystem('heat_convection_A', 
                    HeatConvectionChannel(electrode='anode'), 
                    promotes_inputs=['Nu', 'd_hyd', 'A', 'lambda', 'T_channel_in', 'T_struc_in'], 
                    promotes_outputs=[('Q_conv_ch', 'Q_conv_A')])

prob = om.Problem(model)
prob.setup()


prob.set_val('anode.W_in_A', 1)
prob.set_val('anode.W_out_A', 1)
prob.set_val('anode.h_in_A', 100)
prob.set_val('anode.h_out_A', 200)
prob.set_val('Nu', 4.9)

