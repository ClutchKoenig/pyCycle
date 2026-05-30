import numpy as np
import openmdao.api as om

from pycycle.elements.sofc_segment import SegmentSOFC
from pycycle.thermo.cea.species_data import Properties, janaf

# ================================================================
# Geometry from InputsESC_PaperData.m
# ================================================================
# Cell layers
delta_an  = 35e-6    # m  anode layer thickness
delta_el  = 90e-6    # m  electrolyte thickness
delta_cat = 35e-6    # m  cathode layer thickness
delta_pen = delta_an + delta_el + delta_cat   # 160 um

thickness_IC = 850e-6   # m  interconnect thickness

# Cell planform
length_active  = 0.09    # m  active cell length  (= length_ESC)
length_total   = 0.162   # m
length_passive = length_total - length_active   # 0.072 m
width_total    = 0.142   # m  (= width_total_ESC, used as conduction cross-section height)

N_seg_active  = 10
N_seg_passive = 8

# Channel geometry  (identical for anode and cathode)
N_channels   = 24
w_ch         = 0.0042    # m  channel width
h_ch         = 1e-3      # m  channel height
d_hyd        = 2*w_ch*h_ch / (w_ch + h_ch)   # hydraulic diameter ≈ 1.615e-3 m
A_cross      = w_ch * h_ch * N_channels       # total flow cross-section per electrode

Nu = 4.9   # Nusselt number (both electrodes)

# Per-segment convection areas  [active zone]
seg_len_act      = length_active / N_seg_active                    # 0.009 m
A_pen_act        = N_channels * w_ch * seg_len_act                 # PEN↔channel top wall area
A_ic_side_act    = N_channels * 2 * h_ch * seg_len_act             # IC side-wall area

# Per-segment convection areas  [passive zone]
seg_len_pas      = length_passive / N_seg_passive                  # 0.009 m
A_pen_pas        = N_channels * w_ch * seg_len_pas
A_ic_side_pas    = N_channels * 2 * h_ch * seg_len_pas

# Conduction
lambda_PEN  = 2.0    # W/(m·K)   (lambda_IC assumed equal for now)
lambda_IC   = 2.0    # W/(m·K)   — update if a separate value is available

# --- Compositions (must match what SegmentSOFC is built with) ---
anode_comp   = {'H': 2.0}
cathode_comp = {'N': 0.79, 'O': 0.21}

# Use Properties to get array sizes and b0 element vectors
anode_props   = Properties(janaf, init_elements=anode_comp)
cathode_props = Properties(janaf, init_elements=cathode_comp)

# anode products : ['H', 'H2']          -> n_prod = 2
# cathode products: ['N','NO','NO2','NO3','N2','O','O2'] -> n_prod = 7
print('Anode   products:', anode_props.products)
print('Cathode products:', cathode_props.products)

# --- Build problem ---
prob = om.Problem()
prob.model = SegmentSOFC(
    type='active',
    N_segments=1,
    anode_composition=anode_comp,
    cathode_composition=cathode_comp,
)

prob.setup(force_alloc_complex=True)

# ----------------------------------------------------------------
# Inlet boundary conditions
# ----------------------------------------------------------------
T_op = 1073.0   # K  (800 °C nominal SOFC operating temperature)

prob.set_val('W_in_A', 1e-4, units='kg/s')   # anode  mass flow
prob.set_val('W_in_C', 1e-3, units='kg/s')   # cathode mass flow

prob.set_val('T_A_in', T_op, units='K')
prob.set_val('T_C_in', T_op, units='K')

# Inlet specific enthalpies (set to 0 → solver finds T_out that closes energy balance)
prob.set_val('h_A_in', 0.0, units='J/kg')
prob.set_val('h_C_in', 0.0, units='J/kg')

prob.set_val('p_A_out', 1.01325e5, units='Pa')
prob.set_val('p_C_out', 1.01325e5, units='Pa')

# Element composition vectors (b0 from Properties, units mol/g)
prob.set_val('composition_in_A', anode_props.b0)
prob.set_val('composition_in_C', cathode_props.b0)

# Inlet mole fractions: pure H2 anode, 79/21 N2/O2 cathode
x_in_A = np.zeros(anode_props.num_prod)
x_in_A[anode_props.products.index('H2')] = 1.0          # pure H2

x_in_C = np.zeros(cathode_props.num_prod)
x_in_C[cathode_props.products.index('N2')] = 0.79       # air
x_in_C[cathode_props.products.index('O2')] = 0.21

prob.set_val('x_in_A', x_in_A)
prob.set_val('x_in_C', x_in_C)

# Total molar flows at inlet [mol/s]
prob.set_val('n_in_A', 0.05, units='mol/s')
prob.set_val('n_in_C', 0.05, units='mol/s')

# ----------------------------------------------------------------
# Neighbor solid temperatures (isothermal start: no conduction gradient)
# ----------------------------------------------------------------
prob.set_val('T_PEN_left',  T_op, units='K')
prob.set_val('T_PEN_right', T_op, units='K')
prob.set_val('T_IC_left',   T_op, units='K')
prob.set_val('T_IC_right',  T_op, units='K')

# ----------------------------------------------------------------
# Geometry (per-segment, active zone)
# ----------------------------------------------------------------
# Convection — anode
prob.set_val('Convection.PEN_A.Nu',       Nu)
prob.set_val('Convection.PEN_A.d_hyd',    d_hyd,         units='m')
prob.set_val('Convection.PEN_A.A',        A_pen_act,     units='m**2')
prob.set_val('Convection.PEN_A.A_ic_side',A_ic_side_act, units='m**2')

prob.set_val('Convection.IC_A.Nu',        Nu)
prob.set_val('Convection.IC_A.d_hyd',     d_hyd,         units='m')
prob.set_val('Convection.IC_A.A',         A_pen_act,     units='m**2')
prob.set_val('Convection.IC_A.A_ic_side', A_ic_side_act, units='m**2')

# Convection — cathode
prob.set_val('Convection.PEN_C.Nu',       Nu)
prob.set_val('Convection.PEN_C.d_hyd',    d_hyd,         units='m')
prob.set_val('Convection.PEN_C.A',        A_pen_act,     units='m**2')
prob.set_val('Convection.PEN_C.A_ic_side',A_ic_side_act, units='m**2')

prob.set_val('Convection.IC_C.Nu',        Nu)
prob.set_val('Convection.IC_C.d_hyd',     d_hyd,         units='m')
prob.set_val('Convection.IC_C.A',         A_pen_act,     units='m**2')
prob.set_val('Convection.IC_C.A_ic_side', A_ic_side_act, units='m**2')

# Conduction — PEN  (length_cell = total active length; component divides by N_seg internally)
prob.set_val('Conduction.PEN_Conduction.length_cell',      length_active,  units='m')
prob.set_val('Conduction.PEN_Conduction.height_cell',      width_total,    units='m')
prob.set_val('Conduction.PEN_Conduction.thickness_struct', delta_pen,      units='m')
prob.set_val('Conduction.PEN_Conduction.lambda_struct',    lambda_PEN,     units='W/m/K')

# Conduction — IC
prob.set_val('Conduction.IC_Conduction.length_cell',      length_active,  units='m')
prob.set_val('Conduction.IC_Conduction.height_cell',      width_total,    units='m')
prob.set_val('Conduction.IC_Conduction.thickness_struct', thickness_IC,   units='m')
prob.set_val('Conduction.IC_Conduction.lambda_struct',    lambda_IC,      units='W/m/K')
prob.set_val('Conduction.IC_Conduction.width_IC',         w_ch,           units='m')
prob.set_val('Conduction.IC_Conduction.thickness_ch',     h_ch,           units='m')
prob.set_val('Conduction.IC_Conduction.N_IC_walls',       2*N_channels)   # 2 side walls per channel

# ----------------------------------------------------------------
# Run
# ----------------------------------------------------------------
prob.run_model()

print('\n--- Results ---')
print(f'T_A_out : {prob.get_val("T_A_out", units="K")[0]:.2f} K')
print(f'T_C_out : {prob.get_val("T_C_out", units="K")[0]:.2f} K')
print(f'T_PEN   : {prob.get_val("T_PEN",   units="K")[0]:.2f} K')
print(f'T_IC    : {prob.get_val("T_IC",    units="K")[0]:.2f} K')
print(f'W_out_A : {prob.get_val("W_out_A", units="kg/s")[0]:.6f} kg/s')
print(f'W_out_C : {prob.get_val("W_out_C", units="kg/s")[0]:.6f} kg/s')
