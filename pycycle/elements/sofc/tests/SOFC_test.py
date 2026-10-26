"""
Run with:
    python pycycle/elements/SOFC_test.py
"""

import numpy as np

# Import pycycle.api FIRST: pycycle.elements.sofc.sofc_segment and .api
# import each other at module load time (a pre-existing circular import, see
# sofc_changes_log.md item 21's aside) -- doing a normal top-level pyCycle
# import first, exactly like every other pyCycle test/example does, resolves
# the module graph before either side is imported directly and sidesteps it.
import pycycle.api as pyc
import openmdao.api as om
from pycycle.constants import CEA_AIR_COMPOSITION
from pycycle.elements.sofc import SOFC

# ================================================================
# Geometry (same physical values as SegmentSOFC_test.py, re-derived
# per-segment for this test's own, smaller, segment counts)
# ================================================================
delta_an  = 35e-6    # m  anode layer thickness
delta_el  = 90e-6    # m  electrolyte thickness
delta_cat = 35e-6    # m  cathode layer thickness
delta_pen = delta_an + delta_el + delta_cat   # 160 um

thickness_IC = 850e-6   # m  interconnect thickness

length_active        = 0.09     # m  total active-zone length
length_passive_total = 0.072    # m  total passive-zone length (both sides combined)
length_passive_side  = length_passive_total / 2   # m  passive length on ONE side
width_total          = 0.142    # m  (conduction cross-section height)

N_passive = 5     # segments per PASSIVE side (SOFC.options['N_passive'])
N_active  = 20     # segments in the ACTIVE zone (SOFC.options['N_active'])

N_channels = 24
w_ch  = 0.0042    # m  channel width
h_ch  = 1e-3      # m  channel height
d_hyd = 2 * w_ch * h_ch / (w_ch + h_ch)
Nu    = 4.9

seg_len_act   = length_active / N_active
A_pen_act     = N_channels * w_ch * seg_len_act
A_ic_side_act = N_channels * 2 * h_ch * seg_len_act

seg_len_pas   = length_passive_side / N_passive
A_pen_pas     = N_channels * w_ch * seg_len_pas
A_ic_side_pas = N_channels * 2 * h_ch * seg_len_pas

lambda_PEN = 2.0   # W/(m*K)
lambda_IC  = 24.0   # W/(m*K)


class SOFCCycle(pyc.Cycle):
    """Minimal Cycle: FlowStart (anode, cathode) -> SOFC."""

    def setup(self):
        self.add_subsystem('fs_an',  pyc.FlowStart(composition={'H': 0.95, 'O': 0.05}))
        self.add_subsystem('fs_cat', pyc.FlowStart(composition=CEA_AIR_COMPOSITION))   # air

        self.add_subsystem('sofc', SOFC(N_passive=N_passive, N_active=N_active))

        self.pyc_connect_flow('fs_an.Fl_O',  'sofc.Fl_I_an')
        self.pyc_connect_flow('fs_cat.Fl_O', 'sofc.Fl_I_cat')

        super().setup()


prob = om.Problem()
prob.model = SOFCCycle()
prob.setup(force_alloc_complex=True)

sofc = prob.model.sofc

# ----------------------------------------------------------------
# Inlet flow conditions
# ----------------------------------------------------------------
T_op = 1073.0   # K  (800 degC nominal SOFC operating temperature)

# Operating point: ~0.31 A/cm2 on 20 x 6.4 cm2 active area, flows sized for
# the target utilizations at the inlet compositions above (molar H2/O2
# content per kg of inlet gas taken from the model's stack inlet flows).
F        = 96485.33212   # C/mol
I_stack  = 40.0          # A
U_fuel   = 0.80          # target H2 utilization
U_ox     = 0.25          # target O2 utilization
H2_per_kg_an  = 241.8    # mol H2 per kg anode gas  ({'H': 0.95, 'O': 0.05})
O2_per_kg_cat = CEA_AIR_COMPOSITION['O'] / 2 * 1000   # mol O2 per kg air (O atoms per gram -> O2 per kg, ~7.24)
W_an  = I_stack / (2 * F) / U_fuel / H2_per_kg_an
W_cat = I_stack / (4 * F) / U_ox   / O2_per_kg_cat

prob.set_val('fs_an.T',  T_op, units='K')
prob.set_val('fs_an.P',  1.01325e5, units='Pa')
prob.set_val('fs_an.W',  W_an, units='kg/s')
prob.set_val('fs_an.MN', 0.05)

prob.set_val('fs_cat.T',  T_op, units='K')
prob.set_val('fs_cat.P',  1.01325e5, units='Pa')
prob.set_val('fs_cat.W',  W_cat, units='kg/s')
prob.set_val('fs_cat.MN', 0.05)

# ----------------------------------------------------------------
# Stack-level electrical / cell-count inputs
# ----------------------------------------------------------------
prob.set_val('sofc.n_cell', 1.0)
prob.set_val('sofc.I', I_stack, units='A')

# ----------------------------------------------------------------
# Per-segment geometry
# ----------------------------------------------------------------
N_total = N_passive + N_active + N_passive
for i, name in enumerate(sofc.segments):
    is_active = (N_passive - 1) < i <= (N_passive + N_active - 1)
    seg_len   = seg_len_act if is_active else seg_len_pas
    A_pen     = A_pen_act if is_active else A_pen_pas
    A_ic_side = A_ic_side_act if is_active else A_ic_side_pas

    path = f'sofc.{name}'

    # Convection (present on every segment)
    for electrode, A in (('PEN_A', A_pen), ('PEN_C', A_pen)):
        prob.set_val(f'{path}.Convection.{electrode}.Nu', Nu)
        prob.set_val(f'{path}.Convection.{electrode}.d_hyd', d_hyd, units='m')
        prob.set_val(f'{path}.Convection.{electrode}.A', A, units='m**2')
        prob.set_val(f'{path}.Convection.{electrode}.A_ic_side', A_ic_side, units='m**2')
    for electrode, A in (('IC_A', A_pen), ('IC_C', A_pen)):
        prob.set_val(f'{path}.Convection.{electrode}.Nu', Nu)
        prob.set_val(f'{path}.Convection.{electrode}.d_hyd', d_hyd, units='m')
        prob.set_val(f'{path}.Convection.{electrode}.A', A, units='m**2')
        prob.set_val(f'{path}.Convection.{electrode}.A_ic_side', A_ic_side, units='m**2')


# ----------------------------------------------------------------
# Segment size and conduction geometry (shared by every segment's
# Conduction group; dx = length_seg_active / length_seg_passive)
# ----------------------------------------------------------------
prob.set_val('sofc.length_seg_active',  seg_len_act, units='m')
prob.set_val('sofc.length_seg_passive', seg_len_pas, units='m')
prob.set_val('sofc.A_seg_active',       width_total * seg_len_act, units='m**2')   # active area per segment
prob.set_val('sofc.height_cell',    width_total,          units='m')
prob.set_val('sofc.thickness_PEN',  delta_pen,            units='m')
prob.set_val('sofc.thickness_IC',   thickness_IC,         units='m')
prob.set_val('sofc.lambda_PEN',     lambda_PEN,           units='W/m/K')
prob.set_val('sofc.lambda_IC',      lambda_IC,            units='W/m/K')
prob.set_val('sofc.width_IC',       w_ch,                 units='m')
prob.set_val('sofc.height_ch',      h_ch,                 units='m')
prob.set_val('sofc.N_IC_walls',     2 * N_channels)

print('\n--- Running SOFC ---')
prob.run_model()

print('\n--- Results ---')
for name in sofc.segments:
    T_PEN = prob.get_val(f'sofc.{name}.T_PEN', units='K')[0]
    T_IC  = prob.get_val(f'sofc.{name}.T_IC', units='K')[0]
    print(f'{name:>10s}: T_PEN = {T_PEN:8.2f} K   T_IC = {T_IC:8.2f} K')

for name in sofc.segments_active:
    V_cell = prob.get_val(f'sofc.{name}.V_cell', units='V')[0]
    I_seg  = prob.get_val(f'sofc.balances.{name}_I', units='A')[0]
    P_seg  = prob.get_val(f'sofc.{name}.P_seg', units='W')[0]
    H2_util = prob.get_val(f'sofc.{name}.H2_utilization')[0]
    O2_util = prob.get_val(f'sofc.{name}.O2_utilization')[0]
    print(f'{name:>10s}: V_cell = {V_cell:.4f} V   I = {I_seg:.4f} A   P = {P_seg:.3f} W   '
          f'H2_util = {100*H2_util:6.2f} %   O2_util = {100*O2_util:6.2f} %')

print(f"\nU_cell (shared terminal voltage) = {prob.get_val('sofc.U_cell', units='V')[0]:.4f} V")
print(f"P_el (stack power)               = {prob.get_val('sofc.P_el', units='W')[0]:.3f} W")
print(f"H2 utilization (stack)           = {100*prob.get_val('sofc.H2_utilization_stack')[0]:.2f} %")
print(f"O2 utilization (stack)           = {100*prob.get_val('sofc.O2_utilization_stack')[0]:.2f} %")
print(f"Fl_O_an:tot:T  = {prob.get_val('sofc.Fl_O_an:tot:T', units='K')[0]:.2f} K")
print(f"Fl_O_cat:tot:T = {prob.get_val('sofc.Fl_O_cat:tot:T', units='K')[0]:.2f} K")
