# SOFC Stack Model — Architecture & Implementation Summary

This document summarises the key modelling concepts for a Solid Oxide Fuel Cell (SOFC) stack,
intended as a reference for implementing a simulation in OpenMDAO / PyCycle.

---

## 1. Physical Stack Architecture

### Layer Stack (cross-sectional / through-thickness direction)

Each cell in the stack consists of the following layers, repeated N times in series:

```
┌─────────────────────────────────────────┐
│     Cathode Interconnect (bipolar plate) │  distributes air, collects electrons
├─────────────────────────────────────────┤
│  [AIR CHANNEL]   →  O₂-rich air flows   │
├─────────────────────────────────────────┤
│     Cathode (air electrode, porous)      │  O₂ + 4e⁻ → 2O²⁻
├─────────────────────────────────────────┤
│     Electrolyte (dense YSZ ceramic)      │  conducts O²⁻ only; blocks e⁻ and gas
├─────────────────────────────────────────┤
│     Anode (fuel electrode, porous)       │  H₂ + O²⁻ → H₂O + 2e⁻
├─────────────────────────────────────────┤
│  [FUEL CHANNEL]  →  H₂/H₂O flows        │
├─────────────────────────────────────────┤
│     Anode Interconnect (bipolar plate)   │  distributes fuel, collects electrons
└─────────────────────────────────────────┘
```

- **Operating temperature**: 600–1000 °C (high-temperature solid-state ion conductor)
- **Stack voltage**: `V_stack = N_cells × V_cell`  (~0.7 V per cell at operating point)
- **Current capacity**: set by the active electrode area `A` [m²]
- **Gas distribution**: air and fuel are fed in **parallel** to all cells via manifolds;
  cells are connected **electrically in series**

### PEN Structure
The Cathode + Electrolyte + Anode trilayer is referred to as the **PEN** (Positive-Electrolyte-Negative)
assembly. The PEN temperature `T_pen` is a key state variable shared between the energy balances
of the gas channels and the solid.

---

## 2. Segment Model Along the Flow Direction

The model divides the cell along the **gas flow axis** (x-direction) into three zones:

```
Gas inlet                                                     Gas outlet
   │                                                               │
   ▼                                                               ▼
┌──────────────┬──────────────────────────────────┬──────────────┐
│ Passive Left │           Active Zone             │Passive Right │
│   (I = 0)   │         (I = I_vec)               │   (I = 0)   │
└──────────────┴──────────────────────────────────┴──────────────┘
   ←── L_pass ──→←──────────── L_active ──────────→←── L_pass ──→
```

**Important**: All three zones contain the **complete physical layer stack**
(interconnect + cathode gas channel + cathode + electrolyte + anode + anode gas channel + interconnect).
The zones are NOT different structures — they are the same cell extended along the flow direction.
What differs between zones is **which equation terms are active**.

---

## 3. Passive Left Segment

### Physical meaning
The passive left zone represents the inlet region of the cell. It contains the full PEN + IC structure,
but the electrochemical source terms are zeroed out (`I = 0`). Physically this corresponds to:
- The seal/manifold border region where electrode coating may be absent, or
- A numerical entrance zone where gas is assumed to not yet react (boundary condition buffer)

### Boundary conditions (passive left specifically)
- `T_cell_left` ← **hardwired constant** (e.g. 973 K) — thermal boundary at the very left edge
- `Qdot_conduct_right` ← received from active zone (heat conducted in from the right)
- Gas inlet: fresh air (`ndot_air_in`) and fresh fuel (`ndot_fuel_in`) from system boundary

### Sub-models active in passive left

| Sub-model | Active? | Notes |
|---|---|---|
| Energy balance PEN | ✅ Yes | convection from both gas sides + axial conduction |
| Energy balance IC | ✅ Yes | same structure |
| Energy balance cathode gas | ✅ Yes | enthalpy transport, convection with PEN |
| Energy balance anode gas | ✅ Yes | enthalpy transport, convection with PEN |
| Mass balance cathode | ✅ Yes (I=0) | gas flows through, **no O₂ consumed** |
| Mass balance anode | ✅ Yes (I=0) | gas flows through, **no H₂ consumed** |
| Electrochemistry | ❌ No | I = 0 hardwired → no Nernst, no overpotentials |

### Outputs to next zone
- `T_pen_out` → T_cell_left of active zone
- `Qdot_conduct_right` → thermal coupling to active zone
- `ndot_air_out`, `ndot_fuel_out` → gas composition unchanged from inlet (no reaction)

---

## 4. Active Segment

### Physical meaning
The active zone is where electrochemistry occurs. Current `I = I_vec` is non-zero, O₂ is consumed
at the cathode, H₂ is consumed at the anode, H₂O is produced, and electrical power is extracted.

### Boundary conditions (active zone)
- `T_cell_left` ← `T_pen_out` from passive left
- `Qdot_conduct_left` ← from passive left
- `Qdot_conduct_right` ← from passive right (or 0 if passive right is symmetric boundary)
- Gas inlet: pre-heated gas exiting passive left

### Sub-models active in active zone

| Sub-model | Active? | Notes |
|---|---|---|
| Energy balance PEN | ✅ Yes | + electrochemical heat source + Joule heating |
| Energy balance IC | ✅ Yes | + Joule/conduction contributions |
| Energy balance cathode gas | ✅ Yes | O₂ content decreases along x |
| Energy balance anode gas | ✅ Yes | H₂ decreases, H₂O increases along x |
| Mass balance cathode | ✅ Yes (I=I_vec) | `Δn_O2 = -I / (4F)` per segment |
| Mass balance anode | ✅ Yes (I=I_vec) | `Δn_H2 = -I / (2F)`, `Δn_H2O = +I / (2F)` |
| Electrochemistry / elDATA | ✅ Yes | Nernst potential, overpotentials, V_cell, j(x) |

### Key electrochemical equations

```
V_tn    = −ΔH°(T) / (2F)                                         [thermoneutral voltage]
E_OCV   = −ΔG°(T) / (2F)                                         [standard Nernst voltage]
E_Nernst= E_OCV + (R·T/2F)·ln( x_H2·√(x_O2·P/P_ref) / x_H2O ) [concentration correction]

V_cell  = E_Nernst − η_ASR                                        [operating voltage]
W_elec  = V_cell · I                                              [electrical power]
```

Where:
- `ΔH°`, `ΔG°` are computed from NASA polynomials via `Properties.H0(T)` and `S0(T)` — no polynomial fit
- `x_H2`, `x_H2O` = anode mole fractions;  `x_O2` = **cathode** mole fraction (cross-electrode Nernst)
- `η_ASR` = `ASR(T) · i` — area-specific resistance overpotential (Arrhenius model)
- `i` = current density [A/m²], `F` = Faraday constant

### Composition change along x (species balance per infinitesimal segment dx)

```
dn_H2  / dx = -j(x) / (2F · v_fuel)    [H₂ consumed]
dn_H2O / dx = +j(x) / (2F · v_fuel)    [H₂O produced]
dn_O2  / dx = -j(x) / (4F · v_air)     [O₂ consumed]
```

Where `j(x)` [A/m²] is the local current density and `v` is the molar flow velocity.

### Outputs to next zone
- `T_pen_out` → T_cell_left of passive right
- `Qdot_conduct_right` → thermal coupling to passive right
- `ndot_air_out` → depleted air (lower O₂ fraction)
- `ndot_fuel_out` → depleted fuel (lower H₂, higher H₂O fraction)

---

## 5. Passive Right Segment

### Physical meaning
Mirrors passive left structurally. It handles the thermal mass and conduction at the outlet end,
ensuring the energy balance of the whole cell is closed. The hot depleted gas still exchanges heat
with the solid before exiting. No electrochemistry occurs.

### Boundary conditions (passive right specifically)
- `T_cell_left` ← `T_pen_out` from active zone (solid temperature carried in from active)
- `Qdot_conduct_right` ← **hardwired constant 0** — open right boundary, nothing further to the right
- Gas inlet: depleted air and fuel from active zone outlet

### Sub-models active in passive right
**Identical internal structure to passive left** — same four energy balances, same mass balances
with `I = 0`. Only the boundary conditions differ.

### Outputs
- `T_pen_out` → thermal state at stack outlet
- `ndot_air_out`, `ndot_fuel_out` → final exit composition (unchanged from active zone outlet)

---

## 6. Thermal Coupling Between Zones

Axial heat conduction through the solid (PEN + IC) links all three zones:

```
[Passive Left] ←→ Qdot_conduct ←→ [Active] ←→ Qdot_conduct ←→ [Passive Right]
                T_cell propagates →→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→
```

- `Qdot_conduct = λ_solid · A_cross · (T_left - T_right) / Δx`
- `T_cell_left` is the solid temperature at the left face of each segment, passed as input
  from the segment to its left

---

## 7. Composition Flow Summary

### Anode (fuel) side — along flow direction x

| Location | H₂ | H₂O | Notes |
|---|---|---|---|
| Passive left inlet | `y_H2,in` | `y_H2O,in` | system boundary conditions |
| Passive left outlet | unchanged | unchanged | I = 0, no reaction |
| Active zone (x=0) | `y_H2,in` | `y_H2O,in` | starts reacting |
| Active zone (x=L) | `y_H2,in · (1-FU)` | increases | FU = fuel utilization |
| Passive right outlet | unchanged | unchanged | I = 0, no reaction |

### Cathode (air) side — along flow direction x

| Location | O₂ | N₂ | Notes |
|---|---|---|---|
| Passive left inlet | ~0.21 | ~0.79 | standard air |
| Active zone (x=L) | `0.21 · (1-AU)` | ~0.79 | AU = air utilization |
| Passive right outlet | unchanged | unchanged | I = 0, no reaction |

**Fuel utilization**: `FU = I_total / (2F · ndot_H2_in)` — fraction of H₂ converted  
**Air utilization**: `AU = I_total / (4F · ndot_O2_in)` — fraction of O₂ consumed

---

## 8. OpenMDAO Implementation Architecture

### Component responsibilities

```
SOFCThermoAdd              ExplicitComponent   composition_out, Wout from I and inlet b0/W
                                               (Faraday's law, element-basis update)
heat_convection_electrode  ExplicitComponent   Q_conv_ = α·A·(T_channel_out - T_struc_out)
                                               α = Nu·λ/d_hyd; options: electrode, structure
                                               T_bulk output for thermal_conductivity()
thermal_conductivity       ExplicitComponent   λ_mix(T_bulk, composition)
NernstThermo               ExplicitComponent   E_OCV, E_Nernst, V_tn, Qdot_chem
                                               uses NASA polynomials (Properties.H0/S0) — no fit
                                               x_O2 comes from cathode, x_H2/x_H2O from anode
AreaSpecificResistance     ExplicitComponent   ASR(T), eta_asr = ASR·i
VoltageCalc                ExplicitComponent   V_cell = E_Nernst - eta_asr
ChannelEnergyBalance       ImplicitComponent   solves T_channel_out
                                               residual: h_out·W_out - h_in·W_in
                                                       - Q_conv_PEN - Q_conv_IC - Q_loss = 0
PENEnergyBalance           ImplicitComponent   solves T_cell
                                               residual includes Q_conv_PEN (anode+cathode),
                                               Qdot_chem, -V_cell·I, axial conduction
ICEnergyBalance            ImplicitComponent   solves T_IC
```

### Channel energy balance — full term list

```
0 = (h_in·W_in − h_out·W_out)  ← advective enthalpy flux
  + Q_conv_PEN                  ← convection from PEN to channel
  + Q_conv_IC                   ← convection from IC to channel
  + Q_loss                      ← heat loss to environment (= 0 in passive model)
```

`T_channel_out` is the implicit state; `Thermo(mode='total_TP')` evaluates `h_out` at each Newton iteration.  
`Q_loss` and radiation are wired to zero for the passive model but declared as inputs for extensibility.

### PEN energy balance — active zone term list

```
0 = Q_conv_PEN_an              ← convection from anode channel
  + Q_conv_PEN_cat             ← convection from cathode channel
  + Qdot_chem                  ← total reaction enthalpy rate = V_tn · I
  − V_cell · I                 ← electrical power extracted
  + Qdot_conduct_left          ← axial conduction from left segment
  + Qdot_conduct_right         ← axial conduction from right segment
  + Q_loss                     ← environment losses (= 0 passive)
```

Net electrochemical heat deposited in PEN: `(V_tn − V_cell) · I`

### Component type decision guide

| Situation | Type |
|---|---|
| Compute Nernst potential, V_tn, Qdot_chem | `ExplicitComponent` |
| Compute convective heat transfer | `ExplicitComponent` |
| Compute composition update (Faraday) | `ExplicitComponent` |
| Compute overpotentials from i, T | `ExplicitComponent` |
| CEA equilibrium solve (Gibbs min.) | `ImplicitComponent` (inside Thermo) |
| Close channel energy balance → solve T_channel_out | `ImplicitComponent` |
| Close PEN/IC energy balance → solve T_cell/T_IC | `ImplicitComponent` |

### Passive vs Active component reuse

The passive and active zones use the **same component classes**. The difference is a parameter or
input flag `I` (current):
- Passive: connect `I = 0.0` (constant) into the MassBalance and Electrochemistry components
- Active: connect `I = I_operating` (design variable or solver state)

This avoids code duplication — one `SegmentGroup` can be instantiated three times with different
`I` inputs.

```python
class SegmentGroup(om.Group):
    def setup(self):
        self.add_subsystem('mass_bal_cat', MassBalanceCathode())
        self.add_subsystem('mass_bal_an',  MassBalanceAnode())
        self.add_subsystem('eb_pen',       EnergyBalancePEN())
        self.add_subsystem('eb_ic',        EnergyBalanceIC())
        self.add_subsystem('eb_cat_gas',   EnergyBalanceCathodeGas())
        self.add_subsystem('eb_an_gas',    EnergyBalanceAnodeGas())
        # Only in active zone: add electrochemistry
        # self.add_subsystem('electrochem', Electrochemistry())

# Instantiation
prob.model.add_subsystem('passive_left',  SegmentGroup())
prob.model.add_subsystem('active',        SegmentGroup())   # + Electrochemistry
prob.model.add_subsystem('passive_right', SegmentGroup())

# Wire I=0 for passive zones
prob.model.set_val('passive_left.mass_bal_cat.I',  0.0)
prob.model.set_val('passive_left.mass_bal_an.I',   0.0)
prob.model.set_val('passive_right.mass_bal_cat.I', 0.0)
prob.model.set_val('passive_right.mass_bal_an.I',  0.0)
```

### Thermal signal chain (connect calls)

```python
# T_cell_left boundary: passive left gets hardwired constant
prob.model.set_val('passive_left.eb_pen.T_cell_left', 973.0)   # [K]

# Passive left → Active
prob.model.connect('passive_left.eb_pen.T_pen_out',       'active.eb_pen.T_cell_left')
prob.model.connect('passive_left.eb_pen.Qdot_conduct_out','active.eb_pen.Qdot_conduct_left')

# Active → Passive right
prob.model.connect('active.eb_pen.T_pen_out',             'passive_right.eb_pen.T_cell_left')
prob.model.connect('active.eb_pen.Qdot_conduct_out',      'passive_right.eb_pen.Qdot_conduct_left')

# Right boundary: Qdot_conduct_right = 0 (open boundary)
prob.model.set_val('passive_right.eb_pen.Qdot_conduct_right', 0.0)
```

---

## 9. Implementation Status (as of 2026-05-29)

### File overview

| File | Status | Key contents |
|---|---|---|
| `sofc_thermo_add.py` | ✅ Complete | `SOFCThermoAdd`: element-basis composition + Faraday mass update |
| `sofc_reaction.py` | ✅ Functional | `NernstThermo`, `VoltageCalc`, `AreaSpecificResistanceOverpotential`, `CurrentDensityCalc`, `ElectroChemistry` group |
| `sofc_balances.py` | ✅ Functional | `PENEnergyBalance`, `ICEnergyBalance`, `ChannelEnergyBalance`, `HeatConductionStructure`, `HeatConduction`, `SpeciesFlowCalc` |
| `sofc_heat.py` | ✅ Functional | `HeatConvection`, `HeatConvectionChannel`, `HeatConvectionElectrode`, `ThermalConductivityMixture` (coefficients TBD) |
| `sofc_observables.py` | ✅ Functional | `SpeciesFlowCalc`, `BulkComposition`, `SpeciesUtilization`, `MassSanityCheck` (stub) |
| `sofc_segment.py` | ✅ `p.setup()` passes | `SegmentSOFC(om.Group)` — full active/passive segment |
| `sofc_api.py` | ✅ Complete | Re-exports all classes above |

### SegmentSOFC — subsystem structure

```
SegmentSOFC (om.Group)
├── anode_rxn          SOFCThermoAdd(anode)        W_in_A, composition_in_A, I → composition_out_A, W_out_A
├── cathode_rxn        SOFCThermoAdd(cathode)       W_in_C, composition_in_C, I → composition_out_C, W_out_C
├── anode_thermo_out   Thermo(total_TP)             T_A_out, p_A_out → h_A_out, Cp_A_out, S_A_out, Cv_A_out
├── cathode_thermo_out Thermo(total_TP)             T_C_out, p_C_out → h_C_out, Cp_C_out, S_C_out, Cv_C_out
├── anode_species_out  SpeciesFlowCalc(outlet props) W_out_A → x_out_A, n_out_A
├── cathode_species_out SpeciesFlowCalc(outlet props) W_out_C → x_out_C, n_out_C
├── anode_bulk         BulkComposition             x_in_A, x_out_A → x_bulk_A
├── cathode_bulk       BulkComposition             x_in_C, x_out_C → x_bulk_C
├── Convection         HeatConvection              W_in/out_A/C, h/T in/out → Q_conv_*
├── Conduction         HeatConduction(N_seg)       T_PEN/IC_left/right → Qc_left, Qc_right
├── PEN                PENEnergyBalance(seg_type)  → T_PEN  (implicit)
├── IC                 ICEnergyBalance()           → T_IC   (implicit)
├── Cathode            ChannelEnergyBalance        → T_C_out (implicit)
├── Anode              ChannelEnergyBalance        → T_A_out (implicit)
│── [active only]:
├── Electrochemistry   ElectroChemistry(N_seg)     x_H2/O2/H2O, T_PEN, I → V_cell, Qdot_chem
└── Utilization        SpeciesUtilization          n_in_A, x_in_A, I → H2_utilization, O2_utilization
```

### Key internal wiring decisions

**Composition input to outlet Thermo** — reaction output connected explicitly:
```python
self.connect('anode_rxn.composition_out', 'anode_thermo_out.composition')
```

**n / n_moles from Thermo** — not promoted through Thermo boundary, accessed via path:
```python
self.connect('anode_thermo_out.base_thermo.n',       'anode_species_out.n_i')
self.connect('anode_thermo_out.base_thermo.n_moles', 'anode_species_out.n_moles')
```

**Scalar mole fraction extraction for Electrochemistry** — src_indices on bulk vector:
```python
H2_idx = anode_out_props.products.index('H2')
self.connect('x_bulk_A', 'Electrochemistry.x_H2', src_indices=[H2_idx])
```

**Composition basis strategy** — all species calculations use outlet `Properties` basis so `x_in_A` and `x_out_A` always share the same length. `x_in_A`/`x_in_C`/`n_in_A`/`n_in_C` are external inputs; the parent SOFC Element connects the previous segment's `x_out` to the next segment's `x_in`.

**HeatConvectionChannel mass flow** — uses `W_out·h_out − W_in·h_in` (not `W·Δh`) to account for O²⁻ mass transfer changing W along the channel.

### Boundary variables (SegmentSOFC external interface)

**Inputs from parent:**
`W_in_A/C`, `composition_in_A/C`, `h_A/C_in`, `T_A/C_in`, `p_A/C_out`, `x_in_A/C`, `n_in_A/C`, `I` (active), `T_PEN/IC_left/right`, `N_cell`, geometric parameters

**Outputs to parent:**
`W_out_A/C`, `composition_out_A/C`, `h_A/C_out`, `T_A/C_out`, `p_A/C_out`, `x_out_A/C`, `n_out_A/C`, `T_PEN`, `T_IC`, `V_cell` (active), `H2_utilization`, `O2_utilization` (active)

### Still pending

- Add `NewtonSolver` + `DirectSolver` to `SegmentSOFC`
- Wire `x_bulk_A/C` → `ThermalConductivityMixture.x1/x2` in `HeatConvection`
- Fill `ThermalConductivityMixture` polynomial coefficients from literature
- `InletConverter`: first-segment `x_in_A/C` from SOFC inlet Thermo (handles composition-basis mismatch at inlet)
- Full `SOFC` Element: N-segment chain + inter-segment connections + boundary conditions
- Fix `ICEnergyBalance` residual output pattern (move to separate ExplicitComponent)

---

## 10. Extracting Composition in PyCycle

If using the PyCycle/CEA backend:

```python
# Mole fractions (CEA backend)
n       = prob.get_val('DESIGN.<element>.real_flow.base_thermo.n')       # molar amounts
n_moles = prob.get_val('DESIGN.<element>.real_flow.base_thermo.n_moles') # total moles
mole_fracs = n / n_moles

# Species names
base_thermo = prob.model.DESIGN.<element>.real_flow.base_thermo
species = list(base_thermo.options['spec'].products.keys())

for sp, mf in zip(species, mole_fracs):
    if mf > 1e-8:
        print(f"{sp:12s}: {mf:.6f}")
```

If using the TABULAR backend, `Fl_O:tot:composition` is a scalar FAR (fuel-to-air ratio),
not a species vector. Use offline CEA post-processing to recover individual mole fractions.

---

## 10. Key Parameters Reference

| Symbol | Description | Typical value |
|---|---|---|
| `T_op` | Operating temperature | 700–1000 °C |
| `P_op` | Operating pressure | 1–5 bar |
| `FU` | Fuel utilization | 0.70–0.85 |
| `AU` | Air utilization | 0.20–0.30 |
| `i` | Current density | 0.2–0.8 A/cm² |
| `V_cell` | Cell voltage | 0.6–0.8 V |
| `R_ohmic` | Ohmic area-specific resistance | 0.1–0.3 Ω·cm² |
| `T_cell_left` (passive left BC) | Fixed inlet solid temperature | 973 K |
| `Qdot_conduct_right` (passive right BC) | Open right boundary | 0 W |
| `F` | Faraday constant | 96485 C/mol |
