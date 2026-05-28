# SOFC pyCycle Architecture Notes

## 1. Key Terminology

| Term | Meaning |
|------|---------|
| **PEN** | Positive electrode / Electrolyte / Negative electrode — the three-layer cell core |
| **b0** | Element-basis composition vector: moles of each **atomic element** per kg of mixture |
| **n** | Species-basis vector: moles of each **molecule** per kg of mixture (from ChemEq) |
| **n_moles** | Total moles per kg of mixture = `sum(n)` = `1 / MW_mix * 1000` |

---

## 2. pyCycle Thermo Architecture

### Class hierarchy

```
Thermo(om.Group)                        ← method-agnostic wrapper
└── base_thermo: SetTotalTP(om.Group)   ← CEA implementation
        ├── chem_eq: ChemEq             ← ImplicitComponent, solves for n, n_moles
        └── props:   ThermoCalcs        ← Group, computes h, Cp, S, R, gamma
                └── tp2props: PropsCalcs ← ExplicitComponent, actual calculation
```

### Variable access

| Variable | Promoted path | Physical meaning |
|----------|--------------|-----------------|
| `h` | `Fl_X:tot:h` | Mixture enthalpy [cal/g] |
| `R` | `Fl_X:tot:R` | Specific gas constant [J/(kg·K)] |
| `composition` | `Fl_X:tot:composition` | b0 vector [mol/g_mix] |
| `n` | `anode_out.base_thermo.n` | Species moles/kg — **not promoted** |
| `n_moles` | `anode_out.base_thermo.n_moles` | Total mol/kg — **not promoted** |

> **Note:** `n` and `n_moles` were intentionally removed from the promoted outputs
> (see `TODO` in `thermo.py:57`) because they are CEA-specific and break the
> method-agnostic interface when using `method='TABULAR'`.

### Getting MW_mix and mole fractions

```python
# From promoted variables:
n_moles = R / 8.314          # mol/kg  (R_universal = 8.314 J/mol/K)
MW_mix  = 1000.0 / n_moles   # g/mol

# From internal n vector:
x_i = n[i] / n_moles         # mole fraction of species i
```

### Enthalpy calculation path

```
PropsCalcs.compute():
    h = sum(n * H0(T)) * R_universal * T    ← dot product of species moles with NASA polynomials

Promoted through:
    tp2props → (promotes_outputs=['h',...]) → ThermoCalcs level
    props    → (promotes=['*'])             → SetTotalTP level
    base_thermo → (out_vars includes 'h' only for mode='total_TP','total_hP',etc.) → Thermo level
```

### Why SetTotalTP is a separate class

`Thermo` is method-agnostic — it only knows the agreed interface `(T, P, composition) → (h, Cp, R, ...)`.
Each thermo method (`CEA`, `TABULAR`) owns its internal structure behind `SetTotalTP`.
Adding a new method never requires modifying `Thermo` itself.

---

## 3. Composition Vector (b0)

`b0` tracks **atomic elements**, not molecules:

```python
CEA_AIR_COMPOSITION = {'N': 5.39e-02, 'O': 1.45e-02, 'Ar': 3.23e-04, 'C': 1.10e-05}
```

For an H2/H2O/N2 anode gas the elements are `['H', 'O', 'N']` (alphabetical).

**You cannot recover molecular mole fractions from b0 alone** — ChemEq is needed to
distribute atoms into molecules. The one exception: for a H/O/N system at SOFC
temperatures with H excess, the equilibrium assignment is deterministic:

```
n_H2O = b0_O             (all O → H2O)
n_H2  = (b0_H - 2*b0_O) / 2
n_N2  = b0_N / 2
```

This is only valid because ChemEq gives the same result as Faraday's law for typical
SOFC conditions (FU < 100%, T = 700–1000°C, pure H2 fuel). See Section 5.

---

## 4. output_port_data() Pattern

Called **before** `setup()` in `pyc_setup_output_ports()` to inform the framework
of the outlet element set, so that downstream array sizes can be determined.

| Element type changes? | Method used |
|-----------------------|-------------|
| No (inlet, duct, compressor) | `copy_flow('Fl_I', 'Fl_O')` |
| Yes (combustor, SOFC anode) | `output_port_data()` on ThermoAdd/SOFCThermoAdd |

`SOFCThermoAdd.output_port_data()` simply adds `'O'` to the anode outlet dict
(O²⁻ ions arrive through the electrolyte). The cathode element set is unchanged.

---

## 5. OpenMDAO Component Design Rules

### Options vs Inputs

| Use `options.declare` | Use `add_input` |
|---|---|
| Controls `setup()` structure (array sizes, number of subsystems, which equations run) | Used only in `compute()` / `apply_nonlinear()` |
| `electrode`, `structure`, `mixture` (switches equations) | All geometry: `A`, `d_hyd`, `L`, `thickness_*` |
| `N_segments` (loop count in `setup`) | All material properties: `lambda`, `k` |
| `N_cells` when instantiating N cell subsystems | `n_cell` as scalar multiplier (e.g. `P_elec = I·V·n_cell`) |
| `thermo`, `spec` objects (set array sizes) | Anything you might want to sweep or optimise |

### ImplicitComponent rules

```python
def apply_nonlinear(self, inputs, outputs, residuals):
    # ONLY write to residuals — never to outputs
    # outputs here is the solver's current state guess (read-only)
    residuals['T'] = inputs['Q_in'] - inputs['Q_out']   # ✓
    outputs['T']   = ...                                  # ✗ corrupts solver state

def linearize(self, inputs, outputs, J):
    # outputs IS a parameter here (unlike ExplicitComponent)
    J['T', 'Q_in']  =  1.0
    J['T', 'T']     =  0.0   # must declare even when zero

def compute_partials(self, inputs, J):
    # ExplicitComponent: NO outputs parameter
    pass
```

- Always `declare_partials(state, state, val=0.0)` even when the local partial is zero
- Diagnostic outputs (residual monitors) belong in a **separate ExplicitComponent outside the Newton group** — at convergence the residual is ~0 anyway
- `ICEnergyBalance` still has the wrong pattern (`outputs['IC_residuum']` in `apply_nonlinear`) — needs fixing

### Newton group structure (per segment)

All ImplicitComponents whose states are coupled **and all ExplicitComponents in the same algebraic loop** must share one Group with a Newton solver:

```
sofc_segment (Group)
  nonlinear_solver = NewtonSolver(solve_subsystems=False)
  linear_solver    = DirectSolver()
  ├── anode_thermo_out / cathode_thermo_out   (Thermo, total_TP) — T_out→h_out
  ├── heat_conv_PEN_an/cat, heat_conv_IC_an/cat  (heat_convection_electrode)
  ├── heat_conv_ch_an/cat                     (heat_convection_channel)
  ├── PEN_balance                             → state: T_cell
  ├── IC_balance                              → state: T_IC
  ├── channel_balance_an                      → state: T_an_out
  └── channel_balance_cat                     → state: T_cat_out

[outside Newton group]
  └── ic_residual_monitor  (ExplicitComponent, diagnostic only)
```

The Thermo components **must** be inside the Newton group because they carry the
`T_out → h_out` sensitivity (`dh/dT = cp`) that makes the channel balance Jacobian non-zero.

For a full stack: outer `NonlinearBlockGS` sweeps over segment groups; inter-segment
conduction (`heat_conduction_struct`) carries heat between them.

---

## 5. SOFCThermoAdd

### Purpose

Updates `b0` and `W` after the electrochemical reaction for use by the downstream
`Thermo(mode='total_TP')` component. Mirrors the `ThermoAdd` interface for
pyCycle port-propagation compatibility.

### Reaction physics

| Side | Reaction | Faraday rate |
|------|----------|-------------|
| Anode | H2 + O²⁻ → H2O + 2e⁻ | `ndot_H2 = I / (2F)` mol/s consumed |
| Cathode | O2 + 4e⁻ → 2O²⁻ | `ndot_O2 = I / (4F)` mol/s consumed |

### b0 update logic

```
Anode:   n_out[O_idx] += ndot_H2       (O²⁻ arrives, adds O atoms)
         W_out = W_in + ndot_H2 * MW_O  (mass increase from O²⁻)

Cathode: n_out[O_idx] -= 2 * ndot_O2   (O²⁻ departs, removes O atoms)
         W_out = W_in - ndot_O2 * MW_O2 (mass decrease from O2 consumed)
```

H atoms are conserved — H2 and H2O both contain 2 H atoms, so b0_H does not change.

### Why ChemEq re-equilibration does not lose accuracy

For H2/H2O/N2 at SOFC conditions (H excess, T < 1200°C):
- ChemEq and Faraday's law give **identical** H2/H2O splits
- Both give: `n_H2O_out = n_H2O_in + Δ`, `n_H2_out = n_H2_in - Δ`
- Minor dissociation species (H, OH, O) are negligible below 1200°C
- Accuracy only degrades near 100% fuel utilization or with hydrocarbon fuels

### Known limitations

| Issue | Detail |
|-------|--------|
| No `N_cell` scaling | Current `I` assumed per-cell, no stack scaling |
| CS partials | `declare_partials('*','*', method='cs')` — correct but slow |
| No `Qdot_chem` output | Reaction enthalpy needed by `PENEnergyBalance` not computed here |
| b0 space only | Cannot directly output molecular mole fractions for multi-segment case |

---

## 6. Full Anode Data Flow (per segment)

```
Fl_I_an (W, composition, T, P)
        │
        ▼
┌─────────────────┐
│  SOFCThermoAdd  │  inputs:  Fl_I:stat:W, Fl_I:tot:composition, I
│                 │  outputs: composition_out (b0 updated)
│  Faraday's law  │           Wout
│  b0 update      │           H2_consumed, H2O_produced
└────────┬────────┘
         │ composition_out, Wout
         ▼
┌─────────────────────────────────┐
│  Thermo(mode='total_TP')        │  also inputs: T_cell (from PENEnergyBalance)
│  └── SetTotalTP                 │               P (from inlet)
│       ├── ChemEq                │
│       │    └── n, n_moles       │  ← species distribution (not promoted)
│       └── ThermoCalcs           │
│            └── h, Cp, R, S, ... │  ← promoted to Fl_O_an:tot:*
└──────────┬──────────────────────┘
           │
     ┌─────┴──────────────────┐
     │                        │
     ▼                        ▼
Fl_O_an:tot:*          base_thermo.n
(h, Cp, R, ...)        base_thermo.n_moles
(to energy balance)          │
                             ▼
                   ┌──────────────────────┐
                   │  ChannelMassBalance  │  also input: Wout
                   │                      │
                   │  x_i = n[i]/n_moles  │
                   │  ndot = n_moles * W  │
                   └──────────┬───────────┘
                              │
                              ▼
                   x_H2, x_H2O, x_N2, ndot_out
                   (→ ConvectionHTC, AnodeEnergyBalance)
```

---

## 7. ChannelMassBalance

Takes `n` and `n_moles` from `base_thermo` after ChemEq runs.
Species indices are looked up once in `setup()` from the `products` list.

```python
products = anode_out.base_thermo.thermo.products  # e.g. ['H','H2','H2O','N2',...]
H2_idx   = products.index('H2')
H2O_idx  = products.index('H2O')
N2_idx   = products.index('N2')

x_H2  = n[H2_idx]  / n_moles
x_H2O = n[H2O_idx] / n_moles
x_N2  = n[N2_idx]  / n_moles
ndot  = n_moles * W_kg_s
```

**Connections in segment setup:**
```python
self.connect('anode_out.base_thermo.n',       'anode_mass.n')
self.connect('anode_out.base_thermo.n_moles', 'anode_mass.n_moles')
self.connect('anode_rxn.Wout',                'anode_mass.W')
```

> **Why not compute x_i inside SOFCThermoAdd?**
> For segment 1 (pure H2 inlet) it is possible analytically. For segments 2…N
> the inlet already contains H2O, so b0_H contains H atoms from both H2 and H2O —
> you cannot split them without ChemEq. `ChannelMassBalance` after `Thermo` is
> the general solution.

---

## 8. NernstThermo — Electrochemical Potentials from CEA

`NernstThermo` (`sofc_reaction.py`) replaces the old `NernstPotential` polynomial fit.
It uses `Properties.H0(T)` and `Properties.S0(T)` directly from janaf — no fitting parameters.

### What it computes

```
H0(T), S0(T)   ← NASA polynomial arrays [dimensionless: H/(RT), S/R] for all H/O species

ΔH° = (H0_H2O - H0_H2 - 0.5·H0_O2) · R · T          [J/mol]
ΔS° = (S0_H2O - S0_H2 - 0.5·S0_O2) · R               [J/mol/K]
ΔG° = ΔH° - T·ΔS°                                     [J/mol]

V_tn     = -ΔH° / (2F)                                 [V]  thermoneutral voltage
E_OCV    = -ΔG° / (2F)                                 [V]  standard Nernst
E_Nernst = E_OCV + (R·T/2F)·ln(x_H2·√(x_O2·P/P_ref) / x_H2O)  [V]
Qdot_chem = V_tn · I                                   [W]  total reaction enthalpy rate
```

### Cross-electrode note

`x_O2` is **cathode** mole fraction — O₂ never appears in the anode gas after CEA equilibrium.  
`x_H2`, `x_H2O` come from the anode stream via `ChannelMassBalance`.

### PEN heat from electrochemistry

```
Qdot_chem − V_cell · I  =  (V_tn − V_cell) · I        [W → PEN]
```

`V_tn · I` is total enthalpy released; `V_cell · I` leaves as electrical work.
The remainder heats the PEN. At `V_cell = V_tn` no heat is generated.

### Partials

`declare_partials('*', '*', method='cs')` — NASA polynomial evaluation is CS-safe;
`H0_applyJ` / `S0_applyJ` exist for analytic dH/dT if needed later.

### Old NernstPotential is superseded

`NernstPotential` (quadratic ΔG fit via `GR_ecr_a/b/c`) is fully replaced.
`VoltageCalc` and `AreaSpecificResistanceOverpotential` remain valid as separate components.
