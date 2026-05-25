# pyCycle Thermodynamic Property Chain

Reference for how `Thermo(om.Group)` works internally and how to access property data.

---

## Entry point: `Thermo(om.Group)` — `pycycle/thermo/thermo.py`

`Thermo` is a thin orchestrator. Its `setup()` does three things depending on `mode`:

### Step 1 — Create `base_thermo`

For `method='CEA'`:
```python
base_thermo = cea_thermo.SetTotalTP(spec=..., composition=...)
```

`SetTotalTP` is itself a Group containing two subsystems wired with `promotes=['*']`:

```
SetTotalTP
├── ChemEq      (ImplicitComponent)   — pycycle/thermo/cea/chem_eq.py
└── ThermoCalcs (Group)               — pycycle/thermo/cea/chem_eq.py
```

---

## Step 2 — Inside `ChemEq` (the core solver)

An `ImplicitComponent` with its own Newton solver.

**Inputs:**
- `composition` — b0 vector [mol/g_mix], the elemental makeup
- `T` [degK], `P` [bar]

**State variables solved:**
- `n` — mole amounts of each **species** (H2, H2O, O2, OH, N2, ...), shape `(num_species,)`
- `n_moles` — total moles = sum(n), shape `(1,)`
- `pi` — Lagrange multipliers (internal only)

**What it solves:** Gibbs free energy minimization subject to element conservation:

```
residual_n[j]    = μ_j - Σ(π_i × a_ij) = 0    for each species j
residual_π[i]    = Σ(a_ij × n_j) - b0_i = 0    for each element i
residual_n_moles = Σ(n_j) - n_moles = 0
```

Where the chemical potential:
```
μ_j = H0_j(T)/RT - S0_j(T)/R + ln(n_j) + ln(P) - ln(n_moles)
```
is evaluated from NASA polynomial coefficients for each species.

The output `n` is the **equilibrium species composition** — the molecule distribution
that minimises Gibbs energy at the given T, P and element ratios.

---

## Step 3 — Inside `ThermoCalcs`

Takes `n`, `n_moles`, `T`, `P` (promoted from `ChemEq`) and computes bulk mixture properties.

```
ThermoCalcs
├── PropsRHS              — builds LHS/RHS matrices for property derivatives
├── LinearSystemComp × 2  — solves for dT and dP sensitivities
└── PropsCalcs            — computes h, S, gamma, Cp, Cv, rho, R
                             from n, n_moles, T, P via NASA polynomials
```

**Outputs (all SI):**

| Variable | Units      | Description                        |
|----------|------------|------------------------------------|
| `h`      | cal/g      | Specific enthalpy                  |
| `S`      | cal/(g·K)  | Specific entropy                   |
| `gamma`  | —          | Ratio of specific heats            |
| `Cp`     | cal/(g·K)  | Specific heat at constant pressure |
| `Cv`     | cal/(g·K)  | Specific heat at constant volume   |
| `rho`    | g/cm³      | Density                            |
| `R`      | cal/(g·K)  | Specific gas constant              |

---

## Step 4 — Mode-specific balance (non-TP modes)

| Mode         | What is specified | What is solved | Mechanism |
|--------------|-------------------|----------------|-----------|
| `total_TP`   | T, P              | — (direct)     | none      |
| `total_hP`   | h, P              | T              | `BalanceComp`: drives T until `base_thermo.h == h_input` |
| `total_SP`   | S, P              | T              | `BalanceComp`: drives T until `base_thermo.S == S_input` |
| `static_MN`  | S, MN, W          | T, Ps          | `BalanceComp` + `PsResid` |
| `static_A`   | S, area, W        | T, Ps          | `BalanceComp` + `PsResid` |

For `total_hP`:
```python
bal.add_balance('T', eq_units='cal/g')
self.connect('base_thermo.h', 'balance.lhs:T')
self.promotes('balance', inputs=[('rhs:T', 'h')])
# → T is adjusted until base_thermo.h matches the specified h input
```

---

## Step 5 — `EngUnitProps` (the `flow` subsystem)

Passthrough component that converts SI → English units and promotes outputs as `fl_name:tot:*`.

```
h      [cal/g]        →  fl_name:tot:h      [Btu/lbm]
T      [degK]         →  fl_name:tot:T      [degR]
P      [bar]          →  fl_name:tot:P      [lbf/inch²]
S      [cal/(g·K)]    →  fl_name:tot:S      [Btu/(lbm·degR)]
Cp     [cal/(g·K)]    →  fl_name:tot:Cp     [Btu/(lbm·degR)]
gamma  [—]            →  fl_name:tot:gamma
rho    [g/cm³]        →  fl_name:tot:rho    [lbm/ft³]
R      [cal/(g·K)]    →  fl_name:tot:R      [Btu/(lbm·degR)]
composition [mol/g]   →  fl_name:tot:composition
```

> **Note:** `setup_io()` on `EngUnitProps` is called inside `configure()` (not `setup()`),
> because the composition shape is not known until after `ChemEq` has been set up.
> `configure()` runs after the full system tree is assembled.

---

## Full chain diagram

```
Inputs: composition (b0) [mol/g_mix],  T [degK],  P [bar]
              │
              ▼
┌──── SetTotalTP ──────────────────────────────────────────┐
│  ┌── ChemEq (Newton solver, ImplicitComponent) ────────┐ │
│  │  composition + T + P                                │ │
│  │  → minimize Gibbs free energy                       │ │
│  │  outputs: n  [shape: num_species]  (mole amounts)   │ │
│  │           n_moles  [scalar]                         │ │
│  └─────────────────────────────────────────────────────┘ │
│  ┌── ThermoCalcs ──────────────────────────────────────┐ │
│  │  n + n_moles + T + P (promoted from ChemEq)         │ │
│  │  → h, S, gamma, Cp, Cv, rho, R  (all SI)            │ │
│  └─────────────────────────────────────────────────────┘ │
└──────────────────────────────────────────────────────────┘
              │
              │  (non-TP modes: BalanceComp adjusts T here)
              ▼
┌──── EngUnitProps  ───────────────────────────────────────┐
│  SI → English unit conversion (passthrough, J=I)         │
│  promotes outputs as  fl_name:tot:*                      │
└──────────────────────────────────────────────────────────┘
              │
              ▼
  fl_name:tot:T, :P, :h, :S, :Cp, :gamma, :rho, :composition, ...
```

---

## Accessing flow station variables

After `p.run_model()`:

```python
T   = p.get_val('element.Fl_O:tot:T',    units='degK')   # or degR
P   = p.get_val('element.Fl_O:tot:P',    units='bar')    # or lbf/inch**2
h   = p.get_val('element.Fl_O:tot:h',    units='Btu/lbm')
Cp  = p.get_val('element.Fl_O:tot:Cp',   units='Btu/(lbm*degR)')
b0  = p.get_val('element.Fl_O:tot:composition')          # mol/g_mix
W   = p.get_val('element.Fl_O:stat:W',   units='lbm/s')
```

---

## Accessing mole fractions / activities (NOT on flow station)

`n` and `n_moles` are **not promoted** to the flow station. Access via subsystem path:

```python
n       = p.get_val('element.thermo_comp.base_thermo.n')        # shape: (num_species,)
n_moles = p.get_val('element.thermo_comp.base_thermo.n_moles')  # shape: (1,)
x_i     = n / n_moles                                           # mole fractions

# Species name list (same order as n):
from pycycle.thermo.cea.species_data import Properties, janaf
thermo_obj = Properties(janaf, init_elements={'H': ..., 'O': ...})
species = thermo_obj.products   # e.g. ['H2', 'H2O', 'O2', 'OH', ...]

x_H2  = n[species.index('H2')]  / n_moles
x_H2O = n[species.index('H2O')] / n_moles
x_O2  = n[species.index('O2')]  / n_moles
```

**To connect `n` into a downstream component** (e.g. `NernstPotential`) inside a Group:

```python
self.connect('thermo_an.base_thermo.n',       'nernst.n_species')
self.connect('thermo_an.base_thermo.n_moles', 'nernst.n_moles')
```

---

## Unit summary

| Quantity | Inside `base_thermo` | At `fl_name:tot:*` |
|----------|---------------------|---------------------|
| T        | degK                | degR                |
| P        | bar                 | lbf/inch²           |
| h        | cal/g               | Btu/lbm             |
| S        | cal/(g·K)           | Btu/(lbm·degR)      |
| Cp, Cv   | cal/(g·K)           | Btu/(lbm·degR)      |
| rho      | g/cm³               | lbm/ft³             |
| W        | kg/s (static input) | lbm/s               |

> When connecting SOFC component outputs (typically SI) to `Thermo` inputs, connect to
> the internal names (`thermo_an.T`, `thermo_an.P`) rather than the flow station names —
> these accept degK and bar directly, avoiding double conversion.