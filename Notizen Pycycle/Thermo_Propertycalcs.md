Entry point: Thermo(om.Group)

  Thermo is a thin orchestrator. Its setup() does three things depending on mode:

  Step 1 — Create base_thermo

  For method='CEA', it instantiates:
  base_thermo = cea_thermo.SetTotalTP(spec=..., composition=...)

  SetTotalTP is itself a Group containing two subsystems wired together with promotes=['*']:

  SetTotalTP
  ├── ChemEq      (ImplicitComponent)
  └── ThermoCalcs (Group)

  ---
  Step 2 — Inside ChemEq (the core solver)

  This is where the actual physics happens. It is an ImplicitComponent — it has a Newton solver running inside it.

  Inputs:   
  - composition — the b0 vector [mol/g_mix], i.e. the elemental makeup
  - T [degK], P [bar] 
  
  State variables being solved:
  - n — mole amounts of each individual species (H2, H2O, O2, OH, N2, ...), shape (num_species,)
  - n_moles — total moles = sum(n)
  - pi — Lagrange multipliers (internal, not needed outside)

  What it solves: Gibbs free energy minimization subject to element conservation:

  residual_n[j] = μ_j - Σ(π_i × a_ij) = 0    for each species j
  residual_π[i] = Σ(a_ij × n_j) - b0_i = 0   for each element i
  residual_n_moles = Σ(n_j) - n_moles = 0

  Where μ_j = H0_j/RT - S0_j/R + ln(n_j) + ln(P) - ln(n_moles) is the chemical potential from NASA polynomials.

  The output n is the equilibrium species composition — the distribution of molecules that minimises Gibbs energy at the given T, P and element ratios.
  
  ---
  Step 3 — Inside ThermoCalcs
  
  Takes n, n_moles, T, P (all promoted from ChemEq) and computes bulk mixture properties:

  ThermoCalcs
  ├── PropsRHS     — builds the LHS/RHS matrices for property derivatives
  ├── LinearSystemComp × 2 — solves for dT and dP sensitivities
  └── PropsCalcs   — computes h, S, gamma, Cp, Cv, rho, R
                     from n, n_moles, T, P via NASA polynomials
  
  PropsCalcs outputs (all SI units):
  - h [cal/g], S [cal/(g·K)], gamma [—], Cp [cal/(g·K)], Cv, rho, R
 
  These get promoted up through SetTotalTP and then through Thermo.

  ---
  Step 4 — Mode-specific balance (for non-TP modes)
  
  mode='total_TP' → no extra solver needed, T and P are direct inputs.

  mode='total_hP' → you specify h as input and need T as output. A BalanceComp is added:

  bal.add_balance('T', eq_units='cal/g')
  # drives T until: base_thermo.h == input_h
  self.connect('base_thermo.h', 'balance.lhs:T')
  self.promotes('balance', inputs=[('rhs:T', 'h')])

  mode='static_MN' → additionally solves for static pressure Ps given Mach number.

  ---
  Step 5 — EngUnitProps (the flow subsystem)
  
  After all properties are computed in SI, EngUnitProps passthrough-converts them to English units and promotes them with the fl_name prefix:

  h [cal/g]         →  Fl_O:tot:h     [Btu/lbm]
  T [degK]          →  Fl_O:tot:T     [degR]
  P [bar]           →  Fl_O:tot:P     [lbf/inch²]
  S [cal/(g·K)]     →  Fl_O:tot:S     [Btu/(lbm·degR)]
  Cp [cal/(g·K)]    →  Fl_O:tot:Cp    [Btu/(lbm·degR)]
  composition [b0]  →  Fl_O:tot:composition
  ...

  EngUnitProps.setup_io() is called in configure() (not setup()) because the composition shape isn't known until after ChemEq sets up — configure() runs after the full system is assembled.

  ---
  Full chain diagram

  Input: composition (b0), T, P
             │
             ▼
  ┌─── SetTotalTP ──────────────────────────────────┐
  │  ┌── ChemEq (Newton solver) ─────────────────┐  │
  │  │  b0 + T + P                               │  │
  │  │  → minimize Gibbs energy                  │  │
  │  │  → n [species moles]                      │  │
  │  │  → n_moles                                │  │
  │  └───────────────────────────────────────────┘  │
  │  ┌── ThermoCalcs ─────────────────────────────┐  │
  │  │  n + n_moles + T + P                       │  │
  │  │  → h, S, gamma, Cp, Cv, rho, R  (SI)       │  │
  │  └───────────────────────────────────────────┘  │
  └─────────────────────────────────────────────────┘
             │  (for hP/SP/static modes: BalanceComp solves T here)
             ▼
  ┌─── EngUnitProps ───────────────────────────────┐
  │  SI → English unit conversion                  │
  │  promotes as  fl_name:tot:*                    │
  └────────────────────────────────────────────────┘
             │
             ▼
    Fl_O:tot:T, :P, :h, :S, :Cp, :gamma, :composition, ...
             ▼
    Fl_O:tot:T, :P, :h, :S, :Cp, :gamma, :composition, ...

  ---
  Key implication for your SOFC

  n and n_moles never reach EngUnitProps — they are internal to base_thermo. To use them (for mole fractions → activities) inside a segment you must connect explicitly:

  self.connect('thermo_an.base_thermo.n',       'nernst.n_species')
  self.connect('thermo_an.base_thermo.n_moles', 'nernst.n_moles')

  And the species ordering in n matches Properties(spec, init_elements=...).products — that list tells you which index is H2, H2O, O2, etc.