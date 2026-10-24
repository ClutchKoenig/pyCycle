# SOFC Implementation — Code Changes Log

Changes made to get `SegmentSOFC` to a functional `p.setup()` state.

---

## `sofc_balances.py`

### 1. `ICEnergyBalance` — typo in `initialize`

**Before:**
```python
def initalize(self):
    self.options.declare('seg_type', default='passive', values=['active', 'passive'])
```
**After:**
```python
def initialize(self):
    self.options.declare('seg_type', default='passive', values=['active', 'passive'])
```
**Why:** `initalize` was never called by OpenMDAO, so the option was never declared and `setup()` crashed on `self.options['seg_type']`.

---

### 2. `ICEnergyBalance` — renamed input variables to match promoted names

**Before:**
```python
self.add_input('Qdot_conv_IC_an',  val=0, units='W')
self.add_input('Qdot_conv_IC_cat', val=0, units='W')
self.add_input('Qdot_conduct_IC_left',  val=0, units='W')
self.add_input('Qdot_conduct_IC_right', val=0, units='W')
```
**After:**
```python
self.add_input('Q_conv_IC_A',   val=0, units='W', desc='Convective heat from anode to IC')
self.add_input('Q_conv_IC_C',   val=0, units='W', desc='Convective heat from cathode to IC')
self.add_input('Q_cond_IC_left',  val=0, units='W')
self.add_input('Q_cond_IC_right', val=0, units='W')
```
**Why:** Names now match what `HeatConvection` and `HeatConduction` promote out, so OpenMDAO can resolve connections via `promotes=['*']` in `SegmentSOFC`.

---

### 3. `PENEnergyBalance` — renamed conduction input variables

**Before:**
```python
self.add_input('Qdot_conduct_PEN_left',  val=0.0, units='W')
self.add_input('Qdot_conduct_PEN_right', val=0.0, units='W')
```
**After:**
```python
self.add_input('Qdot_cond_PEN_left',  val=0.0, units='W')
self.add_input('Qdot_cond_PEN_right', val=0.0, units='W')
```
**Why:** Consistent shorter naming convention; referenced consistently in `apply_nonlinear` and `linearize`.

---

## `sofc_heat.py`

### 4. `HeatConvectionChannel` — split `W_channel` into `W_in` / `W_out`

**Before:**
```python
self.add_input('W_channel', units='kg/s')
self.declare_partials('Q_conv_ch', ['W_channel', 'h_in', 'h_out'])

def compute(self, inputs, outputs):
    outputs['Q_conv_ch'] = inputs['W_channel'] * (inputs['h_out'] - inputs['h_in'])

def compute_partials(self, inputs, J):
    J['Q_conv_ch', 'W_channel'] = inputs['h_out'] - inputs['h_in']
    J['Q_conv_ch', 'h_in']      = -inputs['W_channel']
    J['Q_conv_ch', 'h_out']     =  inputs['W_channel']
```
**After:**
```python
self.add_input('W_in',  units='kg/s')
self.add_input('W_out', units='kg/s')
self.declare_partials('Q_conv_ch', ['W_in', 'W_out', 'h_in', 'h_out'])

def compute(self, inputs, outputs):
    outputs['Q_conv_ch'] = inputs['W_out'] * inputs['h_out'] - inputs['W_in'] * inputs['h_in']

def compute_partials(self, inputs, J):
    J['Q_conv_ch', 'W_out'] =  inputs['h_out']
    J['Q_conv_ch', 'h_out'] =  inputs['W_out']
    J['Q_conv_ch', 'W_in']  = -inputs['h_in']
    J['Q_conv_ch', 'h_in']  = -inputs['W_in']
```
**Why:** Anode gains mass from O²⁻ transfer, cathode loses it — `W_in ≠ W_out`. The correct enthalpy flux for a control volume with changing mass flow is `W_out·h_out − W_in·h_in`.

---

### 5. `HeatConvection` — updated promotes for `W_in`/`W_out` split

**Before:**
```python
self.add_subsystem('anode', HeatConvectionChannel(electrode='anode'),
                   promotes_inputs=[('W_channel', 'W_A'), ('h_in', 'h_A_in'), ('h_out', 'h_A_out')],
                   promotes_outputs=[('Q_conv_ch', 'Q_conv_A')])
self.add_subsystem('cathode', HeatConvectionChannel(electrode='cathode'),
                   promotes_inputs=[('W_channel', 'W_C'), ('h_in', 'h_C_in'), ('h_out', 'h_C_out')],
                   promotes_outputs=[('Q_conv_ch', 'Q_conv_C')])
```
**After:**
```python
self.add_subsystem('anode', HeatConvectionChannel(electrode='anode'),
                   promotes_inputs=[('W_in', 'W_in_A'), ('W_out', 'W_out_A'),
                                    ('h_in', 'h_A_in'), ('h_out', 'h_A_out')],
                   promotes_outputs=[('Q_conv_ch', 'Q_conv_A')])
self.add_subsystem('cathode', HeatConvectionChannel(electrode='cathode'),
                   promotes_inputs=[('W_in', 'W_in_C'), ('W_out', 'W_out_C'),
                                    ('h_in', 'h_C_in'), ('h_out', 'h_C_out')],
                   promotes_outputs=[('Q_conv_ch', 'Q_conv_C')])
```

---

### 6. `HeatConduction` — removed invalid `('', '')` placeholder promotes, fixed variable names

**Before:**
```python
self.add_subsystem('IC_Conduction', ...,
                   promotes_inputs=[('T_cell_left', 'T_IC_left'),
                                    ('T_cell_right', 'T_IC_right'),
                                    ('','')],          # crash: no match for empty string
                   promotes_outputs=[('Q_cond_struc_right', 'Q_cond_right_IC'),
                                     ('Q_cond_struc_left',  'Q_cond_left_IC')])
self.add_subsystem('PEN_Conduction', ...,
                   promotes_inputs=[()]                # syntax error / crash
                   promotes_outputs=[...])             # missing comma
```
**After:**
```python
self.add_subsystem('IC_Conduction', ...,
                   promotes_inputs=[('T_cell_left',  'T_IC_left'),
                                    ('T_cell_right', 'T_IC_right')],
                   promotes_outputs=[('Q_cond_struc_right', 'Q_cond_IC_right'),
                                     ('Q_cond_struc_left',  'Q_cond_IC_left')])
self.add_subsystem('PEN_Conduction', ...,
                   promotes_inputs=[('T_cell_left',  'T_PEN_left'),
                                    ('T_cell_right', 'T_PEN_right')],
                   promotes_outputs=[('Q_cond_struc_right', 'Q_cond_PEN_right'),
                                     ('Q_cond_struc_left',  'Q_cond_PEN_left')])
```
Also updated `ExecComp` expressions and `self.connect` calls to use the new `Q_cond_IC_left/right` / `Q_cond_PEN_left/right` naming convention (was `Q_cond_left_IC` etc.).

---

## `sofc_reaction.py`

### 7. `AreaSpecificResistanceOverpotential` — filled `val=,` placeholders

**Before:**
```python
self.add_input('R_ref_asr',  val=,  units='ohm',      desc='Reference asr resistance')
self.add_input('R_c_asr',    val=,  units='ohm',      desc='Reference asr resistance')
self.add_input('ASR_nom_750',val=,  units='ohm*m**2', desc=)
self.add_input('i',          val=,  units='A/m/m',    desc='Current density')
self.add_output('ASR',       val=,  units='ohm*m**2', desc=)
self.add_output('eta_asr',   val=,  units='V',        desc='')
```
**After:**
```python
self.add_input('R_ref_asr',  val=1.0,  units='ohm',      desc='Reference asr resistance')
self.add_input('R_c_asr',    val=1.0,  units='ohm',      desc='Reference asr resistance')
self.add_input('ASR_nom_750',val=1e-4, units='ohm*m**2', desc='Nominal ASR at 750 C')
self.add_input('i',          val=1.0,  units='A/m**2',   desc='Current density')
self.add_output('ASR',       val=1e-4, units='ohm*m**2', desc='Area specific resistance')
self.add_output('eta_asr',   val=0.1,  units='V',        desc='ASR overpotential')
```
**Why:** `val=,` is a Python syntax error. Placeholder defaults added so the file can import.

---

### 8. `CurrentDensityCalc` — wrong base class

**Before:**
```python
class CurrentDensityCalc(Element):
```
**After:**
```python
class CurrentDensityCalc(om.ExplicitComponent):
```
**Why:** `Element` is a pyCycle class for flow-station components. It has a different `setup()` mechanism and doesn't have `add_input`/`add_output` at instantiation time — causing `AttributeError: 'CurrentDensityCalc' object has no attribute 'add_input'`.

---

### 9. `CurrentDensityCalc` — filled `val=,` placeholders

**Before:**
```python
self.add_input('A', val=,   units='m**2',  desc='...')
self.add_output('i', val=,  units='A/m/m', desc='Current density in segment')
```
**After:**
```python
self.add_input('A', val=0.01, units='m**2',   desc='...')
self.add_output('i', val=1.0, units='A/m**2', desc='Current density in segment')
```

---

### 10. Removed incomplete `segment_calc` stub

**Before:**
```python
class segment_calc

class ElectroChemistry(om.Group):
```
**After:**
```python
class ElectroChemistry(om.Group):
```
**Why:** Bare `class segment_calc` with no body is a syntax error that prevented the module from importing.

---

## `sofc_segment.py`

Complete rewrite. The old file had an incomplete skeleton (`sofc_segment`) with placeholder `promotes_inputs=['']`, missing commas, undefined variables, and no internal connections. The new `SegmentSOFC` is a fully wired `om.Group` that passes `p.setup()`.

Key structural additions over the old file:

| What | Old | New |
|---|---|---|
| Class name | `sofc_segment` | `SegmentSOFC` |
| `options.declare` typo | `disc=` | `desc=` |
| `Thermo` mode | `'Total_TP'` | `'total_TP'` |
| Composition options | missing | `anode_composition`, `cathode_composition`, `spec` added |
| `Properties` objects | undefined | built from `output_port_data()` |
| Species indices | missing | computed from `products.index()` for `src_indices` |
| `anode_rxn`/`cathode_rxn` promotes | incomplete, missing commas | fully specified |
| Thermo `promotes_inputs` | missing `composition` input | added |
| Thermo `promotes_outputs` | included `Fl_O:tot:T` (wrong — T is input) | removed |
| `anode_species_out`, `cathode_species_out` | missing | added (`SpeciesFlowCalc`) |
| `anode_bulk`, `cathode_bulk` | missing | added (`BulkComposition`) |
| `Convection` promotes | wrong names, missing commas | `W_in/out_A/C`, `h/T_A/C_in/out` |
| `Conduction` | had empty `promotes_inputs=[]` | `promotes=['*']` |
| `PEN`, `IC` | no `seg_type` option passed | `PENEnergyBalance(seg_type=segment_type)` |
| `ChannelEnergyBalance` | no electrode option, no promotes | `electrode='anode'/'cathode'`, full promotes |
| `Electrochemistry` | only in active, manual `self.connect` (now promoted) | `promotes_inputs/outputs` with `P_cat` rename |
| `Utilization` | missing | added (`SpeciesUtilization`) |
| Internal connections | none | `composition_out→Thermo`, `base_thermo.n/n_moles→SpeciesFlowCalc`, `src_indices` for `x_bulk→Electrochemistry` |

---

## `sofc.py` — 2026-09-02

### 11. `default_des_od_conns` — syntax error

**Before:**
```python
self.default_des_od_conns = 
[
    #tbd
]
```
**After:**
```python
self.default_des_od_conns = [
    #tbd
]
```
**Why:** `=` followed directly by a newline with no bracket open is a `SyntaxError` — the module couldn't be imported at all.

---

### 12. Standalone/root-element fallback for inlet composition

**Added:**
```python
def _seed_default_inflow_composition(self):
    """
    Fall back to the in_composition_* options when nothing upstream has
    populated Fl_I_data (e.g. SOFC used as a standalone/root element with
    no element connected to its Fl_I_an/Fl_I_cat ports). No-op once
    Fl_I_data has already been set, e.g. by the pyCycle graph wiring.
    """
    self.Fl_I_data.setdefault('Fl_I_an', self.options['in_composition_anode'])
    self.Fl_I_data.setdefault('Fl_I_cat', self.options['in_composition_cathode'])
```
Called at the top of both `pyc_setup_output_ports()` and `setup()`.

**Why:** `Fl_I_data` is normally populated by the pyCycle graph-wiring machinery (`mp_cycle.py:151`) from whatever element is connected upstream. With no upstream element (standalone use / unit tests), `self.Fl_I_data['Fl_I_an']` raised `KeyError`. The two `in_composition_anode`/`in_composition_cathode` options existed but were previously never read anywhere — dead code.

---

### 13. Per-segment composition kwarg was static instead of chained

**Before:**
```python
for i, name in enumerate(self.segments):
    if i < (N_passive - 1):
        self.add_subsystem(name, SegmentSOFC(type='passive', ...,
                            anode_composition   = self.Fl_I_data['Fl_I_an'],
                            cathode_composition = self.Fl_I_data['Fl_I_cat'], ...))
    elif (N_passive - 1) < i <= (N_passive + N_active - 1):
        self.add_subsystem(name, SegmentSOFC(type='active', ..., ...))   # same composition kwargs, duplicated block
    elif i > (N_passive + N_active - 1):
        self.add_subsystem(name, SegmentSOFC(type='passive', ..., ...)) # same composition kwargs, duplicated block
```
**After:**
```python
anode_comp = dict(self.Fl_I_data['Fl_I_an'])
cathode_comp = dict(self.Fl_I_data['Fl_I_cat'])
...
for i, name in enumerate(self.segments):
    seg_type = 'active' if (N_passive - 1) < i <= (N_passive + N_active - 1) else 'passive'
    seg = SegmentSOFC(type=seg_type, ...,
                       anode_composition=anode_comp, cathode_composition=cathode_comp, ...)
    self.add_subsystem(name, seg, promotes_input=['n_cell']) if seg_type == 'active' else self.add_subsystem(name, seg)
    anode_comp = SOFCThermoAdd(reaction_type='anode', inflow_composition=anode_comp).output_port_data()
    cathode_comp = SOFCThermoAdd(reaction_type='cathode', inflow_composition=cathode_comp).output_port_data()
```
**Why:** Every segment was built from the *original* SOFC inlet composition dict instead of its own upstream segment's outlet composition. `SOFCThermoAdd.output_port_data()` always adds `'O'` to the anode outlet element set (oxygen ions arriving from the cathode), so every segment's `composition_out_A` was sized 1 element larger than every segment's `composition_in_A`. Chaining `segment[i].composition_out_A → segment[i+1].composition_in_A` then hit a shape mismatch for every anode connection past segment 0. Also collapsed three duplicate `add_subsystem` branches (identical except `type=`) into one.

---

### 14. Segment0 inlet-properties block computed but never connected

**Before:** `anode_thermo_in`/`cathode_thermo_in`/`anode_species_in`/`cathode_species_in` were added with outputs promoted to group-level `h_A_in`/`x_in_A`/`n_in_A`/etc., but nothing connected them to `segment0`. Meanwhile `segment0` was fed directly from the raw flowstation with a name-swapped typo:
```python
self.connect('Fl_I_an:tot:h', 'segment0.h_in_A')   # SegmentSOFC's real port is h_A_in, not h_in_A
self.connect('Fl_I_cat:tot:h', 'segment0.h_in_C')  # real port is h_C_in
```
`x_in_A`/`n_in_A` (species-basis mole fraction / molar flow, needed because segment0 has no upstream segment to inherit them from) weren't connected to `segment0` by either path at all.

**After:**
```python
self.connect('h_A_in',  'segment0.h_A_in')
self.connect('x_in_A',  'segment0.x_in_A')
self.connect('n_in_A',  'segment0.n_in_A')
...
self.connect('h_C_in',  'segment0.h_C_in')
self.connect('x_in_C',  'segment0.x_in_C')
self.connect('n_in_C',  'segment0.n_in_C')
```
**Why:** `SegmentSOFC` needs species-basis `x_in_A`/`n_in_A`, not just elemental composition — segments 1+ get these from the previous segment's `x_out_A`/`n_out_A`, but segment0 has nothing upstream, hence the dedicated `Properties`/`Thermo`/`SpeciesFlowCalc` block. It computed the right values but never delivered them; also fixed the `h_A_in`/`h_in_A` name swap (confirmed against `sofc_segment.py:172`, where `Convection` promotes `h_A_in`/`h_C_in`).

---

### 15. `cathode_thermo_in` promoted T/P names collided with anode's

**Before:**
```python
self.add_subsystem('cathode_thermo_in', cathode_thermo_in,
           promotes_inputs=[('T', 'T_A_in'), ('P', 'p_A_in')],   # anode's names, copy-paste bug
           ...)
```
**After:**
```python
self.add_subsystem('cathode_thermo_in', cathode_thermo_in,
           promotes_inputs=[('T', 'T_C_in'), ('P', 'p_C_in')],
           ...)
```
**Why:** Promoting both components' T/P inputs to the same group-level name would unify them into one variable — cathode's inlet temperature/pressure would have silently followed anode's.

---

### 16. Missing W/T/P connections into the segment0-properties block

**Added:**
```python
self.connect('Fl_I_an:stat:W',   'W_in_A')
self.connect('Fl_I_cat:stat:W',  'W_in_C')
self.connect('Fl_I_an:tot:T',    'T_A_in')
self.connect('Fl_I_an:tot:P',    'p_A_in')
self.connect('Fl_I_cat:tot:T',   'T_C_in')
self.connect('Fl_I_cat:tot:P',   'p_C_in')
```
**Why:** `anode_species_in`/`cathode_species_in` promoted their own `W` input as group-level `W_in_A`/`W_in_C`, and `anode_thermo_in`/`cathode_thermo_in` promoted `T`/`P` similarly — none of these were ever connected to the real flowstation, so `x_in_A`/`n_in_A`/`h_A_in`/etc. (wired to segment0 in #14) were being computed from unconnected default values, not the actual inlet conditions.

---

### 17. Duplicate `W_in_A`/`W_in_C` connections

**Before:** #16's connections plus a second, redundant copy further down in the flowstation-connections block:
```python
self.connect('Fl_I_an:stat:W', 'W_in_A')          # duplicate of #16
self.connect('Fl_I_an:stat:W', 'segment0.W_in_A')
```
**After:** duplicate removed, only the `segment0.W_in_A` connect remains in that block (`W_in_A` is already wired once, in #16).
**Why:** Two `connect()` calls to the same target raises an OpenMDAO "already connected" error at setup.

---

### 18. Segment-chain and segment0 pressure connections targeted nonexistent ports

**Before:**
```python
self.connect(f'{name}.p_A_out', f'{self.segments[i + 1]}.P_A_in')   # P_A_in doesn't exist
self.connect(f'{name}.p_C_out', f'{self.segments[i + 1]}.P_C_in')   # P_C_in doesn't exist
...
self.connect('Fl_I_an:tot:P', 'segment0.')     # empty target
self.connect('Fl_I_cat:tot:P', 'segment0.')    # empty target
```
**After:**
```python
self.connect(f'{name}.p_A_out', f'{self.segments[i + 1]}.p_A_out')
self.connect(f'{name}.p_C_out', f'{self.segments[i + 1]}.p_C_out')
...
self.connect('Fl_I_an:tot:P',  'segment0.p_A_out')
self.connect('Fl_I_cat:tot:P', 'segment0.p_C_out')
```
**Why:** `SegmentSOFC` only exposes a single `p_A_out`/`p_C_out` port per electrode (channel pressure, no drop assumed across the segment) — there's no separate `P_A_in`/`P_C_in` port to target.

---

### 19. `i_sum` — malformed `ExecComp`, wrong segment set, undefined `U_cell`/`I`

**Before:**
```python
self.add_subsystem('i_sum', om.ExecComp(f'i_sum = {i}' for ))   # SyntaxError

self.add_subsystem('Ucell', om.ExecComp('Ucell=Vcell'))          # unconnected, unrelated leftover

balances = self.add_subsystem('balances', om.BalanceComp())
for i, segment in enumerate(self.segments):        # includes passive segments
    balances.add_balance(f'{segment}_I', val=1, units='A', eq_units='V')
    self.connect(f'balances.{segment}_I', f'{segment}.I')       # I doesn't exist on passive segments
    self.connect(f'{segment}.V_cell', f'balances.lhs:{segment}_I')  # V_cell doesn't exist on passive segments
    self.connect('U_cell', f'balances.rhs:{segment}_I')         # 'U_cell' is not produced by anything

balances.add_balance('U_cell', val=1, units='V', eq_units='A', lower=0.3, upper=1.2, use_mult=False)
self.connect('balances.U_cell', 'U_cell')     # 'U_cell' is not consumed by anything
self.connect('i_sum.I_total', 'balances.lhs:U_cell')
self.connect('I', 'balances.rhs:U_cell')      # 'I' is not produced by anything
```
**After:**
```python
segment_current_inputs = {f'I_{segment}': {'units': 'A'} for segment in self.segments_active}
i_sum_expr = 'I_total = ' + ' + '.join(segment_current_inputs)
self.add_subsystem('i_sum', om.ExecComp(i_sum_expr, I_total={'units': 'A'}, **segment_current_inputs))

balances = self.add_subsystem('balances', om.BalanceComp(),
                              promotes_inputs=[('rhs:U_cell', 'I')])
for segment in self.segments_active:      # active segments only
    balances.add_balance(f'{segment}_I', val=1, units='A', eq_units='V')
    self.connect(f'balances.{segment}_I', f'{segment}.I')
    self.connect(f'{segment}.V_cell', f'balances.lhs:{segment}_I')
    self.connect('balances.U_cell', f'balances.rhs:{segment}_I')   # state fed back into its own component's other inputs
    self.connect(f'balances.{segment}_I', f'i_sum.I_{segment}')

balances.add_balance('U_cell', val=1, units='V', eq_units='A', lower=0.3, upper=1.2, use_mult=False)
self.connect('i_sum.I_total', 'balances.lhs:U_cell')
```
Also removed the unfinished module-level stub:

```python
class i_sum(om.ExplicitComponent):
    def initialize(self):
        self.options.declare('N_active', default=10, desc='Number of active Segments')
    def setup(self):
        #TODO: finish
```
**Why:** Three compounding bugs. (1) `i_sum`'s `ExecComp` call was a syntax error, never built. (2) The per-segment balance loop iterated *all* segments, but `I`/`V_cell` only exist on active-type `SegmentSOFC` instances (`sofc_segment.py:204` gates `Electrochemistry` behind `if segment_type == 'active'`) — would have crashed connecting nonexistent ports on passive segments. (3) `'U_cell'` and `'I'` were used as connect sources without any subsystem producing them as outputs — `connect()` requires a real output source, not just a label. Fixed by looping only `self.segments_active`, wiring the shared voltage as `balances.U_cell` fed back into each segment's own `rhs:{segment}_I` input (a self-loop on the same `BalanceComp` instance — the standard way to broadcast one state to several inputs of the same component), and promoting `I` directly via `promotes_inputs=[('rhs:U_cell', 'I')]` so it's a settable top-level input. The result implements two physical conditions: **equipotential** — every active segment's current is solved so its own `V_cell` matches the shared `U_cell`; **sum-of-currents** — `U_cell` itself is solved so the sum of all active-segment currents (`i_sum.I_total`) matches the total stack current `I`. The unused, unrelated `Ucell` `ExecComp` fragment and the dead module-level `i_sum` class stub (superseded by the `ExecComp`, never instantiated) were removed.

---

## `sofc_heat.py` / `sofc_segment.py` / `sofc.py` — 2026-09-03

### 20. Lateral PEN/IC conduction computed net-zero heat into every segment, always

**The bug:** `HeatConduction` (`sofc_heat.py`) built one `HeatConductionStructure` per
segment per structure (PEN, IC). That single instance took `T_cell_left =
T_PEN_left` (the **left neighbor's** own `T_PEN`) and `T_cell_right = T_PEN_right`
(the **right neighbor's** own `T_PEN`) and computed one flux directly between those
two neighbor states — this segment's own `T_PEN` never entered the calculation at
all. Both signed halves of that one flux (`Q_cond_struc_left`/`_right`) were then
fed back into *this same segment's own* `PENEnergyBalance`/`ICEnergyBalance`
residual as `Q_cond_PEN_left`/`Q_cond_PEN_right`. Since both numbers came from the
same `Q`, they were always `-Q` and `+Q` and canceled exactly:

```python
# Before, inside HeatConduction.setup() — one instance per segment per structure:
self.add_subsystem('PEN_Conduction', HeatConductionStructure(structure='PEN', ...),
                   promotes_inputs=[('T_cell_left',  'T_PEN_left'),   # left neighbor's T
                                    ('T_cell_right', 'T_PEN_right')], # right neighbor's T
                   promotes_outputs=[('Q_cond_struc_right', 'Q_cond_PEN_right'),
                                     ('Q_cond_struc_left',  'Q_cond_PEN_left')])
# dT = T_PEN_left - T_PEN_right; Q_cond_PEN_left = -Q, Q_cond_PEN_right = +Q
# => Q_cond_PEN_left + Q_cond_PEN_right == 0, always, regardless of this
#    segment's own T_PEN (which never appears anywhere in this component).
```

Concretely: a segment sitting at 1000 K between two 900 K neighbors (e.g. an
active-zone hot spot) reported **zero** net conductive heat loss, because the
formula only compared the two neighbors to *each other*, never to this segment.
Lateral conduction never actually affected any segment's PEN/IC temperature.
`Qc_left`/`Qc_right` (small diagnostic `ExecComp`s) did show the real nonzero
magnitude, but were confirmed dead ends — never consumed anywhere else in the
model.

**The fix:** compute each interface's conduction resistor exactly once, anchored
to *this segment's own* state on one side, and cross-wire the result to the
neighbor instead of independently recomputing it there — mirroring the pattern
already used for `T_PEN` itself (already cross-wired segment-to-segment in
`sofc.py`). Each segment now builds **at most one link per side**, gated by four
new options on `HeatConduction`/`SegmentSOFC`: `pen_left_bc`, `pen_right_bc`,
`ic_left_bc`, `ic_right_bc` (values `'T'`/`'Q'`):

```python
# After — left link: T_cell_right is now this segment's OWN T_PEN.
if self.options['pen_left_bc'] == 'T':      # default for every segment
    self.add_subsystem('PEN_Conduction_left', HeatConductionStructure(structure='PEN', ...),
                       promotes_inputs=[('T_cell_left',  'T_PEN_left'),  # neighbor or boundary constant
                                        ('T_cell_right', 'T_PEN')],      # this segment's own state
                       promotes_outputs=[('Q_cond_struc_right', 'Q_cond_PEN_left'),      # entering ME
                                         ('Q_cond_struc_left',  'Q_cond_PEN_left_out')])  # entering my left neighbor
# Right link (only ever needed by the last segment in a stack, pen_right_bc='T'):
if self.options['pen_right_bc'] == 'T':     # default 'Q' everywhere else
    self.add_subsystem('PEN_Conduction_right', HeatConductionStructure(structure='PEN', ...),
                       promotes_inputs=[('T_cell_left',  'T_PEN'),         # this segment's own state
                                        ('T_cell_right', 'T_PEN_right')],  # boundary constant
                       promotes_outputs=[('Q_cond_struc_left', 'Q_cond_PEN_right')])
```

Interior segments now rely on their right neighbor's left link for
`Q_cond_PEN_right`/`Q_cond_IC_right`, cross-wired once in `sofc.py`:

```python
self.connect(f'{segments[i + 1]}.Q_cond_PEN_left_out', f'{name}.Q_cond_PEN_right')
self.connect(f'{segments[i + 1]}.Q_cond_IC_left_out',  f'{name}.Q_cond_IC_right')
```

instead of the old (removed) reverse `T_PEN`/`T_IC` cross-wire
(`self.connect(f'{name}.T_PEN', f'{segments[i+1]}.T_PEN_right'... )` and back) that
fed the now-deleted single-shared-component design. Interface conservation is now
automatic by construction (one value, connected both ways) rather than relying on
two independently-computed sides agreeing.

**No changes needed to `PENEnergyBalance`/`ICEnergyBalance`** (`sofc_balances.py`)
— `Q_cond_PEN_left`/`_right` and `Q_cond_IC_left`/`_right` were already plain
inputs with defaults; only what connects to them changed.

### 21. Temperature-vs-flux boundary type switch on the open stack faces

Once `Q_cond_PEN_left`/`_right` were genuinely free inputs (not always produced
internally — see #20), the open faces at each end of the stack (segment0's left
face, `segments[-1]`'s right face) could support a flux-specified ('Q') boundary
in addition to the existing temperature-specified ('T') one, by simply not
building the local link and hardwiring the flux instead — the same way the
temperature is hardwired today via `set_input_defaults`.

Added to `SOFC` (`sofc.py`): `PEN_BC_type_left`/`_right`, `IC_BC_type_left`/`_right`
(`'T'`/`'Q'`, default `'T'` — preserves prior default behavior of both boundaries
hardwired to 973 K) and `Q_PEN_BC_left`/`_right`, `Q_IC_BC_left`/`_right` (hardwired
flux values, default 0.0 W, used only when the corresponding type is `'Q'`).
`setup()` passes the relevant `*_bc` option only to `segment0` (left) and the last
segment (right); every other segment keeps `HeatConduction`'s defaults
(`pen_left_bc='T'`, `pen_right_bc='Q'`), which is exactly correct for interior
segments regardless of what the two true boundaries are doing.

**Test fallout:** `SegmentSOFC_test.py` builds one standalone segment and sets
both `T_PEN_right`/`T_IC_right` directly — under the new default
(`pen_right_bc`/`ic_right_bc == 'Q'`), those inputs no longer exist by default, so
the test now passes `pen_right_bc='T', ic_right_bc='T'` explicitly, and the
absolute-path geometry sets (`Conduction.PEN_Conduction.*` /
`Conduction.IC_Conduction.*`) were updated to the new
`Conduction.PEN_Conduction_left/right.*` / `Conduction.IC_Conduction_left/right.*`
subsystem names (both sides set identically, since this standalone test now
builds both links).

**Verification performed:** `SegmentSOFC_test.py` converges (T_PEN/T_IC ≈ 1073 K,
isothermal start, as before). An isolated `HeatConduction` check confirmed a
1000 K node between two 900 K neighbors now reports `Q_cond_PEN_left = -4.0 W`
(previously exactly 0 regardless of temperature) with the two link halves being
exact opposites. `SOFC(...).setup()` (with `pyc_setup_output_ports()` called
manually, since driving it through the full pyCycle graph-wiring machinery is
still WIP — see the "composition propagation" discussion in git history) resolves
cleanly for the default `'T'/'T'` boundaries as well as `PEN/IC_BC_type_left='Q'`
and `PEN/IC_BC_type_right='Q'`.

**Aside (not fixed, out of scope):** found a pre-existing circular import between
`sofc_segment.py` and `sofc_api.py` (each imports the other at module load time)
that breaks `from pycycle.elements.sofc_segment import SegmentSOFC` as a script's
first pyCycle import — reproducible on the pre-#20 code too (confirmed via `git
stash`), so unrelated to this change. Workaround used for verification here:
`import pycycle.api; import pycycle.sofc_api` before importing `sofc_segment`
directly.

### 22. Collapsed the PEN/IC boundary-type options into one per face

Items #20–#21 as first implemented gave `HeatConduction`/`SegmentSOFC` four
separate options (`pen_left_bc`, `pen_right_bc`, `ic_left_bc`, `ic_right_bc`) and
`SOFC` eight (`PEN_BC_type_left/right`, `IC_BC_type_left/right`, plus
`Q_PEN_BC_left/right`, `Q_IC_BC_left/right` for the hardwired flux values) — noted
as feeling like option clutter on review. Simplified:

- `HeatConduction`/`SegmentSOFC`: `pen_*_bc`/`ic_*_bc` collapsed to a single
  `left_bc`/`right_bc` per segment (PEN and IC always share one boundary *type*
  on a given face — they can still take different geometry/temperature values,
  just not a different `'T'`-vs-`'Q'` choice on the same face).
- `SOFC`: `PEN_BC_type_left/right`/`IC_BC_type_left/right` collapsed to
  `BC_type_left`/`BC_type_right`. The four `Q_PEN_BC_*`/`Q_IC_BC_*` value options
  were dropped entirely — no functionality lost, since
  `PENEnergyBalance`/`ICEnergyBalance` already default `Q_cond_PEN_left`/`right`
  and `Q_cond_IC_left`/`right` to `0.0 W` (a sensible adiabatic default for a
  `'Q'`-type boundary), and a nonzero flux can be set post-`setup()` with
  `prob.set_val('<path>.segment0.Q_cond_PEN_left', val, units='W')` — the same
  pattern `SegmentSOFC_test.py` already uses for the temperature-type boundary.

Net: `SegmentSOFC`/`HeatConduction` options 4→2; `SOFC`'s new boundary-type
surface 8→2 (on top of the 4 pre-existing `T_*_BC_*` value options, unchanged).
Re-verified: `SegmentSOFC_test.py` still converges to the same result;
`SOFC(...).setup()` (manual `pyc_setup_output_ports()`, as in #21) still resolves
cleanly for default, `BC_type_left='Q'`, `BC_type_right='Q'`, and both at once.

---

## `sofc.py` — 2026-09-03 — new `SOFC_test.py`, and `SOFC.setup()` had never actually been run

### 23. `Problem.setup()` alone does not validate connections — `SOFC` had multiple broken ones, undetected until now

Building `pycycle/elements/SOFC_test.py` (a small `pyc.Cycle` wiring
`pyc.FlowStart` (anode, H2) + `pyc.FlowStart` (cathode, air) into `SOFC` via
`pyc_connect_flow`, the idiomatic pattern also used by
`pycycle/elements/test/test_combustor.py`) revealed that `SOFC.setup()`'s wiring
had never actually been exercised by `run_model()` before. The earlier
verification in items #20–#22 only ever called `Problem.setup()`, which does
**not** validate connections — that check only happens in `final_setup()`,
triggered automatically inside `run_model()`. So `p.setup()` reporting "OK"
proved the model *built*, not that it could *run*.

`SOFC_test.py` surfaced four distinct pre-existing bugs, unrelated to #20–#22:

**(a) Multiple `self.connect(<input>, <target>)` calls where the source was
itself an unconnected external input**, not an output — invalid in OpenMDAO
(`connect()` only goes output→input). `Fl_I_an:tot:T/P`, `Fl_I_an:stat:W`,
`Fl_I_an:tot:composition` (and cathode equivalents) are `FlowIn`'s own promoted
*inputs* (meant to be filled externally, by `pyc_connect_flow` or `set_val`), and
`p_A`/`p_C` (the shared stack-wide pressure, promoted from every segment's
`p_A_out`/`p_C_out`) were likewise plain inputs with nothing set as their source.
Fixed by promoting each downstream alias (`T_A_in`, `W_in_A`,
`composition_in_A`, segment `p_A_out`/`p_C_out`, `anode_out`/`cathode_out`'s `P`)
directly to the *same name* as the real port (e.g. `Fl_I_an:tot:P`) instead of
connecting between two inputs — this also eliminates the separate `p_A`/`p_C`
names entirely (every segment's channel pressure, and the outlet Thermo's `P`,
is now literally `Fl_I_an:tot:P`/`Fl_I_cat:tot:P`, matching "no pressure drop
anywhere" physically). Also added a tiny `ExecComp` pass-through for
`Fl_O_an:stat:W`/`Fl_O_cat:stat:W` (`total_TP` mode has no `W` in/out to promote
from, so there was nothing to connect the last segment's outlet mass flow into).

**(b) `self.connect('balances.U_cell', f'balances.rhs:{segment}_I')`** connected
one output of the `balances` `BalanceComp` to another input on that *same*
component instance — also invalid (`connect()` requires two different
subsystems). Fixed by promoting `rhs:{segment}_I` (for every active segment) and
`U_cell`'s own output to one shared name (`'U_cell'`) at the `SOFC` parent level
instead — the standard OpenMDAO idiom for broadcasting one component's own
solved value back into another input on itself.

**(c) `n_in_A`/`n_in_C` only exist on active segments** (only promoted inside
`SegmentSOFC`'s `if segment_type == 'active':` block, for `H2_Utilization`/
`O2_Utilization`), but the inter-segment loop and the segment0 initial-flowstation
block both connected into them unconditionally — broke the moment a passive
segment received the connection. Fixed by checking
`self.segment_types[target] == 'active'` (a new `self.segment_types` dict, built
alongside the existing per-segment type loop) before making either connection.

**(d) NOT fixed — composition propagation.** Once (a)–(c) were fixed,
`SOFC_test.py` runs all the way to the one remaining, already-known issue:
`anode_species_in.x_i` (shape `(2,)`, sized from the raw inlet composition) feeds
every segment's `anode_bulk.x_in` (shape `(8,)`, sized from the *outlet*-basis
composition), which only happens to work for a single always-active segment
(`SegmentSOFC_test.py`) and breaks for any real passive/active mix — exactly the
issue already flagged in this file's own comment ("I made a big error in the way
compositions are passed... every segment gets the same composition") and the
"discussion on the state of composition propagation" commit. Deliberately left
alone per the user's explicit call — they're already mid-way through that
redesign themselves.

**Also added**, separately: `SOFC.setup()` never assigned a top-level solver, so
even with (a)–(c) fixed the segment-to-segment coupling (chained `T_PEN`/`T_IC`,
cross-wired `Q_cond_*_left_out`, and the equipotential/sum-of-currents
`balances`) would never have converged — added
`self.nonlinear_solver = om.NonlinearBlockGS(...)` /
`self.linear_solver = om.LinearBlockGS()` at the end of `setup()`, matching the
already-documented intent for a multi-segment stack ("outer NonlinearBlockGS
sweeps segments").

---

## `sofc.py` — 2026-09-03 — composition propagation fix, and NLBGS replaced with Newton

### 24. Segment0's inlet composition: sized on the wrong (inlet, not outlet) element basis, and never given real values

Picking up item 23(d) (the composition-propagation bug, previously left alone
on request): traced to exactly one root cause, isolated to segment0's own
scratch composition block in `sofc.py` (`anode_in_props`, `anode_thermo_in`,
`anode_species_in`, and cathode equivalents) — not a general architecture
problem. For segment `i >= 1`, `x_in_A`/`n_in_A` are correctly inherited from
segment `i-1`'s own outlet (already consistently sized, since `sofc.py`'s
segment-creation loop already re-chains each segment's construction-time
`anode_composition` from the previous one's `output_port_data()`). Segment0 has
no upstream segment, so it needs this separate scratch block — which was built
from the **raw** inlet dict (e.g. `{'H': 1.0}`, 2 species), while segment0's own
internal `anode_out_props` (used for `anode_bulk.x_in`, per the "all species
calcs inside a segment use the outlet basis" convention) is sized from
`SOFCThermoAdd(...).output_port_data()`, which always adds an `'O'` key for the
anode (`sofc_thermo_add.py`: oxygen arrives from the cathode even before any
current flows) — 8 species. Fixed by building the scratch block from
`SOFCThermoAdd(reaction_type='anode', inflow_composition=anode_comp).output_port_data()`
instead of the raw dict — this adds `'O'` defaulted to exactly `0.0` (correctly
representing "no reaction yet" rather than a placeholder), matching segment0's
own basis exactly. Cathode's `output_port_data()` doesn't add anything
(cathode doesn't gain a new element, only loses O amount at runtime), so this
part is a no-op for cathode — done anyway for symmetry.

Separately, in the same block: `anode_species_in`/`cathode_species_in`
(`SpeciesFlowCalc`) never had their `n_i`/`n_moles` inputs connected to
anything, so they silently sat at component defaults (`n_i=1.0` uniform,
`n_moles=0.01`) instead of real CEA-equilibrium values — every other use of
`SpeciesFlowCalc` in this codebase wires it from a `Thermo`'s
`base_thermo.n`/`base_thermo.n_moles` (e.g. `sofc_segment.py`'s
`anode_thermo_out.base_thermo.n → anode_species_out.n_i`); these two just never
got that wiring. Added the missing `self.connect(...)` calls.

**Verified:** `SOFC_test.py` runs past the shape-mismatch error and reaches
`NLBGS Converged` end-to-end.

### 25. NonlinearBlockGS silently never solves the equipotential/current `BalanceComp`

Getting far enough to actually inspect results (only possible after #24)
surfaced that item 23's `NonlinearBlockGS` solver choice was wrong:
`i_sum.I_total`/`balances.U_cell`/`balances.{segment}_I` all sat exactly at
their `BalanceComp` initial guesses (`val=1`) after `run_model()` — despite
`NLBGS` printing `Converged` with a tiny residual. Cause: `NonlinearBlockGS` has
no mechanism to invert an `ImplicitComponent`'s residual (it isn't a Newton-type
solver — it just re-executes subsystems in sequence and checks the combined
residual norm); the segments' thermal states genuinely do converge under it,
and their much larger residual dominates the combined norm NLBGS reports,
masking that `balances`' own states never moved at all. This is a worse failure
mode than an outright error — it silently reports success while the electrical
solution is simply wrong (e.g. `i_sum.I_total` stayed at `2.0` while `sofc.I`
requested `0.1`).

Fixed by replacing the top-level solver with a flattened Newton solve, which
can actually linearize and invert implicit states wherever they live in the
tree:
```python
newton = self.nonlinear_solver = om.NewtonSolver(solve_subsystems=True)
newton.options['maxiter'] = 50
newton.options['atol'] = 1e-6
newton.options['rtol'] = 1e-6
newton.linesearch = om.ArmijoGoldsteinLS(iprint=0)
self.linear_solver = om.DirectSolver()
```
`solve_subsystems=True` lets each segment's own inner `NewtonSolver` still
reconverge locally each outer iteration, same as before.

**Verified:** re-ran `SOFC_test.py` — `i_sum.I_total = 0.09999999999999991`
(matches the requested `sofc.I = 0.1 A` to solver tolerance, vs. `2.0` before
this fix), `U_cell` and each active segment's `V_cell`/`I` now genuinely moved
from their initial guesses and satisfy the equipotential condition.
