# SOFC — Port Propagation & Segment Setup: Review Notes

Notes from a review/discussion session covering (1) how pyCycle's flow-port
metadata propagates in general, and (2) bugs and design gaps found in
`SOFC.setup()`'s segment-building loop (`sofc.py`) and how it feeds
composition into `SegmentSOFC` (`sofc_segment.py`).

---

## 1. How `Fl_O_data` / `Fl_I_data` propagate in pyCycle (background)

Two separate, parallel systems exist for every flow port — this distinction is
the key to everything below:

1. **Metadata dict (`Fl_O_data` / `Fl_I_data`)** — plain Python dicts holding
   the *element/species set* (`b0` composition), used only to size and shape
   components correctly **before** `setup()` runs anywhere. Populated via
   `copy_flow()` / `init_output_flow()` inside each element's
   `pyc_setup_output_ports()` (`element_base.py`).
2. **Real OpenMDAO variables (`Fl_O:tot:*`, `Fl_O:stat:*`)** — the actual
   numeric ports, built in `setup()` by adding/promoting `Thermo`/`ThermoAdd`
   subsystems. Wired element-to-element via `pyc_connect_flow()`
   (`mp_cycle.py:156`).

The framework's `Cycle.setup()` (`mp_cycle.py:74-153`) does a breadth-first
walk of a connection graph: for every element node it calls
`pyc_setup_output_ports()`, and for every `pyc_connect_flow` edge it copies
`target_element.Fl_I_data[in_port] = src_element.Fl_O_data[out_port]`
(`mp_cycle.py:151`) — **before** the downstream element's own
`pyc_setup_output_ports()` runs. That's why `Compressor`/`Inlet` can safely
read `self.Fl_I_data['Fl_I']` inside their own `pyc_setup_output_ports()` /
`setup()` without ever assigning it themselves.

The chain bottoms out at whichever element has no upstream connection
(`in_degree == 0` in the flow graph, `mp_cycle.py:87`) — typically
`FlowStart`. That element never reads `Fl_I_data` at all; it builds its
`Fl_O_data` straight from its own `options['composition']` /
`options['reactant']` (`flow_start.py:27-44`), which are plain dicts/strings
supplied by whoever instantiates it (e.g. `CEA_AIR_COMPOSITION` from
`constants.py`).

---

## 2. Bugs found in `SOFC.setup()`'s segment loop (`sofc.py:64-84`)

1. **Parenthesis closes too early.** In all three branches, the
   `SegmentSOFC(...)` call closes right after `N_seg=self.N_seg_tot)`, so
   `anode_composition=`, `cathode_composition=`, `spec=` (and
   `promotes_input=` in the third branch) are passed as kwargs to
   `self.add_subsystem(...)` instead of to `SegmentSOFC(...)`. `add_subsystem`
   doesn't accept those kwargs → `TypeError`.
2. **Wrong option name.** `N_seg=self.N_seg_tot` is passed, but
   `SegmentSOFC.initialize()` declares the option as `N_segments`
   (`sofc_segment.py:46`).
3. **Off-by-one gap.** The three branch conditions are all strict (`<`, `<`,
   `>`), so `i == N_passive - 1` and `i == N_passive + N_active - 1` satisfy
   none of them — those two segments never get built at all, breaking every
   downstream `self.connect(...)` that references them by name.
4. **`promotes_input` (should be `promotes_inputs`)**, and it's attached to
   the third (`passive`) branch — but `n_cell` only exists as a promoted
   variable in `SegmentSOFC` when `type='active'`
   (`sofc_segment.py:217,225`). Likely meant for the active branch.
5. **Dangling connect target**: `sofc.py:120`,
   `self.connect('Fl_I_an:tot:P', 'segment0.')` — truncated target name.
6. **Naming mismatches** spotted against `SegmentSOFC`'s real promoted names
   (`sofc_segment.py` docstring + `Convection`/energy-balance promotes):
   - `h_out_A`/`h_in_A` used in `sofc.py:91` should be `h_A_out`/`h_A_in`.
   - `P_A_in`/`P_C_in` used in `sofc.py:103-104` don't exist — `SegmentSOFC`
     only has `p_A_out`/`p_C_out` (no pressure drop assumed per segment, so
     there's no separate "inlet pressure" state to chain between segments).

---

## 3. Logic gap in how composition is passed to segments

**The bug:** `sofc.py` passes the *same* static dict to every segment:
```python
anode_composition   = self.Fl_I_data['Fl_I_an'],
cathode_composition = self.Fl_I_data['Fl_I_cat'],
```
for segment 0 **and** every downstream segment.

**Why that's wrong:** `SOFCThermoAdd.output_port_data()`
(`sofc_thermo_add.py:106-117`) adds an `'O'` element to the *anode* outlet set
that isn't necessarily in the inlet (O²⁻ arrives from the cathode and shows up
as H2O). With a pure-H2 feed, segment 0's anode inlet composition is
`{'H': 1.0}` (1 element) but its outlet becomes `{'H': ..., 'O': 0.0}`
(2 elements). That widened outlet *is* what's really connected forward at
runtime (`sofc.py:97`, `composition_out_A` → next segment's
`composition_in_A`). But segment 1 is still *constructed* with
`anode_composition = self.Fl_I_data['Fl_I_an']` — the original 1-element
dict — so its internal `SOFCThermoAdd.setup()` (`sofc_thermo_add.py:128,150`)
declares its `Fl_I:tot:composition` input with a fixed shape of 1
(`val=in_thermo.b0`, no `shape_by_conn`), while the actual thing connected
into it has shape 2. **Shape mismatch at `setup()` for every segment past the
first**, on the anode side. (The cathode side is unaffected — `cathode`
reactions never add an element to the set, so reusing the same dict there is
fine.)

**Why "just wire it via connections" doesn't fully fix this:** the real
`composition_out_A`/`composition_in_A` connections already carry the correct
*numeric* values forward — that part is fine. The problem is that
`SOFCThermoAdd.setup()` also needs to know, in Python, *which* elements are
at *which* index (to build its permutation matrix `self._map` and to size
`add_input`/`add_output` correctly, `sofc_thermo_add.py:134-145`). An
OpenMDAO connection only carries numbers, not "index 0 = H, index 1 = O" —
that identity/order metadata has to be propagated separately, in Python,
before `setup()` runs. This is exactly the same reason pyCycle's top-level
`Fl_O_data`/`Fl_I_data` system exists (§1) — the segment loop needs its own
mini version of the same mechanism.

**Proposed fix** — chain the composition dict forward with a plain Python
loop variable, advancing it with a throwaway `SOFCThermoAdd` instance (which
is safe to call `.output_port_data()` on immediately after construction,
*without* adding it as a subsystem or waiting for `.setup()` — OpenMDAO
options are available right after `__init__`):

```python
anode_comp   = self.Fl_I_data['Fl_I_an']
cathode_comp = self.Fl_I_data['Fl_I_cat']

for i, name in enumerate(self.segments):
    self.add_subsystem(name, SegmentSOFC(type=seg_type,
                                          N_segments=self.N_seg_tot,
                                          anode_composition=anode_comp,
                                          cathode_composition=cathode_comp,
                                          spec=thermo_data))

    anode_comp = SOFCThermoAdd(reaction_type='anode', spec=thermo_data,
                                inflow_composition=anode_comp).output_port_data()
    cathode_comp = SOFCThermoAdd(reaction_type='cathode', spec=thermo_data,
                                  inflow_composition=cathode_comp).output_port_data()
```

This mirrors the pattern already used in `SOFC.pyc_setup_output_ports()`
(`sofc.py:40-49`), which likewise builds a throwaway `SOFCThermoAdd` purely to
call `.output_port_data()` on it.

**Important distinction confirmed during review:** this fix only propagates
the *element set* (shape/identity), not real molar quantities — that's
correct and by design; `output_port_data()` never runs the actual mass
balance. Real, reaction-affected composition values already do reach
`Thermo`'s property calculation correctly, via an explicit connection that's
easy to miss:
```python
# sofc_segment.py:233-234
self.connect('composition_out_A', 'anode_thermo_out.composition')
self.connect('composition_out_C', 'cathode_thermo_out.composition')
```
`thermo_kwargs={'composition': anode_out_comp, ...}` passed to `Thermo(...)`
(`sofc_segment.py:91`) only sets the *default/shape* for `Thermo`'s
`composition` input; the connect above overrides it every solver iteration
with the real, dynamically computed value from `anode_rxn`
(`SOFCThermoAdd.compute()`, `sofc_thermo_add.py:173-218`, which does the real
Faraday's-law balance using `I` and `W_in`). So that part of the wiring does
not need to change.
