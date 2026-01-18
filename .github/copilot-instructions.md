# copilot-instructions for pyCycle

This file contains concise, project-specific guidance for automated coding agents working on pyCycle. Keep the focus on actionable code patterns, critical files, and workflows that reduce friction for contributors and agents.

## Quick setup
- Install locally (recommended for development)
  ```bash
  python -m pip install -e .[all]
  ```
- Tests: unit tests are run with `testflo` (optional dependency)
  ```bash
  # from repo root
  testflo pycycle
  # run example scripts (long-running, use local running or CI for full coverage)
  python example_cycles/high_bypass_turbofan.py
  ```

## Big-picture architecture (what the codebase expects)
- Built on OpenMDAO: models extend `openmdao.api.Group` and use `Problem` instances (`prob = om.Problem()`), `prob.setup()` and `prob.run_model()`.
- `Cycle`/`MPCycle` (in `pycycle/mp_cycle.py`) are the custom Group types that orchestrate flow graphs and multi-point cycles.
- `Element` base class (in `pycycle/element_base.py`) must be used for any component that contains flow ports and defines `pyc_setup_output_ports()`.
- `pycycle/api.py` exports the public elements, maps, and helpers — prefer importing via `import pycycle.api as pyc` (examples follow this approach).

## Key patterns & conventions
- Flow ports: use `Fl_I` / `Fl_O` naming with `:tot` and `:stat` sub-keys: e.g. `fan.Fl_O:tot:T`, `inlet.Fl_O:stat:W`.
- Thermo: use `Thermo` and `ThermoAdd` objects in elements; `Thermo(mode='total_TP', fl_name='Fl_O:tot', ...)` is common.
  - When building a new `Element`, define `initialize()` to `options.declare(...)`, implement `pyc_setup_output_ports()` to populate `self.Fl_O_data` and then use `setup()` to add `Thermo`/`ThermoAdd` components and promote inputs/outputs.
- Graph-based flow propagation: prefer `self.pyc_connect_flow('fc.Fl_O', 'inlet.Fl_I')` in `Cycle` subclasses; `connect_flow` from `pycycle.connect_flow` is deprecated.
- Use `add_subsystem` and `promotes` to connect subsystems; `pyc_add_element` is deprecated (wraps add_subsystem).
- MPI/Distributed: not a default pattern here; typical development uses local runs and `OpenMDAO` problem-level API.

## Important files and directories
- `pycycle/mp_cycle.py` - Cycle and MPCycle orchestration (flow graph and setup).
- `pycycle/element_base.py` - `Element` base class; every component with a flow port should inherit it.
- `pycycle/thermo/` - thermo solvers (CEA & Tabular). Check `pycycle/thermo/cea` and `pycycle/thermo/tabular`.
- `pycycle/elements/` - elemental building blocks (compressor, turbine, combustor, nozzle, etc.)
- `pycycle/maps/` - compressor/turbine maps used in element constructors (pass in via map_data).
- `example_cycles/` - runnable engine cases. These are canonical examples and tests use them to validate integration.
- `example_cycles/tests/test_all_examples.py` - shows how the repo tests run example scripts (long-running examples are intentionally skipped by the test harness).

## Deprecated or special cases
- `pyc_add_element` and `connect_flow` are deprecated. Use `add_subsystem` and `pyc_connect_flow` / `Cycle.pyc_connect_flow` instead.
- The `Element` base class requires a `pyc_setup_output_ports` method. If you create a new element, define this method so `Cycle.setup()` can propagate flow port data down the graph.

## MPCycle & Multi-point patterns
- Use `pyc.MPCycle` for design/off-design workflows. Key helpers:
  - `pyc_add_pnt(name, pnt, **kwargs)` to add design/off-design points
  - `pyc_add_cycle_param(name, val, units=None)` to promote shared cycle params to points
  - `pyc_connect_des_od(src, target)` to connect values from design → off-design points
  - `pyc_use_default_des_od_conns()` to auto-propagate element default des→od connections

## Typical debugging & PR guidance
- For interactive debugging, instantiate a `Problem`, build a Cycle or MPCycle (example in `example_cycles/high_bypass_turbofan.py`) and call `prob.setup(); prob.run_model()`.
- To reproduce CI tests locally, use conda/virtualenv to install dependencies (OpenMDAO per README compatibility table); GitHub Actions uses Python 3.12 and specific NumPy / SciPy versions.
- CI test command:
  ```bash
  # runs pyCycle unit tests as CI does (produces deprecation & coverage reports)
  testflo -n 1 pycycle --timeout=240 --show_skipped --deprecations_report=deprecations.txt --coverage --coverpkg pycycle --durations=20
  ```

## Small examples & snippets to copy from
- Example `pyc_connect_flow` usage:
  ```py
  # in Cycle.setup()
  self.pyc_connect_flow('fc.Fl_O', 'inlet.Fl_I')
  self.pyc_connect_flow('inlet.Fl_O', 'fan.Fl_I')
  ```
- Add a subcomponent and promote input:
  ```py
  self.add_subsystem('fan', pyc.Compressor(map_data=pyc.FanMap), promotes_inputs=[('Nmech', 'LP_Nmech')])
  ```

## Test & performance notes
- Many example scripts are expensive to run; CI only runs unit-tests for core `pycycle` and short examples via `testflo`. Use `example_cycles` as integration tests locally if needed.
- The test harness `example_cycles/tests/test_all_examples.py` instantiates separate processes for each example; this makes debugging easier but increases runtime.

If something here is ambiguous or you want more detail for a specific file or pattern (such as adding a new Element that uses Tabular thermo), tell me and I will expand this file with concrete code samples.

If anything here is unclear or you want more details about a specific file, example, or workflow (e.g. how to add a new Element that uses tabular thermo), tell me which area to expand. 
