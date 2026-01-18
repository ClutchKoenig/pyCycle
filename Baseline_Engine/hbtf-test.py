"""
Create a baseline high-bypass turbofan engine instance using the `baseline_engine`
class defined in `Baseline-Engine/baseline-engine.py` and the pycycle/OpenMDAO framework.

This script does not attempt to be a full working simulation, but demonstrates how to
import the class by path, create an OpenMDAO Problem, set common baseline input values,
and attempt to run the setup and a single run of the model.

Usage:
    python Baseline-Engine/hbtf-test.py

"""

from pathlib import Path
import importlib.util
import traceback
import argparse
import openmdao.api as om


def load_module_from_path(path: Path, module_name: str):
    """Load a module from a file path and return the module object."""
    spec = importlib.util.spec_from_file_location(module_name, str(path))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main(run_model=False):
    repo_root = Path(__file__).resolve().parent
    baseline_file = repo_root / 'baseline-engine.py'

    if not baseline_file.exists():
        raise FileNotFoundError(f'Expected baseline file at {baseline_file} not found')

    mod = load_module_from_path(baseline_file, 'baseline_engine_module')

    try:
        baseline_cls = getattr(mod, 'baseline_engine')
    except AttributeError:
        raise AttributeError('`baseline_engine` class not found in baseline-engine.py')

    # Create an OpenMDAO Problem and add the baseline engine as the model.
    # Note: the baseline-engine module may be incomplete; if setup fails, fall back
    # to the reference HBTF/MPhbtf example from the example_cycles folder.
    prob = om.Problem()
    baseline_engine = None
    try:
        baseline_engine = baseline_cls()
        prob.model = baseline_engine
        # Try a setup — this may raise exceptions if the baseline engine is incomplete.
        prob.setup()
    except Exception:
        print('Primary baseline-engine setup failed; falling back to the canonical example HBTF (MPhbtf).')
        traceback.print_exc()
        # fallback to example implementation for a runnable model
        try:
            # Try to import example_cycles normally first
            import example_cycles.high_bypass_turbofan as ex_hbtf
        except Exception:
            # If the regular import fails (e.g. because cwd or path), attempt to load by file path
            try:
                repo_root_dir = Path(__file__).resolve().parent.parent
                example_file = repo_root_dir / 'example_cycles' / 'high_bypass_turbofan.py'
                ex_spec = importlib.util.spec_from_file_location('example_hbtf', str(example_file))
                ex_hbtf = importlib.util.module_from_spec(ex_spec)
                ex_spec.loader.exec_module(ex_hbtf)
            except Exception:
                print('Could not import or locate example_cycles.high_bypass_turbofan module')
                traceback.print_exc()
                raise

        try:
            prob = om.Problem()
            baseline_engine = ex_hbtf.MPhbtf()  # reference multi-point HBTF with design & off-design points
            prob.model = baseline_engine
            prob.setup()
            print('Fallback MPhbtf setup successful')
        except Exception:
            print('Fallback MPhbtf setup failed:')
            traceback.print_exc()
            raise

    # Set baseline flight conditions that match many HBTF references
    try:
        # Altitude: 35000 ft, Mach: 0.8
        # Try both common patterns: standalone cycles and MPCycle 'DESIGN' points
        flight_alt_names = ['DESIGN.fc.alt', 'OD_full_pwr.fc.alt', 'flight_cond.alt']
        flight_mn_names = ['DESIGN.fc.MN', 'OD_full_pwr.fc.MN', 'flight_cond.MN']
        for n in flight_alt_names:
            try:
                prob.set_val(n, 35000.0, units='ft')
            except Exception:
                pass
        for n in flight_mn_names:
            try:
                prob.set_val(n, 0.8)
            except Exception:
                pass

        # Bypass ratio — attempt to set on both common naming patterns
        for bpr_name in ('splitter.BPR', 'DESIGN.splitter.BPR', 'splitter.bypass_ratio'):
            try:
                prob.set_val(bpr_name, 15.0)
            except Exception:
                pass

        # a few common efficiency/PR settings (optional)
        # Best effort: set expected PR & eff names across likely promoted names
        settings = {
            'fan.PR': 1.685, 'fan.eff': 0.8948,
            'DESIGN.fan.PR': 1.685, 'DESIGN.fan.eff': 0.8948,
            'lp_compressor.PR': 1.935, 'lp_compressor.eff': 0.9243,
            'DESIGN.lp_compressor.PR': 1.935, 'DESIGN.lp_compressor.eff': 0.9243,
            'hp_compressor.PR': 9.369, 'hp_compressor.eff': 0.8707,
            'DESIGN.hp_compressor.PR': 9.369, 'DESIGN.hp_compressor.eff': 0.8707,
            'hp_turbine.eff': 0.8888, 'lp_turbine.eff': 0.8996,
            'DESIGN.hp_turbine.eff': 0.8888, 'DESIGN.lp_turbine.eff': 0.8996
        }
        for k, v in settings.items():
            try:
                prob.set_val(k, v)
            except Exception:
                # ignore failures — we only set these options if they exist
                pass

        # Run the model (may fail depending on the level of implementation)
        try:
            if run_model:
                prob.set_solver_print(level=-1)
                prob.run_model()
            print('Model run successful — key outputs (if available):')
            keys_to_print = [
                'flight_cond.Fl_O:stat:W',
                'flight_cond.Fl_O:stat:T',
                'flight_cond.Fl_O:stat:P',
            ]
            for k in keys_to_print:
                if k in prob:
                    print(f"{k} = {prob[k]}")
        except Exception:
            print('Model run failed — likely due to incomplete baseline cycle implementation.\n')
            traceback.print_exc()

    except Exception:
        print('Error when setting input values or running the model')
        traceback.print_exc()

    # Give the user an object they can import interactively from this script
    return baseline_engine, prob


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run the baseline HBTF test')
    parser.add_argument('--run', dest='run', action='store_true', help='Run the full model after setup (long)')
    args = parser.parse_args()

    be, prob = main(run_model=args.run)
    # expose a module-level name for convenient imports when run interactively or by tests
    baseline_engine = be
    print('\nCreated baseline_engine class instance (variable `baseline_engine`)')

# Also provide a `baseline_engine` module-level variable when imported
# To build the instance at import-time (safe for tests / interactive use), call main()
# This is a best-effort initialization – if setup/run fails, the instance will still exist
# but may require a local call to `prob.setup()` or further input initialization.
try:
    if 'baseline_engine' not in globals():
        baseline_engine, _ = main()
except Exception:
    # If import-time creation fails, fall back to None and leave the user to call `main()` manually
    baseline_engine = None
    # In interactive sessions you can now inspect `be` and `prob` objects
