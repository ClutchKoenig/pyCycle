"""
Verdichter Syntax Demo (verdichter_syntax.py)

Kurzes Demo-Skript, das zeigt wie man `FlowStart` und `Compressor` gemeinsam
verwendet, so dass `ideal_ht`, `inlet_ht`, `ht_out` und `Fl_O:tot:h` berechnet
werden — auch wenn bei Instantiation keine Werte übergeben wurden.

Dieses Skript demonstriert:
 - dynamische Port- und Kompositionsinitialisierung
 - Datenfluss: FlowStart -> Compressor
 - Entropie-/Enthalpie-Rise Berechnung via `EnthalpyRise`
 - wie man Variablen per `IndepVarComp` setzt und mit `connect` verbindet

Aufrufen:
    python example_cycles/verdichter_syntax.py

"""

import openmdao.api as om
import pycycle.api as pyc

def build_and_run():
    prob = om.Problem()
    model = prob.model = om.Group()

    # Add a small FlowStart to provide composition and initial flow port
    fs = pyc.FlowStart(thermo_method='CEA', thermo_data=pyc.species_data.janaf)
    comp = pyc.Compressor(map_data=pyc.LPCMap, map_extrap=True, thermo_method='CEA', thermo_data=pyc.species_data.janaf)

    model.add_subsystem('fs', fs)
    model.add_subsystem('comp', comp)

    # We must call pyc_setup_output_ports on the FlowStart so the Flow composition
    # is available for the Compressor: normally Cycle.setup() would do this.
    model.fs.pyc_setup_output_ports()

    # copy the FlowStart's Fl_O composition into the Compressor's Fl_I_data
    # so the Compressor can create its Thermo objects in setup()
    model.comp.Fl_I_data['Fl_I'] = model.fs.Fl_O_data['Fl_O']

    # Connect the flow variables from FlowStart to Compressor to ensure the connection
    # metadata (shape/units) is consistent (we don't connect composition by string, as
    # we already copied the composition dict above to comp.Fl_I_data).
    model.connect('fs.Fl_O:tot:T', 'comp.Fl_I:tot:T')
    model.connect('fs.Fl_O:tot:P', 'comp.Fl_I:tot:P')
    model.connect('fs.Fl_O:stat:W', 'comp.Fl_I:stat:W')
    model.connect('fs.Fl_O:tot:h', 'comp.Fl_I:tot:h')

    # Add independent variables for the inlet total properties and engine inputs
    des = om.IndepVarComp()
    des.add_output('Fl_I:tot:T', 500.0, units='degR')
    des.add_output('Fl_I:tot:P', 14.7, units='psi')
    des.add_output('Fl_I:stat:W', 282.115, units='lbm/s')
    des.add_output('Nmech', 4666.1, units='rpm')
    des.add_output('PR', 2.0)
    des.add_output('eff', 0.924)

    model.add_subsystem('des', des)

    # Connect the independent variables to the FlowStart (these set the inlet flow values)
    # FlowStart expects 'T', 'P', 'W' promoted input names rather than 'Fl_I:...'
    model.connect('des.Fl_I:tot:T', 'fs.T')
    model.connect('des.Fl_I:tot:P', 'fs.P')
    model.connect('des.Fl_I:stat:W', 'fs.W')
    model.connect('des.Nmech', 'comp.Nmech')
    model.connect('des.PR', 'comp.PR')
    model.connect('des.eff', 'comp.eff')

    # Provide default composition / R so metadata & shapes are consistent across
    # promoted input names that reference composition and R.
    # Don't set a default for composition here; it's provided by FlowStart
    model.set_input_defaults('comp.Fl_I:tot:R', val=0.0686)

    prob.setup()
    prob.set_solver_print(level=-1)

    # Run model
    prob.run_model()

    # Print relevant intermediate and output values
    print('--- Verdichter Syntax Demo ---')
    try:
        print('Fl_I:tot:h      =', prob['comp.Fl_I:tot:h'])
    except Exception:
        print('Fl_I:tot:h not available')
    try:
        print('ideal_flow.h     =', prob['comp.ideal_flow.h'])
    except Exception:
        print('ideal_flow.h not available (mode mismatch)')
    try:
        print('enth_rise.ht_out =', prob['comp.enth_rise.ht_out'])
    except Exception:
        print('enth_rise.ht_out not available')
    try:
        print('Fl_O:tot:h       =', prob['comp.Fl_O:tot:h'])
        print('Fl_O:tot:P       =', prob['comp.Fl_O:tot:P'])
    except Exception:
        print('Fl_O:tot:* not available')

    # Derived map values
    try:
        print('Wc               =', prob['comp.Wc'])
        print('Nc               =', prob['comp.Nc'])
        print('PR (map)         =', prob['comp.PR'])
        print('eff              =', prob['comp.eff'])
    except Exception:
        pass

    return prob


if __name__ == '__main__':
    build_and_run()
