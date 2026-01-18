"""
Reference HBTF (High Bypass Ratio Turbofan) model: geared + unmixed exhaust

This file provides a compact reference-class `HBTFRef` implemented as a pyc.Cycle
with the reference parameters from the High BPR turbofan case (OPR=50, BPR=15,
T4=1700 K, FN=27 kN). It also supplies a `MPhbtfRef` MPCycle that wraps DESIGN
and off-design points.

The implementation is intentionally similar to `example_cycles/high_bypass_turbofan.py`
but focuses on the reference configuration and sets the main defaults described
in section 8.2 (High Bypass Ratio Turbofan).
"""

import sys
import argparse
import numpy as np
import openmdao.api as om
import pycycle.api as pyc


class HBTFRef(pyc.Cycle):
    """High-Bypass Turbofan reference cycle (geared, unmixed exhaust).

    This cycle defines the component layout (flight cond, inlet, fan, LPC/HPC,
    combustor, turbines, separate nozzles) and a gearbox between fan and LP
    shafts to represent a geared configuration.
    """

    def initialize(self):
        self.options.declare('throttle_mode', default='T4', values=['T4', 'percent_thrust'])
        super().initialize()

    def setup(self):
        design = self.options['design']

        # Choose thermo method and data
        self.options['thermo_method'] = 'CEA'
        self.options['thermo_data'] = pyc.species_data.janaf

        FUEL_TYPE = 'Jet-A(g)'

        # Add the typical components used in HBTF
        self.add_subsystem('fc', pyc.FlightConditions())
        self.add_subsystem('inlet', pyc.Inlet())
        self.add_subsystem('fan', pyc.Compressor(map_data=pyc.FanMap, map_extrap=True),
                           promotes_inputs=[('Nmech', 'LP_Nmech')])
        self.add_subsystem('splitter', pyc.Splitter())
        self.add_subsystem('duct4', pyc.Duct())
        # LPC implemented as Compressor using LPCMap
        self.add_subsystem('lpc', pyc.Compressor(map_data=pyc.LPCMap, map_extrap=True),
                           promotes_inputs=[('Nmech', 'LP_Nmech')])
        self.add_subsystem('duct6', pyc.Duct())
        self.add_subsystem('hpc', pyc.Compressor(map_data=pyc.HPCMap, map_extrap=True),
                           promotes_inputs=[('Nmech', 'HP_Nmech')])
        self.add_subsystem('bld3', pyc.BleedOut(bleed_names=['cool3','cool4']))
        self.add_subsystem('burner', pyc.Combustor(fuel_type=FUEL_TYPE))
        self.add_subsystem('hpt', pyc.Turbine(map_data=pyc.HPTMap, map_extrap=True),
                           promotes_inputs=[('Nmech', 'HP_Nmech')])
        self.add_subsystem('duct11', pyc.Duct())
        self.add_subsystem('lpt', pyc.Turbine(map_data=pyc.LPTMap, map_extrap=True),
                           promotes_inputs=[('Nmech', 'LP_Nmech')])
        self.add_subsystem('duct13', pyc.Duct())
        # Unmixed (separate core and bypass nozzles)
        self.add_subsystem('core_nozz', pyc.Nozzle(nozzType='CV', lossCoef='Cv'))
        self.add_subsystem('byp_bld', pyc.BleedOut(bleed_names=['bypBld']))
        self.add_subsystem('duct15', pyc.Duct())
        self.add_subsystem('byp_nozz', pyc.Nozzle(nozzType='CV', lossCoef='Cv'))

        # Gearbox and shafts: geared configuration
        self.add_subsystem('lp_shaft', pyc.Shaft(num_ports=3), promotes_inputs=[('Nmech', 'LP_Nmech')])
        self.add_subsystem('fan_shaft', pyc.Shaft(num_ports=2), promotes_inputs=[('Nmech', 'Fan_Nmech')])
        # Geared: link fan shaft to lp_shaft
        self.add_subsystem('gearbox', pyc.Gearbox(design=True), promotes_inputs=[('N_in', 'Fan_Nmech'), ('N_out', 'LP_Nmech')])

        self.add_subsystem('hp_shaft', pyc.Shaft(num_ports=2), promotes_inputs=[('Nmech', 'HP_Nmech')])
        self.add_subsystem('perf', pyc.Performance(num_nozzles=2, num_burners=1))

        # Connections (flow graph): similar to canonical HBTF
        self.pyc_connect_flow('fc.Fl_O', 'inlet.Fl_I')
        self.pyc_connect_flow('inlet.Fl_O', 'fan.Fl_I')
        self.pyc_connect_flow('fan.Fl_O', 'splitter.Fl_I')
        self.pyc_connect_flow('splitter.Fl_O1', 'duct4.Fl_I')
        self.pyc_connect_flow('duct4.Fl_O', 'lpc.Fl_I')
        self.pyc_connect_flow('lpc.Fl_O', 'duct6.Fl_I')
        self.pyc_connect_flow('duct6.Fl_O', 'hpc.Fl_I')
        self.pyc_connect_flow('hpc.Fl_O', 'bld3.Fl_I')
        self.pyc_connect_flow('bld3.Fl_O', 'burner.Fl_I')
        self.pyc_connect_flow('burner.Fl_O', 'hpt.Fl_I')
        self.pyc_connect_flow('hpt.Fl_O', 'duct11.Fl_I')
        self.pyc_connect_flow('duct11.Fl_O', 'lpt.Fl_I')
        self.pyc_connect_flow('lpt.Fl_O', 'duct13.Fl_I')
        self.pyc_connect_flow('duct13.Fl_O', 'core_nozz.Fl_I')
        self.pyc_connect_flow('splitter.Fl_O2', 'byp_bld.Fl_I')
        self.pyc_connect_flow('byp_bld.Fl_O', 'duct15.Fl_I')
        self.pyc_connect_flow('duct15.Fl_O', 'byp_nozz.Fl_I')

        # Torque connections
        self.connect('fan.trq', 'fan_shaft.trq_0')
        self.connect('lpc.trq', 'lp_shaft.trq_1')
        self.connect('lpt.trq', 'lp_shaft.trq_2')
        self.connect('hpc.trq', 'hp_shaft.trq_0')
        self.connect('hpt.trq', 'hp_shaft.trq_1')

        # Perf group inputs
        self.connect('inlet.Fl_O:tot:P', 'perf.Pt2')
        self.connect('hpc.Fl_O:tot:P', 'perf.Pt3')
        self.connect('burner.Wfuel', 'perf.Wfuel_0')
        self.connect('inlet.F_ram', 'perf.ram_drag')
        self.connect('core_nozz.Fg', 'perf.Fg_0')
        self.connect('byp_nozz.Fg', 'perf.Fg_1')

        # Simple balances and solver settings follow the canonical example approach
        balance = self.add_subsystem('balance', om.BalanceComp())
        if design:
            balance.add_balance('W', units='lbm/s', eq_units='lbf')
            self.connect('balance.W', 'fc.W')
            self.connect('perf.Fn', 'balance.lhs:W')
            self.promotes('balance', inputs=[('rhs:W', 'Fn_DES')])

            balance.add_balance('FAR', eq_units='degR', lower=1e-4, val=.017)
            self.connect('balance.FAR', 'burner.Fl_I:FAR')
            self.connect('burner.Fl_O:tot:T', 'balance.lhs:FAR')
            self.promotes('balance', inputs=[('rhs:FAR', 'T4_MAX')])

            balance.add_balance('lpt_PR', val=1.5, lower=1.001, upper=8,
                                eq_units='hp', use_mult=True, mult_val=-1)
            self.connect('balance.lpt_PR', 'lpt.PR')
            self.connect('lp_shaft.pwr_in_real', 'balance.lhs:lpt_PR')
            self.connect('lp_shaft.pwr_out_real', 'balance.rhs:lpt_PR')

            balance.add_balance('hpt_PR', val=1.5, lower=1.001, upper=8,
                                eq_units='hp', use_mult=True, mult_val=-1)
            self.connect('balance.hpt_PR', 'hpt.PR')
            self.connect('hp_shaft.pwr_in_real', 'balance.lhs:hpt_PR')
            self.connect('hp_shaft.pwr_out_real', 'balance.rhs:hpt_PR')
        else:
            # Off-design balances similar to example
            balance.add_balance('W', units='lbm/s', lower=10., upper=1000., eq_units='inch**2')
            self.connect('balance.W', 'fc.W')
            self.connect('core_nozz.Throat:stat:area', 'balance.lhs:W')

            balance.add_balance('BPR', lower=2., upper=10., eq_units='inch**2')
            self.connect('balance.BPR', 'splitter.BPR')
            self.connect('byp_nozz.Throat:stat:area', 'balance.lhs:BPR')

        # Solver settings
        newton = self.nonlinear_solver = om.NewtonSolver()
        newton.options['atol'] = 1e-8
        newton.options['rtol'] = 1e-99
        newton.options['iprint'] = 2
        newton.options['maxiter'] = 50
        newton.options['solve_subsystems'] = True
        ls = newton.linesearch = om.ArmijoGoldsteinLS()
        ls.options['maxiter'] = 3
        ls.options['rho'] = 0.75
        self.linear_solver = om.DirectSolver()

        super().setup()


class MPhbtfRef(pyc.MPCycle):

    def setup(self):
        # Add design point using our HBTFRef
        self.pyc_add_pnt('DESIGN', HBTFRef(thermo_method='CEA'))

        # Set baseline flight condition
        self.set_input_defaults('DESIGN.fc.MN', 0.8)
        self.set_input_defaults('DESIGN.fc.alt', 35000., units='ft')

        # High level reference parameters (copied / set per 8.2)
        # OPR, BPR, T4, FN
        OPR_target = 50.0
        BPR_target = 15.0
        T4_K = 1700.0
        Fn_N = 27000.0

        # We'll set HPC PR to 15 (as specified) and solve for fan.PR given a reasonable LPC PR
        hpc_PR = 15.0
        assumed_lpc_PR = 1.935
        # fan.PR so that product fan*lpc*hpc == OPR_target
        fan_PR = OPR_target / (hpc_PR * assumed_lpc_PR)

        # set default component PRs and efficiency guesses
        self.set_input_defaults('DESIGN.hpc.PR', hpc_PR)
        self.set_input_defaults('DESIGN.lpc.PR', assumed_lpc_PR)
        self.set_input_defaults('DESIGN.fan.PR', fan_PR)

        # Set BPR and T4/Fn
        self.set_input_defaults('DESIGN.splitter.BPR', BPR_target)
        # T4: specify Kelvin units
        self.set_input_defaults('DESIGN.T4_MAX', T4_K, units='K')
        # Fn: specify N
        self.set_input_defaults('DESIGN.Fn_DES', Fn_N, units='N')

        # reasonable efficiency initial guesses
        self.set_input_defaults('DESIGN.fan.eff', 0.895)
        self.set_input_defaults('DESIGN.lpc.eff', 0.924)
        self.set_input_defaults('DESIGN.hpc.eff', 0.87)
        self.set_input_defaults('DESIGN.hpt.eff', 0.889)
        self.set_input_defaults('DESIGN.lpt.eff', 0.90)

        # Gearbox preference: choose a gear ratio (fan N relative to LP_Nmech)
        self.set_input_defaults('DESIGN.gearbox.gear_ratio', 0.3225)

        # Off-design points
        self.pyc_add_pnt('OD_full_pwr', HBTFRef(design=False, thermo_method='CEA', throttle_mode='T4'))
        self.pyc_add_pnt('OD_part_pwr', HBTFRef(design=False, thermo_method='CEA', throttle_mode='percent_thrust'))

        self.set_input_defaults('OD_full_pwr.fc.MN', 0.8)
        self.set_input_defaults('OD_full_pwr.fc.alt', 35000., units='ft')
        self.set_input_defaults('OD_part_pwr.fc.MN', 0.8)
        self.set_input_defaults('OD_part_pwr.fc.alt', 35000., units='ft')

        # Link design->offdesign where applicable
        self.pyc_use_default_des_od_conns()

        super().setup()


def main(run=False):
    prob = om.Problem()
    prob.model = mp = MPhbtfRef()
    prob.setup()

    # Example: run the model if requested (runs may be long)
    if run:
        prob.set_solver_print(level=2)
        prob.run_model()
        # Print a few summary outputs
        print('DESIGN.OPR:', prob['DESIGN.perf.OPR'][0])
        print('DESIGN.splitter.BPR:', prob['DESIGN.splitter.BPR'][0])
        print('DESIGN.perf.Fn (N):', prob['DESIGN.perf.Fn'][0] * 4.44822)

    return mp, prob


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--run', action='store_true', help='Run the model after setup (can be slow)')
    args = parser.parse_args()
    mp, prob = main(run=args.run)
