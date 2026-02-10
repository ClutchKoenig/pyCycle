import sys

import numpy as np
import openmdao.api as om

import pycycle.api as pyc


class baseline_engine(pyc.Cycle):
    """
    Design-only baseline engine model (Kurzke 2025, S. 303)
    No Off-Design points — just a single Cycle for the design point.
    """
    def setup(self):
        design = self.options['design']

        USE_TABULAR = False

        self.options['thermo_method'] ='CEA'
        self.options['thermo_data']=pyc.species_data.janaf
        FUEL_TYPE = 'Jet-A(g)'

        # Subsysteme:
        self.add_subsystem('flight_cond', pyc.FlightConditions())
        self.add_subsystem('inlet', pyc.Inlet())
        # Möglicherweise muss man Fan in einen inner und outer Fan teilen, 
        # damit Vergleich zu Kurzke möglichst Akkurat ist
        # Kurzke gibt inner and outer Pressure Ratios und Effizienz an
        self.add_subsystem('fan', pyc.Compressor(map_data=pyc.FanMap, bleed_names=[], map_extrap=True),
                           promotes_inputs=[('Nmech', 'LP_Nmech')])
        self.add_subsystem('splitter', pyc.Splitter())
        self.add_subsystem('duct_lpc_inlet', pyc.Duct())
        # bypass doesnt use nozzle
        self.add_subsystem('duct_bp', pyc.Duct())
        self.add_subsystem('lp_compressor', pyc.Compressor(map_data = pyc.LPCMap, map_extrap=True), 
                           promotes_inputs=[('Nmech', 'LP_Nmech')])
        

        self.add_subsystem('duct_hpc_inlet', pyc.Duct())
        self.add_subsystem('hp_compressor', pyc.Compressor(map_data=pyc.HPCMap, 
                                                           map_extrap=True),
                                                           promotes_inputs = [('Nmech', 'HP_Nmech')])
        self.add_subsystem('bleed_hpc_exit', pyc.BleedOut(bleed_names=['bleed_hpt_vanes_cool', 'bleed_hpt_blades_cool', 'bleed_lpt_vanes_cool', 'bleed_lpt_blades_cool']))
        # p. 305 fig. 8.2 no bleed in this config - 02.12.25: quatsch! 
        # Bleed wird am Austritt des HPC abgenommen siehe p.306 Table

        self.add_subsystem('burner', pyc.Combustor(fuel_type=FUEL_TYPE))

        self.add_subsystem('hp_turbine', pyc.Turbine(map_data=pyc.HPTMap, map_extrap = True, 
                           bleed_names=['bleed_hpt_vanes_cool', 'bleed_hpt_blades_cool']),
                           promotes_inputs=[('Nmech', 'HP_Nmech')])
        self.add_subsystem('duct_lpt_inlet', pyc.Duct())
        self.add_subsystem('lp_turbine', pyc.Turbine(map_data = pyc.LPTMap, map_extrap = True,
                           bleed_names=['bleed_lpt_vanes_cool','bleed_lpt_blades_cool']),
                           promotes_inputs=[('Nmech','LP_Nmech')])
        # no nozzle behind last turbine stage, 
        # warum?
        self.add_subsystem('duct_core_outlet', pyc.Duct())
        self.add_subsystem('core_nozzle', pyc.Nozzle(nozzType='CV', lossCoef='Cv'))
        self.add_subsystem('bypass_nozzle', pyc.Nozzle(nozzType='CV',lossCoef='Cv'))

        #self.add_subsystem('duct')

        #Create shafts (unfinsihed)
        self.add_subsystem('lp_shaft', pyc.Shaft(num_ports=3))
        self.add_subsystem('fan_gearbox', pyc.Gearbox(),
                           promotes_inputs = [('N_in','LP_Nmech'),('N_out','Fan_Nmech')]) # add details
        
        # No fan Shaft, fan is coupled with LP-Shaft
        # self.add_subsystem('fan_shaft', pyc.Shaft(num_ports=2))
        self.add_subsystem('fan_shaft', pyc.Shaft(num_ports=2))
        self.add_subsystem('hp_shaft', pyc.Shaft(num_ports=2))
        self.add_subsystem('performance', pyc.Performance(num_nozzles=2, num_burners=1))
        
        # am 03.12 weiter mit verbindungen
        self.connect('inlet.Fl_O:tot:P', 'performance.Pt2')
        self.connect('hp_compressor.Fl_O:tot:P', 'performance.Pt3')
        self.connect('burner.Wfuel', 'performance.Wfuel_0')
        self.connect('inlet.F_ram', 'performance.ram_drag')
        self.connect('core_nozzle.Fg', 'performance.Fg_0')
        self.connect('bypass_nozzle.Fg', 'performance.Fg_1')

        #connect physische Verbindung
        self.connect('fan_gearbox.trq_out', 'fan_shaft.trq_1')
        self.connect('fan_gearbox.trq_in', 'lp_shaft.trq_0')

        self.connect('fan.trq','fan_shaft.trq_0')

        self.connect('lp_compressor.trq','lp_shaft.trq_1')
        self.connect('lp_turbine.trq','lp_shaft.trq_2')
        
        self.connect('hp_compressor.trq','hp_shaft.trq_0')
        self.connect('hp_turbine.trq','hp_shaft.trq_1')
        
        # Definition BPR Komponente
        self.add_subsystem('opr_comp',om.ExecComp('OPR=F_PR * LPC_PR * HPC_PR',
                           F_PR={'val':1.2, 'units': None},
                           LPC_PR={'val':3, 'units':None},
                           HPC_PR={'val':15, 'units':None},
                           OPR = {'val':50, 'units':None}))

        self.connect('opr_comp.F_PR','fan.PR')
        self.connect('opr_comp.LPC_PR','lp_turbine.PR')
        self.connect('opr_comp.hPC_PR', 'hp_turbine.PR')


        self.add_subsystem('ideal_jet_velocity_ratio', om.ExecComp('vr_id = v18 / v8', 
                                                                   v18 = {'val':300, 'units': 'ft/s'},
                                                                   v8={'val':50, 'units': 'ft/s'},
                                                                   vr_id={'val':0.8, 'units':None}))
        self.connect('ideal_jet_velocity_ratio.v18', 'bypass_nozzle.Fl_O:stat:V')
        self.connect('ideal_jet_velocity_ratio.v8','core_nozzle.Fl_O:stat:V')
        
        balance = self.add_subsystem('balance', om.BalanceComp())
        if design:
            self.add_subsystem('geometry', om.ExecComp(core_nozz_exit={'val':80, 'units':'inch**2'},
                                                       byp_nozzle_exit={'val':200, 'units':'inch**2'}))
            # Balancing of the 
            balance.add_balance('W', units = 'lbm/s', eq_units = 'inch**2')
            self.connect('balance.W','flight_cond.W')
            self.connect('core_nozzle.Fl_O:stat:area','balance.W:lhs')
            self.connect('balance', inputs=['rhs:W','geometry.core_nozz_exit'])
            
            self.add_balance('BPR', units = None, eq_units = 'inch**2')
            self.connect('balance.BPR','splitter.BPR')
            self.connect('bypass_nozzle.Fl_O:stat:area','balance.BPR:lhs')
            self.connect('balance', inputs=['rhs:BPR', 'geometry.byp_nozzle_exit'])

            self.add_balance('OPR', units = None, eq_units = None)
            self.connect('balance.OPR', 'ideal_jet_velocity_ratio.vr_id')
            self.connect('bypass_nozzle')

            # ================================================= 
            # ======== Balances for energy conservation ========
            balance.add_balance('lpt_PR', eq_units='hp', rhs_val=0., res_ref=1e4)
            self.connect('balance.lpt_PR', 'lp_turbine.PR')
            self.connect('lp_shaft.pwr_net', 'balance.lhs:lpt_PR')

            balance.add_balance('hpt_PR', eq_units='hp', rhs_val = 0., res_ref=1e4)
            self.connect('balance.hpt_PR','hp_turbine.PR')
            self.connect('hp_shaft.pwr_net', 'balance.lhs:hpt_PR')
            # ==================================================

        else: 
            # Design-Punkt Balances
            # Thrust requirement = 27kN
            balance.add_balance('W', units = 'lbm/s' , eq_units = 'lbf')
            self.connect('balance.W', 'flight_cond.W')
            self.connect('performance.Fn', 'balance.lhs:W')
            self.promotes('balance', inputs=[('rhs:W', 'Fn_REQUIREMENT')])

            # Turbine entry temperature rq =1700K
            balance.add_balance('FAR', eq_units = 'degR')
            self.connect('balance.FAR', 'burner.Fl_I:FAR')
            self.connect('burner.Fl_O:tot:T', 'balance.lhs:FAR')
            self.promotes('balance', inputs=[('rhs:FAR', 'T4_REQUIREMENT')])
        
            # ================================================= 
            # ======== Balances for energy conservation ========
            balance.add_balance('lpt_PR', eq_units='hp', rhs_val=0., res_ref=1e4)
            self.connect('balance.lpt_PR', 'lp_turbine.PR')
            self.connect('lp_shaft.pwr_net', 'balance.lhs:lpt_PR')

            balance.add_balance('hpt_PR', eq_units='hp', rhs_val = 0., res_ref=1e4)
            self.connect('balance.hpt_PR','hp_turbine.PR')
            self.connect('hp_shaft.pwr_net', 'balance.lhs:hpt_PR')
            # ==================================================

        balance.add_balance('lpc_PR', val=3, units=None, eq_units=None)
        self.connect('balance.lpc_PR', ['opr_comp.LPC_PR','lp_compressor.PR'])
        self.connect('opr_comp.OPR','balance.lhs:lpc_PR')
        self.promotes('balance', inputs=[('rhs:lpc_PR', 'OPR_REQUIREMENT')])
        


        # ======= Inlet Flow =========
        self.pyc_connect_flow('flight_cond.Fl_O','inlet.Fl_I')
        self.pyc_connect_flow('inlet.Fl_O','fan.Fl_I')
        self.pyc_connect_flow('fan.Fl_O','splitter.Fl_I')
        
        # ======= Core Flow =========
        self.pyc_connect_flow('splitter.Fl_O1', 'duct_lpc_inlet.Fl_I')
        self.pyc_connect_flow('duct_lpc_inlet.Fl_O', 'lp_compressor.Fl_I')
        self.pyc_connect_flow('lp_compressor.Fl_O', 'duct_hpc_inlet.Fl_I')
        self.pyc_connect_flow('duct_hpc_inlet.Fl_O', 'hp_compressor.Fl_I')
        self.pyc_connect_flow('hp_compressor.Fl_O', 'bleed_hpc_exit.Fl_I')
        self.pyc_connect_flow('bleed_hpc_exit.Fl_O', 'burner.Fl_I')
        self.pyc_connect_flow('burner.Fl_O', 'hp_turbine.Fl_I')
        self.pyc_connect_flow('hp_turbine.Fl_O', 'duct_lpt_inlet.Fl_I')
        self.pyc_connect_flow('duct_lpt_inlet.Fl_O','lp_turbine.Fl_I')
        self.pyc_connect_flow('lp_turbine.Fl_O', 'duct_core_outlet.Fl_I')
        self.pyc_connect_flow('duct_core_outlet.Fl_O', 'core_nozzle.Fl_I')
        
        # ==== Bleed Flows ====
        self.pyc_connect_flow('bleed_hpc_exit.bleed_hpt_vanes_cool', 'hp_turbine.bleed_hpt_vanes_cool', connect_stat=False)
        self.pyc_connect_flow('bleed_hpc_exit.bleed_hpt_blades_cool', 'hp_turbine.bleed_hpt_blades_cool', connect_stat=False)

        self.pyc_connect_flow('bleed_hpc_exit.bleed_lpt_vanes_cool', 'lp_turbine.bleed_lpt_vanes_cool', connect_stat=False)
        self.pyc_connect_flow('bleed_hpc_exit.bleed_lpt_blades_cool', 'lp_turbine.bleed_lpt_blades_cool', connect_stat=False)
        
        # ======= Bypass Flow =========
        self.pyc_connect_flow('splitter.Fl_O2', 'duct_bp.Fl_I')
        self.pyc_connect_flow('duct_bp.Fl_O', 'bypass_nozzle.Fl_I')

        # ======= vollst. Expansion durch Schließbedingung =========
        #self.connect('fc.Fl_O:stat:p', 'core_nozzle.Fl_O:stat:P')
        self.connect('flight_cond.Fl_O:stat:P', 'core_nozzle.Ps_exhaust')
        self.connect('flight_cond.Fl_O:stat:P', 'bypass_nozzle.Ps_exhaust')


        # ==================================================
        # ============= Solver Settings ====================
        newton = self.nonlinear_solver = om.NewtonSolver()
        newton.options['atol'] = 1e-8
        newton.options['rtol'] = 1e-99
        newton.options['iprint'] = 2
        newton.options['maxiter'] = 50
        newton.options['solve_subsystems'] = True
        newton.options['max_sub_solves'] = 1000
        newton.options['reraise_child_analysiserror'] = False
        ls = newton.linesearch = om.ArmijoGoldsteinLS()
        ls.options['maxiter'] = 3
        ls.options['rho'] = 0.75

        self.linear_solver = om.DirectSolver()

        super().setup()
def viewer(prob, pt, file=sys.stdout):
    """
    print a report of all the relevant cycle properties
    """

    if pt == 'DESIGN':
        MN = prob['DESIGN.fc.Fl_O:stat:MN']
        LPT_PR = prob['DESIGN.balance.lpt_PR']
        HPT_PR = prob['DESIGN.balance.hpt_PR']
        FAR = prob['DESIGN.balance.FAR']
    else:
        MN = prob[pt+'.fc.Fl_O:stat:MN']
        LPT_PR = prob[pt+'.lpt.PR']
        HPT_PR = prob[pt+'.hpt.PR']
        FAR = prob[pt+'.balance.FAR']

    summary_data = (MN[0], prob[pt+'.fc.alt'][0], prob[pt+'.inlet.Fl_O:stat:W'][0], prob[pt+'.perf.Fn'][0],
                        prob[pt+'.perf.Fg'][0], prob[pt+'.inlet.F_ram'][0], prob[pt+'.perf.OPR'][0],
                        prob[pt+'.perf.TSFC'][0], prob[pt+'.splitter.BPR'][0])

    print(file=file, flush=True)
    print(file=file, flush=True)
    print(file=file, flush=True)
    print("----------------------------------------------------------------------------", file=file, flush=True)
    print("                              POINT:", pt, file=file, flush=True)
    print("----------------------------------------------------------------------------", file=file, flush=True)
    print("                       PERFORMANCE CHARACTERISTICS", file=file, flush=True)
    print("    Mach      Alt       W      Fn      Fg    Fram     OPR     TSFC      BPR ", file=file, flush=True)
    print(" %7.5f  %7.1f %7.3f %7.1f %7.1f %7.1f %7.3f  %7.5f  %7.3f" %summary_data, file=file, flush=True)


    fs_names = ['fc.Fl_O', 'inlet.Fl_O', 'fan.Fl_O', 'splitter.Fl_O1', 'splitter.Fl_O2',
                'duct4.Fl_O', 'lpc.Fl_O', 'duct6.Fl_O', 'hpc.Fl_O', 'bld3.Fl_O', 'burner.Fl_O',
                'hpt.Fl_O', 'duct11.Fl_O', 'lpt.Fl_O', 'duct13.Fl_O', 'core_nozz.Fl_O', 'byp_bld.Fl_O',
                'duct15.Fl_O', 'byp_nozz.Fl_O']
    fs_full_names = [f'{pt}.{fs}' for fs in fs_names]
    pyc.print_flow_station(prob, fs_full_names, file=file)

    comp_names = ['fan', 'lpc', 'hpc']
    comp_full_names = [f'{pt}.{c}' for c in comp_names]
    pyc.print_compressor(prob, comp_full_names, file=file)

    pyc.print_burner(prob, [f'{pt}.burner'], file=file)

    turb_names = ['hpt', 'lpt']
    turb_full_names = [f'{pt}.{t}' for t in turb_names]
    pyc.print_turbine(prob, turb_full_names, file=file)

    noz_names = ['core_nozz', 'byp_nozz']
    noz_full_names = [f'{pt}.{n}' for n in noz_names]
    pyc.print_nozzle(prob, noz_full_names, file=file)

    shaft_names = ['hp_shaft', 'lp_shaft']
    shaft_full_names = [f'{pt}.{s}' for s in shaft_names]
    pyc.print_shaft(prob, shaft_full_names, file=file)

    bleed_names = ['hpc', 'bld3', 'byp_bld']
    bleed_full_names = [f'{pt}.{b}' for b in bleed_names]
    pyc.print_bleed(prob, bleed_full_names, file=file)


class MPbaseline_engine(pyc.MPCycle):
    def setup(self):
        self.pyc_add_pnt('DESIGN', baseline_engine(thermo_method='CEA'))
        # self.pyc_add_cycle_param('core_nozzle.Cv',0.93881)
        # Set Input Defaults here
        self.od_pts = ['OD_Thrust']

        self.pyc_add_pnt('OD_Thrust', baseline_engine(design=False, thermo_method='CEA'))
        self.set_input_defaults('OD_Thrust.flight_cond.alt', 35000, units='ft')
        self.set_input_defaults('OD_Thrust.flight_cond.MN',0.8)

        self.pyc_use_default_des_od_conns()
        self.pyc_connect_des_od('core_nozzle.Throat:stat:area','balance.rhs:W')
        self.pyc_connect_des_od('bypass_nozzle.Throat:stat:area','balance.rhs:BPR')
# =============================================================
# Berechne MA aus Kurzke Tabelle 
# Vergleiche Ergebnisse von PC mit Kurzke
# Modelliere Turboprop mit Daten von Pascal
# Bis 06.01.2026 nach GPT Labor 

def main():
    """
    Design-only runner: single point calculation for baseline engine.
    No MPCycle, no Off-Design points.
    """
    # Hier sind noch die Machzahlen zur Dimensionierung des Triebwerks
    # zu implementieren
    prob = om.Problem()
    prob.model = baseline_engine()
    
    prob.setup()
    
    # ========== Set Design Point Parameters ==========
    # Flight conditions
    prob.set_val('flight_cond.alt', 35000., units='ft')
    prob.set_val('flight_cond.MN', 0.8)
    
    # Compressor/Turbine initial guesses
    prob.set_val('DESIGN.fan.PR', 1.37)     # Weiß nicht ob Sinnvoll
    prob.set_val('fan.eff', 0.9)
    prob.set_val('lp_compressor.PR', 2.586)
    prob.set_val('lp_compressor.eff', 0.88)
    prob.set_val('hp_compressor.PR', 15)    # Requirement
    prob.set_val('hp_compressor.eff', 0.85)
    prob.set_val('hp_turbine.eff', 0.91)
    prob.set_val('lp_turbine.eff', 0.92)
    prob.set_val('hp_turbine.PR', 4.605)
    prob.set_val('lp_turbine.PR', 10.168)

    # Balance RHS (targets)
    prob.set_val('Fn_REQUIREMENT', 27 * 224.80894387096 , units='lbf') # kN zu lbf ist * 224.80894387096
    prob.set_val('T4_REQUIREMENT', 1700 * 9/5 , units='degR') # K zu Rankine ist * 9/5
    prob.set_val('OPR_REQUIREMENT', 50.0)  # Overall Pressure Ratio target
    
    # Initial guesses for balance states
    prob['balance.W'] = 100.0
    prob['balance.FAR'] = 0.025
    prob['balance.lpc_PR'] = 3.0
    prob['balance.lpt_PR'] = 4.0
    prob['balance.hpt_PR'] = 3.0
    
    # ========== Run the Model ==========
    prob.set_solver_print(level=2, depth=1)
    prob.run_model()
    
    # ========== Print Results ==========
    print("\n" + "="*70)
    print("BASELINE ENGINE DESIGN POINT RESULTS")
    print("6069.8 lbf="*70)
    print(f"Inlet mass flow rate (W):     {prob['balance.W'][0]:8.2f} lbm/s")
    print(f"Fuel-air ratio (FAR):         {prob['balance.FAR'][0]:.6f}")
    print(f"LPC pressure ratio:           {prob['balance.lpc_PR'][0]:.3f}")
    print(f"LPT pressure ratio:           {prob['balance.lpt_PR'][0]:.3f}")
    print(f"HPT pressure ratio:           {prob['balance.hpt_PR'][0]:.3f}")
    print(f"\nNet thrust (Fn):              {prob['performance.Fn'][0]:8.1f} lbf")
    print(f"Overall Pressure Ratio (OPR): {prob['opr_comp.OPR'][0]:8.2f}")
    print(f"Burner exit temp (T4):        {prob['burner.Fl_O:tot:T'][0]:8.1f} degR")
    print(f"\nLP shaft net power:           {prob['lp_shaft.pwr_net'][0]:10.1f} hp")
    print(f"HP shaft net power:           {prob['hp_shaft.pwr_net'][0]:10.1f} hp")
    print("="*70)
    
    return prob

if __name__ == '__main__':
    main()



"""
Stat. 2-13 = Outer LPC (Kurzke) Bypassstrom Fan
Stat. 02-21 = Inner LPC (Kurzke) Kernstrom Fan 
Stat. 22-24 = IPC (hier LPC genannt)
Stat. 


"""