"""
Standalone test for SOFCThermoAdd.

Tests:
  1. output_port_data()  -- element dict is correct before setup()
  2. compute()           -- physical correctness (mass balance, O-atom balance)
  3. check_partials()    -- derivative correctness via complex step

Run with:
    python example_cycles/sofc_thermo_add_test.py
"""

import openmdao.api as om
import numpy as np

import pycycle.api as pyc
from pycycle.elements.sofc_thermo_add import SOFCThermoAdd

# ---------------------------------------------------------------------------
# Inlet composition definitions
# ---------------------------------------------------------------------------

# Pure H2 anode fuel feed (mole-ratio element dict)
H2_COMPOSITION = {'H': 1.0}

# Simplified air cathode feed: N2/O2 only, molecular ratio ~3.76:1
# Element ratios: N atoms = 2 × 0.79 = 1.58, O atoms = 2 × 0.21 = 0.42
AIR_COMPOSITION = {'N': 0.76, 'O': 0.24}

FARADAY   = 96485.3321233100184  # C/mol
G_PER_LBM = 453.592              # g/lbm
MW_O      = 16.0                 # g/mol
MW_O2     = 32.0                 # g/mol


# ===========================================================================
# Test 1: output_port_data() before setup
# ===========================================================================

def test_output_port_data():
    print("\n=== Test 1: output_port_data() ===")

    anode_comp = SOFCThermoAdd(
        reaction_type='anode',
        spec=pyc.species_data.janaf,
        inflow_composition=H2_COMPOSITION)

    cathode_comp = SOFCThermoAdd(
        reaction_type='cathode',
        spec=pyc.species_data.janaf,
        inflow_composition=AIR_COMPOSITION)

    an_out = anode_comp.output_port_data()
    cat_out = cathode_comp.output_port_data()

    print(f"  Anode   inlet:  {H2_COMPOSITION}")
    print(f"  Anode   outlet: {an_out}")
    print(f"  Cathode inlet:  {AIR_COMPOSITION}")
    print(f"  Cathode outlet: {cat_out}")

    assert 'O' in an_out,  "Anode outlet must contain 'O' (from O²⁻ transfer)"
    assert 'H' in an_out,  "Anode outlet must contain 'H'"
    assert 'N' in cat_out, "Cathode outlet must contain 'N'"
    assert 'O' in cat_out, "Cathode outlet must contain 'O'"
    # Cathode outlet should NOT have gained any new elements
    assert set(cat_out.keys()) == set(AIR_COMPOSITION.keys()), \
        "Cathode outlet element set should equal cathode inlet element set"

    print("  PASS")


# ===========================================================================
# Test 2: Physical correctness (mass balance + O-atom balance)
# ===========================================================================

def _build_problem(reaction_type, inflow_comp, W_in, I):
    """Build and run a single-component OpenMDAO problem."""
    p = om.Problem()
    ivc = p.model.add_subsystem('ivc', om.IndepVarComp(), promotes=['*'])

    comp = SOFCThermoAdd(
        reaction_type=reaction_type,
        spec=pyc.species_data.janaf,
        inflow_composition=inflow_comp)
    p.model.add_subsystem('rxn', comp, promotes=['*'])

    # Determine b0 from Properties so we feed a valid initial composition
    from pycycle.thermo.cea.species_data import Properties
    thermo = Properties(pyc.species_data.janaf, init_elements=inflow_comp)
    b0_in = thermo.b0

    ivc.add_output('Fl_I:stat:W',        val=W_in, units='lbm/s')
    ivc.add_output('Fl_I:tot:composition', val=b0_in)
    ivc.add_output('I',                  val=I,   units='A')

    p.setup(force_alloc_complex=True)
    p.run_model()
    return p


def test_mass_balance():
    """
    At current I, the O atoms lost by the cathode must equal the O atoms gained
    by the anode.  Also checks the total mass balance:
        W_an_out - W_an_in  =  -(W_cat_out - W_cat_in)
    """
    print("\n=== Test 2: Mass balance and O-atom balance ===")

    W_an  = 0.01    # lbm/s  (pure H2 inlet)
    W_cat = 0.20    # lbm/s  (air inlet)
    I     = 200.0   # A

    p_an  = _build_problem('anode',   H2_COMPOSITION,  W_an,  I)
    p_cat = _build_problem('cathode', AIR_COMPOSITION, W_cat, I)

    W_an_out  = float(p_an ['Wout'].flat[0])
    W_cat_out = float(p_cat['Wout'].flat[0])

    delta_W_an  = W_an_out  - W_an    # should be positive (gains O²⁻ mass)
    delta_W_cat = W_cat_out - W_cat   # should be negative (loses O2 mass)

    print(f"  I = {I} A")
    print(f"  Anode   ΔW = {delta_W_an*G_PER_LBM*1e3:.4f} mg/s  (should be > 0)")
    print(f"  Cathode ΔW = {delta_W_cat*G_PER_LBM*1e3:.4f} mg/s  (should be < 0)")
    print(f"  Sum ΔW     = {(delta_W_an + delta_W_cat)*G_PER_LBM*1e3:.6f} mg/s  (should be 0)")

    # Expected from Faraday's law:
    ndot_H2 = I / (2.0 * FARADAY)
    expected_delta_W_an = ndot_H2 * MW_O / G_PER_LBM    # O²⁻ mass gained
    ndot_O2 = I / (4.0 * FARADAY)
    expected_delta_W_cat = -ndot_O2 * MW_O2 / G_PER_LBM  # O2 mass lost

    np.testing.assert_allclose(delta_W_an,  expected_delta_W_an,  rtol=1e-10,
                                err_msg="Anode mass gain wrong")
    np.testing.assert_allclose(delta_W_cat, expected_delta_W_cat, rtol=1e-10,
                                err_msg="Cathode mass loss wrong")
    # Global mass balance: O leaving cathode = O arriving at anode
    np.testing.assert_allclose(delta_W_an + delta_W_cat, 0.0, atol=1e-14,
                                err_msg="Global mass balance violated")

    # --- Composition check: outlet b0_O on anode should be > 0 ---
    b0_an_out = p_an['composition_out']
    from pycycle.thermo.cea.species_data import Properties
    out_thermo = Properties(pyc.species_data.janaf,
                            init_elements=SOFCThermoAdd(
                                reaction_type='anode',
                                spec=pyc.species_data.janaf,
                                inflow_composition=H2_COMPOSITION
                            ).output_port_data())
    O_idx = out_thermo.elements.index('O')
    H_idx = out_thermo.elements.index('H')
    assert b0_an_out[O_idx] > 0, "Anode outlet must have O (from O²⁻ transfer)"
    print(f"  Anode outlet b0_H = {b0_an_out[H_idx]:.6e}  b0_O = {b0_an_out[O_idx]:.6e}")

    print("  PASS")


# ===========================================================================
# Test 3: Scalar I=0 → no reaction, Wout = W_in, composition unchanged
# ===========================================================================

def test_zero_current():
    print("\n=== Test 3: Zero current (no reaction) ===")

    for rtype, comp in [('anode', H2_COMPOSITION), ('cathode', AIR_COMPOSITION)]:
        p = _build_problem(rtype, comp, W_in=0.1, I=0.0)
        W_out = float(p['Wout'].flat[0])
        np.testing.assert_allclose(W_out, 0.1, rtol=1e-12,
                                    err_msg=f"{rtype}: Wout != W_in at I=0")
        # composition_out should equal initial b0 mapped to outlet basis
        from pycycle.thermo.cea.species_data import Properties
        in_thermo  = Properties(pyc.species_data.janaf, init_elements=comp)
        b0_out = p['composition_out']
        # For cathode (same element set), b0_out should equal b0_in
        if rtype == 'cathode':
            np.testing.assert_allclose(b0_out, in_thermo.b0, rtol=1e-12,
                                        err_msg="Cathode composition changed at I=0")
        print(f"  {rtype}: Wout = {W_out:.6f} lbm/s  (inlet = 0.1)  OK")

    print("  PASS")


# ===========================================================================
# Test 4: check_partials (complex-step vs analytic, both CS so just a sanity run)
# ===========================================================================

def test_check_partials():
    print("\n=== Test 4: check_partials ===")

    for rtype, comp, W_in in [
        ('anode',   H2_COMPOSITION,  0.015),
        ('cathode', AIR_COMPOSITION, 0.30),
    ]:
        p = _build_problem(rtype, comp, W_in=W_in, I=150.0)

        # Use FD for the check so it is independent of the CS implementation
        data = p.check_partials(method='fd', compact_print=True, out_stream=None)

        for comp_name, comp_data in data.items():
            for (out, inp), err in comp_data.items():
                rel_err = err['rel error'].forward
                # FD accuracy is ~1e-6; allow 1e-4 to be safe
                assert rel_err < 1e-4 or err['J_fwd'].size == 0, (
                    f"{rtype} | {comp_name}: {out} wrt {inp} "
                    f"rel_err = {rel_err:.3e} (threshold 1e-4)")
        print(f"  {rtype}: partials OK")

    print("  PASS")


# ===========================================================================
# Main
# ===========================================================================

if __name__ == '__main__':
    test_output_port_data()
    test_mass_balance()
    test_zero_current()
    test_check_partials()
    print("\nAll tests passed.")
