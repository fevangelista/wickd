from pathlib import Path
import wickd as w

# Directory containing *this* file
THIS_DIR = Path(__file__).resolve().parent


def load_expressions(filelist, base_dir: Path = THIS_DIR):
    expressions = []
    for filename in filelist:
        p = Path(filename)
        if not p.is_absolute():
            p = base_dir / p  # read relative to this file's folder
        with p.open("r", encoding="utf-8") as f:
            expressions.append(w.utils.string_to_expr(f.read()))
    return expressions


def get_ref_dsrg2():
    return load_expressions(["H0_dsrg2.ref", "H1_dsrg2.ref", "H2_dsrg2.ref"])


def get_ref_dsrg3():
    return load_expressions(
        ["H0_dsrg3.ref", "H1_dsrg3.ref", "H2_dsrg3.ref", "H3_dsrg3.ref"]
    )


def get_ref_dsrg3p():
    return load_expressions(
        ["H0_dsrg3p.ref", "H1_dsrg3p.ref", "H2_dsrg3p.ref", "H3_dsrg3p.ref"]
    )


def print_comparison(val, val2):
    print(f"\n{'-'*30} Result {'-'*30}\n{val}")
    print(f"{'-'*31} Test {'-'*31}\n{val2}")


def initialize():
    w.reset_space()
    w.add_space("c", "fermion", "occupied", ["i", "j", "k", "l", "m", "n"])
    w.add_space("a", "fermion", "general", ["u", "v", "w", "x", "y", "z"])
    w.add_space("v", "fermion", "unoccupied", ["a", "b", "c", "d", "e", "f"])


def test_dsrg1():
    """Test the DSRG(2) commutator expressions [H,A] with A = A1 + A2"""
    initialize()
    Eref = w.op("E_0", [""])
    F = w.utils.gen_op("f", 1, "cav", "cav")
    V = w.utils.gen_op("v", 2, "cav", "cav")
    H = Eref + F + V

    T1 = w.utils.gen_op("t", 1, "av", "ca", diagonal=False)
    T2 = w.utils.gen_op("t", 2, "av", "ca", diagonal=False)
    A1 = T1 - T1.adjoint()
    A2 = T2 - T2.adjoint()
    A = A1 + A2

    wt = w.WickTheorem()

    H0_test = wt.contract(w.commutator(H, A), 0, 0).canonicalize()
    H1_test = wt.contract(w.commutator(H, A), 2, 2).canonicalize()
    H2_test = wt.contract(w.commutator(H, A), 4, 4).canonicalize()

    H0_dsrg2_ref, H1_dsrg2_ref, H2_dsrg2_ref = get_ref_dsrg2()

    assert H0_test == H0_dsrg2_ref
    assert H1_test == H1_dsrg2_ref
    assert H2_test == H2_dsrg2_ref


def test_dsrg2():
    """Test the DSRG(3) commutator expressions [H,A] with A = A1 + A2 + A3"""
    initialize()
    Eref = w.op("E_0", [""])
    F = w.utils.gen_op("f", 1, "cav", "cav")
    V = w.utils.gen_op("v", 2, "cav", "cav")
    H = Eref + F + V

    T1 = w.utils.gen_op("t", 1, "av", "ca", diagonal=False)
    T2 = w.utils.gen_op("t", 2, "av", "ca", diagonal=False)
    T3 = w.utils.gen_op("t", 3, "av", "ca", diagonal=False)
    A1 = T1 - T1.adjoint()
    A2 = T2 - T2.adjoint()
    A3 = T3 - T3.adjoint()
    A = A1 + A2 + A3

    H0_dsrg3_ref, H1_dsrg3_ref, H2_dsrg3_ref, H3_dsrg3_ref = get_ref_dsrg3()

    wt = w.WickTheorem()

    H0_test = wt.contract(w.commutator(H, A), 0, 0).canonicalize()
    assert H0_test == H0_dsrg3_ref

    H1_test = wt.contract(w.commutator(H, A), 2, 2).canonicalize()
    assert H1_test == H1_dsrg3_ref

    H2_test = wt.contract(w.commutator(H, A), 4, 4).canonicalize()
    assert H2_test == H2_dsrg3_ref

    H3_test = wt.contract(w.commutator(H, A), 6, 6).canonicalize()
    assert H3_test == H3_dsrg3_ref


def test_dsrg3():
    """Test the DSRG(3) commutator expressions [H,A] with
    H = H0 + H1 + H2 + H3 and A = A1 + A2 + A3"""
    initialize()
    Eref = w.op("E_0", [""])
    F = w.utils.gen_op("f", 1, "cav", "cav")
    V = w.utils.gen_op("v", 2, "cav", "cav")
    W = w.utils.gen_op("w", 3, "cav", "cav")
    H = Eref + F + V + W

    T1 = w.utils.gen_op("t", 1, "av", "ca", diagonal=False)
    T2 = w.utils.gen_op("t", 2, "av", "ca", diagonal=False)
    T3 = w.utils.gen_op("t", 3, "av", "ca", diagonal=False)
    A1 = T1 - T1.adjoint()
    A2 = T2 - T2.adjoint()
    A3 = T3 - T3.adjoint()
    A = A1 + A2 + A3

    H0_dsrg3p_ref, H1_dsrg3p_ref, H2_dsrg3p_ref, H3_dsrg3p_ref = get_ref_dsrg3p()

    wt = w.WickTheorem()

    H0_test = wt.contract(w.commutator(H, A), 0, 0).canonicalize()
    assert H0_test == H0_dsrg3p_ref

    H1_test = wt.contract(w.commutator(H, A), 2, 2).canonicalize()
    assert H1_test == H1_dsrg3p_ref

    H2_test = wt.contract(w.commutator(H, A), 4, 4).canonicalize()
    assert H2_test == H2_dsrg3p_ref

    H3_test = wt.contract(w.commutator(H, A), 6, 6).canonicalize()
    assert H3_test == H3_dsrg3p_ref
