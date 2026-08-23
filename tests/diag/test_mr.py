import wickd as w


def print_comparison(val, val2):
    print(f"\n{'-'*30} Result {'-'*30}\n{val}")
    print(f"{'-'*31} Test {'-'*31}\n{val2}")


def initialize():
    w.reset_space()
    w.add_space("o", "fermion", "occupied", ["i", "j", "k", "l", "m", "n"])
    w.add_space("a", "fermion", "general", ["u", "v", "w", "x", "y", "z"])
    w.add_space("v", "fermion", "unoccupied", ["a", "b", "c", "d", "e", "f"])


def initialize_two_general_spaces():
    w.reset_space()
    w.add_space("a", "fermion", "general", ["u", "v", "w", "x"])
    w.add_space("b", "fermion", "general", ["p", "q", "r", "s"])


def initialize_three_general_spaces():
    initialize_two_general_spaces()
    w.add_space("c", "fermion", "general", ["i", "j", "k", "l"])


def test_cumulant_contraction_does_not_mix_general_spaces():
    initialize_two_general_spaces()

    left = w.op("L", ["a+ b+"])
    right = w.op("R", ["b a"])

    wt = w.WickTheorem()
    val = wt.contract(left @ right, minrank=0, maxrank=0)
    ref = w.utils.string_to_expr(
        "L^{}_{a0,b0} R^{a1,b1}_{} gamma1^{b0}_{b1} gamma1^{a0}_{a1}"
    )
    print_comparison(val, ref)
    assert val == ref


def test_cumulant_contraction_does_not_mix_disjoint_general_spaces():
    initialize_two_general_spaces()

    left = w.op("L", ["a+ a"])
    right = w.op("R", ["b+ b"])

    wt = w.WickTheorem()
    val = wt.contract(left @ right, minrank=0, maxrank=0)
    ref = w.utils.string_to_expr("")
    print_comparison(val, ref)
    assert val == ref


def test_cumulant_contraction_does_not_mix_general_spaces_at_odd_rank():
    initialize_two_general_spaces()

    left = w.op("L", ["a b a"])
    right = w.op("R", ["a+ b+"])

    wt = w.WickTheorem()
    val = wt.contract(left @ right, minrank=1, maxrank=1)
    ref = w.utils.string_to_expr(
        "L^{a0,a1,b0}_{} R^{}_{a2,b1} "
        "eta1^{b1}_{b0} eta1^{a2}_{a1} a-(a0)"
    )
    print_comparison(val, ref)
    assert val == ref


def test_cumulant_contraction_mixes_general_spaces_when_enabled():
    initialize_two_general_spaces()

    left = w.op("L", ["a+ b+"])
    right = w.op("R", ["b a"])

    wt = w.WickTheorem()
    wt.enable_mixed_general_contractions(True)
    val = wt.contract(left @ right, minrank=0, maxrank=0)
    ref = w.utils.string_to_expr(
        """-L^{}_{a0,b0} R^{a1,b1}_{} gamma1^{a0}_{b1} gamma1^{b0}_{a1}
+L^{}_{a0,b0} R^{a1,b1}_{} gamma1^{b0}_{b1} gamma1^{a0}_{a1}
+L^{}_{a0,b0} R^{a1,b1}_{} lambda2^{a0,b0}_{a1,b1}"""
    )
    print_comparison(val, ref)
    assert val == ref

    wt.do_canonicalize_graph(False)
    val_without_graph_canonicalization = wt.contract(
        left @ right, minrank=0, maxrank=0
    )
    assert val_without_graph_canonicalization == ref


def test_cumulant_contraction_mixes_disjoint_general_spaces_when_enabled():
    initialize_two_general_spaces()

    left = w.op("L", ["a+ a"])
    right = w.op("R", ["b+ b"])

    wt = w.WickTheorem()
    wt.enable_mixed_general_contractions(True)
    val = wt.contract(left @ right, minrank=0, maxrank=0)
    ref = w.utils.string_to_expr(
        """L^{a1}_{a0} R^{b1}_{b0} eta1^{b0}_{a1} gamma1^{a0}_{b1}
+L^{a1}_{a0} R^{b1}_{b0} lambda2^{a0,b0}_{a1,b1}"""
    )
    print_comparison(val, ref)
    assert val == ref


def test_cumulant_contraction_mixes_general_spaces_at_odd_rank_when_enabled():
    initialize_two_general_spaces()

    left = w.op("L", ["a b a"])
    right = w.op("R", ["a+ b+"])

    wt = w.WickTheorem()
    wt.enable_mixed_general_contractions(True)
    val = wt.contract(left @ right, minrank=1, maxrank=1)
    ref = w.utils.string_to_expr(
        """L^{a0,a1,b0}_{} R^{}_{a2,b1} eta1^{b1}_{a1} eta1^{a2}_{a0} a-(b0)
-L^{a0,a1,b0}_{} R^{}_{a2,b1} eta1^{a2}_{b0} eta1^{b1}_{a1} a-(a0)
+L^{a0,a1,b0}_{} R^{}_{a2,b1} eta1^{b1}_{b0} eta1^{a2}_{a1} a-(a0)
+1/2 L^{a0,a1,b0}_{} R^{}_{a2,b1} lambda2^{a2,b1}_{a0,a1} a-(b0)
+L^{a0,a1,b0}_{} R^{}_{a2,b1} lambda2^{a2,b1}_{a1,b0} a-(a0)"""
    )
    print_comparison(val, ref)
    assert val == ref


def test_cumulant_contraction_spans_three_general_spaces_when_enabled():
    initialize_three_general_spaces()

    left = w.op("L", ["a+ b+"])
    right = w.op("R", ["c a"])

    wt = w.WickTheorem()
    wt.enable_mixed_general_contractions(True)
    val = wt.contract(left @ right, minrank=0, maxrank=0)
    ref = w.utils.string_to_expr(
        """-L^{}_{a0,b0} R^{a1,c0}_{} gamma1^{a0}_{c0} gamma1^{b0}_{a1}
+L^{}_{a0,b0} R^{a1,c0}_{} gamma1^{b0}_{c0} gamma1^{a0}_{a1}
+L^{}_{a0,b0} R^{a1,c0}_{} lambda2^{a0,b0}_{a1,c0}"""
    )
    print_comparison(val, ref)
    assert val == ref


def test_max_cumulant_applies_to_mixed_general_contractions():
    initialize_two_general_spaces()

    left = w.op("L", ["a+ b+"])
    right = w.op("R", ["b a"])

    wt = w.WickTheorem()
    wt.enable_mixed_general_contractions(True)
    wt.set_max_cumulant(1)
    val = wt.contract(left @ right, minrank=0, maxrank=0)
    ref = w.utils.string_to_expr(
        """-L^{}_{a0,b0} R^{a1,b1}_{} gamma1^{a0}_{b1} gamma1^{b0}_{a1}
+L^{}_{a0,b0} R^{a1,b1}_{} gamma1^{b0}_{b1} gamma1^{a0}_{a1}"""
    )
    print_comparison(val, ref)
    assert val == ref


def test_mr1():
    initialize()
    T1aa = w.op("t", ["a+ a"])
    Faa = w.op("f", ["a+ a"])

    wt = w.WickTheorem()
    val = wt.contract(w.rational(1), Faa @ T1aa, 0, 0)
    ref = w.utils.string_to_expr(
        """eta1^{a1}_{a0} f^{a0}_{a2} gamma1^{a2}_{a3} t^{a3}_{a1}"
f^{a1}_{a0} lambda2^{a0,a3}_{a1,a2} t^{a2}_{a3}"""
    )
    print_comparison(val, ref)
    assert val == ref


def test_mr2():
    initialize()
    T1aa = w.op("t", ["a+ a"])
    Faa = w.op("f", ["a+ a"])

    wt = w.WickTheorem()
    wt.set_print(w.PrintLevel.summary)
    val = wt.contract(w.rational(1), Faa @ T1aa @ T1aa, 0, 2)
    ref_expr = """eta1^{a1}_{a0} eta1^{a3}_{a2} f^{a2}_{a4} gamma1^{a4}_{a5} t^{a5}_{a1} t^{a0}_{a3}
+eta1^{a1}_{a0} f^{a3}_{a2} lambda2^{a2,a5}_{a3,a4} t^{a0}_{a5} t^{a4}_{a1}
-eta1^{a1}_{a0} f^{a0}_{a2} gamma1^{a4}_{a3} gamma1^{a2}_{a5} t^{a5}_{a4} t^{a3}_{a1}
+eta1^{a1}_{a0} f^{a0}_{a2} lambda2^{a2,a5}_{a3,a4} t^{a4}_{a5} t^{a3}_{a1}
-eta1^{a1}_{a0} f^{a0}_{a2} lambda2^{a2,a5}_{a3,a4} t^{a3}_{a5} t^{a4}_{a1}
-f^{a1}_{a0} gamma1^{a3}_{a2} lambda2^{a0,a5}_{a1,a4} t^{a2}_{a5} t^{a4}_{a3}
-f^{a1}_{a0} gamma1^{a0}_{a2} lambda2^{a4,a5}_{a1,a3} t^{a3}_{a5} t^{a2}_{a4}
+f^{a1}_{a0} gamma1^{a0}_{a2} lambda2^{a4,a5}_{a1,a3} t^{a3}_{a4} t^{a2}_{a5}
-f^{a1}_{a0} lambda3^{a0,a4,a5}_{a1,a2,a3} t^{a2}_{a5} t^{a3}_{a4}"""
    ref = w.utils.string_to_expr(ref_expr)
    # assert val == ref


def test_mr3():
    initialize()
    T2aaaa = w.op("t", ["a+ a+ a a"])
    Vaaaa = w.op("f", ["a+ a+ a a"])

    wt = w.WickTheorem()
    wt.set_print(w.PrintLevel.summary)
    val = wt.contract(w.rational(1), Vaaaa @ T2aaaa, 0, 0)
    ref_expr = """eta1^{a1}_{a0} eta1^{a3}_{a2} f^{a2}_{a4} gamma1^{a4}_{a5} t^{a5}_{a1} t^{a0}_{a3}
+eta1^{a1}_{a0} f^{a3}_{a2} lambda2^{a2,a5}_{a3,a4} t^{a0}_{a5} t^{a4}_{a1}
-eta1^{a1}_{a0} f^{a0}_{a2} gamma1^{a4}_{a3} gamma1^{a2}_{a5} t^{a5}_{a4} t^{a3}_{a1}
+eta1^{a1}_{a0} f^{a0}_{a2} lambda2^{a2,a5}_{a3,a4} t^{a4}_{a5} t^{a3}_{a1}
-eta1^{a1}_{a0} f^{a0}_{a2} lambda2^{a2,a5}_{a3,a4} t^{a3}_{a5} t^{a4}_{a1}
-f^{a1}_{a0} gamma1^{a3}_{a2} lambda2^{a0,a5}_{a1,a4} t^{a2}_{a5} t^{a4}_{a3}
-f^{a1}_{a0} gamma1^{a0}_{a2} lambda2^{a4,a5}_{a1,a3} t^{a3}_{a5} t^{a2}_{a4}
+f^{a1}_{a0} gamma1^{a0}_{a2} lambda2^{a4,a5}_{a1,a3} t^{a3}_{a4} t^{a2}_{a5}
-f^{a1}_{a0} lambda3^{a0,a4,a5}_{a1,a2,a3} t^{a2}_{a5} t^{a3}_{a4}"""
    ref = w.utils.string_to_expr(ref_expr)
