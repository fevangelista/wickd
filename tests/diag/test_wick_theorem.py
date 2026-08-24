import wickd as w


def print_comparison(test, ref):
    print(f"Result: {test}\n")
    print(f"Test:   {ref}\n")


def compare_expressions(test, ref):
    test_expr = w.Expression()
    ref_expr = w.Expression()
    for s in ref:
        ref_expr += w.string_to_expr(s)
    for eq in test:
        test_expr += eq.rhs_expression()
    print_comparison(test_expr, ref_expr)
    assert test_expr == ref_expr


def initialize():
    w.reset_space()
    w.add_space("o", "fermion", "occupied", ["i", "j", "k", "l", "m", "n"])
    w.add_space("v", "fermion", "unoccupied", ["a", "b", "c", "d", "e", "f"])


def test_wick_theorem_1():
    """Test {o+} {1}"""
    initialize()
    A = w.op("a", ["o+"])
    B = w.op("b", [""])

    wt = w.WickTheorem()
    test = wt.contract(w.rational(1), A @ B, 0, 2)
    ref = w.expression("a^{}_{o_0} b^{}_{} a+(o_0)")
    print_comparison(test, ref.canonicalize())
    assert test == ref


def test_wick_theorem_2():
    """Test {1} {o}"""
    initialize()
    A = w.op("a", [""])
    B = w.op("b", ["o"])

    wt = w.WickTheorem()
    test = wt.contract(w.rational(1), A @ B, 0, 2)
    ref = w.expression("a^{}_{} b^{o_0}_{} a-(o_0)")
    print_comparison(test, ref.canonicalize())
    assert test == ref


def test_wick_theorem_3():
    """Test {o+} {o}"""
    initialize()
    A = w.op("a", ["o+"])
    B = w.op("b", ["o"])

    wt = w.WickTheorem()
    test = wt.contract(w.rational(1), A @ B, 0, 4)
    ref = w.expression("a^{}_{o_0} b^{o_0}_{}")
    ref += w.expression("a^{}_{o_0} b^{o_1}_{} a+(o_0) a-(o_1)")
    print_comparison(test, ref.canonicalize())
    assert test == ref


def test_wick_theorem_minrank_boundary():
    """A solution at minrank is kept before its descendants are pruned."""
    initialize()
    A = w.op("a", ["o+"])
    B = w.op("b", ["o"])

    wt = w.WickTheorem()
    test = wt.contract(w.rational(1), A @ B, 2, 2)
    ref = w.expression("a^{}_{o_0} b^{o_1}_{} a+(o_0) a-(o_1)")
    print_comparison(test, ref.canonicalize())
    assert test == ref


def test_wick_theorem_4():
    """Test {o+ o+} {o}"""
    initialize()
    A = w.op("a", ["o+ o+"])
    B = w.op("b", ["o"])

    wt = w.WickTheorem()
    test = wt.contract(w.rational(1), A @ B, 0, 4)
    ref = w.expression("a^{}_{o_0,o_1} b^{o_1}_{} a+(o_0)")
    ref += w.expression("1/2 a^{}_{o_0,o_1} b^{o_2}_{} a+(o_0) a+(o_1) a-(o_2)")
    print_comparison(test, ref.canonicalize())
    assert test == ref


def test_wick_theorem_5():
    """Test {o+} {o o}"""
    initialize()
    A = w.op("a", ["o+"])
    B = w.op("b", ["o o"])

    wt = w.WickTheorem()
    test = wt.contract(w.rational(1), A @ B, 0, 4)
    ref = w.expression("a^{}_{o_0} b^{o_1,o_0}_{} a-(o_1)").canonicalize()
    ref += w.expression("1/2 a^{}_{o_0} b^{o_1,o_2}_{} a+(o_0) a-(o_2) a-(o_1)")
    print_comparison(test, ref.canonicalize())
    assert test == ref


def test_wick_theorem_6():
    """Test {o+ o+} {o o}"""
    initialize()
    A = w.op("a", ["o+ o+"])
    B = w.op("b", ["o o"])

    wt = w.WickTheorem()
    test = wt.contract(w.rational(1), A @ B, 0, 4)
    ref = w.expression("1/2 a^{}_{o_0,o_1} b^{o_0,o_1}_{}")
    ref += w.expression("a^{}_{o_0,o_2} b^{o_1,o_2}_{} a+(o_0) a-(o_1)")
    ref += w.expression(
        "1/4 a^{}_{o_0,o_1} b^{o_2,o_3}_{} a+(o_0) a+(o_1) a-(o_3) a-(o_2)"
    )
    print_comparison(test, ref.canonicalize())
    assert test == ref


def test_wick_theorem_7():
    """Test {o+ o+ o+} {o}"""
    initialize()
    A = w.op("a", ["o+ o+ o+"])
    B = w.op("b", ["o"])

    wt = w.WickTheorem()
    test = wt.contract(w.rational(1), A @ B, 0, 4)

    ref = w.expression("1/2 a^{}_{o_0,o_1,o_2} b^{o_2}_{} a+(o_0) a+(o_1)")
    ref += w.expression(
        "1/6 a^{}_{o_0,o_1,o_2} b^{o_3}_{} " "a+(o_0) a+(o_1) a+(o_2) a-(o_3)"
    )

    print_comparison(test, ref.canonicalize())
    assert test == ref


def test_wick_theorem_8():
    """Test {o+ o+ o+} {o o}"""
    initialize()
    A = w.op("a", ["o+ o+ o+"])
    B = w.op("b", ["o o"])

    wt = w.WickTheorem()
    test = wt.contract(w.rational(1), A @ B, 0, 5)

    ref = w.expression("1/2 a^{}_{o_0,o_1,o_2} b^{o_1,o_2}_{} a+(o_0)")

    ref += w.expression("1/2 a^{}_{o_0,o_2,o_3} b^{o_1,o_3}_{} a+(o_0) a+(o_2) a-(o_1)")

    ref += w.expression(
        "1/12 a^{}_{o_0,o_1,o_2} b^{o_3,o_4}_{} "
        "a+(o_0) a+(o_1) a+(o_2) a-(o_4) a-(o_3)"
    )

    print_comparison(test, ref.canonicalize())
    assert test == ref


def test_wick_theorem_9():
    """Test {o+} {o o o}"""
    initialize()
    A = w.op("a", ["o+"])
    B = w.op("b", ["o o o"])

    wt = w.WickTheorem()
    test = wt.contract(w.rational(1), A @ B, 0, 4)

    ref = w.expression("1/2 a^{}_{o_0} b^{o_1,o_2,o_0}_{} a-(o_2) a-(o_1)")

    ref += w.expression(
        "1/6 a^{}_{o_0} b^{o_1,o_2,o_3}_{} " "a+(o_0) a-(o_3) a-(o_2) a-(o_1)"
    )

    print_comparison(test, ref.canonicalize())
    assert test == ref


def test_wick_theorem_10():
    """Test {o+ o+} {o o o}"""
    initialize()
    A = w.op("a", ["o+ o+"])
    B = w.op("b", ["o o o"])

    wt = w.WickTheorem()
    test = wt.contract(w.rational(1), A @ B, 0, 5)

    ref = w.expression("1/2 a^{}_{o_0,o_1} b^{o_2,o_0,o_1}_{} a-(o_2)")

    ref += w.expression("1/2 a^{}_{o_0,o_2} b^{o_1,o_3,o_2}_{} a+(o_0) a-(o_3) a-(o_1)")

    ref += w.expression(
        "1/12 a^{}_{o_0,o_1} b^{o_2,o_3,o_4}_{} "
        "a+(o_0) a+(o_1) a-(o_4) a-(o_3) a-(o_2)"
    )

    print_comparison(test, ref.canonicalize())
    assert test == ref


def test_wick_theorem_11():
    """Test {o+ o+ o+} {o o o}"""
    initialize()
    A = w.op("a", ["o+ o+ o+"])
    B = w.op("b", ["o o o"])

    wt = w.WickTheorem()
    test = wt.contract(w.rational(1), A @ B, 0, 6)

    ref = w.expression("1/6 a^{}_{o_0,o_1,o_2} b^{o_0,o_1,o_2}_{}")

    ref += w.expression("1/2 a^{}_{o_0,o_2,o_3} b^{o_1,o_2,o_3}_{} a+(o_0) a-(o_1)")

    ref += w.expression(
        "1/4 a^{}_{o_0,o_1,o_4} b^{o_2,o_3,o_4}_{} a+(o_0) a+(o_1) a-(o_3) a-(o_2)"
    )

    ref += w.expression(
        "1/36 a^{}_{o_0,o_1,o_2} b^{o_3,o_4,o_5}_{} "
        "a+(o_0) a+(o_1) a+(o_2) a-(o_5) a-(o_4) a-(o_3)"
    )

    print_comparison(test, ref.canonicalize())
    assert test == ref


def test_wick_theorem_12():
    """Test {o+ o+ o+ o+} {o o o o}"""
    initialize()
    A = w.op("a", ["o+ o+ o+ o+"])
    B = w.op("b", ["o o o o"])

    wt = w.WickTheorem()
    test = wt.contract(w.rational(1), A @ B, 0, 8)

    ref = w.expression("1/24 a^{}_{o_0,o_1,o_2,o_3} b^{o_0,o_1,o_2,o_3}_{}")

    ref += w.expression(
        "1/6 a^{}_{o_0,o_2,o_3,o_4} b^{o_1,o_2,o_3,o_4}_{} " "a+(o_0) a-(o_1)"
    )

    ref += w.expression(
        "1/8 a^{}_{o_0,o_1,o_4,o_5} b^{o_2,o_3,o_4,o_5}_{} "
        "a+(o_0) a+(o_1) a-(o_3) a-(o_2)"
    )

    ref += w.expression(
        "1/36 a^{}_{o_0,o_1,o_2,o_6} b^{o_3,o_4,o_5,o_6}_{} "
        "a+(o_0) a+(o_1) a+(o_2) a-(o_5) a-(o_4) a-(o_3)"
    )

    ref += w.expression(
        "1/576 a^{}_{o_0,o_1,o_2,o_3} b^{o_4,o_5,o_6,o_7}_{} "
        "a+(o_0) a+(o_1) a+(o_2) a+(o_3) "
        "a-(o_7) a-(o_6) a-(o_5) a-(o_4)"
    )

    print_comparison(test, ref.canonicalize())
    assert test == ref


def test_wick_theorem_13():
    """Test {o+ o+ o+ o} {o+ o o o}"""
    initialize()
    A = w.op("a", ["o+ o+ o+ o"])
    B = w.op("b", ["o+ o o o"])

    wt = w.WickTheorem()
    test = wt.contract(w.rational(1), A @ B, 0, 8)

    ref = w.expression(
        "-1/6 a^{o_3}_{o_0,o_1,o_2} b^{o_0,o_1,o_2}_{o_4} a+(o_4) a-(o_3)"
    )

    ref += w.expression(
        "1/2 a^{o_3}_{o_0,o_1,o_2} b^{o_0,o_1,o_4}_{o_5} "
        "a+(o_2) a+(o_5) a-(o_4) a-(o_3)"
    )

    ref += w.expression(
        "-1/4 a^{o_3}_{o_0,o_1,o_2} b^{o_2,o_4,o_5}_{o_6} "
        "a+(o_0) a+(o_1) a+(o_6) a-(o_5) a-(o_4) a-(o_3)"
    )

    ref += w.expression(
        "1/36 a^{o_3}_{o_0,o_1,o_2} b^{o_4,o_5,o_6}_{o_7} "
        "a+(o_0) a+(o_1) a+(o_2) a+(o_7) "
        "a-(o_6) a-(o_5) a-(o_4) a-(o_3)"
    )

    print_comparison(test, ref.canonicalize())
    assert test == ref


def test_wick_theorem_14():
    """Test {o+ o+ o+ o o} {o+ o o o}"""
    initialize()
    A = w.op("a", ["o+ o+ o+ o o"])
    B = w.op("b", ["o+ o o o"])

    wt = w.WickTheorem()
    test = wt.contract(w.rational(1), A @ B, 0, 9)

    ref = w.expression(
        "-1/12 a^{o_3,o_4}_{o_0,o_1,o_2} b^{o_0,o_1,o_2}_{o_5} "
        "a+(o_5) a-(o_4) a-(o_3)"
    )

    ref += w.expression(
        "1/4 a^{o_3,o_4}_{o_0,o_1,o_2} b^{o_0,o_1,o_6}_{o_5} "
        "a+(o_2) a+(o_5) a-(o_6) a-(o_4) a-(o_3)"
    )

    ref += w.expression(
        "-1/8 a^{o_3,o_4}_{o_0,o_1,o_2} b^{o_2,o_6,o_7}_{o_5} "
        "a+(o_0) a+(o_1) a+(o_5) a-(o_7) a-(o_6) a-(o_4) a-(o_3)"
    )

    ref += w.expression(
        "1/72 a^{o_3,o_4}_{o_0,o_1,o_2} b^{o_5,o_6,o_7}_{o_8} "
        "a+(o_0) a+(o_1) a+(o_2) a+(o_8) "
        "a-(o_7) a-(o_6) a-(o_5) a-(o_4) a-(o_3)"
    )

    print_comparison(test, ref.canonicalize())
    assert test == ref
