import wickd as w


def test_expression():
    """Test the Expression class (initialization, printing)"""
    w.reset_space()
    w.add_space("o", "fermion", "occupied", ["i", "j"])
    w.add_space("a", "fermion", "general", ["u", "v"])
    w.add_space("v", "fermion", "unoccupied", ["a", "b", "c"])

    expr = w.Expression()

    term = w.Term()
    term.add([w.cre("v_0"), w.ann("o_0")])
    term.set_normal_ordered(True)
    expr += term
    assert str(expr) == "{ a+(v0) a-(o0) }"

    term = w.SymbolicTerm()
    term.set_normal_ordered(True)
    term.add([w.cre("a_0")])
    expr += (term, w.rational(1, 2))
    expr_str = """1/2 { a+(a0) }
+{ a+(v0) a-(o0) }"""
    assert str(expr) == expr_str

    expr2 = w.Expression()

    term = w.SymbolicTerm()
    term.set_normal_ordered(True)
    term.add([w.cre("a_0")])
    expr2 += term

    expr_str = """{ a+(a0) }"""
    assert str(expr2) == expr_str

    expr += expr2
    expr_str = """3/2 { a+(a0) }
+{ a+(v0) a-(o0) }"""
    assert str(expr) == expr_str
    for idx, i in enumerate(expr):
        if idx == 0:
            assert str(i[0]) == "{ a+(a0) }"
            assert i[1] == w.rational(3, 2)
        if idx == 1:
            assert str(i[0]) == "{ a+(v0) a-(o0) }"
            assert i[1] == w.rational(1, 1)


def test_expression2():
    """Test the Expression class (initialization from operator expression)"""
    w.reset_space()
    w.add_space("o", "fermion", "occupied", ["i", "j"])
    w.add_space("a", "fermion", "general", ["u", "v"])
    w.add_space("v", "fermion", "unoccupied", ["a", "b", "c"])
    expr = w.operator_expr("T", ["v+ v+ o o"], True)
    assert str(expr) == "T^{o0,o1}_{v0,v1} { a+(v0) a+(v1) a-(o1) a-(o0) }"
    expr = w.operator_expr("T", ["v+ v+ v v"], True)
    assert str(expr) == "T^{v2,v3}_{v0,v1} { a+(v0) a+(v1) a-(v3) a-(v2) }"
    expr = w.operator_expr("T", ["v+ a+ a o"], True)
    assert str(expr) == "T^{o0,a1}_{v0,a0} { a+(v0) a+(a0) a-(a1) a-(o0) }"
    expr = w.operator_expr("T", ["v+ a+ o a"], True)
    assert str(expr) == "T^{a1,o0}_{v0,a0} { a+(v0) a+(a0) a-(o0) a-(a1) }"


def test_expression3():
    """Test the Expression class (initialization from operator expression)"""
    w.reset_space()
    w.add_space("o", "fermion", "occupied", ["i", "j"])
    w.add_space("a", "fermion", "general", ["u", "v"])
    w.add_space("v", "fermion", "unoccupied", ["a", "b", "c"])
    # empty string corresponds to no operator
    expr = w.expression("")
    assert str(expr) == ""

    # identity
    expr = w.expression("1")
    assert str(expr) == "1"

    # minus one
    expr = w.expression("-1")
    assert str(expr) == "-1"

    # tensors and operators
    expr = w.expression("-t^{a_1}_{o_0} a+(a_1) a-(o_0)")
    assert str(expr) == "-t^{a1}_{o0} a+(a1) a-(o0)"

    # tensors and operators
    expr = w.expression("-t^{a_1}_{o_0} {a+(a_1) a-(o_0)}")
    assert str(expr) == "-t^{a1}_{o0} { a+(a1) a-(o0) }"


def test_expression4():
    w.reset_space()
    w.add_space("o", "fermion", "occupied", ["i", "j", "k", "l", "m"])
    w.add_space("v", "fermion", "unoccupied", ["a", "b", "c", "d", "e", "f"])
    C = w.op("C", ["o+ o"])
    D = w.op("D", ["v+ v"])
    CD_bch = w.bch_series(C, D, 2)
    CD_bch.canonicalize()
    assert CD_bch.size() == 1


def test_expression5():
    w.reset_space()
    w.add_space("o", "fermion", "occupied", ["i", "j", "k", "l", "m"])
    w.add_space("v", "fermion", "unoccupied", ["a", "b", "c", "d", "e", "f"])
    # empty string corresponds to no operator
    expr = w.expression("f^{o0}_{}")
    assert str(expr) == "f^{o0}_{}"


def test_expression_dot_and_norm():
    w.reset_space()
    w.add_space("o", "fermion", "occupied", ["i", "j", "k", "l", "m"])
    w.add_space("v", "fermion", "unoccupied", ["a", "b", "c", "d", "e", "f"])

    expr = w.expression("2 a+(v_1) a-(o_0)")
    assert expr.dot(expr) == w.rational(4)
    assert expr.norm() == 2.0

    expr2 = w.expression("-1 a+(v_2) a-(o_0)")
    assert expr.dot(expr2) == w.rational(0)
    assert expr2.dot(expr) == w.rational(0)
    assert expr2.norm() == 1.0

    expr3 = w.expression("a+(v_1) a-(o_0)") - expr2 + expr
    assert expr3.dot(expr) == w.rational(6)
    assert expr3.dot(expr2) == w.rational(-1)


def test_str_roundtrip_positive_integer_coefficient():
    """Regression test for Expression::str() dropping the space between a
    positive single-digit integer coefficient and the following tensor name.

    Root cause (expression.cc:68):
        if ((factor_str.size() > 1) and (factor_str != "-")) {
            factor_str += " ";
        }
    The guard `> 1` fails for a one-character coefficient string such as "2",
    so the serialized form becomes "2f^..." instead of "2 f^...".  The parser
    then reads the leading "2" as part of the tensor label, making round-trips
    through str() / string_to_expr() produce a different Expression.
    """
    w.reset_space()
    w.add_space("o", "fermion", "occupied", ["i", "j", "k", "l"])
    w.add_space("v", "fermion", "unoccupied", ["a", "b", "c", "d"])

    # Build a known expression whose first term has coefficient 2.
    # The string "2 f^{o0}_{o1} t^{o1}_{v0}" is valid input for string_to_expr,
    # so expr1 is well-formed.
    expr1 = w.utils.string_to_expr("2 f^{o0}_{o1} t^{o1}_{v0}")

    # The serialized form must have a space after the coefficient so that
    # string_to_expr can parse it back faithfully.
    assert str(expr1) == "2 f^{o0}_{o1} t^{o1}_{v0}", (
        "str(Expression) must emit a space between coefficient '2' and the "
        "first tensor name; got: " + repr(str(expr1))
    )

    # Round-trip: parse(str(expr)) must equal the original expression.
    expr2 = w.utils.string_to_expr(str(expr1))
    assert expr1 == expr2, (
        "str(Expression) → string_to_expr round-trip failed for coefficient 2"
    )

    # The same must hold for all single-digit positive integer coefficients.
    for n in range(2, 10):
        e1 = w.utils.string_to_expr(f"{n} f^{{o0}}_{{o1}} t^{{o1}}_{{v0}}")
        e2 = w.utils.string_to_expr(str(e1))
        assert e1 == e2, (
            f"Round-trip failed for coefficient {n}: "
            f"str produced {repr(str(e1))}"
        )

    # Coefficients that already serialise correctly must still work:
    # +/- sign present (second-or-later term), fraction, coefficient 1, -1.
    for coeff_str in ["1/2", "-1/3", "3/7"]:
        e1 = w.utils.string_to_expr(
            f"{coeff_str} f^{{o0}}_{{o1}} t^{{o1}}_{{v0}}"
        )
        e2 = w.utils.string_to_expr(str(e1))
        assert e1 == e2, f"Round-trip failed for coefficient {coeff_str!r}"
