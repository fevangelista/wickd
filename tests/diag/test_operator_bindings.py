"""Sanity checks for newly exposed accessors/operators on Operator and Equation."""
import wickd as w


def setup():
    w.reset_space()
    w.add_space("o", "fermion", "occupied", ["i", "j"])   # space index 0
    w.add_space("v", "fermion", "unoccupied", ["a", "b"]) # space index 1


def test_operator_label():
    setup()
    op = w.Operator("T1", [0, 1], [1, 0])
    assert op.label() == "T1"


def test_operator_cre_ann():
    """cre(space) and ann(space) return counts per space"""
    setup()
    # T1: 1 cre in v (space 1), 1 ann in o (space 0)
    op = w.Operator("T1", [0, 1], [1, 0])
    assert op.cre(0) == 0
    assert op.cre(1) == 1
    assert op.ann(0) == 1
    assert op.ann(1) == 0


def test_operator_num_ops():
    setup()
    op = w.Operator("T1", [0, 1], [1, 0])
    assert op.num_ops() == 2  # 1 cre + 1 ann

    op2 = w.Operator("T2", [0, 2], [2, 0])
    assert op2.num_ops() == 4


def test_operator_factor():
    """factor() returns a rational"""
    setup()
    op = w.Operator("T1", [0, 1], [1, 0])
    f = op.factor()
    assert isinstance(f, w.rational)


def test_operator_adjoint():
    """adjoint() swaps cre and ann counts"""
    setup()
    op = w.Operator("T1", [0, 1], [1, 0])
    op_adj = op.adjoint()
    assert op_adj.cre(0) == op.ann(0)
    assert op_adj.cre(1) == op.ann(1)
    assert op_adj.ann(0) == op.cre(0)
    assert op_adj.ann(1) == op.cre(1)
    # adjoint of adjoint = original
    assert op_adj.adjoint() == op


def test_operator_eq_ne_lt():
    """__eq__, __ne__, __lt__ on Operator"""
    setup()
    op1 = w.Operator("A", [0, 1], [1, 0])
    op2 = w.Operator("A", [0, 1], [1, 0])
    op3 = w.Operator("B", [0, 1], [1, 0])

    assert op1 == op2
    assert not (op1 != op2)
    assert op1 != op3
    assert not (op1 == op3)
    # __lt__ defines a total order
    assert (op1 < op3) != (op3 < op1)
    assert not (op1 < op1)


def test_equation_eq():
    """Equation.__eq__ compares lhs, rhs, and factor"""
    setup()

    expr = w.expression("T^{v0}_{o0} a+(v0) a-(o0)")
    eqs = expr.to_manybody_equation("T")
    # flatten all equations from all sectors
    all_eqs = [eq for eqs_list in eqs.values() for eq in eqs_list]
    assert len(all_eqs) >= 1

    eq = all_eqs[0]
    # An equation compared to itself must be equal
    assert eq == eq

    # Two equations built from the same expression are equal
    expr2 = w.expression("T^{v0}_{o0} a+(v0) a-(o0)")
    eqs2 = expr2.to_manybody_equation("T")
    all_eqs2 = [eq for eqs_list in eqs2.values() for eq in eqs_list]
    assert all_eqs[0] == all_eqs2[0]
