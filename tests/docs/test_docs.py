import wickd as w


def test_readme_example():
    ref_latex = r"""+ {f}^{j}_{i} {t}^{k}_{a} \hat{a}^{i a}_{j k} \\ 
+ {f}^{a}_{i} {t}^{i}_{a} \\ 
+ {f}^{a}_{i} {t}^{j}_{a} \hat{a}^{i}_{j} \\ 
- {f}^{a}_{i} {t}^{i}_{b} \hat{a}^{b}_{a} \\ 
- {f}^{a}_{i} {t}^{j}_{b} \hat{a}^{i b}_{j a} \\ 
- {f}^{i}_{j} {t}^{j}_{a} \hat{a}^{a}_{i} \\ 
+ {f}^{i}_{a} {t}^{j}_{b} \hat{a}^{a b}_{i j} \\ 
+ {f}^{b}_{a} {t}^{i}_{b} \hat{a}^{a}_{i} \\ 
- {f}^{b}_{a} {t}^{i}_{c} \hat{a}^{a c}_{i b}"""

    ### Example from the README.md file
    # Define a Slater determinant reference with occupied (o) and virtual (v) spaces
    w.reset_space()
    w.add_space("o", "fermion", "occupied", ["i", "j", "k", "l"])
    w.add_space("v", "fermion", "unoccupied", ["a", "b", "c", "d"])

    # Build the operators F
    F = w.utils.gen_op("f", 1, "ov", "ov")
    # Build the one-body operator T1 (occupied to virtual)
    T1 = w.op("t", ["v+ o"])

    # Build the product of F and T1
    F_T1 = F @ T1

    # Apply Wick's theorem and collect all fully contracted terms
    wt = w.WickTheorem()
    expr = wt.contract(F_T1, minrank=0, maxrank=4)

    print(f"F T_{{1}} = {expr.latex()}")
    ### End of example from the README.md file

    assert expr.latex() == ref_latex
