"""
Heavy benchmark test for canonicalize_contraction_graph.

This test intentionally takes ~10 seconds.  It is excluded from the normal
pytest run and must be invoked explicitly:

    pytest -m slow
    pytest tests/diag/test_canonicalize_graph_heavy.py

Purpose
-------
This test stresses graph canonicalization by running the full BCH + Wick
contraction pipeline at excitation level n = 5, which produces 51 954 operator
products for ``canonicalize_contraction_graph`` to process.

The test also quantifies canonicalization overhead directly by repeating the
same contraction with ``do_canonicalize_graph(False)``. Canonicalization
previously took roughly 12× longer because it enumerated every permutation of
the elementary contractions. The sorted canonical-order implementation should
keep that ratio below 3×.

Timing reference after sorting canonical contractions (Apple Silicon M-series,
single thread):
    BCH n=5 expansion : ~0.9 s
    contract canon=ON : ~0.9 s
    contract canon=OFF: ~0.7 s
    ratio             : ~1.3×
"""

import time

import pytest

import wickd as w


# ---------------------------------------------------------------------------
# helpers shared across tests in this module
# ---------------------------------------------------------------------------

def _setup_ov():
    w.reset_space()
    w.add_space("o", "fermion", "occupied",   ["i", "j", "k", "l", "m", "n"])
    w.add_space("v", "fermion", "unoccupied", ["a", "b", "c", "d", "e", "f"])


def _make_hamiltonian():
    """H = E_0 + F + V  (scalar + 1-body + 2-body)."""
    return (
        w.op("E_0", [""])
        + w.utils.gen_op("f", 1, "ov", "ov")
        + w.utils.gen_op("v", 2, "ov", "ov")
    )


def _make_T(n: int):
    """Cluster operator T_1 + T_2 + … + T_n as a single w.op."""
    components = [f"{'v+' * k} {'o' * k}" for k in range(1, n + 1)]
    return w.op("t", components)


# ---------------------------------------------------------------------------
# benchmark test
# ---------------------------------------------------------------------------

@pytest.mark.slow
def test_bch_n5_canonicalization_cost():
    """BCH to 4th order with T up to T_5 (CC excitation level n=5).

        Produces 51 954 operator products; each is passed to
        canonicalize_contraction_graph.

    Assertions
    ----------
    1. The contracted expression is non-trivial (≥100 distinct terms).
        2. Canonicalization overhead must remain below 3× canon=OFF on the same
           input.
        3. Canonicalization may merge terms but must never create more terms.
    """
    _setup_ov()
    H = _make_hamiltonian()
    T = _make_T(5)

    hbar = w.bch_series(H, T, 4)
    assert len(hbar) == 51_954, (
        f"Unexpected hbar size {len(hbar)}; BCH expansion may have changed"
    )

    # ── canon=ON (the slow path) ──────────────────────────────────────────
    wt_on = w.WickTheorem()
    wt_on.set_single_threaded(True)

    t0 = time.perf_counter()
    result_on = wt_on.contract(w.rational(1), hbar, 0, 10)
    t_on = time.perf_counter() - t0

    # ── canon=OFF (the fast baseline) ────────────────────────────────────
    wt_off = w.WickTheorem()
    wt_off.set_single_threaded(True)
    wt_off.do_canonicalize_graph(False)

    t0 = time.perf_counter()
    result_off = wt_off.contract(w.rational(1), hbar, 0, 10)
    t_off = time.perf_counter() - t0

    # ── assertions ───────────────────────────────────────────────────────
    n_on  = len(result_on)
    n_off = len(result_off)

    # 1. non-trivial result
    assert n_on >= 100, f"Expected ≥100 terms, got {n_on}"

    # 2. canonicalization must no longer dominate the contraction time
    ratio = t_on / t_off if t_off > 0 else float("inf")
    assert ratio <= 3.0, (
        f"canon=ON ({t_on:.2f}s) should be ≤3× slower than canon=OFF "
        f"({t_off:.2f}s); got ratio={ratio:.1f}×.  "
        "Canonical contraction ordering may have regressed."
    )

    # 3. canon never introduces more terms than no-canon
    assert n_on <= n_off, (
        f"canon=ON produced more terms ({n_on}) than canon=OFF ({n_off})"
    )
