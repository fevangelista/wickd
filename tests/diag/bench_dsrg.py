"""
Benchmark for the DSRG double nested commutator [[H, A], A].

Measures the wall-clock cost of contracting [[H, A], A] at various operator
ranks, with canonicalization on and off.  The DSRG level (2 = A1+A2,
3 = A1+A2+A3) determines the amplitude operator A.  Because
canonicalize_contraction_graph dominates contract() runtime, the ratio
t_canon / t_nocanon directly quantifies the canonicalization overhead.

Usage
-----
Save a baseline (run before an optimization):
    python tests/diag/bench_dsrg.py --save

Compare against the baseline (run after an optimization):
    python tests/diag/bench_dsrg.py --compare

Just print timings without saving or comparing:
    python tests/diag/bench_dsrg.py

Select DSRG levels and max operator rank (default: level 2, maxrank 4):
    python tests/diag/bench_dsrg.py --levels 2 3 --maxrank 2
"""

import argparse
import json
import time
from pathlib import Path

import wickd as w

BASELINE_FILE = Path(__file__).resolve().parent / "bench_dsrg_baseline.json"
REPS = 1  # [[H,A],A] contractions are expensive; one rep is sufficient


# ── setup ─────────────────────────────────────────────────────────────────────

def setup_cav():
    w.reset_space()
    w.add_space("c", "fermion", "occupied",   ["i", "j", "k", "l", "m", "n"])
    w.add_space("a", "fermion", "general",    ["u", "v", "w", "x", "y", "z"])
    w.add_space("v", "fermion", "unoccupied", ["a", "b", "c", "d", "e", "f"])


def make_hamiltonian():
    """H = E_0 + F + V  (scalar + 1-body + 2-body)."""
    return (
        w.op("E_0", [""])
        + w.utils.gen_op("f", 1, "cav", "cav")
        + w.utils.gen_op("v", 2, "cav", "cav")
    )


def make_A(level):
    """Anti-Hermitian amplitude operator A = sum_{k=1}^{level} (T_k - T_k†)."""
    A = None
    for k in range(1, level + 1):
        Tk = w.utils.gen_op("t", k, "av", "ca", diagonal=False)
        Ak = Tk - Tk.adjoint()
        A = Ak if A is None else A + Ak
    return A


# ── timing helpers ────────────────────────────────────────────────────────────

def time_contract(haa, maxrank: int, *, canon: bool) -> tuple[float, int]:
    """Return (min wall-time in seconds, output term count)."""
    best = float("inf")
    nterms = 0
    for _ in range(REPS):
        wt = w.WickTheorem()
        wt.set_single_threaded(True)
        if not canon:
            wt.do_canonicalize_graph(False)
        t0 = time.perf_counter()
        result = wt.contract(w.rational(1), haa, 0, maxrank)
        elapsed = time.perf_counter() - t0
        best = min(best, elapsed)
        nterms = len(result)
    return best, nterms


# ── benchmark ─────────────────────────────────────────────────────────────────

def run_benchmark(levels, maxrank):
    setup_cav()
    H = make_hamiltonian()

    results = {}
    header = (
        f"{'level':>7}  {'rank':>5}  {'[[H,A],A]':>10}  "
        f"{'build(s)':>9}  {'canon(s)':>10}  {'nocanon(s)':>11}  "
        f"{'ratio':>6}  {'terms':>7}"
    )
    print(header)
    print("-" * len(header))

    for level in levels:
        A = make_A(level)

        t0 = time.perf_counter()
        HA  = w.commutator(H, A)
        HAA = w.commutator(HA, A)
        t_build = time.perf_counter() - t0

        for rank in range(0, maxrank + 1, 2):
            t_on,  nterms = time_contract(HAA, rank, canon=True)
            t_off, _      = time_contract(HAA, rank, canon=False)

            ratio = t_on / t_off if t_off > 0 else float("inf")
            key = f"level={level},rank={rank}"
            results[key] = {
                "haa_size":  len(HAA),
                "t_build":   t_build,
                "t_canon":   t_on,
                "t_nocanon": t_off,
                "terms":     nterms,
            }
            print(
                f"{'DSRG-' + str(level):>7}  {rank:>5}  {len(HAA):>10}  "
                f"{t_build:>9.2f}  {t_on:>10.2f}  {t_off:>11.2f}  "
                f"{ratio:>5.1f}x  {nterms:>7}"
            )

    return results


# ── comparison ────────────────────────────────────────────────────────────────

def compare(current: dict, baseline: dict, tol: float = 0.20):
    print(
        f"\n{'case':<22}  {'baseline(s)':>12}  {'current(s)':>11}  "
        f"{'ratio':>7}  {'status':>9}"
    )
    print("-" * 68)
    all_ok = True
    for key, cur in current.items():
        base = baseline.get(key)
        if base is None:
            print(f"  {key:<20}  {'(no baseline)':>12}")
            continue
        t_base = base["t_canon"]
        t_cur  = cur["t_canon"]
        ratio  = t_cur / t_base if t_base > 0 else float("inf")
        ok     = ratio <= (1.0 + tol)
        status = "OK" if ok else "REGRESSED"
        if not ok:
            all_ok = False
        print(
            f"  {key:<20}  {t_base:>12.2f}  {t_cur:>11.2f}  "
            f"{ratio:>7.2f}x  {status:>9}"
        )
    return all_ok


# ── main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--save",    action="store_true",
                        help="Save results as baseline")
    parser.add_argument("--compare", action="store_true",
                        help="Compare against baseline")
    parser.add_argument("--tol",     type=float, default=0.20,
                        help="Allowed slowdown fraction before REGRESSED (default 0.20)")
    parser.add_argument("--levels",  type=int, nargs="+", default=[2],
                        metavar="N",
                        help="DSRG levels to benchmark (default: 2)")
    parser.add_argument("--maxrank", type=int, default=4,
                        metavar="R",
                        help="Maximum operator rank to contract (even; default: 4)")
    args = parser.parse_args()

    print(f"Running DSRG [[H,A],A] benchmark for levels: {args.levels}, maxrank: {args.maxrank}")
    print(f"({REPS} rep(s) per case; reporting minimum time)\n")

    results = run_benchmark(args.levels, args.maxrank)

    if args.save:
        BASELINE_FILE.write_text(json.dumps(results, indent=2))
        print(f"\nBaseline saved to {BASELINE_FILE}")

    if args.compare:
        if not BASELINE_FILE.exists():
            print(f"\nNo baseline found at {BASELINE_FILE}. Run with --save first.")
            raise SystemExit(1)
        baseline = json.loads(BASELINE_FILE.read_text())
        print("\n=== Comparison vs baseline ===")
        ok = compare(results, baseline, tol=args.tol)
        if not ok:
            print("\nOne or more cases REGRESSED.")
            raise SystemExit(1)
        else:
            print("\nAll cases within tolerance.")


if __name__ == "__main__":
    main()
