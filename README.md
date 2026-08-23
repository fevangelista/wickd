# Wick&d

Wick&d is a Python library for the algebra of second quantization for fermionic systems.
It implements Wick's contraction theorem for general Fermi vacua, making it well-suited for
deriving equations in quantum chemistry methods such as coupled cluster, MBPT, and multireference theories.

<img src="https://github.com/fevangelista/wickd/raw/main/lib/logo.png" width="300">

[![Pip](https://github.com/fevangelista/wickd/actions/workflows/pip.yml/badge.svg)](https://github.com/fevangelista/wickd/actions/workflows/pip.yml)
[![Wheels](https://github.com/fevangelista/wickd/actions/workflows/wheels.yml/badge.svg)](https://github.com/fevangelista/wickd/actions/workflows/wheels.yml)
[![codecov](https://codecov.io/gh/fevangelista/wickd/branch/main/graph/badge.svg?token=oe5ECK9O1N)](https://codecov.io/gh/fevangelista/wickd)
[![DOI](https://zenodo.org/badge/64144811.svg)](https://zenodo.org/badge/latestdoi/64144811)

## Quick start

```python
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
```

When a reference contains multiple general orbital spaces, contractions are
confined to each individual space by default. Enable density and cumulant
contractions spanning those spaces explicitly:

```python
wt.enable_mixed_general_contractions(True)
```

Wick&d supports at most eight total orbital spaces across all field and space
types. Mixed general contractions may involve any subset of the registered
general spaces within this limit.

## Installation

**From PyPI:**
```bash
pip install wickd
```

**Development version:**
```bash
git clone https://github.com/fevangelista/wickd.git
cd wickd
pip install nanobind scikit-build-core[pyproject]
pip install --no-build-isolation -ve .
```

## Requirements

- Python ≥ 3.9
- macOS or Linux
- [Boost](https://www.boost.org/) *(optional)* — enables 1024-bit integer arithmetic for large prefactors; the library builds and runs without it

## Testing

```bash
pip install pytest
pytest
```

## Learning resources

- [Jupyter tutorials](https://github.com/fevangelista/wickd/tree/main/tutorials) — step-by-step introduction to the library
- [Examples](https://github.com/fevangelista/wickd/tree/main/examples) — worked derivations (e.g. coupled cluster equations)

## Citation

If you use Wick&d in your research, please cite:

```bibtex
@article{Evangelista2022,
  author  = {Evangelista, Francesco A.},
  title   = {Automatic derivation of many-body theories based on general {Fermi} vacua},
  journal = {J. Chem. Phys.},
  volume  = {157},
  pages   = {064111},
  year    = {2022},
  doi     = {10.1063/5.0097858},
  url     = {https://doi.org/10.1063/5.0097858},
}
```

## Contributing

Bug reports and pull requests are welcome via [GitHub Issues](https://github.com/fevangelista/wickd/issues).

## License

Wick&d is released under the [MIT License](LICENSE).
