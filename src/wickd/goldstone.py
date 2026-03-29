"""
Goldstone diagram visualization as TikZ/LaTeX.

Each SymbolicTerm in an Expression is rendered as one Goldstone diagram:
  - Operators appear as horizontal bars stacked in time order.
  - Hole lines (occupied) carry downward arrows; particle lines (unoccupied)
    carry upward arrows; active/general lines are dashed with upward arrows.
  - Density-matrix tensors (γ, η, λ, …) are drawn with dashed bars and no hat.
  - Lines connecting non-adjacent operators are curved outward.

Public API
----------
expression_to_tikz(expr, standalone=False) -> str
display_diagrams(expr, resolution=1200)
"""

from collections import defaultdict
import re
import wickd

__all__ = ["expression_to_tikz", "display_diagrams"]

# ── Layout parameters ─────────────────────────────────────────────────────────
_VERT_SEP = 2.0   # vertical gap between operator bars
_BAR_WIDTH = 0.8  # minimum half-width of a bar
_CONN_SEP  = 0.6  # horizontal spacing between line attachment points
_EXT_LEN   = 0.8  # length of external-leg stubs

# ── Greek letter table ────────────────────────────────────────────────────────
_GREEK = {
    "alpha": r"\alpha", "beta": r"\beta",   "gamma": r"\gamma",
    "delta": r"\delta", "epsilon": r"\epsilon", "zeta": r"\zeta",
    "eta":   r"\eta",   "theta": r"\theta", "iota":  r"\iota",
    "kappa": r"\kappa", "lambda": r"\lambda","mu":    r"\mu",
    "nu":    r"\nu",    "xi":    r"\xi",    "pi":    r"\pi",
    "rho":   r"\rho",   "sigma": r"\sigma", "tau":   r"\tau",
    "phi":   r"\phi",   "chi":   r"\chi",   "psi":   r"\psi",
    "omega": r"\omega",
}

# Tensor-label prefixes that identify density-matrix insertions
_DENSITY_MATRIX_PREFIXES = frozenset(
    ("gamma", "eta", "lambda", "kappa", "sigma", "phi")
)


# ── Label & classification helpers ────────────────────────────────────────────

def _label_latex(raw):
    r"""'gamma1' → r'\gamma_1',  'f' → 'f',  'eta' → r'\eta'."""
    m = re.fullmatch(r"([a-zA-Z]+)(\d*)", raw)
    if not m:
        return raw
    base, num = m.group(1).lower(), m.group(2)
    latex = _GREEK.get(base, m.group(1))
    return f"{latex}_{{{num}}}" if num else latex


def _is_density_matrix(tensor):
    return any(tensor.label().startswith(p) for p in _DENSITY_MATRIX_PREFIXES)


def _space_type(idx):
    return wickd.osi().space_type(idx.space())


def _index_key(idx):
    return (idx.space(), idx.pos())


def _index_label(idx):
    osi = wickd.osi()
    return osi.indices(idx.space())[idx.pos()]


def _coeff_latex(coeff):
    n, d = coeff.numerator(), coeff.denominator()
    sign = "-" if n < 0 else ""
    n = abs(n)
    return rf"{sign}\frac{{{n}}}{{{d}}}" if d != 1 else rf"{sign}{n}"


# ── Term parser ───────────────────────────────────────────────────────────────

def _parse_term(term):
    """
    Decompose a SymbolicTerm into vertices, internal lines, and external legs.

    Returns
    -------
    tensors  : list[Tensor]
    lines    : list of (key, idx, ti, tj)   with ti <= tj
    ext_legs : list of (key, idx, tensor_idx, slot)
    """
    tensors = list(term.tensors())

    locs = defaultdict(list)
    for t_idx, tensor in enumerate(tensors):
        for idx in tensor.upper():
            locs[_index_key(idx)].append((t_idx, "upper", idx))
        for idx in tensor.lower():
            locs[_index_key(idx)].append((t_idx, "lower", idx))
    for sqop in term.ops():
        slot = "creation" if sqop.is_creation() else "annihilation"
        locs[_index_key(sqop.index())].append((-1, slot, sqop.index()))

    lines, ext_legs = [], []
    for key, entries in locs.items():
        if len(entries) == 2 and entries[0][0] >= 0 and entries[1][0] >= 0:
            (ti, _, idx_i), (tj, _, _) = entries
            if ti > tj:
                ti, tj = tj, ti
            lines.append((key, idx_i, ti, tj))
        else:
            for t_idx, slot, idx in entries:
                ext_legs.append((key, idx, t_idx, slot))

    return tensors, lines, ext_legs


# ── TikZ line builder ─────────────────────────────────────────────────────────

def _propagator(x, yi, yj, stype, skip):
    """
    Return a TikZ draw command for one propagator line.

    Conventions
    -----------
    occupied  → solid,  arrow downward  (hole line)
    unoccupied→ solid,  arrow upward    (particle line)
    general   → dashed, arrow upward    (active-space line)
    Lines that skip intermediate operators are curved outward.
    """
    is_hole   = stype == "occupied"
    is_active = stype == "general"

    line_style = "dashed, " if is_active else ""
    draw_opts  = rf"[{line_style}-{{Stealth}}]"

    # For hole lines reverse start/end so the arrow ends at the lower operator
    y_from, y_to = (yj, yi) if is_hole else (yi, yj)

    # Bend lines that skip intermediate operators outward from the diagram centre
    if skip > 0:
        magnitude = min(40, 20 * skip)
        upward    = y_to > y_from           # after possible hole-reversal
        outward   = magnitude if x <= 0 else -magnitude
        bend      = outward if upward else -outward
    else:
        bend = 0

    path = (
        rf"({x:.2f},{y_from:.2f}) to[bend left={bend}] ({x:.2f},{y_to:.2f})"
        if bend else
        rf"({x:.2f},{y_from:.2f}) -- ({x:.2f},{y_to:.2f})"
    )
    return rf"\draw{draw_opts} {path};"


# ── Per-term TikZ picture ─────────────────────────────────────────────────────

def _term_to_tikz(tensors, lines, ext_legs, coeff):
    n = len(tensors)
    y = [i * _VERT_SEP for i in range(n)]

    n_lines  = len(lines)
    line_xs  = [(k - (n_lines - 1) / 2.0) * _CONN_SEP for k in range(n_lines)]
    bar_half = max((abs(x) for x in line_xs), default=0.0) + _CONN_SEP
    bar_half = max(bar_half, _BAR_WIDTH)

    out  = []
    dots = []   # filled-circle positions collected last so they render on top

    # ── Operator bars ──────────────────────────────────────────────────────
    for k, tensor in enumerate(tensors):
        dm        = _is_density_matrix(tensor)
        bar_style = "thick, dashed" if dm else "thick"
        lbl       = _label_latex(tensor.label())
        out.append(
            rf"  \draw[{bar_style}] ({-bar_half:.2f},{y[k]:.2f})"
            rf" -- ({bar_half:.2f},{y[k]:.2f})"
            rf" node[right] {{${lbl}$}};"
        )

    # ── Internal propagator lines ──────────────────────────────────────────
    for line_k, (key, idx, ti, tj) in enumerate(lines):
        x    = line_xs[line_k]
        lbl  = _index_label(idx)
        st   = _space_type(idx)
        skip = tj - ti - 1

        out.append("  " + _propagator(x, y[ti], y[tj], st, skip))

        ym = (y[ti] + y[tj]) / 2.0
        out.append(
            rf"  \node[right, font=\small] at ({x:.2f},{ym:.2f}) {{${lbl}$}};"
        )
        dots.append((x, y[ti]))
        dots.append((x, y[tj]))

    # ── External legs ──────────────────────────────────────────────────────
    n_ext  = sum(1 for e in ext_legs if e[2] >= 0)
    ext_xs = [(k - (n_ext - 1) / 2.0) * _CONN_SEP for k in range(n_ext)]
    ext_k  = 0
    for key, idx, t_idx, slot in ext_legs:
        if t_idx < 0:
            continue
        x       = ext_xs[ext_k]; ext_k += 1
        yi_val  = y[t_idx]
        lbl     = _index_label(idx)
        st      = _space_type(idx)
        is_hole = st == "occupied"
        ls      = "dashed, " if st == "general" else ""

        if slot == "upper":
            y_end = yi_val + _EXT_LEN
            opts  = rf"[{ls}{{Stealth}}-]" if is_hole else rf"[{ls}-{{Stealth}}]"
        else:
            y_end = yi_val - _EXT_LEN
            opts  = rf"[{ls}-{{Stealth}}]" if is_hole else rf"[{ls}{{Stealth}}-]"

        out.append(
            rf"  \draw{opts} ({x:.2f},{yi_val:.2f}) -- ({x:.2f},{y_end:.2f});"
        )
        ym = (yi_val + y_end) / 2.0
        out.append(
            rf"  \node[right, font=\small] at ({x:.2f},{ym:.2f}) {{${lbl}$}};"
        )
        dots.append((x, yi_val))

    # ── Filled dots at every line–bar junction (drawn last = on top) ────────
    for dx, dy in dots:
        out.append(rf"  \filldraw ({dx:.2f},{dy:.2f}) circle (2pt);")

    # ── Coefficient ────────────────────────────────────────────────────────
    out.append(
        rf"  \node[below] at (0.00,{y[0] - 0.4:.2f}) {{${_coeff_latex(coeff)}$}};"
    )

    body = "\n".join(out)
    return f"\\begin{{tikzpicture}}[>=Stealth, scale=1]\n{body}\n\\end{{tikzpicture}}"


# ── Standalone document wrapper ───────────────────────────────────────────────

def _wrap_standalone(tikz_picture):
    return "\n".join([
        r"\documentclass[border=8pt]{standalone}",
        r"\usepackage{tikz}",
        r"\usetikzlibrary{arrows.meta}",
        r"\begin{document}",
        tikz_picture,
        r"\end{document}",
    ])


# ── PNG compiler ──────────────────────────────────────────────────────────────

def _compile_to_png(src, resolution=1200):
    """Compile a standalone LaTeX string to PNG and return the file path."""
    import subprocess
    import tempfile
    from pathlib import Path

    d   = tempfile.mkdtemp()
    tex = Path(d) / "diag.tex"
    tex.write_text(src)
    subprocess.run(
        ["pdflatex", "-interaction=nonstopmode", "diag.tex"],
        cwd=d, capture_output=True, check=True,
    )
    pdf = Path(d) / "diag.pdf"
    png = Path(d) / "diag.png"

    for cmd in [
        # macOS built-in: -Z sets max pixel dimension (upscales the vector PDF)
        ["sips", "-s", "format", "png", "-Z", str(resolution),
         str(pdf), "--out", str(png)],
        # Linux (poppler-utils)
        ["pdftoppm", "-png", "-r", "150", "-singlefile",
         str(pdf), str(Path(d) / "diag")],
    ]:
        try:
            subprocess.run(cmd, capture_output=True, check=True)
        except FileNotFoundError:
            continue
        if not png.exists():
            candidates = sorted(Path(d).glob("diag*.png"))
            if candidates:
                png = candidates[0]
        if png.exists():
            return str(png)

    raise RuntimeError(
        "No PDF→PNG converter found. "
        "On Linux install poppler-utils (pdftoppm); on macOS sips is built in."
    )


# ── Public API ────────────────────────────────────────────────────────────────

def display_diagrams(expr, resolution=1200, scale=0.5):
    """
    Render a Wick&d Expression as Goldstone diagrams inline in a Jupyter notebook.

    Parameters
    ----------
    expr       : wickd.Expression
    resolution : int
        Maximum pixel dimension of the compiled PNG (default 1200).
    scale      : float
        Fraction of ``resolution`` used as the display width (default 0.5).

    Requires ``pdflatex`` and either ``sips`` (macOS) or ``pdftoppm`` (Linux).
    """
    from IPython.display import Image, display as ipy_display

    display_width = int(resolution * scale)
    for term, coeff in expr:
        tensors, lines, ext_legs = _parse_term(term)
        if not tensors:
            continue
        src = _wrap_standalone(_term_to_tikz(tensors, lines, ext_legs, coeff))
        png = _compile_to_png(src, resolution=resolution)
        ipy_display(Image(filename=png, width=display_width))


def expression_to_tikz(expr, standalone=False):
    """
    Render a Wick&d Expression as Goldstone diagrams in TikZ/LaTeX.

    Parameters
    ----------
    expr       : wickd.Expression
    standalone : bool
        If True, wrap in a complete compilable LaTeX document.

    Returns
    -------
    str  TikZ source, or a full LaTeX document if standalone=True.
    """
    pictures = []
    for term, coeff in expr:
        tensors, lines, ext_legs = _parse_term(term)
        if not tensors:
            continue
        pictures.append(_term_to_tikz(tensors, lines, ext_legs, coeff))

    diagrams = "\n\\quad\n".join(pictures)

    if not standalone:
        return diagrams

    return "\n".join([
        r"\documentclass[border=8pt]{standalone}",
        r"\usepackage{tikz}",
        r"\usetikzlibrary{arrows.meta}",
        r"\begin{document}",
        diagrams,
        r"\end{document}",
    ])
