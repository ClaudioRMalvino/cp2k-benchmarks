#!/usr/bin/env python3
"""Render the compiler-flags appendix table as a PDF figure so it does not
count toward thesis word count.  Styling matches plots/maqao_tables.py
(DejaVu Serif, navy header, slate stripes).

Output: plots/thesis_figures/compiler_flags_table.pdf
"""
import os
import textwrap

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT_DIR = "/home/crm98/cp2k-benchmarks/plots/thesis_figures"
OUT_PDF = os.path.join(OUT_DIR, "compiler_flags_table.pdf")

# Palette lifted from maqao_tables.py for visual consistency
HEADER_COLOUR     = "#133844"
HEADER_TEXTCOLOUR = "#FFFFFF"
ROW_COLOURS       = ["#ECEEF1", "#FFFFFF"]
GRID_COLOUR       = "#B5BDC8"
BODY_TEXTCOLOUR   = "#232830"

HEADER_FONTSIZE = 10
BODY_FONTSIZE   = 9

plt.rcParams.update({
    "font.family":      "serif",
    "font.serif":       ["DejaVu Serif", "Times New Roman", "Times"],
    "mathtext.fontset": "cm",
    "axes.labelcolor":  BODY_TEXTCOLOUR,
    "text.color":       BODY_TEXTCOLOUR,
})


HEADERS = ("Flag", "Languages", "Purpose / notes")

# Compiler flags passed by the build scripts (cp2k_CSD3_master_build.sh and
# cp2k_CSD3_opt_build.sh).  Both scripts pass the same set; CMake's Release
# preset and CP2K's CMakeLists.txt add further flags downstream, but those
# are build-system plumbing rather than choices made for this work.
FLAGS = [
    ("-O2",
        "C, C++, Fortran",
        "Base optimisation level: aggressive inlining, scalar replacement "
        "of aggregates, and inner-loop transforms enabled."),
    ("-g",
        "C, C++, Fortran",
        "Emit DWARF debug symbols so MAQAO and perf can resolve samples "
        "back to source lines."),
    ("-xCORE-AVX512",
        "C, C++, Fortran",
        "Target the Intel Ice Lake-server AVX-512 ISA; implies all earlier "
        "vector extensions (SSE 4.2, AVX, AVX2)."),
    ("-qopenmp",
        "C, C++, Fortran",
        "Enable OpenMP directives.  Required for the OMP-threaded NN "
        "evaluation in the feature/nnp-native-spline-omp branch and for "
        "the threaded MKL paths used by CP2K's solver."),
    ("-funroll-loops",
        "Fortran",
        "Unroll loops with statically known trip counts; complements the "
        "vectoriser on the symmetry-function kernels."),
    ("-ftree-vectorize",
        "Fortran",
        "Run the inner-loop vectoriser independently of the optimisation "
        "level."),
]

# Column widths (fractions; must sum to ~1).  Purpose gets most space.
COL_WIDTHS = [0.22, 0.20, 0.58]
WRAP_CHARS = (22, 16, 45)   # per-column character wrap width


def _wrap_flag(text, width):
    """Wrap a flag-style string at width, allowing breaks at `_`, `=`, `,`
    and `/` in addition to spaces/hyphens that textwrap already understands.
    Keeps the break characters at the end of the line that precedes the
    break, so the printed flag still reads naturally."""
    if len(text) <= width:
        return text
    breakables = "_=,/ "
    tokens, current = [], ""
    for ch in text:
        current += ch
        if ch in breakables:
            tokens.append(current)
            current = ""
    if current:
        tokens.append(current)
    lines, line = [], ""
    for t in tokens:
        if line and len(line) + len(t) > width:
            lines.append(line.rstrip())
            line = t
        else:
            line += t
    if line:
        lines.append(line.rstrip())
    return "\n".join(lines)


def _wrap_cell(text, width, is_flag=False):
    """Wrap a cell into multiple lines and return the joined string."""
    if width <= 0:
        return text
    if is_flag:
        return _wrap_flag(text, width)
    return "\n".join(textwrap.wrap(text, width=width)) or text


def build_table_cells():
    """Return the wrapped cell-text 2D list."""
    return [
        [
            _wrap_cell(flag,    WRAP_CHARS[0], is_flag=True),
            _wrap_cell(langs,   WRAP_CHARS[1]),
            _wrap_cell(purpose, WRAP_CHARS[2]),
        ]
        for flag, langs, purpose in FLAGS
    ]


def estimate_row_heights(cell_text):
    """One unit per line in the tallest cell of each row."""
    return [max(c.count("\n") + 1 for c in row) for row in cell_text]


def render():
    os.makedirs(OUT_DIR, exist_ok=True)

    cell_text = build_table_cells()
    n_rows = len(cell_text)
    n_cols = len(HEADERS)

    row_lines = estimate_row_heights(cell_text)
    line_inches   = 0.18
    header_inches = 0.40
    total_h = header_inches + sum(r * line_inches for r in row_lines) + 0.40

    fig_w = 8.5
    fig, ax = plt.subplots(figsize=(fig_w, total_h))
    ax.axis("off")

    table = ax.table(cellText=cell_text,
                     colLabels=list(HEADERS),
                     cellLoc="left",
                     loc="upper center",
                     colWidths=COL_WIDTHS)
    table.auto_set_font_size(False)
    table.set_fontsize(BODY_FONTSIZE)

    base_h = 1.0 / (n_rows + 1)
    for col in range(n_cols):
        table[(0, col)].set_height(base_h * 1.0)
    for i, lines in enumerate(row_lines, start=1):
        h = base_h * (lines * 0.80 + 0.4)
        for col in range(n_cols):
            table[(i, col)].set_height(h)

    # Header style
    for col in range(n_cols):
        cell = table[(0, col)]
        cell.set_facecolor(HEADER_COLOUR)
        cell.set_text_props(weight="bold", color=HEADER_TEXTCOLOUR,
                            size=HEADER_FONTSIZE, ha="left")
        cell.set_edgecolor("white")
        cell.set_linewidth(0.5)
        cell.PAD = 0.04

    # Body styling: simple alternating stripes, no section dividers.
    for i in range(1, n_rows + 1):
        for col in range(n_cols):
            cell = table[(i, col)]
            cell.set_facecolor(ROW_COLOURS[(i - 1) % 2])
            cell.set_text_props(color=BODY_TEXTCOLOUR, ha="left",
                                va="center")
            cell.set_edgecolor(GRID_COLOUR)
            cell.set_linewidth(0.4)
            cell.PAD = 0.04

    fig.savefig(OUT_PDF, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {OUT_PDF}")


if __name__ == "__main__":
    render()
