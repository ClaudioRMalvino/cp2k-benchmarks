#!/usr/bin/env python3
"""Render tabular .csv files under plots/maqao_plots/ as styled PDFs."""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

TABLES_DIR = "/home/crm98/cp2k-benchmarks/plots/maqao_plots"

HEADER_COLOUR     = "#133844"
HEADER_TEXTCOLOUR = "#FFFFFF"
ROW_COLOURS       = ["#ECEEF1", "#FFFFFF"]
GRID_COLOUR       = "#B5BDC8"
BODY_TEXTCOLOUR   = "#232830"
GRID_LINEWIDTH    = 0.4
TITLE_FONTSIZE    = 12
HEADER_FONTSIZE   = 10
BODY_FONTSIZE     = 9
FOOTNOTE_FONTSIZE = 8

plt.rcParams.update({
    "font.family":         "serif",
    "font.serif":          ["DejaVu Serif", "Times New Roman", "Times"],
    "mathtext.fontset":    "cm",
    "axes.labelcolor":     BODY_TEXTCOLOUR,
    "text.color":          BODY_TEXTCOLOUR,
})


def render_table_pdf(df, output_path, title="", footnote="",
                     col_widths=None, fig_size=None, row_height=1.00,
                     body_fontsize=None, header_fontsize=None, title_pad=4,
                     cell_styles=None, stripe_group_size=1):
    """Render a DataFrame as a styled PDF table."""
    n_rows, n_cols = df.shape
    if fig_size is None:
        w = max(8.0, 1.6 * n_cols)
        h = 1.0 + 0.36 * n_rows + (0.5 if footnote else 0.0) + (0.4 if title else 0.0)
        fig_size = (w, h)

    fig, ax = plt.subplots(figsize=fig_size)
    ax.axis("off")
    if title:
        ax.set_title(title, fontsize=TITLE_FONTSIZE, fontweight="bold",
                     loc="center", pad=title_pad)

    body_fs   = body_fontsize   if body_fontsize   is not None else BODY_FONTSIZE
    header_fs = header_fontsize if header_fontsize is not None else HEADER_FONTSIZE
    cell_text = df.astype(str).values.tolist()
    # loc="upper center" anchors the table near the title (avoids floating
    # in the middle of the figure on tall multi-row tables).
    table = ax.table(cellText=cell_text, colLabels=df.columns.tolist(),
                     cellLoc="center", loc="upper center", colWidths=col_widths)
    table.auto_set_font_size(False)
    table.set_fontsize(body_fs)
    table.scale(1, row_height * 2.5)

    for j in range(n_cols):
        cell = table[(0, j)]
        cell.set_facecolor(HEADER_COLOUR)
        cell.set_text_props(weight="bold", color=HEADER_TEXTCOLOUR,
                            size=header_fs)
        cell.set_edgecolor("white"); cell.set_linewidth(GRID_LINEWIDTH)
    for i in range(1, n_rows + 1):
        group_idx = (i - 1) // stripe_group_size
        face = ROW_COLOURS[group_idx % 2]
        for j in range(n_cols):
            cell = table[(i, j)]
            cell.set_facecolor(face)
            cell.set_text_props(color=BODY_TEXTCOLOUR)
            cell.set_edgecolor(GRID_COLOUR); cell.set_linewidth(GRID_LINEWIDTH)

    if cell_styles:
        for (i, j), style in cell_styles.items():
            cell = table[(i, j)]
            facecolor     = style.pop("facecolor", None)
            edgecolor     = style.pop("edgecolor", None)
            visible_edges = style.pop("visible_edges", None)
            if facecolor     is not None: cell.set_facecolor(facecolor)
            if edgecolor     is not None: cell.set_edgecolor(edgecolor)
            if visible_edges is not None: cell.visible_edges = visible_edges
            if style: cell.set_text_props(**style)

    if footnote:
        fig.text(0.5, 0.01, footnote, ha="center", va="bottom",
                 fontsize=FOOTNOTE_FONTSIZE, style="italic", wrap=True)

    fig.tight_layout(rect=[0, 0.04 if footnote else 0, 1, 1])
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  saved {output_path}")


def table_global_metrics():
    print("\n[1/2] maqao_global_metrics_table")
    src = os.path.join(TABLES_DIR, "maqao_global_metrics_table.csv")
    df = pd.read_csv(src)
    df = df.rename(columns={
        "branch":                       "Branch",
        "application_time_s":           "App. time (s)",
        "array_access_efficiency_pct":  "Array-access\nefficiency (%)",
        "vec_headroom_speedup":         "Vec. headroom\nspeedup",
        "runtime_overhead_speedup":     "MPI/OMP\noverhead speedup",
        "user_time_pct":                "User-code\ntime (%)",
    })
    df["App. time (s)"]                  = df["App. time (s)"].apply(lambda v: f"{v:.2f}")
    df["Array-access\nefficiency (%)"]   = df["Array-access\nefficiency (%)"].apply(lambda v: f"{v:.2f}")
    df["Vec. headroom\nspeedup"]         = df["Vec. headroom\nspeedup"].apply(lambda v: f"{v:.2f}x")
    df["MPI/OMP\noverhead speedup"]      = df["MPI/OMP\noverhead speedup"].apply(lambda v: f"{v:.2f}x")
    df["User-code\ntime (%)"]            = df["User-code\ntime (%)"].apply(lambda v: f"{v:.2f}")

    render_table_pdf(
        df,
        os.path.join(TABLES_DIR, "maqao_global_metrics_table.pdf"),
        title=("MAQAO global metrics — N = 64 H$_2$O, 16 cores total "
               "(16×1 MPI/OMP for master & native-spline; 8×2 for omp)"),
        col_widths=[0.30, 0.12, 0.16, 0.14, 0.16, 0.12],
    )


def table_vectorisation():
    print("\n[2/2] maqao_vectorisation_table")
    src = os.path.join(TABLES_DIR, "maqao_vectorisation_table.csv")
    df = pd.read_csv(src)
    df = df.rename(columns={
        "loop":                         "Loop",
        "branch":                       "Branch",
        "coverage_pct":                 "Coverage\n(%)",
        "cqa_cycles_per_iter":          "CQA\ncycles/iter",
        "vec_ratio_pct":                "Vec. ratio\n(%)",
        "speedup_if_fully_vectorised":  "Speedup if\nfully vectorised",
    })
    df["Coverage\n(%)"]      = df["Coverage\n(%)"].apply(lambda v: f"{v:.1f}")
    df["CQA\ncycles/iter"]   = df["CQA\ncycles/iter"].astype(int).astype(str)
    df["Vec. ratio\n(%)"]    = df["Vec. ratio\n(%)"].astype(int).astype(str)

    # Visually merge the Loop column across the 3 branches of each group: the
    # name is placed on the middle row and adjacent cell edges are recoloured
    # to match the stripe so three cells read as one merged block. visible_edges
    # would be simpler but produces triangular corner artifacts in the table
    # backend.
    loop_col_idx = df.columns.get_loc("Loop")
    cell_styles  = {}
    for i in range(len(df)):
        pos = i % 3
        if pos != 1:
            df.iat[i, loop_col_idx] = ""
        group_face = ROW_COLOURS[(i // 3) % 2]
        cell_styles[(i + 1, loop_col_idx)] = dict(edgecolor=group_face)

    render_table_pdf(
        df,
        os.path.join(TABLES_DIR, "maqao_vectorisation_table.pdf"),
        title="MAQAO CQA — per-hot-loop vectorisation analysis (N = 64 H$_2$O, 16 cores)",
        col_widths=[0.28, 0.22, 0.10, 0.13, 0.12, 0.15],
        row_height=0.85,
        body_fontsize=10,
        header_fontsize=10,
        title_pad=10,
        fig_size=(11, 5.0),
        cell_styles=cell_styles,
        stripe_group_size=3,
    )


def main():
    os.makedirs(TABLES_DIR, exist_ok=True)
    table_global_metrics()
    table_vectorisation()
    print(f"\nDone. PDFs written to {TABLES_DIR}")


if __name__ == "__main__":
    main()
