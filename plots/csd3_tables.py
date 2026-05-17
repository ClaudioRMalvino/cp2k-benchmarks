#!/usr/bin/env python3
"""Render tabular .csv files under plots/scaling_csd3/ as thesis-ready PDFs."""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

TABLES_DIR = "/home/crm98/cp2k-benchmarks/plots/scaling_csd3"

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
                     col_widths=None, fig_size=None, row_height=0.55,
                     cell_styles=None, body_fontsize=None,
                     header_fontsize=None, title_pad=4):
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

    cell_text = df.astype(str).values.tolist()
    col_labels = df.columns.tolist()
    table = ax.table(cellText=cell_text, colLabels=col_labels,
                     cellLoc="center", loc="center",
                     colWidths=col_widths)
    body_fs   = body_fontsize   if body_fontsize   is not None else BODY_FONTSIZE
    header_fs = header_fontsize if header_fontsize is not None else HEADER_FONTSIZE
    table.auto_set_font_size(False)
    table.set_fontsize(body_fs)
    table.scale(1, row_height * 2.5)

    for j in range(n_cols):
        cell = table[(0, j)]
        cell.set_facecolor(HEADER_COLOUR)
        cell.set_text_props(weight="bold", color=HEADER_TEXTCOLOUR,
                            size=header_fs)
        cell.set_edgecolor("white")
        cell.set_linewidth(GRID_LINEWIDTH)

    for i in range(1, n_rows + 1):
        for j in range(n_cols):
            cell = table[(i, j)]
            cell.set_facecolor(ROW_COLOURS[(i - 1) % 2])
            cell.set_text_props(color=BODY_TEXTCOLOUR)
            cell.set_edgecolor(GRID_COLOUR)
            cell.set_linewidth(GRID_LINEWIDTH)

    if cell_styles:
        for (i, j), style in cell_styles.items():
            cell = table[(i, j)]
            facecolor = style.pop("facecolor", None)
            if facecolor is not None:
                cell.set_facecolor(facecolor)
            if style:
                cell.set_text_props(**style)

    if footnote:
        fig.text(0.5, 0.01, footnote, ha="center", va="bottom",
                 fontsize=FOOTNOTE_FONTSIZE, style="italic", wrap=True)

    fig.tight_layout(rect=[0, 0.04 if footnote else 0, 1, 1])
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  saved {output_path}")


def table_size_scaling_summary():
    print("\n[1/3] size_scaling_summary_table")
    src = os.path.join(TABLES_DIR, "size_scaling_summary_table.csv")
    df = pd.read_csv(src)
    df = df.rename(columns={
        "branch":                   "Branch",
        "decomposition":            "Decomposition",
        "t_step_mean_s":            "Mean (s/step)",
        "t_step_std_s":             "Std (s)",
        "t_step_ci95_s":            "95% CI",
        "walltime_s_for_100_steps": "Walltime (s)\n100 steps",
        "speedup_over_master":      "Speedup\nvs master",
    })
    render_table_pdf(
        df,
        os.path.join(TABLES_DIR, "size_scaling_summary_table.pdf"),
        title="Headline-workload performance: N = 1024 H$_2$O (3072 atoms), one Peta4-IceLake node",
        col_widths=[0.13, 0.13, 0.13, 0.10, 0.12, 0.18, 0.13],
        row_height=1.0,
    )


def _fmt_sci(v):
    """Format a small float as mathtext scientific notation."""
    if pd.isna(v) or v == 0:
        return "—"
    import math
    exp  = int(math.floor(math.log10(abs(v))))
    mant = v / (10 ** exp)
    return f"${mant:.1f}\\times10^{{{exp}}}$"


def table_statistical_comparison():
    print("\n[2/3] statistical_comparison_N1024")
    src = os.path.join(TABLES_DIR, "statistical_comparison_N1024.csv")
    raw = pd.read_csv(src)

    def _fmt_pair(s):
        a, b = [p.strip() for p in s.replace("upstream master", "master").split("vs")]
        return f"{a}\nvs\n{b}"

    df = pd.DataFrame({
        "Branch pair":      raw["pair"].apply(_fmt_pair),
        "Mean A (s)":       raw["mean_a"].apply(lambda v: f"{v:.4f}"),
        "Mean B (s)":       raw["mean_b"].apply(lambda v: f"{v:.4f}"),
        "Diff (%)":         raw["mean_diff_pct"].apply(lambda v: f"{v:+.1f}"),
        "Speedup\n(A ÷ B)": (raw["mean_a"] / raw["mean_b"]).apply(lambda v: f"{v:.2f}x"),
        "p-value\n(Welch)": raw["welch_p"].apply(_fmt_sci),
    })

    pair_col_idx = df.columns.get_loc("Branch pair")
    cell_styles = {(i + 1, pair_col_idx): dict(size=9) for i in range(len(df))}

    render_table_pdf(
        df,
        os.path.join(TABLES_DIR, "statistical_comparison_N1024.pdf"),
        title=("Pairwise inter-branch significance test at N = 1024 H$_2$O "
               "(5 timed reps per branch)"),
        fig_size=(10, 4.2),
        row_height=2.0,
        cell_styles=cell_styles,
        body_fontsize=12,
        header_fontsize=12,
        title_pad=12,
    )


def table_power_law_fits():
    print("\n[3/3] power_law_fits")
    src = os.path.join(TABLES_DIR, "power_law_fits.csv")
    raw = pd.read_csv(src)
    df = pd.DataFrame({
        "Branch":             raw["branch"],
        "Exponent  b":        [f"{r.b:.3f}  ±  {r.b_err:.3f}" for r in raw.itertuples()],
        "Prefactor  a  (s)":  [f"{r.a:.2e}  ±  {r.a_err:.2e}" for r in raw.itertuples()],
        "n":                  raw["n_points"].astype(int).astype(str),
    })
    render_table_pdf(
        df,
        os.path.join(TABLES_DIR, "power_law_fits.pdf"),
        title=("Power-law fits  t(N) = a · N$^b$  per branch  "
               "(weighted log-log regression; N = number of atoms)"),
        col_widths=[0.22, 0.30, 0.36, 0.12],
    )


def main():
    os.makedirs(TABLES_DIR, exist_ok=True)
    table_size_scaling_summary()
    table_statistical_comparison()
    table_power_law_fits()
    print(f"\nDone. PDFs written to {TABLES_DIR}")


if __name__ == "__main__":
    main()
