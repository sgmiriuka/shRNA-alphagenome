from __future__ import annotations

import argparse
import json
from pathlib import Path
import os
from typing import Dict, List, Tuple

import jinja2
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from .config import AppConfig
from . import transcript_bed as tb


def _load_arrays_npz(path: Path):
    """Load numpy arrays from .npz file for plotting.

    Args:
        path: Path to the .npz file.

    Returns:
        dict or None: Dictionary with 'x', 'ref', 'alt', 'delta', 'quantiles' arrays, or None if file doesn't exist.
    """
    if not path.exists():
        return None
    data = np.load(path, allow_pickle=True)
    return dict(x=data.get('x'), ref=data.get('ref'), alt=data.get('alt'), delta=data.get('delta'), quantiles=data.get('quantiles'))


def _plot_track_panel(values: np.ndarray, x: np.ndarray | None, title: str, color: str, out_path: Path):
    """Plot a single track panel (REF, ALT, or DELTA).

    Args:
        values: Array of values to plot.
        x: X-coordinates (genomic positions). If None, uses range(len(values)).
        title: Plot title.
        color: Line color.
        out_path: Path to save the PNG.
    """
    fig = plt.figure(figsize=(8, 3))
    # Use constrained layout for more reliable spacing
    try:
        fig.set_constrained_layout(True)
    except Exception:
        pass
    if x is None:
        x = np.arange(len(values))
    plt.plot(x, values, color=color, alpha=0.9)
    # mark nominal center of window
    if len(x):
        cx = x[len(x)//2]
        plt.axvline(cx, color='k', linestyle='--', linewidth=0.8)
    plt.title(title)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_path, dpi=150)
    plt.close(fig)


def _plot_track_three(ref: np.ndarray, alt: np.ndarray, x: np.ndarray | None, base_title: str, base_out_prefix: Path):
    # REF
    _plot_track_panel(ref, x, f"{base_title} REF", color="#1f77b4", out_path=base_out_prefix.with_name(base_out_prefix.name + "_ref.png"))
    # ALT
    _plot_track_panel(alt, x, f"{base_title} ALT", color="#ff7f0e", out_path=base_out_prefix.with_name(base_out_prefix.name + "_alt.png"))
    # Δ = ALT - REF
    delta = alt - ref
    _plot_track_panel(delta, x, f"{base_title} Δ (ALT−REF)", color="#2ca02c", out_path=base_out_prefix.with_name(base_out_prefix.name + "_delta.png"))


def _plot_track_three_panel(ref: np.ndarray, alt: np.ndarray, x: np.ndarray | None, base_title: str, panel_out_path: Path):
    delta = alt - ref
    fig, axes = plt.subplots(3, 1, figsize=(10, 8), sharex=True)
    try:
        fig.set_constrained_layout(True)
    except Exception:
        pass
    if x is None:
        x = np.arange(len(ref))
    # REF
    axes[0].plot(x, ref, color="#1f77b4", alpha=0.9)
    axes[0].set_title(f"{base_title} REF")
    # Δ
    axes[1].plot(x, delta, color="#2ca02c", alpha=0.9)
    axes[1].set_title(f"{base_title} Δ (ALT−REF)")
    # ALT
    axes[2].plot(x, alt, color="#ff7f0e", alpha=0.9)
    axes[2].set_title(f"{base_title} ALT")
    # center line
    if len(x):
        cx = x[len(x)//2]
        for ax in axes:
            ax.axvline(cx, color='k', linestyle='--', linewidth=0.8)
    panel_out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(panel_out_path, dpi=150)
    plt.close(fig)


def _plot_empty_track_three(base_title: str, base_out_prefix: Path):
    # Create placeholder zero series so images exist even when data missing
    x = np.arange(100)
    zeros = np.zeros_like(x, dtype=float)
    _plot_track_three(zeros, zeros, x, base_title, base_out_prefix)


def _plot_empty_track_three_panel(base_title: str, panel_out_path: Path):
    x = np.arange(100)
    zeros = np.zeros_like(x, dtype=float)
    _plot_track_three_panel(zeros, zeros, x, base_title, panel_out_path)


def _plot_sashimi_three_from_json(junc_json_path: Path, base_title: str, base_out_prefix: Path):
    try:
        from alphagenome.visualization import plot as ag_plot
        from alphagenome.data import genome
    except Exception:
        return
    data = json.loads(junc_json_path.read_text())
    ref = data.get('ref', [])
    alt = data.get('alt', [])
    interval = data.get('interval', None)

    def to_junctions(seq):
        out = []
        for start, end, strand, k in seq:
            out.append(genome.Junction(chromosome='chr?', start=int(start), end=int(end), strand=strand, k=max(float(k), 0.0)))
        return out

    ref_js = to_junctions(ref)
    alt_js = to_junctions(alt)

    # Prepare delta junctions: |ALT-REF| thickness
    def to_map(seq):
        m = {}
        for start, end, strand, k in seq:
            m[(int(start), int(end), strand)] = float(k)
        return m

    m_ref = to_map(ref)
    m_alt = to_map(alt)
    keys = set(m_ref.keys()) | set(m_alt.keys())
    delta_j = []
    for key in keys:
        start, end, strand = key
        kval = abs(m_alt.get(key, 0.0) - m_ref.get(key, 0.0))
        delta_j.append(genome.Junction(chromosome='chr?', start=start, end=end, strand=strand, k=kval))

    # Plot separate PNGs
    def _save(ax_plot, js, suffix, title):
        fig, ax = plt.subplots(1, 1, figsize=(10, 3))
        try:
            fig.set_constrained_layout(True)
        except Exception:
            pass
        ag_plot.sashimi_plot(js, ax=ax)
        if interval:
            ax.set_xlim([interval[0], interval[1]])
        ax.set_title(title)
        out = base_out_prefix.with_name(base_out_prefix.name + f"_{suffix}.png")
        out.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out, dpi=150)
        plt.close(fig)

    _save(ag_plot.sashimi_plot, ref_js, "ref", f"{base_title} REF")
    _save(ag_plot.sashimi_plot, alt_js, "alt", f"{base_title} ALT")
    _save(ag_plot.sashimi_plot, delta_j, "delta", f"{base_title} Δ (|ALT−REF|)")


def _plot_sashimi_three_panel_from_json(junc_json_path: Path, base_title: str, panel_out_path: Path):
    try:
        from alphagenome.visualization import plot as ag_plot
        from alphagenome.data import genome
    except Exception:
        return
    data = json.loads(junc_json_path.read_text())
    ref = data.get('ref', [])
    alt = data.get('alt', [])
    interval = data.get('interval', None)

    def to_map(seq):
        m = {}
        for start, end, strand, k in seq:
            m[(int(start), int(end), strand)] = float(k)
        return m
    m_ref = to_map(ref)
    m_alt = to_map(alt)
    keys = set(m_ref.keys()) | set(m_alt.keys())

    def to_junctions(seq):
        from alphagenome.data import genome as g
        return [g.Junction(chromosome='chr?', start=int(s), end=int(e), strand=st, k=max(float(k), 0.0)) for (s,e,st,k) in seq]

    ref_js = to_junctions([(s,e,st,m_ref.get((s,e,st),0.0)) for (s,e,st) in keys])
    alt_js = to_junctions([(s,e,st,m_alt.get((s,e,st),0.0)) for (s,e,st) in keys])
    delta_js = to_junctions([(s,e,st,abs(m_alt.get((s,e,st),0.0)-m_ref.get((s,e,st),0.0))) for (s,e,st) in keys])

    fig, axes = plt.subplots(3, 1, figsize=(10, 8), sharex=True)
    try:
        fig.set_constrained_layout(True)
    except Exception:
        pass
    ag_plot.sashimi_plot(ref_js, ax=axes[0])
    axes[0].set_title(f"{base_title} REF")
    ag_plot.sashimi_plot(delta_js, ax=axes[1])
    axes[1].set_title(f"{base_title} Δ (|ALT−REF|)")
    ag_plot.sashimi_plot(alt_js, ax=axes[2])
    axes[2].set_title(f"{base_title} ALT")
    if interval:
        for ax in axes:
            ax.set_xlim([interval[0], interval[1]])
    panel_out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(panel_out_path, dpi=150)
    plt.close(fig)


def _render_html(scores_csv: Path, plots_root: Path, out_html: Path):
    """Generate HTML report with ranked candidates and plots.

    Args:
        scores_csv: Path to ranked scores CSV from Ranker.
        plots_root: Directory containing per-candidate plot images.
        out_html: Path to write the HTML report.
    """
    df = pd.read_csv(scores_csv)
    df_sorted = df.sort_values('rank').copy()

    # Compute which images exist to avoid broken links
    def has(path: Path) -> bool:
        try:
            return path.exists()
        except Exception:
            return False

    exist_cols = {
        'has_splicing_panel': [],
        'has_rna_panel': [],
        'has_cage_panel': [],
        'has_procap_panel': [],
        'has_gene_bar': [],
        'has_tf_panel': [],
        'has_histone_panel': [],
    }
    for _, r in df_sorted.iterrows():
        rank = int(r['rank'])
        cand = r['candidate_id']
        exist_cols['has_splicing_panel'].append(has(plots_root / f"{rank}_{cand}_splicing_panel.png"))
        exist_cols['has_rna_panel'].append(has(plots_root / f"{rank}_{cand}_rna_panel.png"))
        exist_cols['has_cage_panel'].append(has(plots_root / f"{rank}_{cand}_cage_panel.png"))
        exist_cols['has_procap_panel'].append(has(plots_root / f"{rank}_{cand}_procap_panel.png"))
        exist_cols['has_gene_bar'].append(has(plots_root / f"{rank}_{cand}_rna_gene_lfc_bar.png"))
        exist_cols['has_tf_panel'].append(has(plots_root / f"{rank}_{cand}_tf_panel.png"))
        exist_cols['has_histone_panel'].append(has(plots_root / f"{rank}_{cand}_histone_panel.png"))
    for k, v in exist_cols.items():
        df_sorted[k] = v

    # Global gene spikes figure (if present)
    gene_spikes_rel = None
    try:
        p = plots_root / "gene_structure_spikes.png"
        if p.exists():
            gene_spikes_rel = os.path.relpath(str(p), start=str(out_html.parent))
    except Exception:
        gene_spikes_rel = None

    env = jinja2.Environment(autoescape=True)
    tmpl = env.from_string(
        """
<!DOCTYPE html>
<html>
  <head>
    <meta charset="utf-8" />
    <title>AlphaGenome Intron-2 shRNA Candidates</title>
    <style>
      body { font-family: Arial, sans-serif; margin: 1.5rem; }
      .card { border: 1px solid #ddd; padding: 1rem; margin-bottom: 1rem; }
      .row { display: flex; gap: 1rem; flex-wrap: wrap; }
      .thumb { max-width: 380px; border: 1px solid #eee; }
      .caption { font-size: 12px; color: #444; margin-top: .25rem; }
      h2 { margin: 0 0 .5rem 0; }
      small { color: #666; }
    </style>
  </head>
  <body>
    <h1>AlphaGenome Intron-2 shRNA Candidate Report</h1>
    {% if gene_spikes_img %}
    <div class="card">
      <h2>Gene Structure and Candidate Spikes</h2>
      <div class="row">
        <div>
          <img class="thumb" src="{{ gene_spikes_img }}" alt="Gene structure with intron highlight and candidate spikes" />
          <div class="caption">Exons (boxes), introns (lines). Highlighted intron shows candidate spikes colored green (best) → red (worst) by composite rank.</div>
        </div>
      </div>
    </div>
    {% endif %}
    {% for _, r in rows.iterrows() %}
    <div class="card">
      <h2>#{{ '%02d'|format(r['rank']) }} — {{ r['candidate_id'] }}</h2>
      <small>Position: {{ r['pos'] if (r.get('pos') == r.get('pos')) else 'NA' }} | Composite: {{ '%.4f'|format(r['composite']) }} | Splicing: {{ '%.4f'|format(r.get('splicing', 0.0) or 0.0) }} | RNA(local): {{ '%.4f'|format(r.get('rna', 0.0) or 0.0) }} | RNA(gene): {{ '%.4f'|format(r.get('rna_gene', 0.0) or 0.0) }} | TSS: {{ '%.4f'|format(r.get('tss', 0.0) or 0.0) }} | TF: {{ '%.4f'|format(r.get('tf', 0.0) or 0.0) }} | Histone: {{ '%.4f'|format(r.get('histone', 0.0) or 0.0) }}</small>
      <div class="row">
        {% if r['has_splicing_panel'] %}
        <div>
          <img class="thumb" src="{{ plots_root }}/{{ r['rank'] }}_{{ r['candidate_id'] }}_splicing_panel.png" alt="Splicing sashimi (REF / |ALT−REF| / ALT)" />
          <div class="caption">Splicing sashimi (REF / |ALT−REF| / ALT)</div>
        </div>
        {% endif %}
        {% if r['has_rna_panel'] %}
        <div>
          <img class="thumb" src="{{ plots_root }}/{{ r['rank'] }}_{{ r['candidate_id'] }}_rna_panel.png" alt="RNA-seq (REF / Δ / ALT)" />
          <div class="caption">RNA-seq track (REF / Δ / ALT)</div>
        </div>
        {% endif %}
      </div>
      <div class="row">
        {% if r['has_cage_panel'] %}
        <div>
          <img class="thumb" src="{{ plots_root }}/{{ r['rank'] }}_{{ r['candidate_id'] }}_cage_panel.png" alt="CAGE (REF / Δ / ALT)" />
          <div class="caption">CAGE (TSS proxy) (REF / Δ / ALT)</div>
        </div>
        {% endif %}
        {% if r['has_procap_panel'] %}
        <div>
          <img class="thumb" src="{{ plots_root }}/{{ r['rank'] }}_{{ r['candidate_id'] }}_procap_panel.png" alt="PROCAP (REF / Δ / ALT)" />
          <div class="caption">PROCAP (TSS proxy) (REF / Δ / ALT)</div>
        </div>
        {% endif %}
        {% if r['has_gene_bar'] %}
        <div>
          <img class="thumb" src="{{ plots_root }}/{{ r['rank'] }}_{{ r['candidate_id'] }}_rna_gene_lfc_bar.png" alt="RNA gene-level LFC per track" />
          <div class="caption">RNA gene-level LFC (per track)</div>
        </div>
        {% endif %}
        {% if r['has_tf_panel'] %}
        <div>
          <img class="thumb" src="{{ plots_root }}/{{ r['rank'] }}_{{ r['candidate_id'] }}_tf_panel.png" alt="TF binding (REF / Δ / ALT)" />
          <div class="caption">TF binding (REF / Δ / ALT)</div>
        </div>
        {% endif %}
        {% if r['has_histone_panel'] %}
        <div>
          <img class="thumb" src="{{ plots_root }}/{{ r['rank'] }}_{{ r['candidate_id'] }}_histone_panel.png" alt="Histone marks (REF / Δ / ALT)" />
          <div class="caption">Histone marks (REF / Δ / ALT)</div>
        </div>
        {% endif %}
      </div>
    </div>
    {% endfor %}
  </body>
</html>
        """
    )
    try:
        plots_root_rel = os.path.relpath(str(plots_root), start=str(out_html.parent))
    except Exception:
        plots_root_rel = str(plots_root)
    html = tmpl.render(rows=df_sorted, plots_root=plots_root_rel, gene_spikes_img=gene_spikes_rel)
    out_html.parent.mkdir(parents=True, exist_ok=True)
    out_html.write_text(html)


def _plot_gene_structure_with_spikes(
    *,
    exons: List[dict],
    chrom: str,
    strand: int,
    intron_start0: int,
    intron_end0: int,
    ranked: pd.DataFrame,
    out_path: Path,
    transcript_id: str | None = None,
):
    """Plot gene structure with candidate insertion spikes.

    Args:
        exons: List of exon dictionaries with 'start', 'end'.
        chrom: Chromosome name.
        strand: Strand (1 for +, -1 for -).
        intron_start0: 0-based intron start.
        intron_end0: 0-based intron end.
        ranked: DataFrame of ranked candidates.
        out_path: Path to save the plot.
        transcript_id: Optional transcript ID for title.
    """
    # Prepare coordinates (use 1-based inclusive for display)
    exon_intervals = [(int(e.get('start')), int(e.get('end'))) for e in exons]
    # Sort by genomic coordinate for plotting left→right
    exon_intervals.sort(key=lambda t: t[0])
    gene_start = min(s for s, _ in exon_intervals)
    gene_end = max(e for _, e in exon_intervals)

    intron_start1 = int(intron_start0) + 1
    intron_end1 = int(intron_end0)

    # Candidate positions (1-based) limited to the highlighted intron region if available
    poses = []
    ranks = []
    comps = []
    if 'pos' in ranked.columns:
        for _, r in ranked.iterrows():
            try:
                pos = int(r['pos'])
            except Exception:
                continue
            if pos >= intron_start1 and pos <= intron_end1:
                poses.append(pos)
                ranks.append(int(r['rank']))
                try:
                    comps.append(float(r.get('composite', np.nan)))
                except Exception:
                    comps.append(np.nan)

    # Color by rank (green=best → red=worst)
    cmap = plt.get_cmap('RdYlGn')
    if ranks:
        rmin = min(ranks)
        rmax = max(ranks)
    else:
        rmin = 1
        rmax = 1

    def rank_to_color(rk: int):
        if rmax == rmin:
            val = 1.0
        else:
            # Normalize so best (lowest rank) → 1.0 (green), worst → 0.0 (red)
            val = 1.0 - (rk - rmin) / (rmax - rmin)
        return cmap(val)

    # Single-panel zoomed intron with spikes
    fig, ax = plt.subplots(1, 1, figsize=(12, 3.8))
    try:
        fig.set_constrained_layout(True)
    except Exception:
        pass

    # Utilities
    def draw_gene_track(ax_, x0, x1, exons_to_draw, *, y=0.05, exon_h=0.10, color_line="#666", color_exon="#333"):
        ax_.hlines(y, x0, x1, color=color_line, linewidth=2, zorder=1)
        for s, e in exons_to_draw:
            s2 = max(s, x0)
            e2 = min(e, x1)
            if e2 > s2:
                ax_.add_patch(plt.Rectangle((s2, y - exon_h / 2), max(e2 - s2, 1), exon_h,
                                            facecolor=color_exon, edgecolor='none', zorder=2))

    def intersect(a0, a1, b0, b1):
        return max(a0, b0) < min(a1, b1)

    # Zoom window around intron (±25% intron width)
    intr_w = max(intron_end1 - intron_start1, 10)
    pad = max(int(0.25 * intr_w), 50)
    zoom_x0 = intron_start1 - pad
    zoom_x1 = intron_end1 + pad

    draw_gene_track(ax, zoom_x0, zoom_x1, [t for t in exon_intervals if intersect(t[0], t[1], zoom_x0, zoom_x1)],
                    y=0.05, exon_h=0.10, color_line="#666", color_exon="#333")

    # Highlight intron
    ax.add_patch(
        plt.Rectangle((intron_start1, 0.0), max(intron_end1 - intron_start1, 1), 0.15,
                      facecolor='#1f77b4', alpha=0.20, edgecolor='#1f77b4', zorder=0)
    )

    # Candidate spikes (height scaled by composite; color by rank)
    if comps:
        comp_vals = np.array([c if np.isfinite(c) else np.nan for c in comps], dtype=float)
        if np.all(~np.isfinite(comp_vals)):
            comp_vals = np.full(len(poses), np.nan)
        cmin = np.nanmin(comp_vals) if np.any(np.isfinite(comp_vals)) else None
        cmax = np.nanmax(comp_vals) if np.any(np.isfinite(comp_vals)) else None
    else:
        comp_vals = np.array([])
        cmin = cmax = None

    base_y0 = 0.20
    min_h = 0.15  # worst → short spike
    max_h = 0.75  # best → tall spike
    for i, (pos, rk) in enumerate(zip(poses, ranks)):
        col = rank_to_color(rk)
        if cmin is not None and cmax is not None and np.isfinite(comp_vals[i]):
            if cmax == cmin:
                norm = 1.0
            else:
                norm = 1.0 - (comp_vals[i] - cmin) / (cmax - cmin)
        else:
            norm = 1.0 if len(ranks) <= 1 else 1.0 - (rk - min(ranks)) / (max(ranks) - min(ranks) + 1e-9)
        height = min_h + (max_h - min_h) * float(np.clip(norm, 0.0, 1.0))
        ax.vlines(pos, base_y0, base_y0 + height, colors=[col], linewidth=3, zorder=3)

    # Legend-like guide for color scale (min/max rank)
    if ranks:
        rmin, rmax = min(ranks), max(ranks)
        ax.text(0.01, 0.97, f"Best rank {rmin}", color=rank_to_color(rmin), transform=ax.transAxes,
                ha='left', va='top', fontsize=10)
        ax.text(0.25, 0.97, f"Worst rank {rmax}", color=rank_to_color(rmax), transform=ax.transAxes,
                ha='left', va='top', fontsize=10)

    # Aesthetics
    ax.set_xlim(zoom_x0, zoom_x1)
    ax.set_ylim(0.0, 1.05)
    ax.set_xlabel("Genomic position (1-based)")
    ax.set_yticks([])
    title_bits = ["Intron region candidate spikes"]
    if transcript_id:
        title_bits.append(transcript_id)
    title_bits.append(f"{chrom} ({'+' if strand==1 else '-'})")
    ax.set_title(" — ".join(title_bits))

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def main(argv: List[str] | None = None) -> None:
    """Generate plots and HTML report for ranked candidates.

    Reads ranked scores and predictions, generates per-candidate plots,
    and creates an HTML summary report.

    Args:
        argv: Command-line arguments. If None, uses sys.argv.
    """
    ap = argparse.ArgumentParser(
        prog="Reporter",
        description="Plot every prediction and build an HTML summary.",
    )
    ap.add_argument("--config", required=False, type=Path)
    ap.add_argument("--scores", required=False, type=Path)
    ap.add_argument("--pred", required=False, type=Path)
    ap.add_argument("--plots", required=False, type=Path)
    ap.add_argument("--html", required=False, type=Path)
    # Optional: supply transcript (or GTF) to render gene structure
    ap.add_argument("--transcript", required=False, type=str)
    ap.add_argument("--gtf", required=False, type=Path)
    ap.add_argument("--intron-bed", required=False, type=Path)
    args = ap.parse_args(argv)

    cfg = AppConfig.from_yaml(args.config) if args.config else AppConfig()
    scores = Path(args.scores) if args.scores else cfg.io.scores_csv
    pred = Path(args.pred) if args.pred else cfg.io.raw_parquet
    plots = Path(args.plots) if args.plots else cfg.io.plots_dir
    html = Path(args.html) if args.html else cfg.io.report_html

    raw_df = pd.read_parquet(pred)
    ranked = pd.read_csv(scores)

    # Build quick lookup per-candidate
    by_cand = raw_df.groupby('candidate_id')

    for _, row in ranked.iterrows():
        cand = row['candidate_id']
        rank = int(row['rank'])
        plots_dir = Path(args.plots)
        plots_dir.mkdir(parents=True, exist_ok=True)

        # Splicing sashimi (separate REF / ALT / Δ files) if we stored json
        j_json = pred.parent / 'arrays' / f"{cand}_splicing_junctions.json"
        if j_json.exists():
            _plot_sashimi_three_from_json(j_json, base_title=f"{cand} splicing", base_out_prefix=plots_dir / f"{rank}_{cand}_splicing")
            # Panel version
            _plot_sashimi_three_panel_from_json(j_json, base_title=f"{cand} splicing", panel_out_path=plots_dir / f"{rank}_{cand}_splicing_panel.png")

        # RNA track (separate REF / ALT / Δ files)
        rna_npz = pred.parent / 'arrays' / f"{cand}_rna.npz"
        if rna_npz.exists():
            arrays = _load_arrays_npz(rna_npz)
            if arrays and arrays.get('ref') is not None and arrays.get('alt') is not None:
                _plot_track_three(arrays['ref'], arrays['alt'], arrays.get('x'), base_title=f"{cand} RNA", base_out_prefix=plots_dir / f"{rank}_{cand}_rna")
                # Panel version
                _plot_track_three_panel(arrays['ref'], arrays['alt'], arrays.get('x'), base_title=f"{cand} RNA", panel_out_path=plots_dir / f"{rank}_{cand}_rna_panel.png")
            else:
                # Generate placeholders if arrays malformed
                _plot_empty_track_three(base_title=f"{cand} RNA", base_out_prefix=plots_dir / f"{rank}_{cand}_rna")
                _plot_empty_track_three_panel(base_title=f"{cand} RNA", panel_out_path=plots_dir / f"{rank}_{cand}_rna_panel.png")
        else:
            # Generate placeholders if arrays missing
            _plot_empty_track_three(base_title=f"{cand} RNA", base_out_prefix=plots_dir / f"{rank}_{cand}_rna")
            _plot_empty_track_three_panel(base_title=f"{cand} RNA", panel_out_path=plots_dir / f"{rank}_{cand}_rna_panel.png")

        # Transcription initiation (CAGE / PROCAP) tracks
        for label in ("cage", "procap"):
            npz = pred.parent / 'arrays' / f"{cand}_{label}.npz"
            if npz.exists():
                arrays = _load_arrays_npz(npz)
                if arrays and arrays.get('ref') is not None and arrays.get('alt') is not None:
                    base = plots_dir / f"{rank}_{cand}_{label}"
                    title = f"{cand} {label.upper()}"
                    _plot_track_three(arrays['ref'], arrays['alt'], arrays.get('x'), base_title=title, base_out_prefix=base)
                    _plot_track_three_panel(arrays['ref'], arrays['alt'], arrays.get('x'), base_title=title, panel_out_path=base.with_name(base.name + "_panel.png"))

        # TF binding and Histone panels if arrays were saved
        for label in ("tf", "histone"):
            npz = pred.parent / 'arrays' / f"{cand}_{label}.npz"
            if npz.exists():
                arrays = _load_arrays_npz(npz)
                if arrays and arrays.get('ref') is not None and arrays.get('alt') is not None:
                    base = plots_dir / f"{rank}_{cand}_{label}"
                    title = f"{cand} {label.upper()}"
                    _plot_track_three(arrays['ref'], arrays['alt'], arrays.get('x'), base_title=title, base_out_prefix=base)
                    _plot_track_three_panel(arrays['ref'], arrays['alt'], arrays.get('x'), base_title=title, panel_out_path=base.with_name(base.name + "_panel.png"))

        # Gene-level RNA LFC bar
        gene_npz = pred.parent / 'arrays' / f"{cand}_rna_gene_lfc.npz"
        if gene_npz.exists():
            try:
                data = np.load(gene_npz, allow_pickle=True)
                vals = data.get('scores')
                names = data.get('track_names')
                if vals is not None:
                    vals_arr = np.asarray(vals, dtype=float).ravel()
                    # Rank by absolute value; treat NaN as 0 for ranking
                    order_base = np.nan_to_num(np.abs(vals_arr), nan=0.0)
                    if order_base.size == 0:
                        # Nothing to plot
                        continue
                    idx = np.argsort(-order_base)[:12]
                    vals_top = vals_arr[idx]
                    if names is not None and len(names) == len(vals_arr):
                        names_arr = np.asarray(names, dtype=object).ravel()
                        names_top = names_arr[idx]
                    else:
                        names_top = [f'track_{i}' for i in idx]
                    fig, ax = plt.subplots(1, 1, figsize=(10, 3.6))
                    # Use manual spacing to avoid layout engine warnings with rotated labels
                    fig.subplots_adjust(left=0.08, right=0.98, top=0.88, bottom=0.35)
                    ax.bar(range(len(vals_top)), vals_top, color="#9467bd")
                    ax.set_title(f"{cand} Gene-level RNA LFC (per track)")
                    ax.set_xticks(range(len(vals_top)))
                    ax.set_xticklabels([str(n) for n in names_top], rotation=45, ha='right', fontsize=8)
                    outp = plots_dir / f"{rank}_{cand}_rna_gene_lfc_bar.png"
                    fig.savefig(outp, dpi=150)
                    plt.close(fig)
            except Exception:
                pass

    # Optional global gene-structure + spikes plot
    try:
        intron_bed_path = Path(args.intron_bed) if args.intron_bed else (cfg.inputs.intron_bed if cfg.inputs.intron_bed else None)
    except Exception:
        intron_bed_path = None
    transcript_id = args.transcript if getattr(args, 'transcript', None) else None
    gtf_path = Path(args.gtf) if getattr(args, 'gtf', None) else None

    # Try to infer transcript_id from intron_bed_path if not provided
    if transcript_id is None and intron_bed_path is not None:
        import re
        bed_name = intron_bed_path.name
        # Look for ENST pattern in filename
        match = re.search(r'(ENST\d+)', bed_name)
        if match:
            transcript_id = match.group(1)
            print(f"[Reporter] Inferred transcript ID from BED filename: {transcript_id}")

    if transcript_id is not None and intron_bed_path is not None:
        try:
            # Fetch exon structures
            if gtf_path is not None:
                chrom, strand, exons = tb._fetch_exons_from_gtf(gtf_path, transcript_id)
            else:
                chrom, strand, exons = tb._fetch_exons_rest_ensembl(transcript_id)
            # Read intron BED
            with open(intron_bed_path, 'r') as f:
                line = next(l for l in f if l.strip() and not l.startswith('#'))
            parts = line.rstrip().split('\t')
            intron_start0 = int(parts[1])
            intron_end0 = int(parts[2])
            _plot_gene_structure_with_spikes(
                exons=exons,
                chrom=(parts[0] if parts and parts[0] else chrom or ''),
                strand=int(strand),
                intron_start0=intron_start0,
                intron_end0=intron_end0,
                ranked=ranked,
                out_path=plots / "gene_structure_spikes.png",
                transcript_id=transcript_id,
            )
        except Exception as e:
            print(f"[Reporter] Failed to create gene_structure_spikes.png: {e}")
    elif intron_bed_path is not None and transcript_id is None:
        print("[Reporter] To create gene_structure_spikes.png, provide --transcript (e.g., ENST00000325495) or ensure the BED filename contains the transcript ID (e.g., gene_ENST00000325495_intron2_hg38.bed).")

    if html:
        _render_html(scores, plots, html)
        print(f"[Reporter] Wrote HTML report to {html}")


if __name__ == "__main__":
    main()
