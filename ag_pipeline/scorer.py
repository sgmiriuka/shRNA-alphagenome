from __future__ import annotations

import argparse
import json
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import numpy as np
import pandas as pd

from .config import AppConfig
import re


def _load_candidates(candidates_tsv: Path) -> pd.DataFrame:
    df = pd.read_csv(candidates_tsv, sep="\t")
    # Normalize columns
    expect = ["candidate_id", "chrom", "pos", "ref", "alt"]
    missing = [c for c in expect if c not in df.columns]
    if missing:
        raise ValueError(f"Missing columns in candidates TSV: {missing}")
    return df


def _resolve_api_key(cfg: AppConfig) -> str:
    """Resolve API key from environment variables.

    Only treats `cfg.alphagenome.api_key_env` as a variable name if it matches
    a safe ENV var pattern to avoid accidentally echoing secrets.
    """
    candidates = []
    if cfg.alphagenome.api_key_env and re.fullmatch(r"[A-Z][A-Z0-9_]*", cfg.alphagenome.api_key_env):
        candidates.append(cfg.alphagenome.api_key_env)
    # Standard fallbacks
    candidates.extend(["ALPHAGENOME_API_KEY", "ALPHA_GENOME_API_KEY"])
    for name in candidates:
        val = os.getenv(name)
        if val:
            return val
    raise RuntimeError(
        "AlphaGenome API key not found. Export ALPHAGENOME_API_KEY before running, "
        "or set alphagenome.api_key_env in ag.yaml to the name of your env var."
    )


def _build_client(api_key: str, address: Optional[str]):
    # Lazy import to avoid hard dependency unless used
    from alphagenome.models import dna_client

    return dna_client.create(api_key=api_key, address=address or None)


def _ag_output_types(modalities: Iterable[str]):
    from alphagenome.models import dna_output as out
    selected = []
    for m in modalities:
        m = m.strip().lower()
        if m == "splicing":
            selected.extend([out.OutputType.SPLICE_JUNCTIONS, out.OutputType.SPLICE_SITE_USAGE])
        elif m == "rna":
            selected.append(out.OutputType.RNA_SEQ)
        elif m in ("tss", "promoter", "initiation"):
            selected.extend([out.OutputType.CAGE, out.OutputType.PROCAP])
        else:
            raise ValueError(f"Unsupported modality: {m}")
    # Deduplicate preserving order
    return list(dict.fromkeys(selected))


def _coerce_center_mask_width(requested: int) -> int:
    """Return a supported CenterMaskScorer width closest to the request.

    Falls back to the nearest allowed integer width (excluding None).
    """
    from alphagenome.models import variant_scorers as vs
    allowed = list(vs.SUPPORTED_WIDTHS[vs.BaseVariantScorer.CENTER_MASK])
    allowed_ints = sorted([w for w in allowed if isinstance(w, int) and w > 0])
    if requested in allowed_ints:
        return requested
    # Choose nearest width; tie-break by smaller
    best = min(allowed_ints, key=lambda w: (abs(w - requested), w))
    print(f"[AlphaGenomeScorer] Requested variant window {requested} is unsupported; using nearest supported width {best}.")
    return best


def _variant_scorers(modalities: Iterable[str], width: int):
    from alphagenome.models import variant_scorers as vs
    from alphagenome.models import dna_output as out

    width = _coerce_center_mask_width(width)
    scorers = []
    for m in modalities:
        m = m.strip().lower()
        if m == "splicing":
            # CenterMaskScorer does not support SPLICE_JUNCTIONS.
            # Use SPLICE_SITE_USAGE (and optionally SPLICE_SITES) for splicing effects.
            for ot in (out.OutputType.SPLICE_SITE_USAGE, out.OutputType.SPLICE_SITES):
                scorers.append(
                    vs.CenterMaskScorer(
                        requested_output=ot,
                        width=width,
                        aggregation_type=vs.AggregationType.DIFF_MEAN,
                    )
                )
            # Add gene-level splicing scorer (gene mask) and junction-centric scorer.
            scorers.append(
                vs.GeneMaskSplicingScorer(
                    requested_output=out.OutputType.SPLICE_SITE_USAGE,
                    width=None,
                )
            )
            scorers.append(vs.SpliceJunctionScorer())
        elif m == "rna":
            scorers.append(
                vs.CenterMaskScorer(
                    requested_output=out.OutputType.RNA_SEQ,
                    width=width,
                    aggregation_type=vs.AggregationType.DIFF_MEAN,
                )
            )
            # Gene-level RNA log fold change
            scorers.append(
                vs.GeneMaskLFCScorer(
                    requested_output=out.OutputType.RNA_SEQ
                )
            )
        elif m in ("tss", "promoter", "initiation"):
            for ot in (out.OutputType.CAGE, out.OutputType.PROCAP):
                scorers.append(
                    vs.CenterMaskScorer(
                        requested_output=ot,
                        width=width,
                        aggregation_type=vs.AggregationType.DIFF_MEAN,
                    )
                )
        else:
            raise ValueError(f"Unsupported modality for scorer: {m}")
    return scorers


def _interval_around(chrom: str, pos_1based: int, width: int):
    from alphagenome.data import genome
    # Center the interval around the 1-based position.
    # genome.Interval is 0-based half-open
    start0 = max(0, pos_1based - 1 - width // 2)
    end0 = start0 + width
    return genome.Interval(chrom, start0, end0, strand='.')


def _make_variant(chrom: str, pos_1based: int, alt: str):
    from alphagenome.data import genome
    return genome.Variant(chromosome=chrom, position=int(pos_1based), reference_bases="", alternate_bases=alt)


def _save_arrays_npz(base_dir: Path, candidate_id: str, label: str, *,
                     x: Optional[np.ndarray] = None,
                     ref_vals: Optional[np.ndarray] = None,
                     alt_vals: Optional[np.ndarray] = None,
                     delta_vals: Optional[np.ndarray] = None,
                     quantiles: Optional[np.ndarray] = None) -> Path:
    base_dir.mkdir(parents=True, exist_ok=True)
    fpath = base_dir / f"{candidate_id}_{label}.npz"
    np.savez_compressed(fpath, x=x, ref=ref_vals, alt=alt_vals, delta=delta_vals, quantiles=quantiles)
    return fpath


def _trackdata_to_arrays(td, aggregate_axis: int = -1) -> Tuple[np.ndarray, np.ndarray]:
    # Returns (x_coords, aggregated_values)
    # td.values shape: [bins, n_tracks] (or multi-d)
    values = td.values
    bins = values.shape[0]
    # Aggregate across tracks if present
    if values.ndim == 1:
        agg = values
    else:
        if values.shape[aggregate_axis] == 0:
            # No tracks available; return x with None to signal skip
            x = np.arange(td.interval.start, td.interval.end, td.resolution) if td.interval else np.arange(bins)
            return x, None
        agg = np.nanmean(values, axis=aggregate_axis)
    x = np.arange(td.interval.start, td.interval.end, td.resolution) if td.interval else np.arange(bins)
    return x, agg


def _junctiondata_to_list(jd) -> List[Tuple]:
    # Returns list of (start, end, strand, k) with k aggregated across tracks
    # jd.values shape: [n_junctions, n_tracks]
    vals = jd.values
    if vals.shape[1] == 0:
        return []
    k = np.nanmean(vals, axis=1)
    out = []
    for junc, kval in zip(jd.junctions, k):
        # Guard NaNs -> 0.0 for storage/plotting
        kv = float(kval)
        if not np.isfinite(kv):
            kv = 0.0
        out.append((junc.start, junc.end, junc.strand, kv))
    return out


def main(argv: List[str] | None = None) -> None:
    ap = argparse.ArgumentParser(
        prog="AlphaGenomeScorer",
        description="Batch-score candidates via AlphaGenome (splicing ± RNA) with a variant-centered scorer.",
    )
    ap.add_argument("--config", required=False, type=Path)
    ap.add_argument("--candidates", required=False, type=Path)
    ap.add_argument("--modalities", nargs='*', default=None)
    ap.add_argument("--variant-window", type=int, default=None)
    ap.add_argument("--out", required=True, type=Path)

    args = ap.parse_args(argv)

    cfg = AppConfig.from_yaml(args.config) if args.config else AppConfig()
    cfg.ensure_out_dirs()
    modalities = args.modalities or cfg.scoring.modalities
    variant_window = args.variant_window or cfg.scoring.variant_window_nt
    outputs = _ag_output_types(modalities)
    variant_scorers = _variant_scorers(modalities, variant_window)

    api_key = _resolve_api_key(cfg)
    client = _build_client(api_key, cfg.alphagenome.address)
    # Import the module for enum access
    from alphagenome.models import dna_client as ag_dna_client

    cand_path = Path(args.candidates) if args.candidates else cfg.io.candidates_tsv
    df = _load_candidates(cand_path)

    records: List[Dict] = []
    arrays_dir = cfg.io.out_dir / "arrays"
    arrays_dir.mkdir(parents=True, exist_ok=True)

    # Optional ontology terms
    ontology_terms = cfg.scoring.tissues if cfg.scoring.tissues else None

    from alphagenome.models import dna_output as out

    for _, row in df.iterrows():
        cand_id = row["candidate_id"]
        chrom = row["chrom"]
        pos = int(row["pos"])  # 1-based
        alt = str(row["alt"]).strip().upper()

        interval = _interval_around(chrom, pos, cfg.alphagenome.sequence_length)
        variant = _make_variant(chrom, pos, alt)

        # Full predictions for REF/ALT
        var_pred = client.predict_variant(
            interval=interval,
            variant=variant,
            requested_outputs=outputs,
            organism=ag_dna_client.Organism.HOMO_SAPIENS,
            ontology_terms=ontology_terms,
        )

        # Variant-centered scoring for effect metrics
        scores = client.score_variant(
            interval=interval,
            variant=variant,
            variant_scorers=variant_scorers,
            organism=ag_dna_client.Organism.HOMO_SAPIENS,
        )

        # Persist arrays for each modality
        # Splicing: junctions and site usage
        if out.OutputType.SPLICE_JUNCTIONS in outputs and var_pred.reference.splice_junctions is not None:
            ref_j = var_pred.reference.splice_junctions
            alt_j = var_pred.alternate.splice_junctions
            # Aggregate k across tracks for storage; reporter can re-plot per track later if desired
            ref_list = _junctiondata_to_list(ref_j)
            alt_list = _junctiondata_to_list(alt_j)
            # Store JSON for junctions to keep shapes small
            j_path = arrays_dir / f"{cand_id}_splicing_junctions.json"
            with open(j_path, 'w') as f:
                json.dump({"ref": ref_list, "alt": alt_list, "interval": [interval.start, interval.end]}, f)

        if out.OutputType.SPLICE_SITE_USAGE in outputs and var_pred.reference.splice_site_usage is not None:
            ref_su = var_pred.reference.splice_site_usage
            alt_su = var_pred.alternate.splice_site_usage
            x_su, ref_su_vals = _trackdata_to_arrays(ref_su)
            _, alt_su_vals = _trackdata_to_arrays(alt_su)
            if ref_su_vals is not None and alt_su_vals is not None:
                delta_su = alt_su_vals - ref_su_vals
                _save_arrays_npz(arrays_dir, cand_id, "splice_site_usage", x=x_su, ref_vals=ref_su_vals, alt_vals=alt_su_vals, delta_vals=delta_su)

        if out.OutputType.RNA_SEQ in outputs and var_pred.reference.rna_seq is not None:
            ref_rna = var_pred.reference.rna_seq
            alt_rna = var_pred.alternate.rna_seq
            x_rna, ref_rna_vals = _trackdata_to_arrays(ref_rna)
            _, alt_rna_vals = _trackdata_to_arrays(alt_rna)
            if ref_rna_vals is not None and alt_rna_vals is not None:
                delta_rna = alt_rna_vals - ref_rna_vals
                _save_arrays_npz(arrays_dir, cand_id, "rna", x=x_rna, ref_vals=ref_rna_vals, alt_vals=alt_rna_vals, delta_vals=delta_rna)

        # Transcription initiation tracks (CAGE, PROCAP)
        if hasattr(var_pred.reference, 'cage') and var_pred.reference.cage is not None and (hasattr(out.OutputType, 'CAGE') and out.OutputType.CAGE in outputs):
            ref_cage = var_pred.reference.cage
            alt_cage = var_pred.alternate.cage
            x_cage, ref_cage_vals = _trackdata_to_arrays(ref_cage)
            _, alt_cage_vals = _trackdata_to_arrays(alt_cage)
            if ref_cage_vals is not None and alt_cage_vals is not None:
                delta_cage = alt_cage_vals - ref_cage_vals
                _save_arrays_npz(arrays_dir, cand_id, "cage", x=x_cage, ref_vals=ref_cage_vals, alt_vals=alt_cage_vals, delta_vals=delta_cage)

        if hasattr(var_pred.reference, 'procap') and var_pred.reference.procap is not None and (hasattr(out.OutputType, 'PROCAP') and out.OutputType.PROCAP in outputs):
            ref_pro = var_pred.reference.procap
            alt_pro = var_pred.alternate.procap
            x_pro, ref_pro_vals = _trackdata_to_arrays(ref_pro)
            _, alt_pro_vals = _trackdata_to_arrays(alt_pro)
            if ref_pro_vals is not None and alt_pro_vals is not None:
                delta_pro = alt_pro_vals - ref_pro_vals
                _save_arrays_npz(arrays_dir, cand_id, "procap", x=x_pro, ref_vals=ref_pro_vals, alt_vals=alt_pro_vals, delta_vals=delta_pro)

        # Summarize score_variant AnnData results for ranking
        for score in scores:
            scorer = score.uns.get('variant_scorer')
            name = str(scorer)
            # score.X is an array (tracks,) possibly, layers['quantiles'] present
            xvals = np.ravel(score.X)
            qvals = None
            if 'quantiles' in score.layers:
                qvals = np.ravel(score.layers['quantiles'])

            # Aggregate metrics
            primary = float(np.nanmean(np.abs(xvals))) if xvals.size else np.nan
            q_mean = float(np.nanmean(qvals)) if qvals is not None and qvals.size else np.nan

            records.append(
                dict(
                    candidate_id=cand_id,
                    chrom=chrom,
                    pos=pos,
                    interval_start=interval.start,
                    interval_end=interval.end,
                    scorer=name,
                    primary=primary,
                    quantile_mean=q_mean,
                )
            )

            # Persist per-track GeneMask RNA LFC for plotting
            try:
                from alphagenome.models import variant_scorers as vs
                import pandas as _pd
                if isinstance(scorer, vs.GeneMaskLFCScorer):
                    X = score.X
                    per_track = None
                    if X.ndim == 2:
                        # Avoid RuntimeWarnings: compute nanmean manually via sums/counts.
                        counts = np.sum(np.isfinite(X), axis=0)
                        if np.any(counts > 0):
                            sums = np.nansum(X, axis=0)
                            per_track = np.full_like(sums, np.nan, dtype=float)
                            mask = counts > 0
                            per_track[mask] = sums[mask] / counts[mask]
                    else:
                        flat = np.ravel(X)
                        if flat.size > 0 and np.any(np.isfinite(flat)):
                            per_track = flat

                    if per_track is not None:
                        track_names = None
                        if isinstance(score.var, _pd.DataFrame) and 'name' in score.var.columns:
                            track_names = [str(n) for n in score.var['name'].tolist()]
                        np.savez_compressed(
                            arrays_dir / f"{cand_id}_rna_gene_lfc.npz",
                            scores=per_track,
                            track_names=np.array(track_names, dtype=object),
                        )
            except Exception:
                pass

    raw_df = pd.DataFrame.from_records(records)
    out_path = Path(args.out) if args.out else cfg.io.raw_parquet
    out_path.parent.mkdir(parents=True, exist_ok=True)
    raw_df.to_parquet(out_path)
    print(f"[AlphaGenomeScorer] Wrote predictions/metrics to {out_path}")


if __name__ == "__main__":
    main()
