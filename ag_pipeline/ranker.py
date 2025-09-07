from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, List

import numpy as np
import pandas as pd
from .config import AppConfig


def compute_composite(scores_df: pd.DataFrame) -> pd.DataFrame:
    # scores_df has rows per (candidate_id, scorer) with columns primary, quantile_mean
    # We will split scorers by modality keywords and aggregate.
    # Heuristic mapping
    def modality_of(s: str) -> str:
        s = s.lower()
        if "splice" in s:
            return "splicing"
        if "genemasklfcscorer" in s and "rna_seq" in s:
            return "rna_gene"
        if "rna" in s:
            return "rna"
        if ("cage" in s) or ("procap" in s):
            return "tss"
        # TF/ChIP: common tokens in scorer repr or output names
        if ("tf_binding" in s) or ("chip_tf" in s) or (" tf " in f" {s} ") or ("tf(" in s) or ("tf)" in s):
            return "tf"
        # Histone marks: look for generic keyword or common H3 tokens
        if ("histone" in s) or ("chip_histone" in s) or ("h3k" in s):
            return "histone"
        return "other"

    scores_df = scores_df.copy()
    scores_df["modality"] = scores_df["scorer"].map(modality_of)

    agg = (
        scores_df.groupby(["candidate_id", "modality"], as_index=False)
        .agg(primary_mean=("primary", "mean"), quantile_mean=("quantile_mean", "mean"))
    )
    # Pivot to wide
    wide = agg.pivot(index="candidate_id", columns="modality", values="primary_mean").reset_index()
    wide.columns.name = None
    for col in ("splicing", "rna", "rna_gene", "tss", "tf", "histone"):
        if col not in wide.columns:
            wide[col] = np.nan

    # Bring through insertion position (1-based) from original scores_df
    if "pos" in scores_df.columns:
        pos_df = scores_df.groupby("candidate_id", as_index=False).agg(pos=("pos", "first"))
        wide = wide.merge(pos_df, on="candidate_id", how="left")
        # Ensure integer if possible
        with np.errstate(all='ignore'):
            wide["pos"] = wide["pos"].astype("Int64")

    # Composite (lower is better): splicing + small contributions from other modalities.
    # Weights chosen heuristically to emphasize splicing, with supporting evidence from RNA, TSS, TF, histone.
    wide["composite"] = (
        wide["splicing"].fillna(0.0)
        + 0.2 * wide["rna"].fillna(0.0)
        + 0.2 * wide["rna_gene"].fillna(0.0)
        + 0.1 * wide["tss"].fillna(0.0)
        + 0.1 * wide["tf"].fillna(0.0)
        + 0.1 * wide["histone"].fillna(0.0)
    )

    # Optional soft penalty: nothing here (buffer already enforced upstream)
    wide = wide.sort_values(by=["composite", "candidate_id"]).reset_index(drop=True)
    wide["rank"] = np.arange(1, len(wide) + 1)
    # Output order: rank, candidate_id, pos, splicing, rna, rna_gene, tss, tf, histone, composite
    cols = ["rank", "candidate_id"]
    if "pos" in wide.columns:
        cols.append("pos")
    cols += ["splicing", "rna", "rna_gene", "tss", "tf", "histone", "composite"]
    # Keep only columns that exist
    cols = [c for c in cols if c in wide.columns]
    return wide[cols]


def main(argv: List[str] | None = None) -> None:
    ap = argparse.ArgumentParser(
        prog="Ranker",
        description="Rank sites by minimal predicted disruption (splicing-driven).",
    )
    ap.add_argument("--config", required=False, type=Path)
    ap.add_argument("--in", dest="inp", required=False, type=Path)
    ap.add_argument("--out", required=False, type=Path)
    args = ap.parse_args(argv)

    cfg = AppConfig.from_yaml(args.config) if args.config else AppConfig()
    inp = Path(args.inp) if args.inp else cfg.io.raw_parquet
    outp = Path(args.out) if args.out else cfg.io.scores_csv

    df = pd.read_parquet(inp)
    ranked = compute_composite(df)
    outp.parent.mkdir(parents=True, exist_ok=True)
    ranked.to_csv(outp, index=False)
    print(f"[Ranker] Wrote ranked candidates to {outp}")


if __name__ == "__main__":
    main()
