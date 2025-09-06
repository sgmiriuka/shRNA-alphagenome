AlphaGenome Variant Insertion Scoring Pipeline
=============================================

Overview
--------

This repository provides a command‑line pipeline for evaluating potential genomic insertion sites for a user‑provided DNA cassette using AlphaGenome. Given a genomic interval (BED) and an insertion sequence (FASTA), the pipeline:

- enumerates candidate insertion positions (optionally respecting splice‑signal buffers),
- predicts REF vs ALT effects for multiple modalities (splicing, local RNA; optional CAGE/PROCAP),
- aggregates variant‑centred and gene‑level scores,
- ranks sites by minimal predicted disruption, and
- produces per‑candidate plots plus an HTML report.

Originally, the purpose was to insert an shRNA in an intron. However, the tools are generic: they accept any interval (not only introns) and any cassette length. AlphaGenome runs remotely and requires an API key; all requests are variant‑centred with a configurable sequence context and aggregation window.

Components
----------

- VariantBuilder: emits insertion candidates respecting splice-signal buffers
- AlphaGenomeScorer: runs AlphaGenome predictions (REF/ALT) and variant-centered scores
- Ranker: ranks by minimal predicted disruption (splicing-driven)
- Reporter: plots per-candidate figures and builds an HTML report

Quick Start (Conda)
------------------

1) Create and activate the Conda env (installs AlphaGenome via pip inside the env):

```
conda env create -f environment.yml
conda activate ag-shrna
```

Ensure you have an API key in env var `ALPHAGENOME_API_KEY` (or `ALPHA_GENOME_API_KEY`).

2) Prepare inputs in `data/`:
- `data/gene_intron2_hg38.bed` (0-based interval; any region works, not just introns)
- `data/shrna_cassette.fa` (FASTA with one sequence; any length supported)

3) Configure `ag.yaml` as needed (modalities, variant window, sequence_length, etc.).

4) Run the full pipeline (standard Python CLI):

```
# (Optional) Create intron-2 BED from a transcript ID (hg38)
python -m ag_pipeline.cli TranscriptToBED \
  --transcript ENST00000325495 \
  --intron-index 2 \
  --out data/gene_intron2_hg38.bed

One-shot pipeline
-----------------

Run all steps (VariantBuilder → Scorer → Ranker → Reporter) in a single command:

```
python -m ag_pipeline.cli Full \
  --config ag.yaml \
  --intron-bed data/gene_intron2_hg38.bed \
  --cassette data/shrna_cassette.fa \
  --modalities splicing rna tss \
  --variant-window 501
```

Notes:
- If you omit `--buffers/--stride/--max`, they default to values from `ag.yaml`.
- Outputs default to `data/candidates.tsv`, `ag_out/raw.parquet`, `ag_out/candidates.csv`, and `ag_out/plots` + `ag_out/report.html`.

Individual steps (optional)
---------------------------

You can still run each step independently (generic to any interval and cassette length):

```
# Candidates (generic interval)
python -m ag_pipeline.cli VariantBuilder \
  --intron-bed data/region.bed \
  --cassette data/insert.fa \
  --buffers 120 18 40 120 \   # use 0 0 0 0 for non‑intronic intervals
  --stride 25 \
  --max 80 \
  --out data/candidates.tsv

# AlphaGenome predictions + scores (now with gene-level RNA and TSS tracks)
python -m ag_pipeline.cli AlphaGenomeScorer \
  --config ag.yaml \
  --candidates data/candidates.tsv \
  --modalities splicing rna tss \
  --variant-window 501 \
  --out ag_out/raw.parquet

# Rank
python -m ag_pipeline.cli Ranker \
  --in ag_out/raw.parquet \
  --out ag_out/candidates.csv

# Plots + HTML
python -m ag_pipeline.cli Reporter \
  --scores ag_out/candidates.csv \
  --pred ag_out/raw.parquet \
  --plots ag_out/plots \
  --html ag_out/report.html
```

Configuration and CLI mapping
-----------------------------

All knobs can be provided either via the config file (`ag.yaml`) or via CLI flags. The Full command reads from `ag.yaml` by default when flags are omitted; individual commands also accept `--config` to use the same defaults.

Config schema (ag.yaml):
- `alphagenome`:
  - `api_key_env`: Name of env var holding your API key (do not put the secret here).
  - `address`: Optional service address override.
  - `sequence_length`: Model context (e.g., 2048, 16384, 131072, ...).
  - `retries`, `batch_size`: Client retry and batching knobs.
- `scoring`:
  - `modalities`: List of modalities to request/pipeline (e.g., ["splicing", "rna", "tss"]).
    - `splicing` → splice_site_usage, splice_sites (CenterMask) + gene-level splicing + splice junction scorer.
    - `rna` → local RNA (CenterMask) + gene-level RNA LFC (GeneMaskLFC).
    - `tss` (alias: `promoter`, `initiation`) → CAGE and PROCAP (CenterMask).
  - `variant_window_nt`: CenterMask width. Must be supported (e.g., 501, 2001, 10001); unsupported values are coerced to nearest.
  - `tissues`: Optional list of ontology CURIEs (e.g., ["UBERON:0002048"]) to filter tracks.
- `buffers`:
  - `donor_min_nt`: Min distance from 5′ donor.
  - `branchpoint_window_nt.start/.end`: Branchpoint exclusion window upstream of 3′ acceptor.
  - `acceptor_min_nt`: Min distance from 3′ acceptor.
- `scan`:
  - `stride_nt`: Step size in VariantBuilder.
  - `max_candidates`: Max candidates to emit.
- `inputs` (optional):
  - `intron_bed`: Path to the interval BED to scan (any region; not restricted to introns).
  - `cassette`: Path to the insertion (e.g., shRNA cassette) FASTA (any length).
- `io`:
  - `out_dir`: Base directory for outputs.
  - `make_html`: Whether to build HTML.
  - `candidates_tsv`, `raw_parquet`, `scores_csv`, `plots_dir`, `report_html`: Default paths for each step.

CLI flags and their config fallbacks:
- VariantBuilder:
  - `--config`: Use ag.yaml defaults.
  - `--intron-bed` ← `inputs.intron_bed` (required via one of them)
  - `--cassette` ← `inputs.cassette` (required via one of them)
  - `--buffers DONOR BP_START BP_END ACCEPTOR` ← `buffers.*` (for non‑intronic intervals, set 0 0 0 0)
  - `--stride` ← `scan.stride_nt`
  - `--max` ← `scan.max_candidates`
  - `--out` ← `io.candidates_tsv`
  - `--focus`/`--window`: optional local refinement; keep candidates with |pos−focus| ≤ window

- AlphaGenomeScorer:
  - `--config`: Required for defaults (API key env, sequence_length, etc.).
  - `--candidates` ← `io.candidates_tsv`
  - `--modalities` ← `scoring.modalities`
  - `--variant-window` ← `scoring.variant_window_nt`
  - `--out` ← `io.raw_parquet`
  - Uses `alphagenome.*` for client config and `scoring.tissues` for ontology filtering.

- Ranker:
  - `--config`: Use paths from ag.yaml.
  - `--in` ← `io.raw_parquet`
  - `--out` ← `io.scores_csv`
  - Output columns: `rank, candidate_id, pos, splicing, rna, rna_gene, tss, composite`
  - Composite: `splicing + 0.2·rna + 0.2·rna_gene`

- Reporter:
  - `--config`: Use paths from ag.yaml.
  - `--scores` ← `io.scores_csv`
  - `--pred` ← `io.raw_parquet`
  - `--plots` ← `io.plots_dir`
  - `--html` ← `io.report_html`
  - Produces per-candidate panels (splicing, RNA, CAGE/PROCAP if enabled) and RNA gene-level LFC bar.

- Full:
  - `--config`: Required for defaults.
  - `--intron-bed` / `--cassette`: Required via CLI or `inputs.*`.
  - Optional `--buffers/--stride/--max`: fall back to `buffers.*`/`scan.*`.
  - Optional `--modalities`/`--variant-window`: fall back to `scoring.*`.
  - Optional `--raw-out/--scores-out/--plots/--html`: fall back to `io.*`.

Environment
- Export your API key in the env var named by `alphagenome.api_key_env` (default `ALPHAGENOME_API_KEY`).
- The code validates CenterMask widths and will coerce unsupported values to nearest supported.
```

Notes
-----

- The AlphaGenome client uses a 2,048 bp sequence context by default (config: `alphagenome.sequence_length`). The variant-centered `CenterMaskScorer` aggregates over a window defined by `scoring.variant_window_nt`. Use a supported width (e.g., 501, 2001, 10001); the CLI will automatically coerce unsupported values to the nearest allowed width.
- Modalities:
  - `splicing`: Splice site usage and sites (CenterMask), plus GeneMask splicing and junction-centric scoring.
  - `rna`: Local RNA (CenterMask) and gene-level RNA log fold change (GeneMaskLFC).
  - `tss`: Transcription initiation proxies (CAGE, PROCAP) with CenterMask.
- Ranking composite now includes gene-level RNA separately: composite = splicing + 0.2·RNA(local) + 0.2·RNA(gene).
- Reporter shows vertical panel images for splicing and RNA, panels for CAGE/PROCAP (if requested), and a small bar chart of gene-level RNA LFC by track.
- `ALPHAGENOME_API_KEY` is required. You can also set `ALPHA_GENOME_API_KEY` or change `alphagenome.api_key_env` in `ag.yaml`.
- Arrays for plotting are persisted to `ag_out/arrays/*.npz` and junctions as JSON. `ag_out/raw.parquet` stores summarised effect metrics for ranking.
- Treat model results as hypothesis-generating; follow up with minigene / RT-PCR validation.

Tip: If you prefer not to activate the env, you can prefix commands with Conda:

```
conda run -n ag-shrna python -m ag_pipeline.cli Ranker --in ag_out/raw.parquet --out ag_out/candidates.csv
```

Alternative direct module calls (skip unified CLI):

```
# VariantBuilder
python -m ag_pipeline.variant_builder --intron-bed ... --cassette ... --buffers 120 18 40 120 --stride 25 --max 80 --out data/candidates.tsv

# Scorer
python -m ag_pipeline.scorer --config ag.yaml --candidates data/candidates.tsv --modalities splicing rna --variant-window 400 --out ag_out/raw.parquet

# Ranker
python -m ag_pipeline.ranker --in ag_out/raw.parquet --out ag_out/candidates.csv

# Reporter
python -m ag_pipeline.reporter --scores ag_out/candidates.csv --pred ag_out/raw.parquet --plots ag_out/plots --html ag_out/report.html
```
