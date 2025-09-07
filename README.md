AlphaGenome Variant Insertion Scoring Pipeline
=============================================

## Overview
--------

This repository provides a command‑line pipeline for evaluating potential genomic insertion sites for a user‑provided DNA cassette using AlphaGenome. Given a genomic interval (BED) and an insertion sequence (FASTA), the pipeline:

- enumerates candidate insertion positions (optionally respecting splice‑signal buffers),
- predicts REF vs ALT effects for multiple modalities (splicing, local RNA; optional CAGE/PROCAP; optional TF binding and histone marks),
- aggregates variant‑centred and gene‑level scores,
- ranks sites by minimal predicted disruption, and
- produces per‑candidate plots plus an HTML report.

Originally, the purpose was to insert an shRNA in an intron. However, the tools are generic: they accept any interval (not only introns) and any cassette length. AlphaGenome runs remotely and requires an API key; all requests are variant‑centred with a configurable sequence context and aggregation window.

## Components
----------

- VariantBuilder: emits insertion candidates respecting splice-signal buffers
- AlphaGenomeScorer: runs AlphaGenome predictions (REF/ALT) and variant-centered scores (supports optional TF/histone)
- Ranker: ranks by minimal predicted disruption (splicing-driven)
- Reporter: plots per-candidate figures and builds an HTML report

## Quick Start (Conda)
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
```

## One-shot pipeline
-----------------

Run all steps (VariantBuilder → Scorer → Ranker → Reporter) in a single command:

```
python -m ag_pipeline.cli Full \
  --config ag.yaml \
  --intron-bed data/gene_intron2_hg38.bed \
  --cassette data/shrna_cassette.fa \
  --modalities splicing rna tss tf histone \
  --variant-window 501
```

Notes:
- If you omit `--buffers/--stride/--max`, they default to values from `ag.yaml`.
- Outputs default to `data/candidates.tsv`, `ag_out/raw.parquet`, `ag_out/candidates.csv`, and `ag_out/plots` + `ag_out/report.html`.

## Individual steps (optional)
---------------------------

You can still run each step independently (generic to any interval and cassette length):


### Candidates (generic interval)

```
python -m ag_pipeline.cli VariantBuilder \
  --intron-bed data/region.bed \
  --cassette data/insert.fa \
  --buffers 120 18 40 120 \   # use 0 0 0 0 for non‑intronic intervals
  --stride 25 \
  --max 80 \
  --out data/candidates.tsv
```

### AlphaGenome predictions + scores (now with gene-level RNA and TSS tracks)

```
python -m ag_pipeline.cli AlphaGenomeScorer \
  --config ag.yaml \
  --candidates data/candidates.tsv \
  --modalities splicing rna tss tf histone \
  --variant-window 501 \
  --out ag_out/raw.parquet
```

### Rank

```
python -m ag_pipeline.cli Ranker \
  --in ag_out/raw.parquet \
  --out ag_out/candidates.csv
```

### Plots + HTML

```
python -m ag_pipeline.cli Reporter \
  --scores ag_out/candidates.csv \
  --pred ag_out/raw.parquet \
  --plots ag_out/plots \
  --html ag_out/report.html
```

## Configuration and CLI mapping
-----------------------------

All knobs can be provided either via the config file (`ag.yaml`) or via CLI flags. The Full command reads from `ag.yaml` by default when flags are omitted; individual commands also accept `--config` to use the same defaults.

Config schema (ag.yaml):
- `alphagenome`:
  - `api_key_env`: Name of env var holding your API key (do not put the secret here).
  - `address`: Optional service address override.
  - `sequence_length`: Model context (e.g., 2048, 16384, 131072, ...).
  - `retries`, `batch_size`: Client retry and batching knobs.
- `scoring`:
  - `modalities`: List of modalities to request/pipeline (e.g., ["splicing", "rna", "tss", "tf", "histone"]).
    - `splicing` → splice_site_usage, splice_sites (CenterMask) + gene-level splicing + splice junction scorer.
    - `rna` → local RNA (CenterMask) + gene-level RNA LFC (GeneMaskLFC).
    - `tss` (alias: `promoter`, `initiation`) → CAGE and PROCAP (CenterMask).
    - `tf` (aliases: `chip_tf`, `tf_binding`) → TF binding tracks (CenterMask) where supported by the AlphaGenome build.
    - `histone` (alias: `chip_histone`) → aggregate histone marks or common H3K* marks (CenterMask) where supported.
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
  - Output columns: `rank, candidate_id, pos, splicing, rna, rna_gene, tss, tf, histone, composite`
  - Composite: `splicing + 0.2·rna + 0.2·rna_gene + 0.1·tss + 0.1·tf + 0.1·histone`

- Reporter:
  - `--config`: Use paths from ag.yaml.
  - `--scores` ← `io.scores_csv`
  - `--pred` ← `io.raw_parquet`
  - `--plots` ← `io.plots_dir`
  - `--html` ← `io.report_html`
- Produces per-candidate panels (splicing, RNA, CAGE/PROCAP if enabled, plus TF/histone if enabled) and RNA gene-level LFC bar.

New: Gene structure figure with candidate spikes
- Reporter can render a gene-wide exon/intron track with the intron-of-interest highlighted and one spike per candidate colored green→red by composite rank (green = best). Provide a transcript ID (or a local GTF) and the intron BED:

  ag Reporter \
    --scores ag_out/candidates.csv \
    --pred ag_out/raw.parquet \
    --plots ag_out/plots \
    --html ag_out/report.html \
    --transcript ENST00000325495 \
    --intron-bed data/gene_intron2_hg38.bed

  - Optional: add --gtf path/to/your.gtf to avoid Ensembl REST.
  - The figure is saved to ag_out/plots/gene_structure_spikes.png and embedded near the top of the HTML report.
  - Spike height is scaled by composite score: lower/better composite → taller spike. If composite is missing, height falls back to rank.

- Full:
  - `--config`: Required for defaults.
  - `--intron-bed` / `--cassette`: Required via CLI or `inputs.*`.
  - Optional `--buffers/--stride/--max`: fall back to `buffers.*`/`scan.*`.
  - Optional `--modalities`/`--variant-window`: fall back to `scoring.*`.
  - Optional `--raw-out/--scores-out/--plots/--html`: fall back to `io.*`.

Environment
- Export your API key in the env var named by `alphagenome.api_key_env` (default `ALPHAGENOME_API_KEY`).
- The code validates CenterMask widths and will coerce unsupported values to nearest supported.


## Notes
-----

- The AlphaGenome client uses a 2,048 bp sequence context by default (config: `alphagenome.sequence_length`). The variant-centered `CenterMaskScorer` aggregates over a window defined by `scoring.variant_window_nt`. Use a supported width (e.g., 501, 2001, 10001); the CLI will automatically coerce unsupported values to the nearest allowed width.
- Modalities:
  - `splicing`: Splice site usage and sites (CenterMask), plus GeneMask splicing and junction-centric scoring.
  - `rna`: Local RNA (CenterMask) and gene-level RNA log fold change (GeneMaskLFC).
  - `tss`: Transcription initiation proxies (CAGE, PROCAP) with CenterMask.
- `tf`: Transcription factor binding (if supported by the AlphaGenome build), scored with CenterMask.
- `histone`: Histone marks (aggregate or common H3K* marks), scored with CenterMask.
- Ranking composite now includes TSS/TF/histone contributions: composite = splicing + 0.2·RNA(local) + 0.2·RNA(gene) + 0.1·TSS + 0.1·TF + 0.1·Histone.
- Reporter shows vertical panel images for splicing and RNA, panels for CAGE/PROCAP (if requested), and a small bar chart of gene-level RNA LFC by track.
- `ALPHAGENOME_API_KEY` is required. You can also set `ALPHA_GENOME_API_KEY` or change `alphagenome.api_key_env` in `ag.yaml`.
- Arrays for plotting are persisted to `ag_out/arrays/*.npz` and junctions as JSON. `ag_out/raw.parquet` stores summarised effect metrics for ranking.
- Treat model results as hypothesis-generating; follow up with minigene / RT-PCR validation.

Tip: If you prefer not to activate the env, you can prefix commands with Conda:

```
conda run -n ag-shrna python -m ag_pipeline.cli Ranker --in ag_out/raw.parquet --out ag_out/candidates.csv
```

Direct module entry points (advanced)
------------------------------------

If you prefer to call the underlying modules directly (skipping the unified `ag_pipeline.cli` wrapper), you can invoke each stage as shown below. These accept the same arguments as their CLI counterparts; use `--config ag.yaml` to inherit defaults where applicable.

VariantBuilder (module)
```
python -m ag_pipeline.variant_builder \
  --intron-bed data/region.bed \
  --cassette data/insert.fa \
  --buffers 120 18 40 120 \   # for non‑intronic regions use: 0 0 0 0
  --stride 25 \
  --max 80 \
  --out data/candidates.tsv
```
- --intron-bed: BED interval to scan (can be any region, not just introns).
- --cassette: FASTA with insertion sequence.
- --buffers DONOR BP_START BP_END ACCEPTOR: Splice-signal exclusion buffers; set to `0 0 0 0` for generic regions.
- --stride: Step in nucleotides when enumerating candidate breakpoints.
- --max: Maximum number of candidates to emit.
- --out: TSV written by this step.
  Optional: `--focus POS --window N` to keep only candidates within N nt of POS.

Scorer (module; with optional TF and histone)
```
python -m ag_pipeline.scorer \
  --config ag.yaml \
  --candidates data/candidates.tsv \
  --modalities splicing rna tss tf histone \
  --variant-window 501 \
  --out ag_out/raw.parquet
```
- --config: Provides AlphaGenome client defaults (API key env, sequence_length, tissues).
- --candidates: Candidates TSV produced by VariantBuilder.
- --modalities: Any subset of `splicing rna tss tf histone` supported by your AlphaGenome build.
- --variant-window: CenterMask aggregation width (use a supported value, e.g., 501, 2001, 10001).
- --out: Long/tidy Parquet with per-scorer aggregates used downstream.

Ranker (module)
```
python -m ag_pipeline.ranker \
  --in ag_out/raw.parquet \
  --out ag_out/candidates.csv
```
- --in: Parquet from Scorer.
- --out: Ranked CSV with modality scores and composite.

Reporter (module)
```
python -m ag_pipeline.reporter \
  --scores ag_out/candidates.csv \
  --pred ag_out/raw.parquet \
  --plots ag_out/plots \
  --html ag_out/report.html
```
- --scores: Ranked CSV from Ranker.
- --pred: Parquet from Scorer (for plots).
- --plots/--html: Output locations for images and the HTML report.

## Outputs
-------

- `data/candidates.tsv`: Raw candidates from VariantBuilder (one row per candidate site). Columns: `candidate_id, chrom, pos, ref, alt`. See notes on coordinate conventions below.
- `ag_out/raw.parquet`: Long/tidy table with per-candidate, per-scorer aggregated metrics used for ranking.
- `ag_out/candidates.csv`: Final ranked summary (one row per candidate) with scores per modality and the composite. See “Candidates CSV (ranked scores)” below.
- `ag_out/plots/`: Per-candidate panels (splicing, RNA; CAGE/PROCAP if enabled; TF/histone if enabled) and gene structure figure when configured.
- `ag_out/report.html`: HTML report that embeds the ranked table and plots.
- `ag_out/arrays/*.npz` and `ag_out/arrays/*_splicing_junctions.json`: Arrays used for plotting and splice junction metadata.
- Optional: `ag_out/candidates_scored.tsv` written by AlphaGenomeScorer merges `data/candidates.tsv` with the wide per‑candidate scores for convenience.

### Candidates CSV (ranked scores)
------------------------------

This is the main output produced by Ranker at `ag_out/candidates.csv`. It contains one row per candidate with modality scores and a composite used for ordering.

- rank: 1-based rank; 1 = best (lowest predicted disruption).
- candidate_id: Identifier matching the VariantBuilder output (e.g., `cand0001`).
- pos: 1‑based genomic insertion breakpoint used for scoring.
- splicing: Mean absolute predicted splicing disruption around the variant (CenterMask aggregate). Lower is better.
- rna: Mean absolute predicted local RNA expression change around the variant (CenterMask). Lower is better.
- rna_gene: Mean absolute gene‑level RNA log fold change (GeneMaskLFC across RNA‑seq tracks). Lower magnitude is better.
- tss: Mean absolute predicted change in transcription initiation proxies (CAGE/PROCAP; CenterMask). Lower is better.
- tf: Mean absolute predicted change across TF binding tracks if available (CenterMask). Lower is better.
- histone: Mean absolute predicted change across histone marks if available (CenterMask or selected H3K* aggregates). Lower is better.
- composite: Splicing‑weighted sum used for ranking: `splicing + 0.2·rna + 0.2·rna_gene + 0.1·tss + 0.1·tf + 0.1·histone`. Lower composite → better rank.

##### Coordinate conventions (applies to `data/candidates.tsv` and `ag_out/candidates.csv`)
- 1‑based position: `pos` follows the AlphaGenome variant schema as a 1‑based coordinate.
- Pure insertion: `ref` is empty and `alt` contains the full inserted sequence in the TSV. Conceptually, the cassette is inserted at the breakpoint without deleting reference bases.
- Model context: Scoring uses the model context (`alphagenome.sequence_length`) centered on `pos`; variant‑centered scorers aggregate over `scoring.variant_window_nt`.
