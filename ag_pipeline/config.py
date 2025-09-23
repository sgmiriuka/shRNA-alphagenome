from __future__ import annotations

import dataclasses
import os
from pathlib import Path
from typing import Any, Dict, List, Optional

import yaml


@dataclasses.dataclass
class AlphaGenomeConfig:
    """Configuration for AlphaGenome API client settings.

    Attributes:
        api_key_env: Environment variable name containing the AlphaGenome API key.
        address: Optional custom service address override.
        sequence_length: Model context width in nucleotides (must be supported by the service).
        retries: Number of retries for RPC calls.
        batch_size: Batch size for internal iteration (does not affect server batching).
    """
    api_key_env: str = "ALPHAGENOME_API_KEY"
    address: Optional[str] = None
    sequence_length: int = 2048
    retries: int = 3
    batch_size: int = 25


@dataclasses.dataclass
class ScoringConfig:
    """Configuration for scoring modalities and variant-centered aggregation.

    Attributes:
        modalities: List of modalities to request from AlphaGenome (e.g., splicing, rna, tss, tf, histone).
        variant_window_nt: CenterMask aggregation window width in nucleotides.
        tissues: Optional list of ontology CURIEs to filter tracks.
    """
    modalities: List[str] = dataclasses.field(default_factory=lambda: ["splicing", "rna"])  # noqa: E501
    variant_window_nt: int = 400
    tissues: List[str] = dataclasses.field(default_factory=list)


@dataclasses.dataclass
class BuffersConfig:
    """Configuration for splice-signal buffer distances in introns.

    Attributes:
        donor_min_nt: Minimum distance from 5' donor splice site.
        branchpoint_window_nt_start: Start of branchpoint exclusion window upstream of 3' acceptor.
        branchpoint_window_nt_end: End of branchpoint exclusion window upstream of 3' acceptor.
        acceptor_min_nt: Minimum distance from 3' acceptor splice site.
    """
    donor_min_nt: int = 120
    branchpoint_window_nt_start: int = 18
    branchpoint_window_nt_end: int = 40
    acceptor_min_nt: int = 120


@dataclasses.dataclass
class ScanConfig:
    """Configuration for scanning parameters in VariantBuilder.

    Attributes:
        stride_nt: Step size in nucleotides between candidate insertion positions.
        max_candidates: Maximum number of candidate positions to emit.
    """
    stride_nt: int = 25
    max_candidates: int = 80


@dataclasses.dataclass
class IOConfig:
    """Configuration for input/output file paths.

    Attributes:
        out_dir: Base directory for all outputs.
        make_html: Whether to generate HTML reports.
        candidates_tsv: Path to output candidates TSV from VariantBuilder.
        raw_parquet: Path to output raw predictions/scores Parquet from AlphaGenomeScorer.
        scores_csv: Path to output ranked scores CSV from Ranker.
        plots_dir: Directory for per-candidate plot images from Reporter.
        report_html: Path to output HTML report from Reporter.
    """
    out_dir: Path = Path("ag_out")
    make_html: bool = True
    # Common paths (can be overridden per CLI):
    candidates_tsv: Path = Path("data/candidates.tsv")
    raw_parquet: Path = Path("ag_out/raw.parquet")
    scores_csv: Path = Path("ag_out/candidates.csv")
    plots_dir: Path = Path("ag_out/plots")
    report_html: Path = Path("ag_out/report.html")


@dataclasses.dataclass
class InputsConfig:
    """Configuration for input file paths.

    Attributes:
        intron_bed: Path to the BED file defining the genomic interval to scan.
        cassette: Path to the FASTA file containing the insertion sequence.
        cassettes: Optional list of FASTA files for multi-cassette runs.
        transcript: Optional transcript ID associated with the intron BED (used for plotting gene structure).
        gtf: Optional GTF file for transcript exon structure lookup.
    """
    intron_bed: Optional[Path] = None
    cassette: Optional[Path] = None
    cassettes: List[Path] = dataclasses.field(default_factory=list)
    transcript: Optional[str] = None
    gtf: Optional[Path] = None


@dataclasses.dataclass
class AppConfig:
    """Main application configuration aggregating all sub-configurations.

    Attributes:
        alphagenome: AlphaGenome API client settings.
        scoring: Scoring modalities and variant window settings.
        buffers: Splice-signal buffer distances for introns.
        scan: Scanning parameters for candidate generation.
        io: Input/output file path configurations.
        inputs: Input file path configurations.
    """
    alphagenome: AlphaGenomeConfig = dataclasses.field(default_factory=AlphaGenomeConfig)
    scoring: ScoringConfig = dataclasses.field(default_factory=ScoringConfig)
    buffers: BuffersConfig = dataclasses.field(default_factory=BuffersConfig)
    scan: ScanConfig = dataclasses.field(default_factory=ScanConfig)
    io: IOConfig = dataclasses.field(default_factory=IOConfig)
    inputs: InputsConfig = dataclasses.field(default_factory=InputsConfig)

    @staticmethod
    def from_yaml(path: str | os.PathLike[str]) -> "AppConfig":
        """Load application configuration from a YAML file.

        Args:
            path: Path to the YAML configuration file.

        Returns:
            AppConfig: Parsed configuration object.

        Raises:
            FileNotFoundError: If the YAML file does not exist.
            yaml.YAMLError: If the YAML is malformed.
        """
        with open(path, "r") as f:
            raw = yaml.safe_load(f) or {}

        # Flatten helpers
        alpha = raw.get("alphagenome", {}) or {}
        scoring = raw.get("scoring", {}) or {}
        buffers = raw.get("buffers", {}) or {}
        scan = raw.get("scan", {}) or {}
        io = raw.get("io", {}) or {}

        alphagenome_cfg = AlphaGenomeConfig(
            api_key_env=alpha.get("api_key_env", "ALPHAGENOME_API_KEY"),
            address=alpha.get("address"),
            sequence_length=int(alpha.get("sequence_length", 2048)),
            retries=int(alpha.get("retries", 3)),
            batch_size=int(alpha.get("batch_size", 25)),
        )

        scoring_cfg = ScoringConfig(
            modalities=list(scoring.get("modalities", ["splicing", "rna"])),
            variant_window_nt=int(scoring.get("variant_window_nt", 400)),
            tissues=list(scoring.get("tissues", []) or []),
        )

        buffers_cfg = BuffersConfig(
            donor_min_nt=int(buffers.get("donor_min_nt", 120)),
            branchpoint_window_nt_start=int(
                (buffers.get("branchpoint_window_nt", {}) or {}).get("start", 18)
            ),
            branchpoint_window_nt_end=int(
                (buffers.get("branchpoint_window_nt", {}) or {}).get("end", 40)
            ),
            acceptor_min_nt=int(buffers.get("acceptor_min_nt", 120)),
        )

        scan_cfg = ScanConfig(
            stride_nt=int(scan.get("stride_nt", 25)),
            max_candidates=int(scan.get("max_candidates", 80)),
        )

        io_cfg = IOConfig(
            out_dir=Path(io.get("out_dir", "ag_out")),
            make_html=bool(io.get("make_html", True)),
            candidates_tsv=Path(io.get("candidates_tsv", "data/candidates.tsv")),
            raw_parquet=Path(io.get("raw_parquet", "ag_out/raw.parquet")),
            scores_csv=Path(io.get("scores_csv", "ag_out/candidates.csv")),
            plots_dir=Path(io.get("plots_dir", "ag_out/plots")),
            report_html=Path(io.get("report_html", "ag_out/report.html")),
        )

        inputs_raw = raw.get("inputs", {}) or {}
        multi = inputs_raw.get("cassettes") or []
        inputs_cfg = InputsConfig(
            intron_bed=Path(inputs_raw["intron_bed"]) if inputs_raw.get("intron_bed") else None,
            cassette=Path(inputs_raw["cassette"]) if inputs_raw.get("cassette") else None,
            cassettes=[Path(p) for p in multi],
            transcript=str(inputs_raw["transcript"]) if inputs_raw.get("transcript") else None,
            gtf=Path(inputs_raw["gtf"]) if inputs_raw.get("gtf") else None,
        )

        return AppConfig(
            alphagenome=alphagenome_cfg,
            scoring=scoring_cfg,
            buffers=buffers_cfg,
            scan=scan_cfg,
            io=io_cfg,
            inputs=inputs_cfg,
        )

    def ensure_out_dirs(self) -> None:
        """Ensure that output directories exist, creating them if necessary.

        Creates the main output directory and plots subdirectory.
        """
        self.io.out_dir.mkdir(parents=True, exist_ok=True)
        # Ensure plots directory
        plots_dir = self.io.plots_dir if self.io.plots_dir else (self.io.out_dir / "plots")
        Path(plots_dir).mkdir(parents=True, exist_ok=True)
