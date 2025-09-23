from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path
from typing import List, Set

from . import variant_builder, scorer, ranker, reporter
from . import transcript_bed
from . import variant_builder_extracting as vb_ext
from .config import AppConfig


def main(argv=None):
    """Main entry point for the AlphaGenome shRNA intron-2 pipeline CLI.

    Parses command-line arguments and dispatches to the appropriate subcommand
    (VariantBuilder, AlphaGenomeScorer, Ranker, Reporter, TranscriptToBED, or Full pipeline).

    Args:
        argv: List of command-line arguments. If None, uses sys.argv.

    Returns:
        int: Exit code (0 for success, 1 for error).
    """
    parser = argparse.ArgumentParser(prog="ag", description="AlphaGenome shRNA intron-2 pipeline")
    sub = parser.add_subparsers(dest="cmd", required=True)

    # VariantBuilder
    vb = sub.add_parser("VariantBuilder", help="Emit insertion or deletion variants within intron 2")
    vb.add_argument("--intron-bed", required=True)
    vb.add_argument("--cassette", required=True)
    vb.add_argument("--buffers", nargs=4, required=True)
    vb.add_argument("--stride", required=True)
    vb.add_argument("--max", dest="max_candidates", required=True)
    vb.add_argument("--out", required=True)
    vb.add_argument("--variant-type", choices=["insertion", "extraction"], default="insertion", help="Type of variant: insertion (add cassette) or extraction (delete cassette-length bases)")
    vb.add_argument("--genome-fasta", required=False, help="Genome FASTA file (required for extraction)")

    # AlphaGenomeScorer
    sc = sub.add_parser("AlphaGenomeScorer", help="Score candidates via AlphaGenome")
    sc.add_argument("--config", required=True)
    sc.add_argument("--candidates", required=True)
    sc.add_argument("--modalities", nargs='*', default=None)
    sc.add_argument("--variant-window", type=int, default=None)
    sc.add_argument("--out", required=True)

    # Ranker
    rk = sub.add_parser("Ranker", help="Rank sites by disruption")
    rk.add_argument("--in", dest="inp", required=True)
    rk.add_argument("--out", required=True)

    # Reporter
    rp = sub.add_parser("Reporter", help="Generate plots and HTML")
    rp.add_argument("--scores", required=True)
    rp.add_argument("--pred", required=True)
    rp.add_argument("--plots", required=True)
    rp.add_argument("--html", required=True)
    rp.add_argument("--transcript", required=False, help="Optional transcript ID to draw gene structure")
    rp.add_argument("--gtf", required=False, help="Optional local GTF for exon structure")
    rp.add_argument("--intron-bed", required=False, help="BED with highlighted intron interval")

    # TranscriptToBED
    tb = sub.add_parser("TranscriptToBED", help="Create intron BED from a transcript ID")
    tb.add_argument("--transcript", required=True)
    tb.add_argument("--intron-index", type=int, default=2)
    tb.add_argument("--gtf", default=None)
    tb.add_argument("--out", required=True)

    # Full pipeline: VariantBuilder -> Scorer -> Ranker -> Reporter
    full = sub.add_parser("Full", help="Run VariantBuilder → Scorer → Ranker → Reporter")
    full.add_argument("--config", default="ag.yaml")
    full.add_argument("--intron-bed", required=True)
    full.add_argument("--cassette", required=False, default=None, help="Single cassette FASTA path (can be combined with --cassette-list/--cassette-dir)")
    full.add_argument("--cassette-list", nargs='+', default=None, help="One or more cassette FASTA paths to run sequentially")
    full.add_argument("--cassette-dir", default=None, help="Directory containing cassette FASTA files to run sequentially")
    full.add_argument("--multi-out-root", default=None, help="Base output directory for multi-cassette runs (defaults to io.out_dir)")
    full.add_argument("--buffers", nargs=4, type=int, default=None, help="DONOR BP_START BP_END ACCEPTOR (optional, falls back to ag.yaml)")
    full.add_argument("--stride", type=int, default=None)
    full.add_argument("--max", dest="max_candidates", type=int, default=None)
    full.add_argument("--candidates", default=None)
    full.add_argument("--modalities", nargs='*', default=None)
    full.add_argument("--variant-window", type=int, default=None)
    full.add_argument("--raw-out", default=None)
    full.add_argument("--scores-out", default=None)
    full.add_argument("--plots", default=None)
    full.add_argument("--html", default=None)
    full.add_argument("--transcript", required=False, help="Optional transcript ID to draw gene structure")
    full.add_argument("--gtf", required=False, help="Optional local GTF for exon structure")
    full.add_argument("--variant-type", choices=["insertion", "extraction"], default="insertion", help="Type of variant: insertion (add cassette) or extraction (delete cassette-length bases)")
    full.add_argument("--genome-fasta", required=False, help="Genome FASTA file (required for extraction)")

    args, rest = parser.parse_known_args(argv)

    if args.cmd == "VariantBuilder":
        if args.variant_type == "extraction":
            if not args.genome_fasta:
                raise SystemExit("--genome-fasta is required for extraction variant type")
            vb_ext.main(
                [
                    "--intron-bed", args.intron_bed,
                    "--cassette", args.cassette,
                    "--genome-fasta", args.genome_fasta,
                    "--buffers", *[str(x) for x in args.buffers],
                    "--stride", str(args.stride),
                    "--max", str(args.max_candidates),
                    "--out", args.out,
                ]
            )
        else:
            variant_builder.main(
                [
                    "--intron-bed", args.intron_bed,
                    "--cassette", args.cassette,
                    "--buffers", *[str(x) for x in args.buffers],
                    "--stride", str(args.stride),
                    "--max", str(args.max_candidates),
                    "--out", args.out,
                ]
            )
    elif args.cmd == "AlphaGenomeScorer":
        scorer.main(
            [
                "--config", args.config,
                "--candidates", args.candidates,
                *( ["--modalities", *args.modalities] if args.modalities else [] ),
                *( ["--variant-window", str(args.variant_window)] if args.variant_window else [] ),
                "--out", args.out,
            ]
        )
    elif args.cmd == "Ranker":
        ranker.main(["--in", args.inp, "--out", args.out])
    elif args.cmd == "Reporter":
        call = ["--scores", args.scores, "--pred", args.pred, "--plots", args.plots, "--html", args.html]
        if args.transcript:
            call.extend(["--transcript", args.transcript])
        if args.gtf:
            call.extend(["--gtf", args.gtf])
        if args.intron_bed:
            call.extend(["--intron-bed", args.intron_bed])
        reporter.main(call)
    elif args.cmd == "TranscriptToBED":
        call = ["--transcript", args.transcript, "--intron-index", str(args.intron_index), "--out", args.out]
        if args.gtf:
            call.extend(["--gtf", args.gtf])
        transcript_bed.main(call)
    elif args.cmd == "Full":
        cfg = AppConfig.from_yaml(args.config)

        intron_bed_path = Path(args.intron_bed).expanduser()
        if not intron_bed_path.exists():
            raise SystemExit(f"--intron-bed not found: {intron_bed_path}")
        intron_bed_str = str(intron_bed_path)

        if args.buffers is None:
            args.buffers = [
                cfg.buffers.donor_min_nt,
                cfg.buffers.branchpoint_window_nt_start,
                cfg.buffers.branchpoint_window_nt_end,
                cfg.buffers.acceptor_min_nt,
            ]
        if args.stride is None:
            args.stride = cfg.scan.stride_nt
        if args.max_candidates is None:
            args.max_candidates = cfg.scan.max_candidates

        genome_fasta_str = None
        if getattr(args, 'variant_type', 'insertion') == "extraction":
            genome_fasta_value = getattr(args, 'genome_fasta', None)
            if not genome_fasta_value:
                raise SystemExit("--genome-fasta is required for extraction variant type")
            genome_fasta_path = Path(genome_fasta_value).expanduser()
            if not genome_fasta_path.exists():
                raise SystemExit(f"--genome-fasta not found: {genome_fasta_path}")
            genome_fasta_str = str(genome_fasta_path)

        buffers_args = [str(x) for x in args.buffers]
        stride_str = str(args.stride)
        max_str = str(args.max_candidates)

        cassette_paths: List[Path] = []
        seen: Set[str] = set()

        def _add_cassette(path_like, source: str) -> None:
            if not path_like:
                return
            p = Path(path_like).expanduser()
            if not p.exists():
                raise SystemExit(f"{source}: cassette FASTA not found: {p}")
            if p.is_dir():
                raise SystemExit(f"{source}: expected FASTA file but found directory: {p}")
            key = str(p.resolve())
            if key in seen:
                return
            seen.add(key)
            cassette_paths.append(p)

        def _looks_like_fasta(name: str) -> bool:
            lowered = name.lower()
            return lowered.endswith((".fa", ".fasta", ".fna", ".fas"))

        def _add_from_dir(directory: str) -> None:
            dir_path = Path(directory).expanduser()
            if not dir_path.exists():
                raise SystemExit(f"--cassette-dir: directory not found: {dir_path}")
            if not dir_path.is_dir():
                raise SystemExit(f"--cassette-dir: expected directory but found file: {dir_path}")
            matched = sorted(p for p in dir_path.iterdir() if p.is_file() and _looks_like_fasta(p.name))
            if not matched:
                raise SystemExit(f"--cassette-dir: no FASTA files found in {dir_path}")
            for child in matched:
                _add_cassette(child, f"--cassette-dir {dir_path}")

        _add_cassette(args.cassette, "--cassette")
        if args.cassette_list:
            for idx, item in enumerate(args.cassette_list, start=1):
                _add_cassette(item, f"--cassette-list[{idx}]")
        if args.cassette_dir:
            _add_from_dir(args.cassette_dir)

        if not cassette_paths:
            for idx, item in enumerate(cfg.inputs.cassettes, start=1):
                _add_cassette(item, f"ag.yaml inputs.cassettes[{idx}]")
            if cfg.inputs.cassette:
                _add_cassette(cfg.inputs.cassette, "ag.yaml inputs.cassette")

        if not cassette_paths:
            raise SystemExit("No cassette FASTA provided: supply --cassette/--cassette-list/--cassette-dir or configure inputs.cassette(s) in ag.yaml")

        multi_mode = len(cassette_paths) > 1

        if multi_mode and any(val is not None for val in (args.candidates, args.raw_out, args.scores_out, args.plots, args.html)):
            raise SystemExit("Per-output overrides (--candidates/--raw-out/--scores-out/--plots/--html) are not supported with multi-cassette mode. Use --multi-out-root instead.")

        def _sanitize(name: str) -> str:
            sanitized = re.sub(r"[^A-Za-z0-9._-]+", "_", name).strip("_.")
            return sanitized or "cassette"

        transcript_arg = getattr(args, 'transcript', None) or (cfg.inputs.transcript if cfg.inputs.transcript else None)
        gtf_arg = getattr(args, 'gtf', None) or (str(cfg.inputs.gtf) if getattr(cfg.inputs, 'gtf', None) else None)

        def _run_single(*, cassette_path: Path, candidates: Path, raw: Path, scores: Path, plots: Path, html: Path) -> None:
            candidates.parent.mkdir(parents=True, exist_ok=True)
            raw.parent.mkdir(parents=True, exist_ok=True)
            scores.parent.mkdir(parents=True, exist_ok=True)
            html.parent.mkdir(parents=True, exist_ok=True)
            plots.mkdir(parents=True, exist_ok=True)

            if getattr(args, 'variant_type', 'insertion') == "extraction":
                vb_ext.main([
                    "--intron-bed", intron_bed_str,
                    "--cassette", str(cassette_path),
                    "--genome-fasta", genome_fasta_str,
                    "--buffers", *buffers_args,
                    "--stride", stride_str,
                    "--max", max_str,
                    "--out", str(candidates),
                ])
            else:
                variant_builder.main([
                    "--intron-bed", intron_bed_str,
                    "--cassette", str(cassette_path),
                    "--buffers", *buffers_args,
                    "--stride", stride_str,
                    "--max", max_str,
                    "--out", str(candidates),
                ])

            scorer_call = [
                "--config", args.config,
                "--candidates", str(candidates),
                "--out", str(raw),
            ]
            if args.modalities:
                scorer_call.extend(["--modalities", *args.modalities])
            if args.variant_window:
                scorer_call.extend(["--variant-window", str(args.variant_window)])
            scorer.main(scorer_call)

            ranker.main(["--in", str(raw), "--out", str(scores)])

            reporter_call = [
                "--scores", str(scores),
                "--pred", str(raw),
                "--plots", str(plots),
                "--html", str(html),
                "--intron-bed", intron_bed_str,
            ]
            if transcript_arg:
                reporter_call.extend(["--transcript", transcript_arg])
            if gtf_arg:
                reporter_call.extend(["--gtf", gtf_arg])
            reporter.main(reporter_call)

        if multi_mode:
            base_root = Path(args.multi_out_root).expanduser() if args.multi_out_root else Path(cfg.io.out_dir)
            base_root.mkdir(parents=True, exist_ok=True)
            total = len(cassette_paths)
            summaries = []
            for idx, cassette_path in enumerate(cassette_paths, start=1):
                label = _sanitize(cassette_path.stem or cassette_path.name)
                run_root = base_root / f"{idx:02d}_{label}"
                plots_dir = run_root / "plots"
                print(f"[Full] ({idx}/{total}) Running cassette {cassette_path} → {run_root}")
                _run_single(
                    cassette_path=cassette_path,
                    candidates=run_root / "candidates.tsv",
                    raw=run_root / "raw.parquet",
                    scores=run_root / "candidates.csv",
                    plots=plots_dir,
                    html=run_root / "report.html",
                )
                summaries.append((cassette_path, run_root))
            print("[Full] Completed multi-cassette run. Outputs:")
            for cassette_path, run_root in summaries:
                print(f"  - {cassette_path} → {run_root}")
        else:
            cassette_path = cassette_paths[0]
            candidates_path = Path(args.candidates).expanduser() if args.candidates else Path(cfg.io.candidates_tsv)
            raw_path = Path(args.raw_out).expanduser() if args.raw_out else Path(cfg.io.raw_parquet)
            scores_path = Path(args.scores_out).expanduser() if args.scores_out else Path(cfg.io.scores_csv)
            plots_path = Path(args.plots).expanduser() if args.plots else Path(cfg.io.plots_dir)
            html_path = Path(args.html).expanduser() if args.html else Path(cfg.io.report_html)
            print(f"[Full] Running cassette {cassette_path}")
            _run_single(
                cassette_path=cassette_path,
                candidates=candidates_path,
                raw=raw_path,
                scores=scores_path,
                plots=plots_path,
                html=html_path,
            )
    else:
        parser.print_help()
        return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
