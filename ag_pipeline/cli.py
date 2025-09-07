from __future__ import annotations

import argparse
import sys

from . import variant_builder, scorer, ranker, reporter
from . import transcript_bed
from .config import AppConfig


def main(argv=None):
    parser = argparse.ArgumentParser(prog="ag", description="AlphaGenome shRNA intron-2 pipeline")
    sub = parser.add_subparsers(dest="cmd", required=True)

    # VariantBuilder
    vb = sub.add_parser("VariantBuilder", help="Emit insertion variants within intron 2")
    vb.add_argument("--intron-bed", required=True)
    vb.add_argument("--cassette", required=True)
    vb.add_argument("--buffers", nargs=4, required=True)
    vb.add_argument("--stride", required=True)
    vb.add_argument("--max", dest="max_candidates", required=True)
    vb.add_argument("--out", required=True)

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
    full.add_argument("--cassette", required=True)
    full.add_argument("--buffers", nargs=4, type=int, default=None, help="DONOR BP_START BP_END ACCEPTOR (optional, falls back to ag.yaml)")
    full.add_argument("--stride", type=int, default=None)
    full.add_argument("--max", dest="max_candidates", type=int, default=None)
    full.add_argument("--candidates", default="data/candidates.tsv")
    full.add_argument("--modalities", nargs='*', default=None)
    full.add_argument("--variant-window", type=int, default=None)
    full.add_argument("--raw-out", default="ag_out/raw.parquet")
    full.add_argument("--scores-out", default="ag_out/candidates.csv")
    full.add_argument("--plots", default="ag_out/plots")
    full.add_argument("--html", default="ag_out/report.html")
    full.add_argument("--transcript", required=False, help="Optional transcript ID to draw gene structure")
    full.add_argument("--gtf", required=False, help="Optional local GTF for exon structure")

    args, rest = parser.parse_known_args(argv)

    if args.cmd == "VariantBuilder":
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
        # Load config to fill in defaults
        cfg = AppConfig.from_yaml(args.config)
        # Fill buffers / stride / max from config if not provided
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
        # Run VariantBuilder
        variant_builder.main([
            "--intron-bed", args.intron_bed,
            "--cassette", args.cassette,
            "--buffers", *[str(x) for x in args.buffers],
            "--stride", str(args.stride),
            "--max", str(args.max_candidates),
            "--out", args.candidates,
        ])
        # Run Scorer
        scorer_call = [
            "--config", args.config,
            "--candidates", args.candidates,
            "--out", args.raw_out,
        ]
        if args.modalities:
            scorer_call.extend(["--modalities", *args.modalities])
        if args.variant_window:
            scorer_call.extend(["--variant-window", str(args.variant_window)])
        scorer.main(scorer_call)
        # Run Ranker
        ranker.main(["--in", args.raw_out, "--out", args.scores_out])
        # Run Reporter
        reporter.main([
            "--scores", args.scores_out,
            "--pred", args.raw_out,
            "--plots", args.plots,
            "--html", args.html,
            *( ["--transcript", args.transcript] if getattr(args, 'transcript', None) else [] ),
            *( ["--gtf", args.gtf] if getattr(args, 'gtf', None) else [] ),
            "--intron-bed", args.intron_bed,
        ])
    else:
        parser.print_help()
        return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
