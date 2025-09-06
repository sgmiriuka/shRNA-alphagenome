from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import List, Tuple

from .config import AppConfig


@dataclass
class Intron:
    chrom: str
    start0: int  # 0-based, inclusive
    end0: int    # 0-based, exclusive
    strand: str = "."

    @property
    def width(self) -> int:
        return self.end0 - self.start0


def read_bed_first_interval(path: Path) -> Intron:
    with open(path, "r") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip().split("\t")
            chrom = parts[0]
            start0 = int(parts[1])
            end0 = int(parts[2])
            strand = parts[5] if len(parts) >= 6 else "."
            return Intron(chrom=chrom, start0=start0, end0=end0, strand=strand)
    raise ValueError(f"No intervals found in BED: {path}")


def read_fasta_single(path: Path) -> str:
    seq = []
    with open(path, "r") as f:
        for line in f:
            if not line.strip():
                continue
            if line.startswith(">") or line.lstrip().startswith("#"):
                continue
            seq.append(line.strip())
    sequence = "".join(seq).upper()
    if not sequence:
        raise ValueError(f"Empty sequence in FASTA: {path}")
    return sequence


def generate_candidate_positions(
    intron: Intron,
    donor_min_nt: int,
    bp_start: int,
    bp_end: int,
    acceptor_min_nt: int,
    stride: int,
) -> List[int]:
    # Allowed region is intron.start + donor_min .. intron.end - acceptor_min
    start_allowed = intron.start0 + donor_min_nt
    end_allowed = intron.end0 - acceptor_min_nt
    if start_allowed >= end_allowed:
        return []

    # Exclude the branchpoint window (end-bp_end .. end-bp_start)
    bp_excl_start = intron.end0 - bp_end
    bp_excl_end = intron.end0 - bp_start

    # Build allowed sub-intervals
    allowed_intervals: List[Tuple[int, int]] = []
    if start_allowed < bp_excl_start:
        allowed_intervals.append((start_allowed, bp_excl_start))
    if bp_excl_end < end_allowed:
        allowed_intervals.append((bp_excl_end, end_allowed))
    if not allowed_intervals:
        allowed_intervals.append((start_allowed, end_allowed))

    positions_1based: List[int] = []
    for (a, b) in allowed_intervals:
        # iterate 0-based half-open [a, b)
        pos0 = a
        while pos0 < b:
            # Convert to 1-based position for insertion site
            positions_1based.append(pos0 + 1)
            pos0 += stride
    return positions_1based


def write_candidates_tsv(
    out_path: Path,
    chrom: str,
    positions_1based: List[int],
    alt_seq: str,
    max_candidates: int,
) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as out:
        out.write("candidate_id\tchrom\tpos\tref\talt\n")
        for i, pos in enumerate(positions_1based[:max_candidates], start=1):
            cand_id = f"cand{i:04d}"
            out.write(f"{cand_id}\t{chrom}\t{pos}\t\t{alt_seq}\n")


def main(argv: List[str] | None = None) -> None:
    p = argparse.ArgumentParser(
        prog="VariantBuilder",
        description=(
            "Emit insertion variants within a user-provided interval (e.g., an intron) "
            "that respect configurable splice-signal buffers."
        ),
    )
    p.add_argument("--config", default=None, type=Path, help="Path to ag.yaml for defaults")
    p.add_argument("--intron-bed", required=False, type=Path, help="BED with the interval to scan (not restricted to introns)")
    p.add_argument("--cassette", required=False, type=Path, help="FASTA containing the insertion sequence (any length)")
    p.add_argument(
        "--buffers",
        nargs=4,
        metavar=("DONOR", "BP_START", "BP_END", "ACCEPTOR"),
        required=False,
        type=int,
        help="Buffers in nt: min distance to donor, branchpoint window start/end upstream of acceptor, min distance to acceptor. For non-intronic intervals set 0 0 0 0.",
    )
    p.add_argument("--stride", required=False, type=int, help="Scan stride in nt between candidate positions")
    p.add_argument("--max", dest="max_candidates", required=False, type=int, help="Max number of candidates to emit")
    p.add_argument("--out", required=False, type=Path, help="Output TSV path (defaults to io.candidates_tsv in ag.yaml)")
    p.add_argument("--focus", required=False, type=int, help="Optional 1-based genomic position to refine around")
    p.add_argument("--window", required=False, type=int, help="Half-window (nt) around --focus to keep (positions with |pos-focus| <= window)")

    args = p.parse_args(argv)

    cfg = AppConfig.from_yaml(args.config) if args.config else AppConfig()
    # Fill defaults from config if not provided
    intron_bed = args.intron_bed or (cfg.inputs.intron_bed if cfg.inputs.intron_bed else None)
    cassette = args.cassette or (cfg.inputs.cassette if cfg.inputs.cassette else None)
    if intron_bed is None or cassette is None:
        raise SystemExit("--intron-bed and --cassette are required (via CLI or ag.yaml inputs.*)")

    if args.buffers is not None:
        donor_min_nt, bp_start, bp_end, acceptor_min_nt = args.buffers
    else:
        donor_min_nt = cfg.buffers.donor_min_nt
        bp_start = cfg.buffers.branchpoint_window_nt_start
        bp_end = cfg.buffers.branchpoint_window_nt_end
        acceptor_min_nt = cfg.buffers.acceptor_min_nt

    stride = args.stride or cfg.scan.stride_nt
    max_candidates = args.max_candidates or cfg.scan.max_candidates
    out_path = args.out or cfg.io.candidates_tsv

    intron = read_bed_first_interval(Path(intron_bed))
    alt_seq = read_fasta_single(Path(cassette))

    positions = generate_candidate_positions(
        intron,
        donor_min_nt=donor_min_nt,
        bp_start=bp_start,
        bp_end=bp_end,
        acceptor_min_nt=acceptor_min_nt,
        stride=stride,
    )
    # Optional local refinement around a focal position
    if args.focus is not None and args.window is not None:
        positions = [p for p in positions if abs(int(p) - int(args.focus)) <= int(args.window)]
    if not positions:
        raise SystemExit("No candidate positions satisfy buffer constraints.")

    write_candidates_tsv(
        Path(out_path), chrom=intron.chrom, positions_1based=positions, alt_seq=alt_seq, max_candidates=max_candidates
    )

    print(f"[VariantBuilder] Wrote {min(len(positions), max_candidates)} candidates to {out_path}")


if __name__ == "__main__":
    main()
