from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import List, Tuple, Dict

from .config import AppConfig


@dataclass
class Intron:
    """Represents a genomic intron interval.

    Attributes:
        chrom: Chromosome name.
        start0: 0-based inclusive start position.
        end0: 0-based exclusive end position.
        strand: Strand (+, -, or .).
    """
    chrom: str
    start0: int  # 0-based, inclusive
    end0: int    # 0-based, exclusive
    strand: str = "."

    @property
    def width(self) -> int:
        """Calculate intron width in nucleotides.

        Returns:
            int: Width of the intron.
        """
        return self.end0 - self.start0


def read_fasta_dict(path: Path) -> Dict[str, str]:
    """Read FASTA file into a dictionary of chrom: sequence.

    Args:
        path: Path to FASTA file.

    Returns:
        Dict[str, str]: Dictionary with chromosome names as keys and sequences as values.
    """
    fasta_dict = {}
    current_chrom = None
    seq_lines = []
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_chrom:
                    fasta_dict[current_chrom] = "".join(seq_lines).upper()
                current_chrom = line[1:].split()[0]  # Take first word as chrom name
                seq_lines = []
            else:
                seq_lines.append(line)
        if current_chrom:
            fasta_dict[current_chrom] = "".join(seq_lines).upper()
    return fasta_dict


def read_bed_first_interval(path: Path) -> Intron:
    """Read the first interval from a BED file.

    Args:
        path: Path to BED file.

    Returns:
        Intron: Parsed intron object.

    Raises:
        ValueError: If no valid intervals found.
    """
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
    """Read sequence from a single-entry FASTA file.

    Args:
        path: Path to FASTA file.

    Returns:
        str: Uppercase sequence string.

    Raises:
        ValueError: If sequence is empty.
    """
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
    cassette_len: int,
) -> List[int]:
    """Generate candidate deletion positions within an intron, respecting splice buffers.

    Args:
        intron: Intron interval.
        donor_min_nt: Minimum distance from donor.
        bp_start: Branchpoint window start.
        bp_end: Branchpoint window end.
        acceptor_min_nt: Minimum distance from acceptor.
        stride: Step size between positions.
        cassette_len: Length of the cassette to delete.

    Returns:
        List of 1-based candidate positions.
    """
    # Allowed region is intron.start + donor_min .. intron.end - acceptor_min - cassette_len + 1
    # To ensure the deletion fits within the intron
    start_allowed = intron.start0 + donor_min_nt
    end_allowed = intron.end0 - acceptor_min_nt - cassette_len + 1
    if start_allowed >= end_allowed:
        return []

    # Exclude the branchpoint window (end-bp_end .. end-bp_start)
    bp_excl_start = intron.end0 - bp_end
    bp_excl_end = intron.end0 - bp_start

    # Build allowed sub-intervals
    allowed_intervals: List[Tuple[int, int]] = []
    if start_allowed < bp_excl_start:
        allowed_intervals.append((start_allowed, min(bp_excl_start, end_allowed)))
    if bp_excl_end < end_allowed:
        allowed_intervals.append((max(bp_excl_end, start_allowed), end_allowed))
    if not allowed_intervals:
        allowed_intervals.append((start_allowed, end_allowed))

    positions_1based: List[int] = []
    for (a, b) in allowed_intervals:
        # iterate 0-based half-open [a, b)
        pos0 = a
        while pos0 < b:
            # Convert to 1-based position for deletion start site
            positions_1based.append(pos0 + 1)
            pos0 += stride
    return positions_1based


def write_candidates_tsv(
    out_path: Path,
    chrom: str,
    positions_1based: List[int],
    genome_seq: str,
    intron_start0: int,
    cassette_len: int,
    max_candidates: int,
) -> None:
    """Write candidate positions to TSV file for deletions.

    Args:
        out_path: Output TSV path.
        chrom: Chromosome name.
        positions_1based: List of 1-based positions.
        genome_seq: Genome sequence for the chromosome.
        intron_start0: 0-based intron start.
        cassette_len: Length of cassette (deletion length).
        max_candidates: Maximum number to write.
    """
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as out:
        out.write("candidate_id\tchrom\tpos\tref\talt\n")
        for i, pos in enumerate(positions_1based[:max_candidates], start=1):
            cand_id = f"cand{i:04d}"
            # 0-based start of deletion
            del_start0 = pos - 1
            # Relative to intron start
            rel_start = del_start0 - intron_start0
            rel_end = rel_start + cassette_len
            if rel_end > len(genome_seq):
                continue  # Skip if deletion extends beyond available sequence
            ref_seq = genome_seq[rel_start:rel_end]
            out.write(f"{cand_id}\t{chrom}\t{pos}\t{ref_seq}\t\n")


def main(argv: List[str] | None = None) -> None:
    """Generate candidate deletion positions within a genomic interval.

    Reads intron BED, cassette FASTA, and genome FASTA, generates positions respecting buffers,
    and writes candidates to TSV with extracted ref sequences.

    Args:
        argv: Command-line arguments. If None, uses sys.argv.
    """
    p = argparse.ArgumentParser(
        prog="VariantBuilderExtracting",
        description=(
            "Emit deletion variants within a user-provided interval (e.g., an intron) "
            "that respect configurable splice-signal buffers, extracting bases equal to cassette length."
        ),
    )
    p.add_argument("--config", default=None, type=Path, help="Path to ag.yaml for defaults")
    p.add_argument("--intron-bed", required=False, type=Path, help="BED with the interval to scan (not restricted to introns)")
    p.add_argument("--cassette", required=False, type=Path, help="FASTA containing the cassette sequence (to determine deletion length)")
    p.add_argument("--genome-fasta", required=False, type=Path, help="Genome FASTA file")
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
    genome_fasta = args.genome_fasta or (cfg.inputs.genome_fasta if hasattr(cfg.inputs, 'genome_fasta') and cfg.inputs.genome_fasta else None)
    if intron_bed is None or cassette is None or genome_fasta is None:
        raise SystemExit("--intron-bed, --cassette, and --genome-fasta are required (via CLI or ag.yaml inputs.*)")

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
    cassette_len = len(alt_seq)

    genome_dict = read_fasta_dict(Path(genome_fasta))
    if intron.chrom not in genome_dict:
        raise SystemExit(f"Chromosome {intron.chrom} not found in genome FASTA")
    genome_seq = genome_dict[intron.chrom][intron.start0:intron.end0]  # Extract intron sequence

    positions = generate_candidate_positions(
        intron,
        donor_min_nt=donor_min_nt,
        bp_start=bp_start,
        bp_end=bp_end,
        acceptor_min_nt=acceptor_min_nt,
        stride=stride,
        cassette_len=cassette_len,
    )
    # Optional local refinement around a focal position
    if args.focus is not None and args.window is not None:
        positions = [p for p in positions if abs(int(p) - int(args.focus)) <= int(args.window)]
    if not positions:
        raise SystemExit("No candidate positions satisfy buffer constraints.")

    write_candidates_tsv(
        Path(out_path), chrom=intron.chrom, positions_1based=positions, genome_seq=genome_seq,
        intron_start0=intron.start0, cassette_len=cassette_len, max_candidates=max_candidates
    )

    print(f"[VariantBuilderExtracting] Wrote {min(len(positions), max_candidates)} candidates to {out_path}")


if __name__ == "__main__":
    main()