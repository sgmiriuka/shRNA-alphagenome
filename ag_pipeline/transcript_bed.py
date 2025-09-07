from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path
from typing import List, Tuple


def _fetch_exons_rest_ensembl(transcript_id: str) -> tuple[str, int, List[dict]]:
    """Fetch exon structures from Ensembl REST (GRCh38) for a transcript.

    Returns tuple (chrom, strand, exons_list) where exons have 'start','end','rank' (if available).
    """
    import urllib.request

    url = f"https://rest.ensembl.org/lookup/id/{transcript_id}?expand=1"
    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=30) as r:
        data = json.loads(r.read().decode("utf-8"))
    chrom = data.get("seq_region_name")
    strand = int(data.get("strand", 1))
    exons = data.get("Exon", [])
    if not exons:
        raise ValueError(f"No exons found via Ensembl REST for {transcript_id}")
    return chrom, strand, exons


def _fetch_exons_from_gtf(gtf_path: Path, transcript_id: str) -> tuple[str, int, List[dict]]:
    """Parse a GTF file and extract exons for transcript id.

    Returns (chrom, strand(+1/-1), exons as dicts with start/end and optional rank).
    """
    chrom = None
    strand_char = '+'
    exons: List[dict] = []
    attr_re = re.compile(r'(\w+) "([^"]+)"')
    with open(gtf_path, "r") as f:
        for line in f:
            if not line or line.startswith('#'):
                continue
            chrom_f, src, feat, start, end, score, strand, frame, attrs = line.rstrip().split('\t')
            if feat != 'exon':
                continue
            d = {m.group(1): m.group(2) for m in attr_re.finditer(attrs)}
            if d.get('transcript_id') != transcript_id:
                continue
            chrom = chrom_f
            strand_char = strand
            # Convert to 1-based inclusive coordinates from GTF
            start0 = int(start)  # GTF start is 1-based; keep 1-based for now
            end0 = int(end)
            # rank sometimes present as exon_number
            rank = None
            if 'exon_number' in d:
                try:
                    rank = int(d['exon_number'])
                except Exception:
                    rank = None
            exons.append({'start': start0, 'end': end0, 'rank': rank})
    if not exons:
        raise ValueError(f"Transcript {transcript_id} not found in {gtf_path}")
    strand = 1 if strand_char == '+' else -1
    return chrom or '', strand, exons


def _sort_exons(exons: List[dict], strand: int) -> List[dict]:
    """Sort exons by rank or genomic position.

    Args:
        exons: List of exon dictionaries with 'start', 'end', 'rank'.
        strand: Strand (1 for +, -1 for -).

    Returns:
        Sorted list of exons.
    """
    if all('rank' in e and isinstance(e['rank'], int) for e in exons):
        return sorted(exons, key=lambda e: e['rank'])
    # fallback by genomic order, strand-aware
    return sorted(exons, key=lambda e: e['start'], reverse=(strand != 1))


def compute_intron_bed_line(transcript_id: str, intron_index: int = 2, *, gtf: Path | None = None) -> str:
    """Compute a BED line for the given transcript's intron index (default 2).

    Uses Ensembl REST if `gtf` is not provided; otherwise parses the GTF.
    Returns a single-line BED string.
    """
    if gtf:
        chrom, strand, exons = _fetch_exons_from_gtf(gtf, transcript_id)
    else:
        chrom, strand, exons = _fetch_exons_rest_ensembl(transcript_id)

    exons_sorted = _sort_exons(exons, strand)
    if len(exons_sorted) < intron_index + 0:  # ensure enough exons for requested intron
        # Need at least intron_index+1 exons (intron between exon i and i+1). For intron 2, >=3 exons
        if len(exons_sorted) < (intron_index + 1):
            raise ValueError(
                f"Transcript has only {len(exons_sorted)} exons; cannot get intron {intron_index}."
            )

    e_up = exons_sorted[intron_index - 1]
    e_dn = exons_sorted[intron_index]

    # BED coordinates are 0-based start, end exclusive.
    # Intron spans the gap between exon_up and exon_dn in genomic coordinates.
    # For any strand, the intron range in reference orientation is:
    #   start0 = min(e_up.end, e_dn.end)
    #   end0   = max(e_up.start, e_dn.start)
    start0 = min(int(e_up['end']), int(e_dn['end']))
    end0 = max(int(e_up['start']), int(e_dn['start']))
    if end0 <= start0:
        raise ValueError("Computed intron is empty or negative width.")

    # Normalize chromosome
    if chrom and not chrom.startswith('chr'):
        chrom = 'chr' + chrom
    strand_char = '+' if strand == 1 else '-'
    name = f"intron{intron_index}"
    return f"{chrom}\t{start0}\t{end0}\t{name}\t0\t{strand_char}"


def main(argv: List[str] | None = None) -> None:
    """Create BED line for a transcript's intron.

    Fetches exon data and computes intron coordinates, writing to BED file.

    Args:
        argv: Command-line arguments. If None, uses sys.argv.
    """
    ap = argparse.ArgumentParser(
        prog="TranscriptToBED",
        description="Create intron BED line for a transcript (default intron 2).",
    )
    ap.add_argument("--transcript", required=True, help="Ensembl transcript ID, e.g., ENST00000335295")
    ap.add_argument("--intron-index", type=int, default=2, help="Intron index (1-based; intron 2 is between exon 2 and 3)")
    ap.add_argument("--gtf", type=Path, default=None, help="Optional local GTF to parse instead of Ensembl REST")
    # If not provided, we'll auto-build an output filename from the feature name
    # contained in the BED line (e.g., intron2, exon3, promoter), to stay agnostic
    # to class and index.
    ap.add_argument(
        "--out",
        type=Path,
        default=None,
        help="Output BED path. If omitted, auto-named as data/gene_<feature>_hg38.bed",
    )

    args = ap.parse_args(argv)
    line = compute_intron_bed_line(args.transcript, args.intron_index, gtf=args.gtf)

    # Build default output name if not provided, using the feature name from the BED line
    # (4th column). This keeps naming agnostic to class (e.g., intron, exon, promoter)
    # and index (e.g., intron2, exon5).
    if args.out is None:
        parts = line.split("\t")
        feature = parts[3] if len(parts) >= 4 and parts[3] else "feature"
        out_path = Path("data") / f"gene_{feature}_hg38.bed"
    else:
        out_path = args.out

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as f:
        f.write(line + "\n")
    print(f"[TranscriptToBED] Wrote {out_path}:\n{line}")


if __name__ == "__main__":
    main()
