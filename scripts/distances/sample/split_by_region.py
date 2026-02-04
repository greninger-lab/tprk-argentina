#!/usr/bin/env python3
"""
split_by_region.py

Split a multi-FASTA file whose sequence headers are of the form {region}-{code}
into per-region FASTA files.

Each output directory (named after the region) will contain:
    alleles.fasta   -> sequences belonging to that region
"""

import os
import sys
from collections import defaultdict

def parse_fasta(filepath):
    """Generator that yields (header, sequence) tuples from a FASTA file."""
    header = None
    seq_lines = []
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header:
                    yield (header, "".join(seq_lines))
                header = line[1:].strip()
                seq_lines = []
            else:
                seq_lines.append(line)
        if header:
            yield (header, "".join(seq_lines))

def main(infile):
    # Collect sequences by region prefix
    regions = defaultdict(list)

    for header, seq in parse_fasta(infile):
        # Split header like "V1-000126" -> region="V1", code="000126"
        if "-" not in header:
            sys.stderr.write(f"Warning: header '{header}' missing '-', skipping.\n")
            continue
        region = header.split("-", 1)[0]
        regions[region].append((header, seq))

    # Create per-region folders and write alleles.fasta
    for region, entries in regions.items():
        os.makedirs(region, exist_ok=True)
        outpath = os.path.join(region, "alleles.fasta")
        with open(outpath, "w") as out:
            for header, seq in entries:
                out.write(f">{header}\n")
                # Write sequences wrapped at 60 nt per line
                # for i in range(0, len(seq), 60):
                out.write(seq + "\n")
        print(f"Wrote {len(entries)} sequences to {outpath}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit("Usage: ./split_by_region.py input.fasta")
    main(sys.argv[1])
