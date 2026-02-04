#!/usr/bin/env python3
# dereplicate.py
# Find unique reads and create FASTA with copy counts
# Usage: ./dereplicate.py input.fastq output.fasta

import sys
from collections import Counter
import math

def round_sig(x, sig=2):
    return round(x, sig - int(math.floor(math.log10(abs(x))))) if x != 0 else 0

def main():
    if len(sys.argv) < 3:
        print("Usage: ./dereplicate.py input.fastq output.fasta qual")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    qual = sys.argv[3]
    
    # print("Dereplicating sequences...")
    # print(f"Input: {input_file}")
    # print(f"Output: {output_file}")
    
    # Read all sequences from FASTQ
    sequences = []
    with open(input_file, 'r') as f:
        line_num = 0
        for line in f:
            line_num += 1
            if line_num % 4 == 2:  # Sequence line in FASTQ (every 4th line starting at 2)
                sequences.append(line.strip())
    
    # Count occurrences of each unique sequence
    seq_counts = Counter(sequences)
    
    # Sort by count (descending) for better organization
    sorted_seqs = sorted(seq_counts.items(), key=lambda x: x[1], reverse=True)
    
    # Write FASTA with sequential names and counts
    with open(output_file, 'w') as f:
        for i, (seq, count) in enumerate(sorted_seqs, 1):
            f.write(f">seq{i}_{count}\n")
            f.write(f"{seq}\n")
    
    # Report statistics
    total_reads = len(sequences)
    unique_seqs = len(seq_counts)
    compression_ratio = total_reads / unique_seqs if unique_seqs > 0 else 0
    compression_ratio = round_sig(compression_ratio, 2)

    print(f"Q{qual}_total\tQ{qual}_unique\tQ{qual}_compression")
    print(f"{total_reads}\t{unique_seqs}\t{compression_ratio}")

    # print()
    # print("=== Results ===")
    # print(f"Total input reads: {total_reads}")
    # print(f"Unique sequences: {unique_seqs}")
    # print(f"Compression ratio: {compression_ratio:.1f}x")
    
    # # Show top 10 most abundant sequences
    # print()
    # print("=== Top 10 most abundant sequences ===")
    # for i, (seq, count) in enumerate(sorted_seqs[:10], 1):
    #     print(f"seq{i}: {count} copies")
    
    # print()
    # print(f"Output written to: {output_file}")
    # print("Done.")

if __name__ == "__main__":
    main()
    