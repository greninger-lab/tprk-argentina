#!/usr/bin/env python3
"""
GenBank Kmer Analysis Script
Analyzes kmers in fasta sequences.
"""

import sys
from Bio import SeqIO
from Bio.Seq import Seq


def generate_kmers_and_analyze(fasta_file, kmer_length):
    """
    Main function to process GenBank file and generate kmer analysis.
    """
    # Parse GenBank file
    records = list(SeqIO.parse(fasta_file, "fasta"))
    
    for record in records:
        # Get sequence name from LOCUS (record.name in BioPython)
        sequence_name = record.name
        sequence = str(record.seq)
        original_seq_length = len(sequence)
        
        working_sequence = sequence
        
        # Generate output filenames
        fasta_file = f"{sequence_name}_k{kmer_length}_kmers.fasta"
        
        # Open output files
        with open(fasta_file, 'w') as fasta:
            
            # Generate kmers
            for i in range(len(working_sequence) - kmer_length + 1):
                kmer = working_sequence[i:i + kmer_length]
                
                # Write to FASTA file (kmer counter is i+1 for 1-based numbering)
                fasta.write(f">{sequence_name}_kmer_{i+1}\n{kmer}\n")
        
        print(f"Analysis complete for {sequence_name}.")
        print(f"  FASTA output: {fasta_file}")


def main():
    """Main entry point."""
    if len(sys.argv) != 3:
        print("Usage: python3 generate_variable_region_kmers.py <fasta_file> <kmer_length>")
        sys.exit(1)
    
    fasta_file = sys.argv[1]
    try:
        kmer_length = int(sys.argv[2])
        if kmer_length < 1:
            raise ValueError("Kmer length must be positive")
    except ValueError as e:
        print(f"Error: Invalid kmer length - {e}")
        sys.exit(1)
    
    try:
        generate_kmers_and_analyze(fasta_file, kmer_length)
    except Exception as e:
        print(f"Error processing file: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
    