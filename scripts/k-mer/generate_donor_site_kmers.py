#!/usr/bin/env python3
"""
GenBank Kmer Analysis Script
Analyzes kmers in GenBank sequences and checks if they fall within variable regions (V1-V7).
Handles reverse complementation for tprKds1 sequences.
"""

import sys
from Bio import SeqIO
from Bio.Seq import Seq


def reverse_complement(seq):
    """Return reverse complement of a sequence."""
    return str(Seq(seq).reverse_complement())


def convert_coordinates_for_rc(start, end, seq_length):
    """
    Convert feature coordinates from complement notation to reverse complement sequence coordinates.
    GenBank complement(start..end) means the feature is on positions start to end on the reverse strand.
    After reverse complementing, these positions map to new locations.
    """
    # In GenBank, positions are 1-based, but Python is 0-based
    # For a complement feature at positions start..end (1-based):
    # After RC, start maps to (seq_length - end + 1) and end maps to (seq_length - start + 1)
    new_start = seq_length - end + 1
    new_end = seq_length - start + 1
    return new_start, new_end


def extract_variable_regions(record, is_reverse_complemented, seq_length):
    """
    Extract V1-V7 variable region coordinates from GenBank record.
    Returns a list of tuples: (label, start, end) in 1-based coordinates.
    """
    variable_regions = []
    
    for feature in record.features:
        # Check if feature type is one of V1-V7
        if feature.type in ['V1', 'V2', 'V3', 'V4', 'V5', 'V6', 'V7']:
            # Get the label from qualifiers
            label = feature.qualifiers.get('label', [f"{feature.type}-unknown"])[0]
            
            # GenBank uses 1-based coordinates, BioPython converts to 0-based
            # We need to work with 1-based coordinates for output
            start = int(feature.location.start) + 1  # Convert back to 1-based
            end = int(feature.location.end)  # End is already exclusive, so this is correct
            
            # If this is a tprKds1 sequence and features are on complement strand
            # convert coordinates
            if is_reverse_complemented:
                start, end = convert_coordinates_for_rc(start, end, seq_length)
            
            variable_regions.append((label, start, end))
    
    return variable_regions


def is_kmer_in_region(kmer_start, kmer_end, regions):
    """
    Check if kmer (1-based coordinates) is completely within any variable region.
    Returns (True/False, list of matching region labels).
    """
    matching_labels = []
    
    for label, region_start, region_end in regions:
        # Check if kmer is completely within the region
        if kmer_start >= region_start and kmer_end <= region_end:
            matching_labels.append(label)
    
    return len(matching_labels) > 0, matching_labels


def generate_kmers_and_analyze(genbank_file, kmer_length):
    """
    Main function to process GenBank file and generate kmer analysis.
    """
    # Parse GenBank file
    records = list(SeqIO.parse(genbank_file, "genbank"))
    
    for record in records:
        # Get sequence name from LOCUS (record.name in BioPython)
        sequence_name = record.name
        sequence = str(record.seq)
        original_seq_length = len(sequence)
        
        # Determine if sequence should be reverse complemented
        is_tprKds1 = "tprKds1" in sequence_name
        
        # Process sequence
        if is_tprKds1:
            working_sequence = reverse_complement(sequence)
        else:
            working_sequence = sequence
        
        # Extract variable regions
        variable_regions = extract_variable_regions(record, is_tprKds1, original_seq_length)
        
        # Generate output filenames
        output_file = f"{sequence_name}_k{kmer_length}_analysis.tsv"
        fasta_file = f"{sequence_name}_k{kmer_length}_kmers.fasta"
        
        # Open output files
        with open(output_file, 'w') as out, open(fasta_file, 'w') as fasta:
            # Write header for TSV
            out.write("kmer\tstart\tend\tin_ds\tds\n")
            
            # Generate kmers
            for i in range(len(working_sequence) - kmer_length + 1):
                kmer = working_sequence[i:i + kmer_length]
                
                # Calculate coordinates in the working sequence (1-based)
                kmer_start_working = i + 1
                kmer_end_working = i + kmer_length
                
                # Convert back to original coordinates if sequence was reverse complemented
                if is_tprKds1:
                    # Map back to original coordinates
                    kmer_end_original = original_seq_length - kmer_end_working + 1
                    kmer_start_original = original_seq_length - kmer_start_working + 1
                else:
                    kmer_start_original = kmer_start_working
                    kmer_end_original = kmer_end_working
                
                # Check if kmer is in any variable region (using working coordinates)
                in_region, matching_labels = is_kmer_in_region(
                    kmer_start_working, kmer_end_working, variable_regions
                )
                
                # Format output for TSV
                in_ds = "Y" if in_region else "N"
                ds_labels = ",".join(matching_labels) if matching_labels else ""
                
                out.write(f"{kmer}\t{kmer_start_original}\t{kmer_end_original}\t{in_ds}\t{ds_labels}\n")
                
                # Write to FASTA file (kmer counter is i+1 for 1-based numbering)
                fasta.write(f">{sequence_name}_kmer_{i+1}\n{kmer}\n")
        
        print(f"Analysis complete for {sequence_name}.")
        print(f"  TSV output: {output_file}")
        print(f"  FASTA output: {fasta_file}")


def main():
    """Main entry point."""
    if len(sys.argv) != 3:
        print("Usage: python3 generate_donor_site_kmers.py <genbank_file> <kmer_length>")
        sys.exit(1)
    
    genbank_file = sys.argv[1]
    try:
        kmer_length = int(sys.argv[2])
        if kmer_length < 1:
            raise ValueError("Kmer length must be positive")
    except ValueError as e:
        print(f"Error: Invalid kmer length - {e}")
        sys.exit(1)
    
    try:
        generate_kmers_and_analyze(genbank_file, kmer_length)
    except Exception as e:
        print(f"Error processing file: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
    