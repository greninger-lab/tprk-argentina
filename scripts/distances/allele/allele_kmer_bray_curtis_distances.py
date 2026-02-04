import pandas as pd
import numpy as np
from collections import Counter
from itertools import combinations
import os

def generate_kmers(sequence, k):
    """Generate k-mers from a sequence (forward direction only)"""
    kmers = []
    for i in range(len(sequence) - k + 1):
        kmers.append(sequence[i:i+k])
    return kmers

def calculate_bray_curtis(counts1, counts2):
    """Calculate Bray-Curtis dissimilarity between two k-mer count dictionaries"""
    all_kmers = set(counts1.keys()) | set(counts2.keys())
    
    sum_min = 0
    sum_total = 0
    
    for kmer in all_kmers:
        c1 = counts1.get(kmer, 0)
        c2 = counts2.get(kmer, 0)
        sum_min += min(c1, c2)
        sum_total += c1 + c2
    
    if sum_total == 0:
        return 0.0
    
    return 1 - (2 * sum_min / sum_total)

def process_sequences(csv_file, output_dir='kmer_analysis', save_kmers=True):
    """Process sequences and generate distance matrices"""
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Read the CSV file
    df = pd.read_csv(csv_file)
    
    # Define k-mer lengths
    k_lengths = [9, 11, 13, 15, 17, 19, 21, 23, 25, 27]
    
    # Define regions
    regions = ['V1', 'V2', 'V3', 'V4', 'V5', 'V6', 'V7']
    
    # Process each k-mer length
    for k in k_lengths:
        print(f"Processing k-mer length: {k}")
        
        # Create a subdirectory for this k-mer length
        k_dir = os.path.join(output_dir, f'k{k}')
        os.makedirs(k_dir, exist_ok=True)
        
        # Process each region
        for region in regions:
            print(f"  Processing region: {region}")
            
            # Generate k-mers for each sequence
            kmer_counts = {}
            for idx, row in df.iterrows():
                seq_code = row['seq_code']
                sequence = row[region]
                kmers = generate_kmers(sequence, k)
                kmer_counts[seq_code] = Counter(kmers)
            
            # Save k-mer counts if requested
            if save_kmers:
                kmer_output = []
                for seq_code, counts in kmer_counts.items():
                    for kmer, count in sorted(counts.items()):
                        kmer_output.append({
                            'seq_code': seq_code,
                            'kmer': kmer,
                            'count': count
                        })
                kmer_df = pd.DataFrame(kmer_output)
                kmer_file = os.path.join(k_dir, f'{region}_kmers.csv')
                kmer_df.to_csv(kmer_file, index=False)
                print(f"    Saved k-mers: {kmer_file}")
            
            # Calculate distance matrix
            seq_codes = list(kmer_counts.keys())
            n = len(seq_codes)
            dist_matrix = np.zeros((n, n))
            
            for i, j in combinations(range(n), 2):
                dist = calculate_bray_curtis(kmer_counts[seq_codes[i]], 
                                            kmer_counts[seq_codes[j]])
                dist_matrix[i, j] = dist
                dist_matrix[j, i] = dist
            
            # Create DataFrame and save
            dist_df = pd.DataFrame(dist_matrix, 
                                  index=seq_codes, 
                                  columns=seq_codes)
            output_file = os.path.join(k_dir, f'{region}_bray_curtis.csv')
            dist_df.to_csv(output_file)
            print(f"    Saved: {output_file}")
        
        # Process combined regions (all V1-V7 concatenated)
        print(f"  Processing combined regions")
        
        kmer_counts_combined = {}
        for idx, row in df.iterrows():
            seq_code = row['seq_code']
            # Concatenate all regions
            combined_sequence = ''.join([row[region] for region in regions])
            kmers = generate_kmers(combined_sequence, k)
            kmer_counts_combined[seq_code] = Counter(kmers)
        
        # Save k-mer counts for combined if requested
        if save_kmers:
            kmer_output_combined = []
            for seq_code, counts in kmer_counts_combined.items():
                for kmer, count in sorted(counts.items()):
                    kmer_output_combined.append({
                        'seq_code': seq_code,
                        'kmer': kmer,
                        'count': count
                    })
            kmer_df_combined = pd.DataFrame(kmer_output_combined)
            kmer_file_combined = os.path.join(k_dir, 'combined_kmers.csv')
            kmer_df_combined.to_csv(kmer_file_combined, index=False)
            print(f"    Saved k-mers: {kmer_file_combined}")
        
        # Calculate distance matrix for combined
        seq_codes = list(kmer_counts_combined.keys())
        n = len(seq_codes)
        dist_matrix_combined = np.zeros((n, n))
        
        for i, j in combinations(range(n), 2):
            dist = calculate_bray_curtis(kmer_counts_combined[seq_codes[i]], 
                                        kmer_counts_combined[seq_codes[j]])
            dist_matrix_combined[i, j] = dist
            dist_matrix_combined[j, i] = dist
        
        # Create DataFrame and save
        dist_df_combined = pd.DataFrame(dist_matrix_combined, 
                                       index=seq_codes, 
                                       columns=seq_codes)
        output_file = os.path.join(k_dir, 'combined_bray_curtis.csv')
        dist_df_combined.to_csv(output_file)
        print(f"    Saved: {output_file}")
    
    print(f"\nAll analysis complete! Results saved in '{output_dir}' directory")
    print(f"Generated {len(k_lengths)} k-mer lengths × {len(regions) + 1} matrices = {len(k_lengths) * (len(regions) + 1)} total distance matrices")

# Run the analysis
if __name__ == "__main__":
    # Input file should be named unique_alleles.csv
    csv_file = 'unique_alleles.csv'
    
    # Set save_kmers=True to save k-mer counts (default)
    # Set save_kmers=False to skip saving k-mers and only generate distance matrices
    process_sequences(csv_file, save_kmers=True)
    