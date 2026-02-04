#!/usr/bin/env bash

KMER_LEN="$1"

cp ../alleles.fasta .
cp ../abundances.csv .
../count_kmers.sh "$KMER_LEN"
python3 ../kmer_weighted_profiles.py --kmers_dir allele_kmers --abundances_csv "../abundances.csv" --normalize --clr
python3 ../plot_kmer_distances.py --dist_files dist_brays.tsv dist_jsd.tsv dist_clr_euclidean.tsv \
				--pcoa_csv pcoa_brays_coords.csv \
				--grouping_files ../Stage.csv ../Location.csv ../Stage_Location.csv \
				--output_dir kmer_plots \
				--eigenvalues_csv pcoa_brays_eigenvalues.csv \
				--region "tprK_k$KMER_LEN" \
				--ypos ../ypositions.csv

python3 ../split_by_region.py ./alleles.fasta
python3 ../split_abundances_by_region.py ./abundances.csv

# --- NEW LOOP FOR V1–V7 SUBDIRECTORIES ---
for region in V1 V2 V3 V4 V5 V6 V7; do
  if [ -d "$region" ]; then
    echo "Entering $region and running kmer_analysis_region.sh..."
    (
      cd "$region"
      ../../kmer_analysis_region.sh "$region" "$KMER_LEN"
    )
  else
    echo "Warning: directory $region not found — skipping"
  fi
done
