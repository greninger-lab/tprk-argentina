#!/usr/bin/env bash

ALLELE_SET="$1"  # optional argument
KMER_LEN="$2"

../../count_kmers.sh "$KMER_LEN"
python3 ../../kmer_weighted_profiles.py --kmers_dir allele_kmers --abundances_csv "abundances.csv" --normalize --clr
PLOT_CMD=(python3 ../../plot_kmer_distances.py
  --dist_files dist_brays.tsv dist_jsd.tsv dist_clr_euclidean.tsv
  --pcoa_csv pcoa_brays_coords.csv
  --grouping_files ../../Stage.csv ../../Location.csv ../../Stage_Location.csv
  --output_dir kmer_plots
  --eigenvalues_csv pcoa_brays_eigenvalues.csv
  --ypos ../../ypositions.csv
)

# If ALLELE_SET is provided, add the --region flag
if [ -n "$ALLELE_SET" ]; then
  PLOT_CMD+=(--region "${ALLELE_SET}_k${KMER_LEN}")
fi

# Run plotting
echo "Running: ${PLOT_CMD[*]}"
"${PLOT_CMD[@]}"
