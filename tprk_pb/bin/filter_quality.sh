#!/bin/bash
# filter_quality.sh
# Remove reads containing any base with quality < threshold
# Keep reads where ALL bases have quality >= threshold
# Usage: ./filter_quality.sh input.fastq output.fastq [min_quality]
set -euo pipefail

INPUT=$1
OUTPUT=$2
MIN_QUAL=${3:-10}

REJECTED="${OUTPUT%.fastq}.rejected.fastq"

# echo "Filtering reads with minimum quality threshold: Q${MIN_QUAL}"
# echo "Keeping reads where ALL bases have quality ≥ Q${MIN_QUAL}"
# echo "Rejecting reads with ANY base < Q${MIN_QUAL}"

# Filter reads - keep only those where ALL bases have quality >= MIN_QUAL
seqkit seq --min-qual "$MIN_QUAL" "$INPUT" > "$OUTPUT"

# Get IDs of passed reads
seqkit seq -n -i "$OUTPUT" > "${OUTPUT}.passed.ids"

# Get rejected reads by inverse matching
seqkit grep -v -f "${OUTPUT}.passed.ids" "$INPUT" > "$REJECTED"

# Count and report
TOTAL=$(seqkit stats -T "$INPUT" 2>/dev/null | awk 'NR==2{print $4}')
PASSED=$(seqkit stats -T "$OUTPUT" 2>/dev/null | awk 'NR==2{print $4}')
REJECTED_COUNT=$(seqkit stats -T "$REJECTED" 2>/dev/null | awk 'NR==2{print $4}')

# echo ""
# echo "=== Filtering Results ==="
# echo "Total input reads: ${TOTAL:-0}"
# echo "Passed (all bases ≥ Q${MIN_QUAL}): ${PASSED:-0} ($(awk -v p=${PASSED:-0} -v t=${TOTAL:-1} 'BEGIN{printf "%.1f", 100*p/t}')%)"
# echo "Rejected (≥1 base < Q${MIN_QUAL}): ${REJECTED_COUNT:-0} ($(awk -v r=${REJECTED_COUNT:-0} -v t=${TOTAL:-1} 'BEGIN{printf "%.1f", 100*r/t}')%)"
# echo "Sanity check: ${PASSED:-0} + ${REJECTED_COUNT:-0} = $((${PASSED:-0} + ${REJECTED_COUNT:-0})) (should equal ${TOTAL:-0})"

# echo ""
# echo "Passed reads: $OUTPUT"
# echo "Rejected reads: $REJECTED"

echo -e "q${MIN_QUAL}_pass\tq${MIN_QUAL}_fail"
echo -e "${PASSED:-0}\t${REJECTED_COUNT:-0}"

# Clean up
rm -f "${OUTPUT}.passed.ids"

# echo "Done."
