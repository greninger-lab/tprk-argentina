#!/bin/bash
# flip_to_forward.sh
# Flip orientation of PacBio HiFi reads so all are in reverse-strand orientation
# Only keeps reads with exactly one FP and one RP in consistent orientation
# Usage: ./flip_to_forward.sh input.fastq output.fw.fastq [0.2] [1400] [1700]
set -euo pipefail

INPUT=$1
OUTPUT=$2
ERROR_RATE=${3:-0.2}
MIN_LENGTH=${4:-1400}
MAX_LENGTH=${5:-1700}

FP="GGAAAGAAAAGAACCATACATCC"
RP="TCAGAATCCGGAACTGCG"

FP_MM=$(awk -v l=${#FP} -v e=$ERROR_RATE 'BEGIN{printf "%d", l*e}')
RP_MM=$(awk -v l=${#RP} -v e=$ERROR_RATE 'BEGIN{printf "%d", l*e}')

# Create rejected output filename
REJECTED="${OUTPUT%.fastq}.rejected.fastq"

echo "Normalizing orientation to reverse-strand with error rate $ERROR_RATE..."
echo "FP max mismatches: $FP_MM, RP max mismatches: $RP_MM"
echo "Length range: $MIN_LENGTH - $MAX_LENGTH bp (inclusive)"

# Locate primers in all reads (both strands)
echo "Locating primers..."

seqkit locate -i -p "$FP" -m "$FP_MM" "$INPUT" 2>/dev/null | awk 'NR>1 {print $1, $4, $5, "FP"}' > "${OUTPUT}.fp_loc.txt"
seqkit locate -i -p "$RP" -m "$RP_MM" "$INPUT" 2>/dev/null | awk 'NR>1 {print $1, $4, $5, "RP"}' > "${OUTPUT}.rp_loc.txt"

# Combine and process
cat "${OUTPUT}.fp_loc.txt" "${OUTPUT}.rp_loc.txt" | awk 'BEGIN {OFS="\t"} 
{ 
    seqid = $1
    strand = $2
    start = $3
    primer = $4
    
    # Store all positions for this seqid/primer/strand combination
    key = seqid SUBSEP primer SUBSEP strand
    count[key]++
    pos[key, count[key]] = start
    
    seqids[seqid] = 1
    primers_found[seqid, primer] = 1
} 
END { 
    # For each sequence, find the innermost primers (closest to center)
    for (seqid in seqids) { 
        # Check FP on plus strand (rightmost = largest position = innermost)
        fp_plus_key = seqid SUBSEP "FP" SUBSEP "+"
        if (count[fp_plus_key] > 0) {
            max_fp_plus = pos[fp_plus_key, 1]
            for (i = 2; i <= count[fp_plus_key]; i++) {
                if (pos[fp_plus_key, i] > max_fp_plus) {
                    max_fp_plus = pos[fp_plus_key, i]
                }
            }
            terminal_fp_plus[seqid] = max_fp_plus
            has_fp_plus[seqid] = 1
        }
        
        # Check RP on plus strand (leftmost = smallest position = innermost)
        rp_plus_key = seqid SUBSEP "RP" SUBSEP "+"
        if (count[rp_plus_key] > 0) {
            min_rp_plus = pos[rp_plus_key, 1]
            for (i = 2; i <= count[rp_plus_key]; i++) {
                if (pos[rp_plus_key, i] < min_rp_plus) {
                    min_rp_plus = pos[rp_plus_key, i]
                }
            }
            terminal_rp_plus[seqid] = min_rp_plus
            has_rp_plus[seqid] = 1
        }
        
        # Check FP on minus strand (leftmost = smallest position = innermost)
        fp_minus_key = seqid SUBSEP "FP" SUBSEP "-"
        if (count[fp_minus_key] > 0) {
            min_fp_minus = pos[fp_minus_key, 1]
            for (i = 2; i <= count[fp_minus_key]; i++) {
                if (pos[fp_minus_key, i] < min_fp_minus) {
                    min_fp_minus = pos[fp_minus_key, i]
                }
            }
            terminal_fp_minus[seqid] = min_fp_minus
            has_fp_minus[seqid] = 1
        }
        
        # Check RP on minus strand (rightmost = largest position = innermost)
        rp_minus_key = seqid SUBSEP "RP" SUBSEP "-"
        if (count[rp_minus_key] > 0) {
            max_rp_minus = pos[rp_minus_key, 1]
            for (i = 2; i <= count[rp_minus_key]; i++) {
                if (pos[rp_minus_key, i] > max_rp_minus) {
                    max_rp_minus = pos[rp_minus_key, i]
                }
            }
            terminal_rp_minus[seqid] = max_rp_minus
            has_rp_minus[seqid] = 1
        }
    }
    
    # Now determine orientation based on terminal primers
    for (seqid in seqids) {
        # Must have both FP and RP somewhere
        if (!((seqid, "FP") in primers_found) || !((seqid, "RP") in primers_found)) {
            exclude[seqid] = 1
            if (!((seqid, "FP") in primers_found)) {
                exclude_reason[seqid] = "missing_FP"
            } else {
                exclude_reason[seqid] = "missing_RP"
            }
            continue
        }
        
        # Plus strand: FP (left) before RP (right) = forward orientation, needs flipping
        if (has_fp_plus[seqid] && has_rp_plus[seqid]) {
            if (terminal_fp_plus[seqid] < terminal_rp_plus[seqid]) {
                orientation[seqid] = "forward"
            } else {
                exclude[seqid] = 1
                exclude_reason[seqid] = "wrong_order_plus"
            }
        }
        
        # Minus strand: RP (left) before FP (right) = reverse orientation, keep as-is
        if (has_fp_minus[seqid] && has_rp_minus[seqid]) {
            if (terminal_fp_minus[seqid] > terminal_rp_minus[seqid]) {
                if (orientation[seqid] == "forward") {
                    exclude[seqid] = 1
                    exclude_reason[seqid] = "both_orientations"
                    delete orientation[seqid]
                } else {
                    orientation[seqid] = "reverse"
                }
            } else {
                exclude[seqid] = 1
                exclude_reason[seqid] = "wrong_order_minus"
            }
        }
        
        # Mixed strands (one primer on + and one on -)
        if ((has_fp_plus[seqid] && has_rp_minus[seqid] && !has_rp_plus[seqid] && !has_fp_minus[seqid]) || 
            (has_fp_minus[seqid] && has_rp_plus[seqid] && !has_fp_plus[seqid] && !has_rp_minus[seqid])) {
            exclude[seqid] = 1
            exclude_reason[seqid] = "mixed_strands"
            delete orientation[seqid]
        }
    }
    
    for (seqid in orientation) { 
        if (!exclude[seqid]) { 
            print seqid, orientation[seqid]
        } 
    } 
    
    for (seqid in exclude) { 
        if (exclude[seqid]) { 
            reason = exclude_reason[seqid] 
            if (reason == "") reason = "unknown" 
            print seqid, "excluded", reason 
        } 
    } 
}' > "${OUTPUT}.orientations.txt"


# Split into forward, reverse, and problematic reads
awk '$2=="forward" {print $1}' "${OUTPUT}.orientations.txt" > "${OUTPUT}.forward.ids"
awk '$2=="reverse" {print $1}' "${OUTPUT}.orientations.txt" > "${OUTPUT}.reverse.ids"
awk '$2=="excluded" {print $1, $3}' "${OUTPUT}.orientations.txt" > "${OUTPUT}.excluded.txt"

# Extract reads that are already in reverse orientation
echo "Extracting reverse-strand orientation reads..."
if [ -s "${OUTPUT}.reverse.ids" ]; then
    seqkit grep -f "${OUTPUT}.reverse.ids" "$INPUT" > "${OUTPUT}.reverse.tmp"
else
    touch "${OUTPUT}.reverse.tmp"
fi

# Extract and flip forward orientation reads to reverse-strand
echo "Extracting and flipping forward orientation reads to reverse-strand..."
if [ -s "${OUTPUT}.forward.ids" ]; then
    seqkit grep -f "${OUTPUT}.forward.ids" "$INPUT" \
        | seqkit seq -t dna -r -p > "${OUTPUT}.flipped.tmp"
else
    touch "${OUTPUT}.flipped.tmp"
fi

# Merge
echo "Merging files..."
cat "${OUTPUT}.reverse.tmp" "${OUTPUT}.flipped.tmp" | seqkit rmdup -n > "${OUTPUT}.merged.tmp"

# Filter by length range and split into passed and rejected
echo "Filtering by length range ($MIN_LENGTH - $MAX_LENGTH bp)..."
seqkit seq -m "$MIN_LENGTH" -M "$MAX_LENGTH" "${OUTPUT}.merged.tmp" > "$OUTPUT"

# Get sequences outside length range
seqkit seq -M $((MIN_LENGTH - 1)) "${OUTPUT}.merged.tmp" > "${OUTPUT}.tooshort.tmp"
seqkit seq -m $((MAX_LENGTH + 1)) "${OUTPUT}.merged.tmp" > "${OUTPUT}.toolong.tmp"

# Combine all rejected reads
echo "Collecting rejected reads..."
cat "${OUTPUT}.tooshort.tmp" "${OUTPUT}.toolong.tmp" > "$REJECTED"
if [ -s "${OUTPUT}.excluded.txt" ]; then
    seqkit grep -f <(awk '{print $1}' "${OUTPUT}.excluded.txt") "$INPUT" >> "$REJECTED"
fi

# Count results
TOTAL=$(seqkit stats -T "$INPUT" 2>/dev/null | awk 'NR==2{print $4}')
FWD_COUNT=$(wc -l < "${OUTPUT}.forward.ids")
REV_COUNT=$(wc -l < "${OUTPUT}.reverse.ids")
EXCLUDED_COUNT=$(wc -l < "${OUTPUT}.excluded.txt")
TOO_SHORT_COUNT=$(seqkit stats -T "${OUTPUT}.tooshort.tmp" 2>/dev/null | awk 'NR==2{print $4}')
TOO_SHORT_COUNT=${TOO_SHORT_COUNT:-0}
TOO_LONG_COUNT=$(seqkit stats -T "${OUTPUT}.toolong.tmp" 2>/dev/null | awk 'NR==2{print $4}')
TOO_LONG_COUNT=${TOO_LONG_COUNT:-0}
FINAL_COUNT=$(seqkit stats -T "$OUTPUT" 2>/dev/null | awk 'NR==2{print $4}')
REJECTED_COUNT=$(seqkit stats -T "$REJECTED" 2>/dev/null | awk 'NR==2{print $4}')

echo ""
echo "=== Results ==="
echo "Total input reads: ${TOTAL:-0}"
echo "Already reverse-strand: $REV_COUNT"
echo "Forward-strand (flipped to reverse): $FWD_COUNT"
echo "Excluded (primers): $EXCLUDED_COUNT"
echo "Excluded (too short <$MIN_LENGTH): $TOO_SHORT_COUNT"
echo "Excluded (too long >$MAX_LENGTH): $TOO_LONG_COUNT"
echo "Final output reads: ${FINAL_COUNT:-0}"
echo "Total rejected reads: ${REJECTED_COUNT:-0}"

# Show breakdown of exclusions
if [ -s "${OUTPUT}.excluded.txt" ]; then
    echo ""
    echo "=== Primer exclusion reasons ==="
    awk '{reason[$2]++} END {for (r in reason) print r ": " reason[r]}' "${OUTPUT}.excluded.txt" | sort -t: -k2 -rn
fi

echo ""
echo "Output written to: $OUTPUT"
echo "Rejected reads written to: $REJECTED"

# Clean up
rm -f "${OUTPUT}.fp_loc.txt" "${OUTPUT}.rp_loc.txt" "${OUTPUT}.orientations.txt" \
      "${OUTPUT}.forward.ids" "${OUTPUT}.reverse.ids" \
      "${OUTPUT}.reverse.tmp" "${OUTPUT}.flipped.tmp" "${OUTPUT}.excluded.txt" \
      "${OUTPUT}.merged.tmp" "${OUTPUT}.tooshort.tmp" "${OUTPUT}.toolong.tmp"

echo "Done."
