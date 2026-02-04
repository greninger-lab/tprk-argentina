#!/bin/bash
# check_primers.sh
# Summarize primer counts in input, flipped, and trimmed FASTQ files
# Usage: ./check_primers.sh input.fastq flipped.fastq trimmed.fastq 0.2
set -euo pipefail

INPUT=$1
FLIPPED=$2
TRIMMED=$3
ERROR_RATE=${4:-0.2}

FP="GGAAAGAAAAGAACCATACATCC"
RP="TCAGAATCCGGAACTGCG"

# Compute max mismatches
FP_MM=$(awk -v l=${#FP} -v e=$ERROR_RATE 'BEGIN{printf "%d", l*e}')
RP_MM=$(awk -v l=${#RP} -v e=$ERROR_RATE 'BEGIN{printf "%d", l*e}')

count_primers() {
    local FILE=$1
    
    # Total reads - skip header, get num_seqs column (field 4)
    local TOTAL=$(seqkit stats -T "$FILE" 2>/dev/null | awk 'NR==2{print $4}')
    TOTAL=${TOTAL:-0}
    
    # FP matches (search sequence, both strands)
    local FP_COUNT=$(seqkit grep -s -i --by-seq -p "$FP" -m "$FP_MM" --quiet "$FILE" 2>/dev/null | seqkit stats -T 2>/dev/null | awk 'NR==2{print $4}')
    FP_COUNT=${FP_COUNT:-0}
    
    # RP matches (search sequence, both strands)
    local RP_COUNT=$(seqkit grep -s -i --by-seq -p "$RP" -m "$RP_MM" --quiet "$FILE" 2>/dev/null | seqkit stats -T 2>/dev/null | awk 'NR==2{print $4}')
    RP_COUNT=${RP_COUNT:-0}
    
    # Reads containing both FP and RP
    local BOTH_COUNT=$(seqkit grep -s -i --by-seq -p "$FP" -m "$FP_MM" --quiet "$FILE" 2>/dev/null \
                        | seqkit grep -s -i --by-seq -p "$RP" -m "$RP_MM" --quiet 2>/dev/null \
                        | seqkit stats -T 2>/dev/null | awk 'NR==2{print $4}')
    BOTH_COUNT=${BOTH_COUNT:-0}
    
    echo "$TOTAL $FP_COUNT $RP_COUNT $BOTH_COUNT"
}

read TOTAL_IN FP_IN RP_IN BOTH_IN <<< $(count_primers "$INPUT")
read TOTAL_FW FP_FW RP_FW BOTH_FW <<< $(count_primers "$FLIPPED")
read TOTAL_TRIM FP_TRIM RP_TRIM BOTH_TRIM <<< $(count_primers "$TRIMMED")

echo -e "total\ttotal_fp\ttotal_rp\ttotal_fprp\ttotal_fw\ttotal_fw_fp\ttotal_fw_rp\ttotal_fw_fprp\ttotal_trim\ttotal_trim_fp\ttotal_trim_rp\ttotal_trim_fprp"
echo -e "${TOTAL_IN}\t${FP_IN}\t${RP_IN}\t${BOTH_IN}\t${TOTAL_FW}\t${FP_FW}\t${RP_FW}\t${BOTH_FW}\t${TOTAL_TRIM}\t${FP_TRIM}\t${RP_TRIM}\t${BOTH_TRIM}"
