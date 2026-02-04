#!/usr/bin/env bash
# ======================================================================
#  High-performance k-mer counting with Jellyfish
#  Usage: ./count_kmers.sh <KMER_LEN>
# ======================================================================

KMER_LEN="$1"

mkdir -p alleles_split
seqkit split "alleles.fasta" --by-id -O alleles_split

mkdir -p allele_jf allele_kmers
for f in alleles_split/*.fasta; do
  id=$(basename "$f" .fasta)
  jellyfish count -m "${KMER_LEN}" -s 10M -t 50 "$f" -o allele_jf/"$id".jf
  jellyfish dump -c allele_jf/"$id".jf > allele_kmers/"$id".kmers.tsv
done

# set -euo pipefail
#
# KMER_LEN="${1:-}"
# if [[ -z "$KMER_LEN" ]]; then
#   echo "Usage: $0 <KMER_LEN>"
#   exit 1
# fi
#
# INPUT_FASTA="alleles.fasta"
# SPLIT_DIR="alleles_split"
# JF_DIR="allele_jf"
# KMERS_DIR="allele_kmers"
#
# # ──────────────────────────────────────────────────────────────────────
# # 🧠 Auto-detect system resources
# # ──────────────────────────────────────────────────────────────────────
# TOTAL_CPUS=$(sysctl -n hw.logicalcpu)
# TOTAL_MEM_GB=$(($(sysctl -n hw.memsize) / 1024 / 1024 / 1024))
#
# # Reserve 10–20% headroom for macOS background processes
# RESERVED_CPUS=2
# RESERVED_MEM_GB=8
# USABLE_CPUS=$(( TOTAL_CPUS - RESERVED_CPUS ))
# USABLE_MEM_GB=$(( TOTAL_MEM_GB - RESERVED_MEM_GB ))
#
# # Decide concurrency
# THREADS_PER_JOB=10
# JOBS=$(( USABLE_CPUS / THREADS_PER_JOB ))
# HASH_SIZE_GB=$(( USABLE_MEM_GB / JOBS * 8 / 10 ))  # 80% per-job hash
# HASH_SIZE_OPT="${HASH_SIZE_GB}G"
#
# # ──────────────────────────────────────────────────────────────────────
# # ⚙️ Display configuration
# # ──────────────────────────────────────────────────────────────────────
# cat <<EOF
# ⚙️  Configuration
# ─────────────────────────────
# • K-mer length:       $KMER_LEN
# • CPU cores (total):  $TOTAL_CPUS
# • RAM (total):        ${TOTAL_MEM_GB} GB
# • Parallel jobs:      $JOBS
# • Threads per job:    $THREADS_PER_JOB
# • Hash size (-s):     $HASH_SIZE_OPT
# ─────────────────────────────
# EOF
#
# # ──────────────────────────────────────────────────────────────────────
# # 📂 Prepare directories
# # ──────────────────────────────────────────────────────────────────────
# mkdir -p "$SPLIT_DIR" "$JF_DIR" "$KMERS_DIR"
#
# # Only split if needed
# if [ -z "$(find "$SPLIT_DIR" -maxdepth 1 -name '*.fasta' -print -quit)" ]; then
#   echo "📂 Splitting $INPUT_FASTA into individual FASTA files..."
#   seqkit split "$INPUT_FASTA" --by-id -O "$SPLIT_DIR"
# else
#   echo "📂 Using existing split files in $SPLIT_DIR"
# fi
#
# # ──────────────────────────────────────────────────────────────────────
# # 🚀 Run Jellyfish in parallel
# # ──────────────────────────────────────────────────────────────────────
# echo "🧬 Counting k-mers with Jellyfish using $JOBS concurrent jobs..."
# export KMER_LEN THREADS_PER_JOB HASH_SIZE_OPT JF_DIR KMERS_DIR
#
# find "$SPLIT_DIR" -type f -name "*.fasta" | parallel -j "$JOBS" --halt soon,fail=1 '
#   f="{}"
#   id=$(basename "$f" .fasta)
#   echo "🧩 Processing $id ..."
#   jellyfish count -m "$KMER_LEN" -s "$HASH_SIZE_OPT" -t "$THREADS_PER_JOB" -C "$f" -o "$JF_DIR/$id.jf"
#   jellyfish dump -c "$JF_DIR/$id.jf" > "$KMERS_DIR/$id.kmers.tsv"
#   echo "✅ Finished: $id"
# '
#
# echo
# echo "🎉 All k-mer counting complete!"
# echo "Results saved to: $KMERS_DIR"
