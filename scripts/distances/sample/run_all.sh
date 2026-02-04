#!/usr/bin/env bash

for kmer_len in 9 11 13 15 17 19 21 23 25 27; do
	mkdir "tprK_k${kmer_len}"
	cd "tprK_k${kmer_len}"
	../kmer_analysis.sh "$kmer_len"
	cd ../
done
