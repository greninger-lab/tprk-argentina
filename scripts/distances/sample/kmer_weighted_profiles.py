#!/usr/bin/env python3
"""
kmer_weighted_profiles.py (region-aware version)

Inputs:
 - allele_kmers_dir: dir containing files named <allele_id>.kmers.tsv  (kmer[TAB]count)
 - abundances.csv: columns like sample,region,v_code,av (others ignored)
   * av can be 0–1 or 0–100 (%)
   * script normalizes within each region per sample
   * all regions combined to a single per-sample profile

Outputs:
 - dist_brays.tsv
 - dist_jsd.tsv
 - (optional) dist_clr_euclidean.tsv
 - pcoa_brays_coords.csv
 - pcoa_brays_eigenvalues.csv
"""

import os, sys
import argparse
from collections import defaultdict
import pandas as pd
import numpy as np
from scipy.spatial.distance import pdist, squareform

# --------------------------------------------------------------------
def read_allele_kmers(kdir):
    allele_kmers = {}
    for fname in os.listdir(kdir):
        if not (fname.endswith(".kmers.tsv") or fname.endswith(".kmers.txt") or fname.endswith(".kmers")):
            continue
        fname_base = fname.rsplit(".kmers", 1)[0]
        if "part_" in fname_base:
            allele_id = fname_base.split("part_")[-1]
        else:
            allele_id = fname_base
        d = {}
        path = os.path.join(kdir, fname)
        with open(path) as fh:
            for line in fh:
                parts = line.strip().split()
                if len(parts) < 2:
                    continue
                kmer = parts[0]
                try:
                    cnt = float(parts[1])
                except:
                    cnt = 0.0
                if cnt > 0:
                    d[kmer] = cnt
        allele_kmers[allele_id] = d
    if not allele_kmers:
        print("Warning: No allele kmer files found in directory", file=sys.stderr)
    return allele_kmers

# --------------------------------------------------------------------
def build_sample_vectors(allele_kmers, abund_df):
    # Flexible column handling
    colmap = {}
    for c in abund_df.columns:
        lc = c.lower()
        if lc in ("sample", "sample_id"):
            colmap["sample"] = c
        elif lc in ("allele_id", "v_code", "vcode", "allele"):
            colmap["allele_id"] = c
        elif lc in ("av", "abundance", "freq", "fraction"):
            colmap["abundance"] = c
        elif lc in ("region", "segment"):
            colmap["region"] = c
    missing = {"sample", "allele_id", "abundance"} - set(colmap.keys())
    if missing:
        raise ValueError(f"Missing required columns in abundance file: {missing}")

    df = abund_df.copy()
    df.rename(columns={colmap["sample"]: "sample",
                       colmap["allele_id"]: "allele_id",
                       colmap["abundance"]: "abundance"}, inplace=True)
    if "region" in colmap:
        df.rename(columns={colmap["region"]: "region"}, inplace=True)
    else:
        df["region"] = "all"

    # Handle percentage abundances
    if df["abundance"].max() > 1.0:
        print("Detected abundances > 1.0 — treating as percentages and dividing by 100.")
        df["abundance"] = df["abundance"] / 100.0

    samples = sorted(df["sample"].unique())
    out = {}

    for s in samples:
        sub = df[df["sample"] == s].copy()
        sub["region_total"] = sub.groupby("region")["abundance"].transform("sum")
        sub["norm_abund"] = sub["abundance"] / sub["region_total"]

        vec = defaultdict(float)
        for _, row in sub.iterrows():
            allele = row["allele_id"]
            w = float(row["norm_abund"])
            if allele not in allele_kmers:
                print(f"Warning: allele {allele} not found for sample {s}, skipping", file=sys.stderr)
                continue
            for kmer, cnt in allele_kmers[allele].items():
                vec[kmer] += cnt * w

        if not vec:
            print(f"Warning: sample {s} produced empty kmer vector", file=sys.stderr)
        out[s] = dict(vec)

    return out

# --------------------------------------------------------------------
def make_dense_matrix(sample_kmers, kmer_whitelist=None):
    samples = list(sample_kmers.keys())
    if kmer_whitelist is None:
        kset = set()
        for s in samples:
            kset.update(sample_kmers[s].keys())
        kmers = sorted(kset)
    else:
        kmers = list(kmer_whitelist)
    df = pd.DataFrame(0.0, index=samples, columns=kmers, dtype=float)
    for s in samples:
        for k, v in sample_kmers[s].items():
            if k in df.columns:
                df.at[s, k] = v
    if df.shape[1] == 0:
        print("Warning: no kmers found across all samples — distance matrices will be zero", file=sys.stderr)
    if (df.nunique(axis=0) == 1).all():
        print("Warning: all sample vectors are identical — distances will be zero", file=sys.stderr)
    return df

# --------------------------------------------------------------------
def bray_curtis_from_matrix(mat):
    arr = mat.values
    if (mat.sum(axis=1) == 0).any():
        print("Warning: some samples have zero total counts; distances may be unreliable", file=sys.stderr)
    d = pdist(arr, metric='braycurtis')
    return pd.DataFrame(squareform(d), index=mat.index, columns=mat.index)

# --------------------------------------------------------------------
def safe_jensen_shannon(P, base=2.0):
    n = P.shape[0]
    d = np.zeros((n*(n-1))//2)
    idx = 0
    for i in range(n):
        pi = P[i]
        for j in range(i+1, n):
            pj = P[j]
            m = 0.5 * (pi + pj)
            kl1 = np.sum(np.where(pi != 0, pi * np.log(pi/m), 0))
            kl2 = np.sum(np.where(pj != 0, pj * np.log(pj/m), 0))
            js = 0.5 * (kl1 + kl2) / np.log(base)
            js = max(js, 0.0)
            d[idx] = np.sqrt(js / 2.0)
            idx += 1
    return d

# --------------------------------------------------------------------
def jensen_shannon_from_matrix(mat, pseudocount=1e-12):
    arr = mat.values.astype(float) + pseudocount
    row_sums = arr.sum(axis=1, keepdims=True)
    zero_rows = (row_sums.flatten() == 0)
    if zero_rows.any():
        print("Warning: some samples are zero, replacing with uniform distribution", file=sys.stderr)
        arr[zero_rows, :] = 1.0 / arr.shape[1]
        row_sums = arr.sum(axis=1, keepdims=True)
    P = arr / row_sums
    n = P.shape[0]
    d = safe_jensen_shannon(P)
    return pd.DataFrame(squareform(d), index=mat.index, columns=mat.index)

# --------------------------------------------------------------------
def clr_transform(mat, pseudocount=1e-6):
    arr = mat.values + pseudocount
    loga = np.log(arr)
    means = loga.mean(axis=1, keepdims=True)
    clr = loga - means
    if np.allclose(clr.var(axis=1), 0):
        print("Warning: CLR transform produced zero variance for all samples", file=sys.stderr)
    return pd.DataFrame(clr, index=mat.index, columns=mat.columns)

# --------------------------------------------------------------------
def pcoa_from_distance(dm, n_components=2):
    arr = dm.values
    n = arr.shape[0]
    H = np.eye(n) - np.ones((n, n))/n
    B = -0.5 * H.dot(arr**2).dot(H)
    eigvals, eigvecs = np.linalg.eigh(B)
    idx = np.argsort(eigvals)[::-1]
    eigvals = eigvals[idx]
    eigvecs = eigvecs[:, idx]
    if not np.any(eigvals > 0):
        print("Warning: no positive eigenvalues in PCoA; check distance matrix", file=sys.stderr)
    pos = eigvals > 0
    L = np.diag(np.sqrt(eigvals[pos]))
    coords = eigvecs[:, pos].dot(L)[:, :n_components]
    return pd.DataFrame(coords, index=dm.index, columns=[f'PC{i+1}' for i in range(coords.shape[1])]), eigvals

# --------------------------------------------------------------------
def main(args):
    abund = pd.read_csv(args.abundances_csv)
    print(f"Reading allele kmer files from: {args.kmers_dir}")
    allele_kmers = read_allele_kmers(args.kmers_dir)
    print(f"Found {len(allele_kmers)} allele kmer files")
    sample_kmers = build_sample_vectors(allele_kmers, abund)
    print(f"Built weighted kmer vectors for {len(sample_kmers)} samples")

    if args.top_kmers is not None:
        totals = defaultdict(float)
        for s, vec in sample_kmers.items():
            for k, v in vec.items():
                totals[k] += v
        top = sorted(totals.items(), key=lambda x: x[1], reverse=True)[:args.top_kmers]
        whitelist = set([k for k,_ in top])
        print(f"Using top {len(whitelist)} kmers across all samples (memory reduction)")
    else:
        whitelist = None

    mat = make_dense_matrix(sample_kmers, kmer_whitelist=whitelist)
    if args.normalize:
        mat = mat.div(mat.sum(axis=1), axis=0).fillna(0.0)

    print("Computing Bray-Curtis distance...")
    dm_bray = bray_curtis_from_matrix(mat)
    dm_bray.to_csv("dist_brays.tsv", sep="\t")

    print("Computing Jensen-Shannon distance...")
    dm_jsd = jensen_shannon_from_matrix(mat)
    dm_jsd.to_csv("dist_jsd.tsv", sep="\t")

    if args.clr:
        print("Computing CLR + Euclidean...")
        clr = clr_transform(mat, pseudocount=args.pseudocount)
        dm_clr = pd.DataFrame(squareform(pdist(clr.values, metric='euclidean')), index=mat.index, columns=mat.index)
        dm_clr.to_csv("dist_clr_euclidean.tsv", sep="\t")

    print("Running PCoA (Bray)...")
    coords, eigvals = pcoa_from_distance(dm_bray, n_components=2)
    coords.to_csv("pcoa_brays_coords.csv")
    pd.DataFrame(eigvals, columns=["eigenvalue"]).to_csv("pcoa_brays_eigenvalues.csv", index=False)

    print("\nDone. Outputs:")
    print(" - dist_brays.tsv")
    print(" - dist_jsd.tsv")
    if args.clr:
        print(" - dist_clr_euclidean.tsv")
    print(" - pcoa_brays_coords.csv")
    print(" - pcoa_brays_eigenvalues.csv")

# --------------------------------------------------------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--kmers_dir", required=True, help="directory with <allele_id>.kmers.tsv files")
    parser.add_argument("--abundances_csv", required=True, help="CSV with sample,region,v_code,av,...")
    parser.add_argument("--top_kmers", type=int, default=None, help="limit to top N global kmers to reduce memory")
    parser.add_argument("--normalize", action='store_true', help="normalize sample vectors to relative abundances before distances")
    parser.add_argument("--clr", action='store_true', help="compute CLR+Euclidean distance (writes dist_clr_euclidean.tsv)")
    parser.add_argument("--pseudocount", type=float, default=1e-6)
    args = parser.parse_args()
    main(args)
