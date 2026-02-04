#!/usr/bin/env python3
"""
split_abundances_by_region.py

Splits a CSV file named abundances.csv by the 'region' column.
Each region gets its own subfolder (e.g., V1/, V2/, V3/) containing:
    abundances.csv   -> all rows with that region
"""

import csv
import os
import sys
from collections import defaultdict

def main(infile):
    # Read all rows grouped by region
    with open(infile, newline='') as f:
        reader = csv.DictReader(f)
        if "region" not in reader.fieldnames:
            sys.exit("Error: input CSV must have a 'region' column.")
        
        regions = defaultdict(list)
        for row in reader:
            region = row["region"].strip()
            regions[region].append(row)

        print(f"Found {len(regions)} regions: {', '.join(sorted(regions.keys()))}")

    # Write each region’s CSV into its folder
    for region, rows in regions.items():
        os.makedirs(region, exist_ok=True)
        outpath = os.path.join(region, "abundances.csv")
        with open(outpath, "w", newline='') as out:
            writer = csv.DictWriter(out, fieldnames=reader.fieldnames)
            writer.writeheader()
            writer.writerows(rows)
        print(f"Wrote {len(rows)} rows to {outpath}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit("Usage: ./split_abundances_by_region.py abundances.csv")
    main(sys.argv[1])
