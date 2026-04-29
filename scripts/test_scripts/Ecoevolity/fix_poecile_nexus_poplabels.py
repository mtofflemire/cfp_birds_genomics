#!/usr/bin/env python3
import sys
import csv
import re

if len(sys.argv) != 4:
    sys.exit(
        "Usage: fix_poecile_nexus_poplabels.py <input.nexus> <metadata.csv> <output.nexus>"
    )

nexus_path = sys.argv[1]
meta_path  = sys.argv[2]
out_path   = sys.argv[3]

# -----------------------------
# Load metadata (CSV + BOM safe)
# -----------------------------
lookup = {}

with open(meta_path, newline="", encoding="utf-8-sig") as csvfile:
    reader = csv.DictReader(csvfile, delimiter=",")

    required = {"sampleID", "pop"}
    if not required.issubset(reader.fieldnames):
        sys.exit(
            f"Metadata must contain columns: {required}\n"
            f"Found: {reader.fieldnames}"
        )

    for row in reader:
        sample = row["sampleID"].strip()
        pop = row["pop"].strip()

        # Extract museum catalog ID
        m = re.search(r"(UWBM\d+|MVZ\d+|SDNHM\d+|UMMZ\d+)", sample)
        if not m:
            sys.exit(f"Could not extract catalog ID from sampleID: {sample}")

        catalog = m.group(1)
        lookup[catalog] = f"{pop}_{catalog}"

# -----------------------------
# Rewrite NEXUS
# -----------------------------
with open(nexus_path) as fin, open(out_path, "w") as fout:
    in_matrix = False

    for line in fin:
        stripped = line.strip()

        if stripped.upper().startswith("MATRIX"):
            in_matrix = True
            fout.write(line)
            continue

        if in_matrix and stripped.startswith(";"):
            in_matrix = False
            fout.write(line)
            continue

        if in_matrix and stripped:
            parts = stripped.split(None, 1)
            if len(parts) != 2:
                sys.exit(f"Cannot parse MATRIX row:\n{line}")

            old_name, seq = parts

            m = re.search(r"(UWBM\d+|MVZ\d+|SDNHM\d+|UMMZ\d+)", old_name)
            if not m:
                sys.exit(f"Could not find catalog ID in taxon name: {old_name}")

            catalog = m.group(1)
            if catalog not in lookup:
                sys.exit(f"{catalog} not found in metadata")

            fout.write(f"{lookup[catalog]} {seq}\n")
        else:
            fout.write(line)

print(f"✔ Rewrote taxon labels → {out_path}")
