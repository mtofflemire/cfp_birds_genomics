#!/usr/bin/env python3
import random
import sys
import re
from collections import defaultdict

def replace_ntax(line: str, new_ntax: int) -> str:
    return re.sub(r'(?i)\bNTAX\s*=\s*\d+\b', f'NTAX={new_ntax}', line)

def main():
    if len(sys.argv) < 4:
        print("Usage: thin_nexus_by_population.py <input.nexus> <output.nexus> <n_per_pop> [seed]", file=sys.stderr)
        sys.exit(1)

    infile = sys.argv[1]
    outfile = sys.argv[2]
    n_per_pop = int(sys.argv[3])
    seed = int(sys.argv[4]) if len(sys.argv) >= 5 else None
    if seed is not None:
        random.seed(seed)

    pre = []
    matrix_line = None
    taxa = []
    seqs = []
    matrix_end = None
    post = []

    state = "pre"

    with open(infile) as f:
        for raw in f:
            line = raw.rstrip("\n")

            if state == "pre":
                if line.strip().upper().startswith("MATRIX"):
                    matrix_line = raw
                    state = "matrix"
                else:
                    pre.append(raw)

            elif state == "matrix":
                if line.strip().startswith(";"):
                    matrix_end = raw
                    state = "post"
                    continue
                if not line.strip():
                    continue

                name, seq = line.strip().split(None, 1)
                taxa.append(name)
                seqs.append(seq)

            else:
                post.append(raw)

    if matrix_line is None or matrix_end is None:
        raise ValueError("Could not find MATRIX block")

    # group taxa by population
    pop_indices = defaultdict(list)
    for i, taxon in enumerate(taxa):
        pop = taxon.split("_", 1)[0]
        pop_indices[pop].append(i)

    # sample taxa
    keep_indices = []
    for pop, inds in pop_indices.items():
        if len(inds) < n_per_pop:
            raise ValueError(f"Population {pop} has only {len(inds)} samples (< {n_per_pop})")
        keep_indices.extend(random.sample(inds, n_per_pop))

    keep_indices = sorted(keep_indices)

    new_taxa = [taxa[i] for i in keep_indices]
    new_seqs = [seqs[i] for i in keep_indices]

    # update NTAX
    new_pre = []
    for raw in pre:
        if "NTAX" in raw.upper() and "DIMENSIONS" in raw.upper():
            new_pre.append(replace_ntax(raw, len(new_taxa)))
        else:
            new_pre.append(raw)

    with open(outfile, "w") as out:
        out.writelines(new_pre)
        out.write(matrix_line)
        for name, seq in zip(new_taxa, new_seqs):
            out.write(f"{name} {seq}\n")
        out.write(matrix_end)
        out.writelines(post)

if __name__ == "__main__":
    main()
