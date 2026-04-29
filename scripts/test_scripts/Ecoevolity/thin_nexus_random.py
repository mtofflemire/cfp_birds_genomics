#!/usr/bin/env python3
import random
import sys
import re

def replace_nchar_in_dimensions(line: str, new_nchar: int) -> str:
    # Handles NCHAR=123 or NCHAR = 123 (any spacing/case)
    return re.sub(r'(?i)\bNCHAR\s*=\s*\d+\b', f'NCHAR={new_nchar}', line)

def main():
    if len(sys.argv) < 4:
        print("Usage: thin_nexus_random.py <input.nexus> <output.nexus> <nsites> [seed]", file=sys.stderr)
        sys.exit(1)

    infile  = sys.argv[1]
    outfile = sys.argv[2]
    nsites  = int(sys.argv[3])
    seed    = int(sys.argv[4]) if len(sys.argv) >= 5 else None
    if seed is not None:
        random.seed(seed)

    pre = []          # lines before MATRIX
    matrix_line = None
    taxa = []
    seqs = []
    matrix_end = None # the line containing ';' that ends MATRIX
    post = []         # lines after MATRIX end

    state = "pre"

    with open(infile, "r") as f:
        for raw in f:
            line = raw.rstrip("\n")

            if state == "pre":
                if line.strip().upper().startswith("MATRIX"):
                    matrix_line = raw  # keep original newline
                    state = "matrix"
                else:
                    pre.append(raw)

            elif state == "matrix":
                # end of MATRIX block
                if line.strip().startswith(";"):
                    matrix_end = raw
                    state = "post"
                    continue

                # skip blank lines inside matrix
                if not line.strip():
                    continue

                # parse "taxon sequence"
                parts = line.strip().split(None, 1)
                if len(parts) != 2:
                    raise ValueError(f"Could not parse MATRIX row: {line}")
                name, seq = parts
                taxa.append(name)
                seqs.append(seq)

            else:  # post
                post.append(raw)

    if matrix_line is None or matrix_end is None:
        raise ValueError("Could not find a MATRIX block terminated by ';'")

    if not seqs:
        raise ValueError("No sequences found in MATRIX block")

    nchar = len(seqs[0])
    if any(len(s) != nchar for s in seqs):
        raise ValueError("Not all sequences have the same length (NEXUS matrix is ragged)")

    if nsites > nchar:
        raise ValueError(f"Requested {nsites} sites but alignment has only {nchar}")

    cols = sorted(random.sample(range(nchar), nsites))
    new_seqs = [''.join(seq[i] for i in cols) for seq in seqs]

    # Update NCHAR in any DIMENSIONS line(s) in the preamble (keep END; placement intact)
    new_pre = []
    for raw in pre:
        if "NCHAR" in raw.upper() and "DIMENSIONS" in raw.upper():
            new_pre.append(replace_nchar_in_dimensions(raw, nsites))
        else:
            new_pre.append(raw)

    with open(outfile, "w") as out:
        out.writelines(new_pre)
        out.write(matrix_line)
        for name, seq in zip(taxa, new_seqs):
            out.write(f"{name} {seq}\n")
        out.write(matrix_end)
        out.writelines(post)

if __name__ == "__main__":
    main()
