#!/bin/bash
set -euo pipefail

VCF='/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/Cfp_birds_wgs/data/toxostoma/Toxostoma.pruned.vcf.gz'
OUTDIR='/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/Cfp_birds_wgs/out/1_pca'
PREFIX="toxostoma"

mkdir -p "$OUTDIR"

plink \
  --vcf "$VCF" \
  --double-id \
  --allow-extra-chr \
  --set-missing-var-ids @:# \
  --make-bed \
  --out "$OUTDIR/$PREFIX"


plink \
  --bfile "$OUTDIR/$PREFIX" \
  --allow-extra-chr \
  --pca 20 header \
  --out "$OUTDIR/${PREFIX}_pca"


wc -l "$OUTDIR/$PREFIX.bim"

