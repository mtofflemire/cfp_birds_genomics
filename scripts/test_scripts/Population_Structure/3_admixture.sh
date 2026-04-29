#!/bin/bash

VCF='/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/California_birds_wgs/data/picoides/picoides.pruned.vcf.gz'
RESULTS='/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/California_birds_wgs/out/2_structure/admixture'
BED_PREFIX="$RESULTS/picoides"
ADMIX=/Users/michaeltofflemire/softs/dist/admixture_macosx-1.3.0/admixture

# Convert pruned VCF directly to PLINK BED
plink --vcf "$VCF" \
      --make-bed \
      --out "$BED_PREFIX" \
      --set-missing-var-ids @:# \
      --allow-extra-chr \
      --double-id

# Fix chromosome column for ADMIXTURE
awk '{ $1="0"; print }' "${BED_PREFIX}.bim" > "${BED_PREFIX}.bim.tmp"
mv "${BED_PREFIX}.bim.tmp" "${BED_PREFIX}.bim"


for K in 1 2 3 4 5; do
    "$ADMIX" -s 43 --cv=10 -j9 "${BED_PREFIX}.bed" "$K" > "$RESULTS/log${K}.out"
done









