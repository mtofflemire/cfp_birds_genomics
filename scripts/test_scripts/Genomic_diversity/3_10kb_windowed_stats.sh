#!/bin/bash

OUT_DIR='/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/Cfp_birds_wgs/out/pi'

declare -A VCFs=(
  [poecile]='/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/Cfp_birds_wgs/data/poecile/poecile.diversity.autosomes.vcf.gz'
  [picoides]='/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/Cfp_birds_wgs/data/picoides/picoides.diversity.autosomes.vcf.gz'
  [chamaea]='/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/Cfp_birds_wgs/data/chamaea/41-Chamaea_filteredQC.vcf.gz'
  [toxostoma]='/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/Cfp_birds_wgs/data/toxostoma/Toxostoma_filteredQC_2.vcf.gz'
)

for species in "${!VCFs[@]}"; do
  VCF_FILE="${VCFs[$species]}"
  echo "Running site-pi for $species"

  vcftools --gzvcf "$VCF_FILE" \
           --site-pi \
           --out "$OUT_DIR/${species}.sites"
done





#!/bin/bash

OUT_DIR='/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/Cfp_birds_wgs/out/tajima'

declare -A VCFs=(
  [poecile]='/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/Cfp_birds_wgs/data/poecile/poecile.diversity.autosomes.vcf.gz'
  [picoides]='/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/Cfp_birds_wgs/data/picoides/picoides.diversity.autosomes.vcf.gz'
  [chamaea]='/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/Cfp_birds_wgs/data/chamaea/41-Chamaea_filteredQC.vcf.gz'
  [toxostoma]='/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/Cfp_birds_wgs/data/toxostoma/Toxostoma_filteredQC_2.vcf.gz'
)

for species in "${!VCFs[@]}"; do
  VCF_FILE="${VCFs[$species]}"
  echo "Running Tajima's D for $species"

  vcftools --gzvcf "$VCF_FILE" \
           --TajimaD 10000 \
           --out "$OUT_DIR/${species}_tajima_10kb"
done
