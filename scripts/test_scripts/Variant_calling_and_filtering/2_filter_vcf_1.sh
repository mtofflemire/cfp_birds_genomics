#!/bin/bash

#chamaea filter for structure
#this was filter setting for structure analysis
INPUT='/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/Cfp_birds_wgs/data/chamaea/41-Chamaea_filteredQC.vcf.gz'
OUTPUT="/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/Cfp_birds_wgs/data/chamaea.pruned.vcf.gz"

echo "Running pruning pipeline..."

bcftools +setGT "$INPUT" -- -t q -n . -i 'FMT/DP<6' \
| bcftools +fill-tags -- -t MAF,F_MISSING \
| bcftools view -m2 -M2 -v snps -i 'MAF >= 0.05 && F_MISSING <= 0.10' \
| bcftools +prune -m 0.25 -w 2500bp \
| bcftools view -Oz -o "$OUTPUT"


echo "Indexing VCF..."
tabix -f -p vcf "$OUTPUT"

echo "DONE — Output written to:"
echo "$OUTPUT"





#picoides filter for structure
#this was filter setting for structure analysis
INPUT='/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/Cfp_birds_wgs/data/picoides/picoides.filteredQC.vcf.gz'
OUTPUT="/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/Cfp_birds_wgs/data/picoides/picoides.pruned.vcf.gz"

echo "Running pruning pipeline..."

bcftools +setGT "$INPUT" -- -t q -n . -i 'FMT/DP<3' \
| bcftools +fill-tags -- -t MAF,F_MISSING \
| bcftools view -m2 -M2 -v snps -i 'MAF >= 0.05 && F_MISSING <= 0.10' \
| bcftools +prune -m 0.25 -w 2500bp \
| bcftools view -Oz -o "$OUTPUT"


echo "Indexing VCF..."
tabix -f -p vcf "$OUTPUT"

echo "DONE — Output written to:"
echo "$OUTPUT"





#Mountain Chikcadee filter for structure
#this was filter setting for structure analysis
INPUT='/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/Cfp_birds_wgs/data/poecile/poecile.no_outgroup.vcf.gz'
OUTPUT="/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/Cfp_birds_wgs/data/poecile/poecile.pruned.vcf.gz"

echo "Running pruning pipeline..."

bcftools +setGT "$INPUT" -- -t q -n . -i 'FMT/DP<3' \
| bcftools +fill-tags -- -t MAF,F_MISSING \
| bcftools view -m2 -M2 -v snps -i 'MAF >= 0.05 && F_MISSING <= 0.10' \
| bcftools +prune -m 0.25 -w 2500bp \
| bcftools view -Oz -o "$OUTPUT"


echo "Indexing VCF..."
tabix -f -p vcf "$OUTPUT"

echo "DONE — Output written to:"
echo "$OUTPUT"






#Mountain Chikcadee filter for structure
#this was filter setting for structure analysis
INPUT='/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/Cfp_birds_wgs/data/toxostoma/Toxostoma_filteredQC.vcf.gz'
OUTPUT="/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/Cfp_birds_wgs/data/toxostoma/Toxostoma.pruned.vcf.gz"

echo "Running pruning pipeline..."

bcftools +setGT "$INPUT" -- -t q -n . -i 'FMT/DP<3' \
| bcftools +fill-tags -- -t MAF,F_MISSING \
| bcftools view -m2 -M2 -v snps -i 'MAF >= 0.05 && F_MISSING <= 0.10' \
| bcftools +prune -m 0.25 -w 2500bp \
| bcftools view -Oz -o "$OUTPUT"


echo "Indexing VCF..."
tabix -f -p vcf "$OUTPUT"

echo "DONE — Output written to:"
echo "$OUTPUT"







