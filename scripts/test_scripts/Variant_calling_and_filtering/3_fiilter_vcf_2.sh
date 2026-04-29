
#Wrentit diversity filtering
INPUT='/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/Cfp_birds_wgs/data/chamaea/41-Chamaea_filteredQC.vcf.gz'
OUTPUT="/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/Cfp_birds_wgs/data/chamaea/41-Chamaea.diversity.vcf.gz"

bcftools +setGT "$INPUT" -- \
    -t q -n . -i 'FMT/DP<6' \
| bcftools +fill-tags -- -t F_MISSING \
| bcftools view \
    -m2 -M2 -v snps \
    -i 'F_MISSING <= 0.05' \
-Oz -o "$OUTPUT"

tabix -p vcf "$OUTPUT"


#California Thrasher diversity stats filtering
INPUT='/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/Cfp_birds_wgs/data/toxostoma/Toxostoma_filteredQC.vcf.gz'
OUTPUT="/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/Cfp_birds_wgs/data/toxostoma/Toxostoma.diversity.vcf.gz"

bcftools +setGT "$INPUT" -- \
    -t q -n . -i 'FMT/DP<3' \
| bcftools +fill-tags -- -t F_MISSING \
| bcftools view \
    -m2 -M2 -v snps \
    -i 'F_MISSING <= 0.05' \
-Oz -o "$OUTPUT"

tabix -p vcf "$OUTPUT"


#White-headed Woodpecker diversity stats filtering
INPUT='/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/Cfp_birds_wgs/data/picoides/picoides.filteredQC.vcf.gz'
OUTPUT="/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/Cfp_birds_wgs/data/picoides/picoides.diversity.vcf.gz"

bcftools +setGT "$INPUT" -- \
    -t q -n . -i 'FMT/DP<3' \
| bcftools +fill-tags -- -t F_MISSING \
| bcftools view \
    -m2 -M2 -v snps \
    -i 'F_MISSING <= 0.05' \
-Oz -o "$OUTPUT"

tabix -p vcf "$OUTPUT"



#Mountain Chickadee diversity stats filtering
INPUT='/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/Cfp_birds_wgs/data/poecile/poecile.no_outgroup.vcf.gz'
OUTPUT="/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/Cfp_birds_wgs/data/poecile/poecile.diversity.vcf.gz"

bcftools +setGT "$INPUT" -- \
    -t q -n . -i 'FMT/DP<3' \
| bcftools +fill-tags -- -t F_MISSING \
| bcftools view \
    -m2 -M2 -v snps \
    -i 'F_MISSING <= 0.05' \
-Oz -o "$OUTPUT"

tabix -p vcf "$OUTPUT"














