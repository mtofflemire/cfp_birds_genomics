cd '/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/California_birds_wgs/data/chamaea'


bcftools view \
  -R '/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/California_birds_wgs/out/rda_outliers/chamaea_RDA_outliers.pos.txt' \
  '/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/California_birds_wgs/data/chamaea/chamaea.pruned.vcf.gz' \
  -Oz -o "/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/California_birds_wgs/data/chamaea/chamaea.pruned.selective.snps.vcf.gz"



bcftools index \
"/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/California_birds_wgs/data/chamaea/chamaea.pruned.selective.snps.vcf.gz"



