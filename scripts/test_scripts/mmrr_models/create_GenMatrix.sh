#need to make genomic dissimilarity matrix with the 1-ibs command in PLINK1.9 
plink \
  --vcf '/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/California_birds_wgs/data/chamaea/chamaea.pruned.selective.snps.vcf.gz' \
  --double-id \
  --allow-extra-chr \
  --distance square 1-ibs \
  --out '/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/California_birds_wgs/data/chamaea/chamaea.gendist.selective.snp'



