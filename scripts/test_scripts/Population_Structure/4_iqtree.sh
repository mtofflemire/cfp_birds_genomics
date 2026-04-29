
python /Users/michaeltofflemire/softs/vcf2phylip-master/vcf2phylip.py -i '/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/poecile_genotyping/data2/poecile.snps.pruned.vcf.gz' -n






/Users/michaeltofflemire/softs/iqtree-2.4.0-macOS/bin/iqtree2 \
  -s '/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/cfp_comparative_genomics/data/Poecile/iqtree/poecile.pruned.recoded.min4.phy.varsites.phy' \
  -st DNA \
  -nt 8 \
  -m MFP+ASC \
  -bb 1000 
