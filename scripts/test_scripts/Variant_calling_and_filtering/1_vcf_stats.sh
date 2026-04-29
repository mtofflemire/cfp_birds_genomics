#summarize sequencing depth per individual for each species
vcftools \
  --gzvcf '/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/Cfp_birds_wgs/data/chamaea/41-Chamaea_filteredQC.vcf.gz' \
  --depth \
  --out '/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/Cfp_birds_wgs/data/chamaea/chamaea'

  vcftools \
  --gzvcf '/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/Cfp_birds_wgs/data/toxostoma/Toxostoma_filteredQC.vcf.gz' \
  --depth \
  --out '/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/Cfp_birds_wgs/data/toxostoma/toxostoma'

  vcftools \
  --gzvcf '/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/Cfp_birds_wgs/data/poecile/poecile.filteredQC.vcf.gz' \
  --depth \
  --out '/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/Cfp_birds_wgs/data/poecile/poecile'

  vcftools \
  --gzvcf '/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/Cfp_birds_wgs/data/picoides/picoides.filteredQC.vcf.gz' \
  --depth \
  --out '/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/Cfp_birds_wgs/data/picoides/picoides'


#summarize results
#assess mean, min and max coverage values per species
awk 'NR>1 {sum+=$3; if(min=="" || $3<min) min=$3; if($3>max) max=$3}
END {print "Mean:",sum/(NR-1), "Min:",min, "Max:",max}' \
'/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/Cfp_birds_wgs/data/picoides/picoides.idepth'




awk 'NR>1 {x[NR]=$3; sum+=$3}
END {
mean=sum/(NR-1);
for(i in x){ss+=(x[i]-mean)^2}
sd=sqrt(ss/(NR-2));
print "Mean:",mean,"SD:",sd,"Min:",min,"Max:",max
}' \
'/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/Cfp_birds_wgs/data/poecile/poecile.idepth'

