#calculate heterozygosity with plinkv1.9
#chamaea heterozysosity
plink --vcf '/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/Cfp_birds_wgs/data/chamaea/chamaea.diversity.autosomes.vcf.gz' --double-id --allow-extra-chr --make-bed --out '/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/Cfp_birds_wgs/data/chamaea/41-Chamaea.diversity'
plink --bfile '/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/Cfp_birds_wgs/data/chamaea/41-Chamaea.diversity' --allow-extra-chr --het --out '/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/Cfp_birds_wgs/out/het/41-Chamaea.diversity'


#Toxostoma heterozygosity
plink --vcf '/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/Cfp_birds_wgs/data/toxostoma/Toxostoma_filteredQC_2.vcf.gz' --double-id --allow-extra-chr --make-bed --out '/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/Cfp_birds_wgs/out/het/Toxostoma_filteredQC_2'
plink --bfile '/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/Cfp_birds_wgs/out/het/Toxostoma_filteredQC_2' --allow-extra-chr --het --out '/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/Cfp_birds_wgs/out/het/Toxostoma_filteredQC_2'



#Poecile heterozygosity
plink --vcf '/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/Cfp_birds_wgs/data/poecile/poecile.diversity.autosomes.vcf.gz' --double-id --allow-extra-chr --make-bed --out '/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/Cfp_birds_wgs/out/het/poecile.diversity.autosomes'
plink --bfile '/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/Cfp_birds_wgs/out/het/poecile.diversity.autosomes' --allow-extra-chr --het --out '/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/Cfp_birds_wgs/out/het/poecile.diversity.autosomes'



#Picoides heterozygosity
plink --vcf '/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/Cfp_birds_wgs/data/picoides/picoides.diversity.autosomes.vcf.gz' --double-id --allow-extra-chr --make-bed --out '/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/Cfp_birds_wgs/out/het/picoides.diversity'
plink --bfile '/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/Cfp_birds_wgs/out/het/picoides.diversity' --allow-extra-chr --het --out '/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/Cfp_birds_wgs/out/het/picoides.diversity'



