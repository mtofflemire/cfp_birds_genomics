META='/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/California_birds_wgs/data/toxostoma/toxostoma-meta.csv'
VCF='/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/California_birds_wgs/data/toxostoma/Toxostoma_filteredQC_2.vcf.gz'
OUT='/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/California_birds_wgs/data/toxostoma/coords.txt'

VCF_SAMPLES='/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/California_birds_wgs/data/toxostoma/vcf_samples.txt'
META_SUB='/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/California_birds_wgs/data/toxostoma/meta_subset.csv'

# 1️⃣ Extract VCF sample order and remove lane suffix (_S##_L###)
bcftools query -l "$VCF" \
  | sed -E 's/_S[0-9]+_L[0-9]+$//' \
  > "$VCF_SAMPLES"

# 2️⃣ Safely extract needed metadata columns (handles commas correctly)
csvcut -c sampleID,decimallatitude,decimallongitude "$META" > "$META_SUB"

# 3️⃣ Reorder metadata to match VCF order
awk -F',' '
NR==FNR {
  order[++n]=$1
  next
}
NR==1 { next }
{
  coords[$1]=$2 "\t" $3
}
END {
  for (i=1; i<=n; i++)
    print coords[order[i]]
}
' "$VCF_SAMPLES" "$META_SUB" > "$OUT"

