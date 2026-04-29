META='/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/Cfp_birds_wgs/data/chamaea/Cfa-metadata-nsamp159copy.csv'
SAMPLES='/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/Cfp_birds_wgs/data/chamaea/samples.txt'
OUT='/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/Cfp_birds_wgs/data/chamaea/Cfa_coords_ordered.txt'

awk -F',' -v samples="$SAMPLES" '
BEGIN {
  while ((getline line < samples) > 0) {
    order[line] = ++n
  }
}

NR==1 {
  for (i=1; i<=NF; i++) {
    if ($i=="sampleID") id=i
    if ($i=="lat") lat=i
    if ($i=="long") lon=i
  }
  next
}

{
  key = $id
  if (key in order) {
    coords[order[key]] = $lat "\t" $lon
  }
}

END {
  for (i=1; i<=n; i++)
    print coords[i]
}
' "$META" > "$OUT"
