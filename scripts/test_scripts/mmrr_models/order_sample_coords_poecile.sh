META='/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/Cfp_birds_wgs/data/poecile/poecile_metadata.csv'
SAMPLES='/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/Cfp_birds_wgs/data/poecile/poecile.samples.keep.txt'
OUT='/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/Cfp_birds_wgs/data/poecile/poecile_coords_ordered.txt'

awk -F',' -v samples="$SAMPLES" '
BEGIN {
  while ((getline line < samples) > 0) {
    gsub(/\r/, "", line)
    gsub(/^ +| +$/, "", line)
    sub(/_S[0-9]+.*/, "", line)
    order[line] = ++n
  }
}

NR==1 {
  for (i=1; i<=NF; i++) {
    gsub(/"/, "", $i)

    if ($i=="sampleID") id=i
    if ($i=="decimallatitude") lat=i
    if ($i=="decimallongitude") lon=i
  }
  next
}

{
  gsub(/"/, "", $0)
  key = $id

  if (key in order) {

    latitude  = $lat
    longitude = $lon

    # Check numeric and valid ranges
    if (latitude ~ /^-?[0-9.]+$/ &&
        longitude ~ /^-?[0-9.]+$/ &&
        latitude  >= -90  && latitude  <= 90 &&
        longitude >= -180 && longitude <= 180) {

      coords[order[key]] = latitude "\t" longitude
    }
  }
}

END {
  for (i=1; i<=n; i++) {
    if (coords[i] == "")
      print "WARNING: Missing or invalid coordinates for sample #" i > "/dev/stderr"
    else
      print coords[i]
  }
}
' "$META" > "$OUT"
