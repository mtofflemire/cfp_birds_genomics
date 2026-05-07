#!/bin/bash
#PBS -V
#PBS -N markdup.poecile
#PBS -l nodes=1:ppn=20,mem=128gb,walltime=500:00:00
#PBS -q batch
#PBS -j oe
#PBS -o /usr/scratch2/userdata2/mtofflemire/projects/cfp_birds/04_MAPPING/06_markdup.log

set -euo pipefail
shopt -s nullglob
source ~/miniconda3/etc/profile.d/conda.sh
conda activate genomics
cd "${PBS_O_WORKDIR:-$PWD}"

BAM_DIR="/usr/scratch2/userdata2/mtofflemire/projects/cfp_birds/04_MAPPING"

for BAM in "$BAM_DIR"/*.sorted.bam; do
  SAMPLE_NAME="$(basename "$BAM" .sorted.bam)"
  echo "Marking duplicates for $SAMPLE_NAME"
  gatk --java-options "-Xms4g -Xmx16g" MarkDuplicates \
    -I "$BAM" \
    -O "$BAM_DIR/${SAMPLE_NAME}.dedup.bam" \
    -M "$BAM_DIR/${SAMPLE_NAME}.dup_metrics.txt" \
    --CREATE_INDEX true \
    --VALIDATION_STRINGENCY SILENT
  samtools flagstat "$BAM_DIR/${SAMPLE_NAME}.dedup.bam" > "$BAM_DIR/${SAMPLE_NAME}_dedup.flagstat.txt"
  echo "Finished $SAMPLE_NAME at $(date)"
done
