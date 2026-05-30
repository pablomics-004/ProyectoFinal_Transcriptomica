#!/usr/bin/env bash

set -euo pipefail
shopt -s nullglob

# ====================== DIRECTORIOS ======================

DATA="./data"
BAMs="${DATA}/star_pe/BAM"
FC_DIR="${DATA}/featureCounts"

mkdir -p "$FC_DIR"

# ====================== RENOMBRADO ======================

bam_files=("$BAMs"/*.bam)

for bam in "${bam_files[@]}"; do
    
    # Renombrado de archivos BAM
    base=$(basename "$bam")
    new_bam="${base%.Aligned.sortedByCoord.out.bam}.bam"

    echo "$bam -> ${BAMs}/$new_bam"

    if [[ "$base" != "$new_bam" ]]; then
        mv -f "$bam" "${BAMs}/$new_bam"
    fi
done

# ====================== FEATURECOUNTS ======================

wks=9
ANNOTATION="./data/references/gencode.v49.annotation.gtf"

echo "================================================================================================="
current_date_time="$(date "+%Y-%m-%d %H:%M:%S")"
echo "Iniciando FeatureCounts - $current_date_time"
echo ""

bam_files=("$BAMs"/*.bam)
mapfile -t sorted_bams < <(printf "%s\n" "${bam_files[@]}" | sort)

if [[ ${#sorted_bams[@]} -eq 0 ]]; then
    echo "No hay BAMs disponibles, se omite featureCounts"
    exit 1
fi

featureCounts \
    -T $wks \
    -a "$ANNOTATION" \
    -o "${FC_DIR}/counts_table.txt" \
    --largestOverlap \
    -s 0 \
    -p -B -C \
    "${sorted_bams[@]}"

current_date_time="$(date "+%Y-%m-%d %H:%M:%S")"
echo "FeatureCounts finalizado, tabla guardada en $FC_DIR - $current_date_time"
echo "================================================================================================="

current_date_time="$(date "+%Y-%m-%d %H:%M:%S")"
echo "Script finalizado - $current_date_time"
echo "================================================================================================="