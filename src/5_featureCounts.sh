#!/usr/bin/env bash

set -euo pipefail
shopt -s nullglob

# ====================== DIRECTORIOS ======================

DATA="./data"
BAMs="${DATA}/star_pe/BAM"
FILTERED_BAMs="${DATA}/star_pe/BAM_filtered"
FC_DIR="${DATA}/featureCounts"
TMP_PARENT="${DATA}/tmp_samtools"

mkdir -p "$FC_DIR" "$FILTERED_BAMs" "$TMP_PARENT"

# ====================== CONFIGURACIÓN ======================

wks=9
MAPQ=10
ANNOTATION="./data/references/gencode.v49.annotation.gtf"

# ====================== RENOMBRADO ======================

echo "================================================================================================="
echo "Renombrando BAMs de STAR"
echo "================================================================================================="

bam_files=("$BAMs"/*.bam)

for bam in "${bam_files[@]}"; do

    base=$(basename "$bam")
    new_bam="${base%.Aligned.sortedByCoord.out.bam}.bam"

    echo "$bam -> ${BAMs}/$new_bam"

    if [[ "$base" != "$new_bam" ]]; then
        mv -f "$bam" "${BAMs}/$new_bam"
    fi
done

# ====================== FILTRADO, ORDENADO E INDEXADO ======================

echo "================================================================================================="
current_date_time="$(date "+%Y-%m-%d %H:%M:%S")"
echo "Iniciando filtrado, ordenado e indexado con samtools - $current_date_time"
echo "================================================================================================="

bam_files=("$BAMs"/*.bam)

if [[ ${#bam_files[@]} -eq 0 ]]; then
    echo "No hay BAMs disponibles en $BAMs"
    exit 1
fi

for RAW_BAM in "${bam_files[@]}"; do

    base=$(basename "$RAW_BAM" .bam)
    SORTED_BAM="${FILTERED_BAMs}/${base}.filtsort.bam"
    tmp_dir="${TMP_PARENT}/${base}"

    mkdir -p "$tmp_dir"

    echo "-------------------------------------------------------------------------------------------------"
    echo "Procesando BAM: $RAW_BAM"
    echo "Salida filtrada/ordenada: $SORTED_BAM"
    echo "Filtro: mapped reads, MAPQ >= $MAPQ"

    if samtools view -F 4 -q "$MAPQ" -u "$RAW_BAM" | \
       samtools sort -@ "$wks" -T "${tmp_dir}/sort" -o "$SORTED_BAM"; then

        echo "Procesado con éxito: $RAW_BAM"

    else
        echo "[ERROR] Falló al filtrar/ordenar con samtools en $RAW_BAM" >&2
        rm -rf "$tmp_dir"
        exit 1
    fi

    echo "Indexando BAM: $SORTED_BAM"

    if samtools index -@ "$wks" "$SORTED_BAM"; then
        echo "Índice generado: ${SORTED_BAM}.bai"
    else
        echo "[ERROR] Falló el indexado de $SORTED_BAM" >&2
        rm -rf "$tmp_dir"
        exit 1
    fi

    rm -rf "$tmp_dir"

done

# ====================== FEATURECOUNTS ======================

echo "================================================================================================="
current_date_time="$(date "+%Y-%m-%d %H:%M:%S")"
echo "Iniciando FeatureCounts - $current_date_time"
echo "================================================================================================="

bam_files=("$FILTERED_BAMs"/*.filtered.sorted.bam)
mapfile -t sorted_bams < <(printf "%s\n" "${bam_files[@]}" | sort)

if [[ ${#sorted_bams[@]} -eq 0 ]]; then
    echo "No hay BAMs filtrados disponibles, se omite featureCounts"
    exit 1
fi

featureCounts \
    -T "$wks" \
    -a "$ANNOTATION" \
    -o "${FC_DIR}/counts_table.txt" \
    -s 0 \
    -p -B -C \
    "${sorted_bams[@]}"

current_date_time="$(date "+%Y-%m-%d %H:%M:%S")"
echo "FeatureCounts finalizado, tabla guardada en $FC_DIR - $current_date_time"
echo "================================================================================================="

current_date_time="$(date "+%Y-%m-%d %H:%M:%S")"
echo "Script finalizado - $current_date_time"
echo "================================================================================================="