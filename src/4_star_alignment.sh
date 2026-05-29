#!/usr/bin/env bash

set -euo pipefail
shopt -s nullglob

# Directorios
CLEAN="./data/fastp"
STAR_PE="./data/star_pe"
STAR_IDX="./data/references/star_index/"
METADATA="./data/metadata/Runs_metadata.csv"

mkdir -p "$STAR_PE"
[[ ! -d "$STAR_IDX" ]] && echo "[ERROR] Directorio de índice STAR_PE no existente en $STAR_IDX" && exit 1

# Hilos para STAR
threads_star=12

# Archivos R1 limpios
files=("$CLEAN"/*_clean_1.fastq.gz)

# Mapa: SRR -> nombre nuevo
declare -A sample_names

while IFS="," read -r Run age_of_onset tissue BioSample SRR_file_name; do
    [[ "$Run" == "Run" ]] && continue
    sample_names["$Run"]="$SRR_file_name"
done < "$METADATA"

for f in "${files[@]}"; do
    base=$(basename "$f" _clean_1.fastq.gz)
    r2="${CLEAN}/${base}_clean_2.fastq.gz"

    out_name="${sample_names[$base]:-$base}"

    echo "============================================================"
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Iniciando alineamiento STAR PE"
    echo "Muestra SRR: $base"
    echo "Nombre de salida: $out_name"
    echo "R1: $f"
    echo "R2: $r2"
    echo "Salida: ${STAR_PE}/${out_name}"

    if [[ ! -f "$r2" ]]; then
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: falta archivo R2 para $base: $r2"
        continue
    fi

    STAR \
        --runThreadN "$threads_star" \
        --genomeDir "$STAR_IDX" \
        --readFilesIn "$f" "$r2" \
        --outFileNamePrefix "${STAR_PE}/${out_name}." \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMunmapped None \
	--readFilesCommand zcat

    status=$?

    if [[ $status -eq 0 ]]; then
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Alineamiento completado exitosamente: $base → $out_name"
    else
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: STAR falló para $base con código $status"
        exit "$status"
    fi

    echo "============================================================"
done
