#!/usr/bin/env bash

set -euo pipefail
shopt -s nullglob

# ====================== DIRECTORIOS ======================

CLEAN="./data/fastp"
STAR_PE="./data/star_pe"
BAM="${STAR_PE}/BAM"
STAR_IDX="./data/references/star_index/"
METADATA="./data/metadata/Runs_metadata.csv"

mkdir -p "$STAR_PE" "$BAM"

[[ ! -d "$STAR_IDX" ]] && echo "[ERROR] Directorio de índice STAR no existente en $STAR_IDX" && exit 1
[[ ! -f "$METADATA" ]] && echo "[ERROR] No existe metadata: $METADATA" && exit 1

# ====================== FUNCIONES ======================

clean_up() {
    echo "============================================================"
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Liberando el genoma de la memoria"

    STAR \
        --runMode alignReads \
        --genomeDir "$STAR_IDX" \
        --genomeLoad Remove || true
}

trap clean_up EXIT

# ====================== MAIN ======================

threads_star=12

files=("$CLEAN"/*_clean_1.fastq.gz)

if [[ ${#files[@]} -eq 0 ]]; then
    echo "[ERROR] No se encontraron archivos *_clean_1.fastq.gz en $CLEAN"
    exit 1
fi

# Mapa: SRR -> nombre nuevo
declare -A sample_names

while IFS="," read -r Run source_name BioSample SRR_file_name; do
    [[ "$Run" == "Run" ]] && continue

    SRR_file_name="${SRR_file_name//$'\r'/}"
    sample_names["$Run"]="$SRR_file_name"
done < "$METADATA"

# ====================== STAR ======================

echo "============================================================"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Pre-cargando genoma en memoria compartida"

STAR \
    --runMode alignReads \
    --genomeDir "$STAR_IDX" \
    --genomeLoad LoadAndExit

for f in "${files[@]}"; do
    base=$(basename "$f" _clean_1.fastq.gz)
    r2="${CLEAN}/${base}_clean_2.fastq.gz"

    if [[ -v "sample_names[$base]" ]]; then
        out_name="${sample_names[$base]}"
    else
        out_name="$base"
        echo "[WARNING] No se encontró $base en metadata. Se usará nombre original."
    fi

    echo "============================================================"
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Iniciando alineamiento STAR PE"
    echo "Muestra SRR: $base"
    echo "Nombre de salida: $out_name"
    echo "R1: $f"
    echo "R2: $r2"
    echo "Salida: ${BAM}/${out_name}."

    if [[ ! -f "$r2" ]]; then
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: falta archivo R2 para $base: $r2"
        continue
    fi

    if STAR \
        --runThreadN "$threads_star" \
        --genomeDir "$STAR_IDX" \
        --readFilesIn "$f" "$r2" \
        --outFileNamePrefix "${BAM}/${out_name}." \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMunmapped None \
        --readFilesCommand zcat \
        --genomeLoad LoadAndKeep \
	--limitBAMsortRAM 30000000000
    then
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Alineamiento completado exitosamente: $base → $out_name"
    else
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: STAR falló para $base" >&2
        exit 1
    fi

    echo "============================================================"
done

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Todos los alineamientos terminaron."

# ========================== ORDENANDO RESULTADOS ==========================

LOG_FINAL="${STAR_PE}/Logs_final"
LOG_OUT="${STAR_PE}/Logs_out"
SJ="${STAR_PE}/SJ"

mkdir -p "$LOG_FINAL" "$LOG_OUT" "$SJ"

mv -f "$BAM"/*.Log.final.out "$LOG_FINAL"
mv -f "$BAM"/*.out "$LOG_OUT"
mv -f "$BAM"/*.tab "$SJ"
