#!/usr/bin/env bash

set -euo pipefail
shopt -s nullglob

# # Cambio interno de ambientes de conda
# source "$(conda info --base)/etc/profile.d/conda.sh"
# conda activate /home/pablosm/.conda/envs/mi_bioinfo

# Directorios de trabajo
DATA="./data"

FASTQ="$DATA/fastq_raw"
CLEAN="$DATA/fastp"

FASTQC_RAW_1="$DATA/fastqc_1_raw"
FASTQC_RAW_2="$DATA/fastqc_2_raw"

FASTQC_CLEAN_1="$DATA/fastqc_1_clean"
FASTQC_CLEAN_2="$DATA/fastqc_2_clean"

MULTIQC_DIR="$DATA/multiqc"

mkdir -p \
    "$FASTQ" \
    "$CLEAN" \
    "$FASTQC_RAW_1" \
    "$FASTQC_RAW_2" \
    "$FASTQC_CLEAN_1" \
    "$FASTQC_CLEAN_2" \
    "$MULTIQC_DIR" \
    ./logs

# Configuración de hilos
threads_fastp=10
threads_fastqc=4

# Número de bases a remover al inicio
erase_bp=13

echo "========== FASTP =========="

files=("$FASTQ"/SRR*_1.fastq)

for f in "${files[@]}"; do
    base=$(basename "$f" _1.fastq)

    echo "Procesando muestra: $base"

    # Limpieza PE por plataforma Illumina NovaSeq
    fastp \
        -i "${FASTQ}/${base}_1.fastq" \
        -I "${FASTQ}/${base}_2.fastq" \
        -o "$CLEAN/${base}_clean_1.fastq" \
        -O "$CLEAN/${base}_clean_2.fastq" \
        -w "$threads_fastp" \
        --trim_poly_g \
        --trim_front1 "$erase_bp" \
        --trim_front2 "$erase_bp" \
        --detect_adapter_for_pe
done

echo "========== FASTQC SOBRE READS LIMPIOS =========="

clean_files_1=("$CLEAN"/*_clean_1.fastq)

for f1 in "${clean_files_1[@]}"; do
    base=$(basename "$f1" _clean_1.fastq)
    f2="$CLEAN/${base}_clean_2.fastq"

    echo "Evaluando calidad de muestra limpia: $base"

    fastqc \
        -o "$FASTQC_CLEAN_1" \
        -t "$threads_fastqc" \
        "$f1"

    if [[ -f "$f2" ]]; then
        fastqc \
            -o "$FASTQC_CLEAN_2" \
            -t "$threads_fastqc" \
            "$f2"
    else
        echo "Falta archivo par: $f2"
    fi
done

echo "========== MULTIQC =========="

multiqc "$DATA" -o "$MULTIQC_DIR"

echo "Proceso terminado."