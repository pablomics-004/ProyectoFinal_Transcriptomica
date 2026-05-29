#!/usr/bin/env bash

set -euo pipefail
shopt -s nullglob

# Cambio interno de ambientes de conda
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /home/pablosm/.conda/envs/mi_bioinfo

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
threads_fastp=8
threads_fastqc=4

# Número de procesos fastp en paralelo
wks=1

# Número de bases a remover al inicio
erase_bp=13

# Manejo de procesos en paralelo
pids=()
batch_processing() {
    if [[ ${#pids[@]} -eq $wks ]]; then
        wait "${pids[@]}"
        current_date_time="$(date "+%Y-%m-%d %H:%M:%S")"
        echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
        echo "$current_date_time"
        echo "Lote de $wks procesos completados. Continuando con el script..."
        pids=()
    fi
}

echo "========== FASTP =========="

files=("$FASTQ"/SRR*_1.fastq)

for f in "${files[@]}"; do
    base=$(basename "$f" _1.fastq)

    echo "Procesando muestra: $base"

    # Limpieza PE por plataforma Illumina HiSeq 4000
    fastp \
        -i "${FASTQ}/${base}_1.fastq" \
        -I "${FASTQ}/${base}_2.fastq" \
        -o "$CLEAN/${base}_clean_1.fastq.gz" \
        -O "$CLEAN/${base}_clean_2.fastq.gz" \
        -w "$threads_fastp" \
        --trim_front1 "$erase_bp" \
        --trim_front2 "$erase_bp" \
        --detect_adapter_for_pe &

    # Guardado del JOBID para control de procesos
    pids+=("$!")
    batch_processing
done

# Manejo de procesos restantes
if [[ ${#pids[@]} -gt 0 ]]; then
    wait "${pids[@]}"
    pids=()
fi

echo "========== FASTQC SOBRE READS LIMPIOS =========="
conda activate base

clean_files_1=("$CLEAN"/*_clean_1.fastq.gz)

for f1 in "${clean_files_1[@]}"; do
    base=$(basename "$f1" _clean_1.fastq.gz)
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
conda activate /home/pablosm/.conda/envs/mi_bioinfo

multiqc "$DATA" -o "$MULTIQC_DIR"

echo "Proceso terminado."
