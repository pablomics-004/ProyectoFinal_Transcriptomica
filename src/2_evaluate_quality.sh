#!/usr/bin/env bash

echo "Entró al script"

set -euo pipefail # Si un comando falla, el script se detiene
shopt -s nullglob # Evita errores si no se encuentran archivos coincidentes

# Directorios de trabajo
FASTQ="./data/fastq"
FASTQC1="${FASTQ}c_1_raw"
FASTQC2="${FASTQ}c_2_raw"

# Hilos para FASTQC
threads=8

mkdir -p "${FASTQ}_raw" ./logs/ "$FASTQC1" "$FASTQC2"

echo "================================ FASTQC ======================================="

# Evaluación de calidad de las lecturas utilizando fastqc en paralelo

echo "Procesando FASTQ_1"
nohup fastqc "${FASTQ}_raw"/*_1.fastq -t "$threads" -o "$FASTQC1" > ./logs/fastqc_1_raw_explore.log 2>&1 &

echo "Procesando FASTQ_2"
nohup fastqc "${FASTQ}_raw"/*_2.fastq -t "$threads" -o "$FASTQC2" > ./logs/fastqc_2_raw_explore.log 2>&1 &

status=$?

if [[ $status -eq 0 ]]; then
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Script completado"
else
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: FASTQC falló para $base con código $status"
    exit "$status"
fi
