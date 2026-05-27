#!/usr/bin/env bash

set -euo pipefail # Si un comando falla, el script se detiene
shopt -s nullglob # Evita errores si no se encuentran archivos coincidentes

# Directorios de trabajo
FASTQ="./data/fastq"

# Hilos para agilizar el proceso
threads=8

mkdir -p "${FASTQ}_raw" ./logs/

# Evaluación de calidad de las lecturas utilizando fastqc en paralelo
nohup fastqc "${FASTQ}_raw"/*_1.fastq -t "$threads" -o "${FASTQ}c_1_raw" > ./logs/fastqc_1_raw_explore.log 2>&1 &
nohup fastqc "${FASTQ}_raw"/*_2.fastq -t "$threads" -o "${FASTQ}c_2_raw" > ./logs/fastqc_2_raw_explore.log 2>&1 &
