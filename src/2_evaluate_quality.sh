#!/usr/bin/env bash

set -euo pipefail # Si un comando falla, el script se detiene
shopt -s nullglob # Evita errores si no se encuentran archivos coincidentes

# Activando ambiente de conda para fastqc
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate multiqc

# Directorios de trabajo
DATA="./data"
FASTQ="$DATA/fastq_raw"

mkdir -p "$FASTQ" ./logs/

# Evaluación de calidad de las lecturas utilizando fastqc en paralelo
nohup fastqc "$FASTQ"/*_1.fastq -o "${FASTQ}c_1_raw" > ./logs/fastqc_1_explore.log 2>&1 &
nohup fastqc "$FASTQ"/*_2.fastq -o "${FASTQ}c_2_raw" > ./logs/fastqc_2_explore.log 2>&1 &
