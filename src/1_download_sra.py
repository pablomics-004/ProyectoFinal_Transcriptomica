#!/usr/bin/env python3

"""
Ejecutar de la siguiente forma desde la raíz del repositorio:

```shell

conda activate /export/space4/users/yaelmont/TRANSCRIPTOMICA/rna

./src/1_download_sra.py
```
"""

# =============== LIBRERÍAS ===============

from datetime import datetime
from pathlib import Path
import subprocess as sb
import pandas as pd
import os

# =============== FUNCIONES ===============

def srr_file_name(row) -> str:
    treatment = "treated" if "GLPG3970" in row["treatment"] else "placebo"
    sex = "F" if row["sex"] == "female" else "M"

    new_name = (
        f"{treatment}_" +
        f"{row['AGE']}_" +
        f"{sex}_" +
        f"{row['Run']}"
    )

    return new_name

def to_root() -> None:
    dirs = {"src", "data", "results", "docs"}
    
    if (pwd := os.getcwd()).split("/")[-1] in dirs:
        os.chdir("..")
        print(f"Working directory changed to: {os.getcwd()}", flush=True)
    else:
        print(f"Current working directory: {pwd}", flush=True)

    return

# =============== MAIN ===============

def main() -> None:

    to_root()

    # ================== GENERACIÓN DE METADATOS ==================

    sra_table_path = Path("data/metadata/SraRunTable.csv")
    if not (sra_table_path.is_file() or sra_table_path.exists()):
        raise FileExistsError(
            f"No existe el archivo o no se trata de un archivo regular: {str(sra_table_path)}"
        )

    print("-" * 50, flush=True)
    print(f"Adecuación de metadatos: {sra_table_path}", flush=True)

    raw_table = pd.read_csv(
        sra_table_path, sep=",", header=0
    )

    # Nos quedamos solo con muestras Week 6
    srr_table = raw_table[
        raw_table["timepoint"] == "Week 6"
    ].copy()

    # Selección y orden de columnas
    srr_table = srr_table[
        ["Run", "AGE", "treatment", "sex"]
    ].copy()

    # Orden categórico para referencia posterior en DESeq2/inmoose
    srr_table["treatment"] = pd.Categorical(
        srr_table["treatment"],
        categories=["Placebo", "GLPG3970 350 mg q.d."],
        ordered=True
    )

    srr_table["sex"] = pd.Categorical(
        srr_table["sex"],
        categories=["female", "male"],
        ordered=True
    )

    # Orden visual de la tabla
    srr_table = srr_table.sort_values(
        ["treatment", "AGE", "sex", "Run"]
    ).reset_index(drop=True)

    # Nombre útil para renombrar archivos después
    srr_table["SRR_file_name"] = srr_table.apply(srr_file_name, axis=1)

    # Guardado de la tabla
    srr_table.to_csv(
        "data/metadata/Runs_metadata.csv", sep=",",
        index=False
    )
    
    print("-" * 50, flush=True)
    print("Adecuación completada.")

    del raw_table, sra_table_path

    # ================== DESCARGA DE SRRs ==================

    srr_files = list(srr_table["Run"])
    threads = "15"

    print(f"Descarga de los archivos SRR: {srr_files}", flush=True)
    print("-" * 50, flush=True)

    fastq_dir = Path("./data/fastq_raw/")
    tmp_dir = Path("./tmp/")

    for dir_path in (fastq_dir, tmp_dir):
        dir_path.mkdir(parents=True, exist_ok=True)

    for srr in srr_files:
        
        # Nombres de archivos
        fastq_1 = fastq_dir / f"{srr}_1.fastq"
        fastq_2 = fastq_dir / f"{srr}_2.fastq"

        if fastq_1.exists() and fastq_2.exists():
            print(f"[{datetime.now().strftime('%d-%H:%M:%S')}] {srr} ya existe. Saltando.", flush=True)
            continue

        print(f"[{datetime.now().strftime('%d-%H:%M:%S')}] Descargando {srr} con fasterq-dump...", flush=True)

        sb.run([
            "fasterq-dump", srr,
            "--split-files",
            "--threads", threads,
            "--outdir", str(fastq_dir),
            "-t", str(tmp_dir)
        ], check=True)

        print(f"[{datetime.now().strftime('%d-%H:%M:%S')}] Finalizando el procesado de {srr}.", flush=True)
        print("-" * 50, flush=True)

    return

if __name__ == "__main__":
    main()
