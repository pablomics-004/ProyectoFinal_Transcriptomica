#!/usr/bin/env python3

# =============== LIBRERÍAS ===============

from datetime import datetime
from pathlib import Path
import subprocess as sb
import pandas as pd
import os

# =============== FUNCIONES ===============

def srr_file_name(row) -> str:
    tissue = "PS" if "Psoriatic" in row["tissue"] else "NS"

    # XX_tissue_SRR{X}+8
    new_name = (
        f"{row["age_of_onset"]}_" +
        f"{tissue}_" + f"{row["Run"]}"
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

    srr_table = raw_table[["Run", "age_of_onset", "tissue", "BioSample"]]
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
    threads = "5"

    print(f"Descarga de los archivos SRR: {srr_files}", flush=True)
    print("-" * 50, flush=True)

    sra_dir = Path("./data/sra/")
    fastq_dir = Path("./data/fastq_raw/")
    tmp_dir = Path("./tmp/")

    for dir in (sra_dir, fastq_dir, tmp_dir):
        dir.mkdir(exist_ok=True)

    for srr in srr_files:
        
        # Nombres de archivos
        sra_file = sra_dir / srr / f"{srr}.sra"
        fastq_1 = fastq_dir / srr / f"{srr}_1.fastq"
        fastq_2 = fastq_dir / srr / f"{srr}_2.fastq"

        if not sra_file.exists():
            print(f"[{datetime.now().strftime('%d-%H:%M:%S')}] Descargando {srr}...", flush=True)
            sb.run(["prefetch", srr, "--output-directory", sra_dir], check=True)

        print(f"[{datetime.now().strftime('%d-%H:%M:%S')}] Procesando {srr}...", flush=True)

        if not (fastq_1.exists() and fastq_2.exists()):
            sb.run([
                "fasterq-dump", str(sra_file),
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