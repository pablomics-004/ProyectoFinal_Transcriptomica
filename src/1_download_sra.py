#!/usr/bin/env python3

"""
Ejecutar de la siguiente forma desde la raíz del repositorio:

```shell

conda activate /export/space4/users/yaelmont/TRANSCRIPTOMICA/rna

./src/1_download_sra.py
```
"""

# =============== LIBRERÍAS ===============

from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime
from pathlib import Path
import subprocess as sb
import pandas as pd
import os
import shutil

# =============== FUNCIONES ===============

def srr_file_name(row) -> str:
    condition = "PS" if row["source_name"] == "Conventional psoriatic skin" else "NS"
    return f"{condition}_{row['Run']}"


def to_root() -> None:
    dirs = {"src", "data", "results", "docs"}

    pwd = os.getcwd()

    if pwd.split("/")[-1] in dirs:
        os.chdir("..")
        print(f"Working directory changed to: {os.getcwd()}", flush=True)
    else:
        print(f"Current working directory: {pwd}", flush=True)


def download_srr(
    srr: str,
    fastq_dir: Path,
    tmp_parent: Path,
    threads: str
) -> str:

    fastq_1 = fastq_dir / f"{srr}_1.fastq"
    fastq_2 = fastq_dir / f"{srr}_2.fastq"

    if fastq_1.exists() and fastq_2.exists():
        return f"[{datetime.now().strftime('%d-%H:%M:%S')}] {srr} ya existe. Saltando."

    tmp_dir = tmp_parent / srr
    tmp_dir.mkdir(parents=True, exist_ok=True)

    print(
        f"[{datetime.now().strftime('%d-%H:%M:%S')}] "
        f"Descargando {srr} con fasterq-dump...",
        flush=True
    )

    try:
        sb.run([
            "fasterq-dump", srr,
            "--split-files",
            "--threads", threads,
            "--outdir", str(fastq_dir),
            "-t", str(tmp_dir)
        ], check=True)

    except Exception:
        print(f"[ERROR] Falló fasterq-dump para {srr}", flush=True)
        raise

    finally:
        if tmp_dir.exists():
            shutil.rmtree(tmp_dir, ignore_errors=True)

    return f"[{datetime.now().strftime('%d-%H:%M:%S')}] Finalizó {srr}"


# =============== MAIN ===============

def main() -> None:

    to_root()

    # ================== GENERACIÓN DE METADATOS ==================

    sra_table_path = Path("data/metadata/SraRunTable.csv")

    if not sra_table_path.is_file():
        raise FileExistsError(
            f"No existe el archivo o no se trata de un archivo regular: {str(sra_table_path)}"
        )

    print("-" * 50, flush=True)
    print(f"Adecuación de metadatos: {sra_table_path}", flush=True)

    raw_table = pd.read_csv(
        sra_table_path,
        sep=",",
        header=0
    )

    srr_table = raw_table[
        raw_table["source_name"].isin([
            "Conventional psoriatic skin",
            "Control skin"
        ])
    ].copy()

    srr_table = srr_table[
        ["Run", "source_name", "Instrument"]
    ].copy()

    srr_table["SRR_file_name"] = srr_table.apply(
        srr_file_name,
        axis=1
    )

    srr_table = srr_table.sort_values(
        ["source_name", "Run"]
    ).reset_index(drop=True)

    output_metadata = Path("data/metadata/Runs_metadata.csv")

    srr_table.to_csv(
        output_metadata,
        sep=",",
        index=False
    )

    print("Muestras seleccionadas:", flush=True)
    print(srr_table["source_name"].value_counts(), flush=True)
    print(f"Metadata guardada en: {output_metadata}", flush=True)
    print("-" * 50, flush=True)
    print("Adecuación completada.", flush=True)

    # ================== DESCARGA DE SRRs ==================

    srr_files = list(srr_table["Run"])

    threads = "4"
    max_workers = 3

    print(f"Descarga de los archivos SRR: {srr_files}", flush=True)
    print(f"Descargas en paralelo: {max_workers}", flush=True)
    print(f"Threads por descarga: {threads}", flush=True)
    print("-" * 50, flush=True)

    fastq_dir = Path("./data/fastq_raw/")
    tmp_parent = Path("./tmp/")

    for dir_path in (fastq_dir, tmp_parent):
        dir_path.mkdir(parents=True, exist_ok=True)

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {
            executor.submit(
                download_srr,
                srr,
                fastq_dir,
                tmp_parent,
                threads
            ): srr
            for srr in srr_files
        }

        for future in as_completed(futures):
            srr = futures[future]

            try:
                message = future.result()
                print(message, flush=True)
                print("-" * 50, flush=True)

            except Exception as e:
                print(f"[ERROR] Falló la descarga de {srr}: {e}", flush=True)
                raise

    print(
        f"[{datetime.now().strftime('%d-%H:%M:%S')}] Todas las descargas terminaron.",
        flush=True
    )


if __name__ == "__main__":
    main()
