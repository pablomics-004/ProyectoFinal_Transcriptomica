#!/usr/bin/env python

# ===================== LIBRERÍAS =====================

import matplotlib.pylab as plt
from pathlib import Path
import seaborn as sns
import pandas as pd
import numpy as np
import re

# ===================== FUNCIONES =====================

def parse_feature_counts(
    fc_path: Path, gene_length: bool = False
) -> tuple[pd.DataFrame, pd.DataFrame] | pd.DataFrame:

    if not (fc_path.exists() and fc_path.is_file()):
        raise ValueError(f"El archivo {fc_path} no es un archivo regular.")

    # =================== LIMPIEZA DE COLUMNAS ===================
    
    def clean_column_name(col: object) -> list[str]:
        if col == "Geneid":
            return col
        
        return col.split("/")[-1].split(".")[0]
    
    # Cargando tabla
    featcount_df = pd.read_csv(fc_path, sep="\t", comment="#")

    # Retirando las versiones del Geneid
    featcount_df["Geneid"] = featcount_df["Geneid"].str.replace(
        r"\.\d+$", "", regex=True
    )

    # ==================== CUENTAS ====================
    counts = featcount_df.drop(
        columns=["Chr", "Start", "End", "Strand", "Length"]
    )
    ## Cambiando el nombre de las columnas para las cuentas
    counts.columns = [
        clean_column_name(col)
        for col in counts.columns
    ]

    ## Ordenando las columnas por edad
    gene_col = ["Geneid"]
    sample_cols = sorted(col for col in counts.columns if col != "Geneid")

    counts = counts[gene_col + sample_cols]

    # ==================== GeneID - Length ====================
    if gene_length:
        geneid_length = featcount_df.drop(
            columns=[
                e 
                for e in featcount_df.columns 
                if e not in {"Geneid", "Length"}
            ]
        )
    
        return counts, geneid_length

    return counts

def parse_gtf_gene_names(
    gtf_path: Path,
    gene_type: bool = False
) -> pd.DataFrame:

    if not (gtf_path.exists() and gtf_path.is_file()):
        raise ValueError(f"El archivo {gtf_path} no es un archivo regular.")

    # =================== FUNCIONES INTERNAS ===================

    def extract_attribute(attribute: str, key: str) -> str | None:
        match = re.search(fr'{key} "([^"]+)"', attribute)

        if match:
            return match.group(1)

        return None

    # =================== CARGANDO GTF ===================

    gtf_df = pd.read_csv(
        gtf_path,
        sep="\t",
        comment="#",
        header=None,
        names=[
            "Chr",
            "source",
            "feature",
            "Start",
            "End",
            "score",
            "Strand",
            "frame",
            "attribute"
        ]
    )

    # =================== FILTRANDO GENES ===================

    genes_df = gtf_df[
        gtf_df["feature"] == "gene"
    ].copy()

    # =================== EXTRACCIÓN DE ATRIBUTOS ===================

    genes_df["Geneid"] = genes_df["attribute"].apply(
        lambda x: extract_attribute(x, "gene_id")
    )

    genes_df["gene_name"] = genes_df["attribute"].apply(
        lambda x: extract_attribute(x, "gene_name")
    )

    genes_df["gene_type"] = genes_df["attribute"].apply(
        lambda x: extract_attribute(x, "gene_type")
    )

    # Retirando versiones del gene_id
    genes_df["Geneid"] = genes_df["Geneid"].str.replace(
        r"\.\d+$", "",
        regex=True
    )

    # =================== TABLA FINAL ===================

    columns = ["Geneid", "gene_name"]

    if gene_type:
        columns.append("gene_type")

    geneid_name = genes_df[columns].drop_duplicates()

    return geneid_name

def extract_star_info(log_file) -> np.ndarray:
    """
    Extrae el porcentaje de lecturas alineadas y el tiempo de ejecución a partir del archivo de log de STAR.
    """

    def standarize_time(time_str: str) -> float:
        """
        Convierte un tiempo en formato "h:m:s" a segundos.
        """
        h, m, s = map(float, time_str.strip().split(":"))
        return h * 3600 + m * 60 + s

    info = np.empty(4, dtype=np.float64) # [aligned_uniq, >1, overall, time_sec]

    with open(log_file, "r") as f:
        for line in f.readlines():
            line = line.strip()
            if "Started mapping on" in line:
                t0 = standarize_time(line.split()[-1])
            if "Finished on" in line:
                t1 = standarize_time(line.split()[-1])
                info[3] = t1 - t0

            if "Uniquely mapped reads %" in line:
                info[0] = float(line.split()[-1].replace("%", ""))

            if "% of reads mapped to multiple loci" in line:
                info[1] = float(line.split()[-1].replace("%", ""))
                info[2] = info[0] + info[1]

    return info

# ===================== MAIN =====================

def main():

    # ===================== DIRECTORIOS + ARCHIVOS =====================
    
    # Directorios
    data = Path("./data/")
    fc_dir = data / "featureCounts"
    images_dir = Path("./docs/images/")

    # Archivos
    gtf = sorted(data.rglob("*.gtf"))[0]
    star_logs = sorted(data.rglob("*.Log.final.out"))
    counts_table, summary = sorted(data.rglob("counts_table.*"))

    # ===================== PARSEO DE TABLAS =====================
    
    counts, geneid_length = parse_feature_counts(counts_table, True)
    geneid_name = parse_gtf_gene_names(gtf)

    # Guardado de los archivos
    counts.to_csv(fc_dir / "counts_cleaned.csv", sep=",", index=False)
    geneid_length.to_csv(fc_dir / "geneid-length.csv", sep=",", index=False)
    geneid_name.to_csv(fc_dir / "geneid-name.csv", sep=",", index=False)

    del counts, counts_table, geneid_length, geneid_name

    # ===================== GRÁFICOS DEL ALINEAMIENTO =====================

    alignment_cols = [
        "Alineamientos únicos",
        "Múltiples loci",
        "Total de alineamientos"
    ]

    time_col = "Tiempo (seg)"

    columns = alignment_cols + [time_col]

    alignment_data = pd.DataFrame(
        {
            srr.name.split(".")[0]: extract_star_info(srr)
            for srr in star_logs
        },
        index=columns
    ).T

    # ===================== PALETAS =====================

    palette = sns.color_palette("viridis", n_colors=len(alignment_cols))
    time_color = sns.color_palette("viridis", n_colors=6)[-2]

    # Para formato largo, útil en seaborn
    alignment_long = alignment_data[alignment_cols].melt(
        var_name="Tipo de alineamiento",
        value_name="Porcentaje"
    )

    # ===================== BARRAS: TASAS =====================

    plt.close("all")
    plt.figure(figsize=(8, 5))

    mean_rates = alignment_data[alignment_cols].mean().reset_index()
    mean_rates.columns = ["Tipo de alineamiento", "Porcentaje promedio"]

    ax = sns.barplot(
        data=mean_rates,
        x="Tipo de alineamiento",
        y="Porcentaje promedio",
        palette=palette,
        edgecolor="black",
        hue="Tipo de alineamiento",
        legend=False
    )

    ax.set_ylabel("Porcentaje promedio de alineamiento (%)")
    ax.set_xlabel("")
    ax.set_title("Promedio de tasas de alineamiento")
    plt.xticks(rotation=20)
    plt.tight_layout()

    plt.savefig(
        images_dir / "promedio_alineamientos.png",
        format="png",
        dpi=400
    )

    # ===================== VIOLINPLOT: TASAS =====================

    plt.close("all")
    plt.figure(figsize=(8, 5))

    ax = sns.violinplot(
        data=alignment_long,
        x="Tipo de alineamiento",
        y="Porcentaje",
        palette=palette,
        inner=None,
        cut=0,
        linewidth=1.2,
        hue="Tipo de alineamiento",
        legend=False
    )

    # Mediana blanca
    medians = alignment_data[alignment_cols].median()

    for i, median in enumerate(medians):
        ax.hlines(
            y=median,
            xmin=i - 0.18,
            xmax=i + 0.18,
            color="white",
            linewidth=2.5,
            zorder=5
        )

    ax.set_ylabel("Porcentaje de alineamiento (%)")
    ax.set_xlabel("")
    ax.set_title("Tasas de alineamiento")
    plt.xticks(rotation=20)
    plt.tight_layout()

    plt.savefig(
        images_dir / "violinplot_alineamientos.png",
        format="png",
        dpi=400
    )

    # ===================== HISTOGRAMA: TIEMPO =====================

    plt.close("all")
    plt.figure(figsize=(7, 5))

    ax = sns.histplot(
        data=alignment_data,
        x=time_col,
        bins=8,
        color=time_color,
        edgecolor="black"
    )

    ax.set_xlabel("Tiempo (seg)")
    ax.set_ylabel("Frecuencia")
    ax.set_title("Tiempo de alineamiento")
    plt.tight_layout()

    plt.savefig(
        images_dir / "hist_tiempos.png",
        format="png",
        dpi=400
    )

    # ===================== VIOLINPLOT: TIEMPO =====================

    plt.close("all")
    plt.figure(figsize=(5, 5))

    ax = sns.violinplot(
        y=alignment_data[time_col],
        color=time_color,
        inner=None,
        cut=0,
        linewidth=1.2
    )

    # Mediana blanca
    median_time = alignment_data[time_col].median()

    ax.hlines(
        y=median_time,
        xmin=-0.15,
        xmax=0.15,
        color="white",
        linewidth=2.5,
        zorder=5
    )

    ax.set_xticks([0])
    ax.set_xticklabels(["Tiempo"])
    ax.set_ylabel("Tiempo (seg)")
    ax.set_title("Tiempo de alineamiento")
    plt.tight_layout()

    plt.savefig(
        images_dir / "violinplot_tiempos.png",
        format="png",
        dpi=400
    )

    return

if __name__ == "__main__":
    main()