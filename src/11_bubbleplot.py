#!/usr/bin/env python3

# ===================== LIBRERÍAS =====================

import matplotlib.pyplot as plt
from pathlib import Path
import seaborn as sns
import pandas as pd
import numpy as np

# ===================== FUNCIONES =====================

def bubbleplot(
  data: pd.DataFrame,
  ax,
  ratio_col: str | None = None,
  term_col: str = "Term",
  count_col: str = "Count",
  padj_col: str = "FDR",
  top_n: int = 25,
  title: str | None = None,
  cmap: str = "viridis"
) -> None:

  data2plot = data.copy()

  # Reescalando el p-adj
  data2plot[padj_col] = (
    data2plot[padj_col]
    .replace(0, np.nextafter(0, 1))
  )
  data2plot["-log10(FDR)"] = -np.log10(data2plot[padj_col])

  # Selección de los términos más significativos
  data2plot = (
    data2plot.sort_values(padj_col)
    .head(top_n)
    .copy()
  )

  # Ordenar términos para que el más significativo quede arriba
  data2plot[term_col] = pd.Categorical(
    data2plot[term_col],
    categories=data2plot[term_col][::-1],
    ordered=True
  )

  if ratio_col is None:
    ratio_col = "Gene Ratio"
    data2plot[ratio_col] = data2plot[count_col] / data2plot["List Total"]

  elif ratio_col == "Fold Enrichment":
    data2plot["Fold Enrichment / max(Fold Enrichment)"] = (
      data2plot[ratio_col] / data[ratio_col].max()
    )
    ratio_col = "Fold Enrichment / max(Fold Enrichment)"

  sns.scatterplot(
    data=data2plot,
    x=ratio_col,
    y=term_col,
    size=count_col,
    hue="-log10(FDR)",
    sizes=(80, 600),
    palette=cmap,
    edgecolor="black",
    linewidth=0.4,
    alpha=0.85,
    ax=ax
  )

  ax.set_title(title, fontweight="bold")
  ax.set_xlabel(ratio_col)
  ax.set_ylabel("")
  ax.grid(axis="x", linestyle="--", alpha=0.3)

  sns.despine(ax=ax, left=True)

  return ax

# ===================== MAIN =====================

def main() -> None:

    # ===================== DIRECTORIOS =====================

    enrichment_dir = Path("./results/functional_enrichment")
    images_dir = Path("./docs/images")

    images_dir.mkdir(parents=True, exist_ok=True)

    # ===================== ARCHIVOS DAVID =====================

    david_files = sorted(enrichment_dir.glob("DAVID*.csv"))

    if len(david_files) == 0:
        raise FileNotFoundError(
            f"No se encontraron archivos con prefijo DAVID en {enrichment_dir}"
        )

    # ===================== GRAFICADO =====================

    for david_file in david_files:

        print(f"Procesando: {david_file}")

        data = pd.read_csv(david_file)

        # Verificación mínima de columnas necesarias
        required_cols = {"Term", "Count", "FDR"}

        missing_cols = required_cols - set(data.columns)

        if missing_cols:
            print(
                f"[WARNING] Se omite {david_file.name}. "
                f"Faltan columnas: {missing_cols}"
            )
            continue

        # Título a partir del nombre del archivo
        title = david_file.stem.replace("DAVID_", "").replace("_", " ")

        # Nombre de salida
        out_png = images_dir / f"{david_file.stem}.png"

        plt.close("all")

        fig, ax = plt.subplots(figsize=(9, 7))

        bubbleplot(
            data=data,
            ax=ax,
            ratio_col="Fold Enrichment" if "Fold Enrichment" in data.columns else None,
            term_col="Term",
            count_col="Count",
            padj_col="FDR",
            top_n=25,
            title=title,
            cmap="viridis"
        )

        plt.tight_layout()

        plt.savefig(
            out_png,
            dpi=400,
            bbox_inches="tight"
        )

        plt.close("all")

        print(f"Guardado: {out_png}")

    return

if __name__ == "__main__":
    main()