#!/usr/bin/env python3

# ===================== LIBRERÍAS =====================

from pathlib import Path
import matplotlib.pyplot as plt
from gseapy import prerank
from gseapy.plot import gseaplot2
import pandas as pd


# ===================== FUNCIONES =====================

def make_ranking(
    deg_path: Path,
    gene_col: str = "gene_name",
    rank_col: str = "stat"
) -> pd.DataFrame:

    data = pd.read_csv(
        deg_path,
        sep=",",
        header=0,
        index_col=0
    )

    required_cols = {gene_col, rank_col}

    missing = required_cols - set(data.columns)

    if missing:
        raise ValueError(
            f"Faltan columnas en {deg_path.name}: {missing}"
        )

    data = data.dropna(
        subset=[gene_col, rank_col]
    )

    ranking = (
        data[[gene_col, rank_col]]
        .rename(columns={gene_col: "Gene_name", rank_col: "Rank"})
        .assign(Gene_name=lambda x: x["Gene_name"].astype(str).str.upper())
        .drop_duplicates("Gene_name")
        .sort_values("Rank", ascending=False)
    )

    return ranking


def run_gsea_prerank(
    ranking: pd.DataFrame,
    out_dir: Path,
    prefix: str,
    gene_sets: list[str],
    organism: str = "Human",
    min_size: int = 5,
    max_size: int = 1000,
    permutation_num: int = 1000,
    seed: int = 42,
    top_n_terms: int = 8
) -> None:

    out_dir.mkdir(parents=True, exist_ok=True)

    pre_res = prerank(
        rnk=ranking,
        gene_sets=gene_sets,
        organism=organism,
        seed=seed,
        min_size=min_size,
        max_size=max_size,
        permutation_num=permutation_num,
        outdir=None,
        verbose=True
    )

    out_results = pre_res.res2d.sort_values(
        by="FDR q-val",
        ascending=True
    )

    out_results.to_csv(
        out_dir / f"{prefix}_gseapy_results.csv",
        sep=",",
        index=False
    )

    sig_results = out_results[
        out_results["FDR q-val"] < 0.25
    ].copy()

    sig_results.to_csv(
        out_dir / f"{prefix}_gseapy_results_FDR025.csv",
        sep=",",
        index=False
    )

    terms = out_results["Term"].head(top_n_terms).tolist()

    if len(terms) == 0:
        print(f"[WARNING] No hay términos para graficar en {prefix}")
        return

    hits = [
        pre_res.results[t]["hits"]
        for t in terms
    ]

    runes = [
        pre_res.results[t]["RES"]
        for t in terms
    ]

    plt.close("all")

    fig = gseaplot2(
        terms=terms,
        RESs=runes,
        hits=hits,
        rank_metric=pre_res.ranking,
        legend_kws={"loc": (1.02, 0.2)},
        figsize=(6, 6)
    )

    plt.suptitle(
        f"GSEA prerank: {prefix}",
        fontweight="bold"
    )

    plt.tight_layout()

    plt.savefig(
        out_dir / f"{prefix}_gseapy_top_terms.png",
        dpi=400,
        bbox_inches="tight"
    )

    plt.close("all")


# ===================== MAIN =====================

def main() -> None:

    # ===================== DIRECTORIOS =====================

    diffexp_dir = Path("./results/diffExp")
    gsea_dir = Path("./results/functional_enrichment/GSEA")
    images_dir = Path("./docs/images")

    gsea_dir.mkdir(parents=True, exist_ok=True)
    images_dir.mkdir(parents=True, exist_ok=True)

    # ===================== ARCHIVOS =====================

    analyses = {
        "all_hiseq_PS_vs_NS": diffexp_dir / "all_hiseq_PS_vs_NS_all_results.csv",
        "hiseq2000_PS_vs_NS": diffexp_dir / "hiseq2000_PS_vs_NS_all_results.csv"
    }

    # ===================== GENE SETS =====================

    gene_sets = [
        "GO_Biological_Process_2023",
        "KEGG_2021_Human",
        "Reactome_2022"
    ]

    # ===================== GSEA =====================

    for prefix, deg_path in analyses.items():

        if not deg_path.is_file():
            print(f"[WARNING] No existe el archivo: {deg_path}")
            continue

        print("=" * 80)
        print(f"Procesando: {prefix}")
        print(f"Archivo: {deg_path}")

        ranking = make_ranking(
            deg_path=deg_path,
            gene_col="gene_name",
            rank_col="stat"
        )

        ranking.to_csv(
            gsea_dir / f"{prefix}_ranking.csv",
            sep=",",
            index=False
        )

        run_gsea_prerank(
            ranking=ranking,
            out_dir=gsea_dir,
            prefix=prefix,
            gene_sets=gene_sets,
            organism="Human",
            min_size=5,
            max_size=1000,
            permutation_num=1000,
            seed=42,
            top_n_terms=8
        )

        # Copia visual al directorio de imágenes
        src_img = gsea_dir / f"{prefix}_gseapy_top_terms.png"
        dst_img = images_dir / f"{prefix}_gseapy_top_terms.png"

        if src_img.is_file():
            dst_img.write_bytes(src_img.read_bytes())

    print("=" * 80)
    print("GSEA finalizado.")


if __name__ == "__main__":
    main()