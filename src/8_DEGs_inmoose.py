#!/usr/bin/env python3

# ===================== LIBRERÍAS =====================

from pathlib import Path
from inmoose.utils import Factor
import inmoose.deseq2 as ds2
import pandas as pd

# ===================== FUNCIONES =====================

def build_metadata(counts: pd.DataFrame, metadata_path: Path) -> pd.DataFrame:

    metadata_raw = pd.read_csv(metadata_path, sep=",", header=0)

    metadata_raw = metadata_raw.set_index("SRR_file_name")
    metadata_raw = metadata_raw.loc[counts.columns]

    condition = metadata_raw["source_name"].map(
        {
            "Control skin": "NS",
            "conventional psoriatic skin": "PS",
            "Conventional psoriatic skin": "PS"
        }
    )

    metadata = pd.DataFrame(
        {
            "condition": Factor(condition),
            "Instrument": Factor(metadata_raw["Instrument"])
        },
        index=counts.columns
    )

    metadata["condition"] = metadata["condition"].cat.reorder_categories(
        ["NS", "PS"],
        ordered=True
    )

    metadata["Instrument"] = metadata["Instrument"].cat.reorder_categories(
        ["Illumina HiSeq 2000", "Illumina HiSeq 4000"],
        ordered=True
    )

    return metadata

def add_gene_names(df: pd.DataFrame, gene_names: pd.DataFrame) -> pd.DataFrame:

    df_out = df.copy()

    df_out = df_out.merge(
        gene_names,
        left_index=True,
        right_on="Geneid",
        how="left"
    )

    df_out = df_out.set_index("Geneid")

    cols = ["gene_name"] + [
        col for col in df_out.columns
        if col != "gene_name"
    ]

    df_out = df_out[cols]

    return df_out

def run_deseq_model(
    counts: pd.DataFrame,
    metadata: pd.DataFrame,
    design: str,
    contrast: list[str],
    out_prefix: Path,
    gene_names: pd.DataFrame,
    fdr_cutoff: float = 1e-4,
    lfc_cutoff: float = 0.5,
    min_count: int = 10
) -> pd.DataFrame:

    dds = ds2.DESeqDataSet(
        countData=counts.T,
        clinicalData=metadata,
        design=design
    )

    keep = dds.counts().sum(axis=0) >= min_count
    dds = dds[:, keep]

    dds = ds2.DESeq(dds)

    print(f"\nModelo: {design}")
    print(f"Contrastes disponibles: {dds.resultsNames()}")

    res = dds.results(
        contrast=contrast
    )

    res = pd.DataFrame(res)

    if "betaConv" in dds.var.columns:
        beta_conv = dds.var["betaConv"].astype("boolean")
        n_not_converged = beta_conv.eq(False).sum()

        print(f"Genes sin converger: {n_not_converged}")

        beta_conv_res = beta_conv.reindex(res.index).fillna(False)
        res = res.loc[beta_conv_res]

        not_converged = dds.var.index[beta_conv.eq(False)]

        pd.Series(not_converged).to_csv(
            f"{out_prefix}_genes_not_beta_converged.csv",
            index=False,
            header=["Geneid"]
        )

    res = add_gene_names(
        res,
        gene_names
    )

    degs = res[
        (res["adj_pvalue"] < fdr_cutoff) &
        (res["log2FoldChange"].abs() > lfc_cutoff)
    ].copy()

    degs_up = degs[
        degs["log2FoldChange"] > lfc_cutoff
    ].copy()

    degs_down = degs[
        degs["log2FoldChange"] < -lfc_cutoff
    ].copy()

    res.to_csv(
        f"{out_prefix}_all_results.csv",
        index=True
    )

    degs.to_csv(
        f"{out_prefix}_DEGs_FDR_{fdr_cutoff}_LFC_{lfc_cutoff}.csv",
        index=True
    )

    degs_up.to_csv(
        f"{out_prefix}_UP_in_PS_FDR_{fdr_cutoff}_LFC_{lfc_cutoff}.csv",
        index=True
    )

    degs_down.to_csv(
        f"{out_prefix}_DOWN_in_PS_FDR_{fdr_cutoff}_LFC_{lfc_cutoff}.csv",
        index=True
    )

    summary = pd.DataFrame(
        {
            "model": [design],
            "contrast": ["PS_vs_NS"],
            "FDR_cutoff": [fdr_cutoff],
            "LFC_cutoff": [lfc_cutoff],
            "total_tested_genes": [res.shape[0]],
            "total_DEGs": [degs.shape[0]],
            "up_in_PS": [degs_up.shape[0]],
            "down_in_PS": [degs_down.shape[0]]
        }
    )

    summary.to_csv(
        f"{out_prefix}_summary.csv",
        index=False
    )

    print(summary.to_string(index=False))

    return res

# ===================== MAIN =====================

def main() -> None:

    # ===================== RUTAS =====================

    feature_counts = Path("./data/featureCounts")
    metadata_path = Path("./data/metadata/Runs_metadata.csv")
    out_dir = Path("./results/diffExp")

    out_dir.mkdir(parents=True, exist_ok=True)

    counts_path = feature_counts / "counts_cleaned.csv"

    # Nombres de los genes
    gene_names_path = feature_counts / "geneid-name.csv"
    gene_names = pd.read_csv(
        gene_names_path,
        sep=",",
        header=0
    )

    gene_names = gene_names[["Geneid", "gene_name"]].drop_duplicates()

    # ===================== PARÁMETROS =====================

    FDR_cutoff = 1e-4
    LFC_cutoff = 0.5
    min_count = 10

    contrast = ["condition", "PS", "NS"]

    # ===================== CARGA =====================

    counts = pd.read_csv(
        counts_path,
        sep=",",
        header=0,
        index_col=0
    )

    counts.index.name = None

    metadata = build_metadata(
        counts=counts,
        metadata_path=metadata_path
    )

    metadata.to_csv(
        out_dir / "metadata_model.csv",
        index_label="sample"
    )

    # ===================== MODELO PRINCIPAL =====================

    run_deseq_model(
        counts=counts,
        metadata=metadata,
        design="~ Instrument + condition",
        contrast=contrast,
        out_prefix=out_dir / "all_hiseq_PS_vs_NS",
        fdr_cutoff=FDR_cutoff,
        lfc_cutoff=LFC_cutoff,
        min_count=min_count,
        gene_names=gene_names
    )

    # ===================== MODELO DE SENSIBILIDAD =====================

    hiseq2000_samples = metadata[
        metadata["Instrument"] == "Illumina HiSeq 2000"
    ].index

    counts_hiseq2000 = counts[hiseq2000_samples].copy()
    metadata_hiseq2000 = metadata.loc[hiseq2000_samples].copy()

    metadata_hiseq2000["condition"] = Factor(
        metadata_hiseq2000["condition"].astype(str)
    )

    metadata_hiseq2000["condition"] = (
        metadata_hiseq2000["condition"]
        .cat.reorder_categories(
            ["NS", "PS"],
            ordered=True
        )
    )

    metadata_hiseq2000 = metadata_hiseq2000[
        ["condition"]
    ].copy()

    metadata_hiseq2000.to_csv(
        out_dir / "metadata_hiseq2000_only.csv",
        index_label="sample"
    )

    run_deseq_model(
        counts=counts_hiseq2000,
        metadata=metadata_hiseq2000,
        design="~ condition",
        contrast=contrast,
        out_prefix=out_dir / "hiseq2000_PS_vs_NS",
        fdr_cutoff=FDR_cutoff,
        lfc_cutoff=LFC_cutoff,
        min_count=min_count,
        gene_names=gene_names
    )


if __name__ == "__main__":
    main()