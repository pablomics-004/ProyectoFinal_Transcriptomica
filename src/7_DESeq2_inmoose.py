#!/usr/bin/env python3

# ===================== LIBRERÍAS =====================

from sklearn.decomposition import PCA
from inmoose.utils import Factor
import matplotlib.pylab as plt
import inmoose.deseq2 as ds2
from pathlib import Path
from scipy import sparse
import seaborn as sns
import pandas as pd
import numpy as np

# ===================== FUNCIONES =====================

def plot_pca_vst(
    vst_obj,
    ax,
    title=None,
    hue="condition",
    style="Instrument",
    palette=None,
    show_labels=True,
    label_size=8
):

    X = vst_obj.X

    if hasattr(X, "toarray"):
        X = X.toarray()

    X = np.asarray(X)

    pca = PCA(n_components=2)
    coords = pca.fit_transform(X)

    var_exp = pca.explained_variance_ratio_ * 100

    pca_df = pd.DataFrame(
        {
            "PC1": coords[:, 0],
            "PC2": coords[:, 1],
            hue: vst_obj.obs[hue].values,
            "sample": vst_obj.obs_names
        },
        index=vst_obj.obs_names
    )

    if style is not None:
        pca_df[style] = vst_obj.obs[style].values

    plot_kwargs = {
        "data": pca_df,
        "x": "PC1",
        "y": "PC2",
        "hue": hue,
        "palette": palette,
        "ax": ax,
        "s": 90
    }

    if style is not None:
        plot_kwargs["style"] = style

    sns.scatterplot(**plot_kwargs)

    if show_labels:
        for _, row in pca_df.iterrows():
            ax.text(
                row["PC1"],
                row["PC2"],
                row["sample"],
                fontsize=label_size,
                ha="left",
                va="bottom"
            )

    ax.set_xlabel(f"PC1: {var_exp[0]:.0f}% variance")
    ax.set_ylabel(f"PC2: {var_exp[1]:.0f}% variance")

    if title is not None:
        ax.set_title(title, fontweight="bold")

    sns.despine(ax=ax)

    return ax, pca_df

def build_dds_and_vst(
    counts: pd.DataFrame,
    metadata: pd.DataFrame,
    design: str,
    min_count: int = 10
):

    dds = ds2.DESeqDataSet(
        countData=counts.T,
        clinicalData=metadata,
        design=design
    )

    keep = dds.counts().sum(axis=0) >= min_count
    dds = dds[:, keep]

    vst = ds2.varianceStabilizingTransformation(dds)

    vst_values = vst.X

    if sparse.issparse(vst_values):
        vst_values = vst_values.toarray()

    vst_df = pd.DataFrame(
        vst_values,
        index=vst.obs_names,
        columns=vst.var_names
    )

    return dds, vst, vst_df

# ===================== MAIN =====================

def main():

    # ===================== DIRECTORIOS =====================
    
    # Directorios
    results = Path("./results")
    inmoose = results / "diffExp"
    images_dir = Path("./docs/images/")
    featureCounts = Path("./data/featureCounts")

    # Asegurando su generación
    for dir in (results, inmoose, images_dir, featureCounts):
        dir.mkdir(exist_ok=True)

    # ===================== METADATA =====================

    counts = pd.read_csv(
        featureCounts / "counts_cleaned.csv",
        sep=",",
        header=0,
        index_col=0
    )
    counts.index.name = None

    # Cargar metadata del CSV
    metadata_raw = pd.read_csv(
        "./data/metadata/Runs_metadata.csv",
        sep=",",
        header=0
    )

    # Usar SRR_file_name como nombre de muestra, porque debe coincidir con columns de counts
    metadata_raw = metadata_raw.set_index("SRR_file_name")

    # Reordenar metadata para que tenga el mismo orden que counts.columns
    metadata_raw = metadata_raw.loc[counts.columns]

    # Crear condition a partir de source_name
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

    # Contraste explícito: PS vs NS
    metadata["condition"] = metadata["condition"].cat.reorder_categories(
        ["NS", "PS"],
        ordered=True
    )

    metadata["Instrument"] = metadata["Instrument"].cat.reorder_categories(
        ["Illumina HiSeq 2000", "Illumina HiSeq 4000"],
        ordered=True
    )

    metadata.to_csv(
        inmoose / "metadata_model.csv",
        index_label="sample"
    )

    del metadata_raw

        # ===================== ANÁLISIS PRINCIPAL =====================

    dds, vst, vst_df = build_dds_and_vst(
        counts=counts,
        metadata=metadata,
        design="~ Instrument + condition",
        min_count=10
    )

    vst_df.T.to_csv(
        inmoose / "inmoose-vst-matrix_all_samples.csv",
        index_label="Geneid"
    )

    plt.close("all")

    fig, ax = plt.subplots(figsize=(8, 6))

    ax, pca_df = plot_pca_vst(
        vst,
        ax=ax,
        title="PCA con VST: todas las muestras",
        hue="condition",
        style="Instrument",
        palette={"NS": "#4C72B0", "PS": "#DD8452"},
        show_labels=False
    )

    pca_df.to_csv(
        inmoose / "pca_vst_all_samples.csv",
        index_label="sample"
    )

    plt.tight_layout()
    plt.savefig(
        images_dir / "pca_vst_all_samples_condition_instrument.png",
        dpi=400
    )
    plt.close("all")

    # ===================== ANÁLISIS DE SENSIBILIDAD =====================

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

    metadata_hiseq2000.to_csv(
        inmoose / "metadata_sensitivity_hiseq2000.csv",
        index_label="sample"
    )

    dds_hiseq2000, vst_hiseq2000, vst_df_hiseq2000 = build_dds_and_vst(
        counts=counts_hiseq2000,
        metadata=metadata_hiseq2000,
        design="~ condition",
        min_count=10
    )

    vst_df_hiseq2000.T.to_csv(
        inmoose / "inmoose-vst-matrix_hiseq2000_only.csv",
        index_label="Geneid"
    )

    # ===================== PCA =====================

    plt.close("all")

    fig, ax = plt.subplots(figsize=(8, 6))

    ax, pca_df_hiseq2000 = plot_pca_vst(
        vst_hiseq2000,
        ax=ax,
        title="PCA con VST: solo Illumina HiSeq 2000",
        hue="condition",
        style=None,
        palette={"NS": "#4C72B0", "PS": "#DD8452"},
        show_labels=False
    )

    pca_df_hiseq2000.to_csv(
        inmoose / "pca_vst_hiseq2000_only.csv",
        index_label="sample"
    )

    plt.tight_layout()
    plt.savefig(
        images_dir / "pca_vst_hiseq2000_only.png",
        dpi=400
    )
    plt.close("all")

    return

if __name__ == "__main__":
    main()