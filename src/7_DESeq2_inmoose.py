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

    sns.scatterplot(
        data=pca_df,
        x="PC1",
        y="PC2",
        hue=hue,
        palette=palette,
        ax=ax,
        s=70
    )

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

    # Cargado de datos
    counts = pd.read_csv(
        featureCounts / "counts_cleaned.csv", 
        sep=",", header=0, index_col=0
    )
    counts.index.name = None

    # Generación de los metadatos
    sample_names = list(counts.columns)
    condition = [col.split("_")[0] for col in counts.columns]

    metadata = pd.DataFrame(
        {"condition" : Factor(condition)}, 
        index=sample_names
    )

    # Contraste explícito
    metadata.condition = (
        metadata.condition
        .cat.reorder_categories(
            ["NS", "PS"], ordered=True
        )
    )

    # ===================== DESeq2 (Inmoose) =====================

    dds = ds2.DESeqDataSet(
        countData=counts.T,
        clinicalData=metadata,
        design="~ condition"
    )

    # Criterio arbitrario para el filtrado de genes con baja expresión
    keep = dds.counts().sum(axis=0) >= 10
    dds = dds[:, keep]
    del keep

    # Transformación estabilizadora de varianza
    vst = ds2.varianceStabilizingTransformation(dds)
    vst_values = vst.X

    if sparse.issparse(vst_values):
        vst_values = vst_values.toarray()

    # DataFrame: muestras x genes
    vst_df = pd.DataFrame(
        vst_values,
        index=vst.obs_names,
        columns=vst.var_names
    )

    # Guardado
    vst_df.T.to_csv(
        inmoose / "inmoose-vst-matrix.csv",
        index_label="Geneid"
    )

    # ===================== PCA =====================

    plt.close("all")

    fig, ax = plt.subplots(figsize=(8, 6))

    ax, pca_df = plot_pca_vst(
        vst,
        ax=ax,
        title="PCA con VST",
        hue="condition",
        palette={"NS": "#4C72B0", "PS": "#DD8452"},
        show_labels=False
    )

    plt.tight_layout()
    plt.savefig(images_dir / "pca_vst_labeled.png", dpi=400)
    
    plt.close("all")

    return

if __name__ == "__main__":
    main()