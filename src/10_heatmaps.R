#!/usr/bin/env Rscript

# ===================== LIBRERÍAS =====================

packages = c(
  "tidyverse",
  "ComplexHeatmap",
  "circlize"
)

for (p in packages) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p)
  }
  library(p, character.only = TRUE)
}

# ===================== FUNCIONES =====================

make_deg_heatmap = function(
  res,
  expr_mat,
  metadata,
  lfc_col = "log2FoldChange",
  fdr_col = "adj_pvalue",
  gene_col = "gene_name",
  LFC_cutoff = 0.5,
  FDR_cutoff = 1e-4,
  top_n = 80,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  heatmap_title = "Heatmap",
  heatmap_name = "Z score"
) {
  
  res_df = res %>%
    as.data.frame() %>%
    rownames_to_column("Geneid") %>%
    mutate(
      lfc_value = as.numeric(.data[[lfc_col]]),
      fdr_value = as.numeric(.data[[fdr_col]]),
      gene_label = ifelse(
        is.na(.data[[gene_col]]) | .data[[gene_col]] == "",
        Geneid,
        .data[[gene_col]]
      ),
      direction = case_when(
        lfc_value > LFC_cutoff & fdr_value < FDR_cutoff ~ "UP in PS",
        lfc_value < -LFC_cutoff & fdr_value < FDR_cutoff ~ "DOWN in PS",
        TRUE ~ "NO"
      )
    )
  
  significant = res_df %>%
    filter(direction != "NO") %>%
    arrange(fdr_value, desc(abs(lfc_value)))
  
  if (nrow(significant) == 0) {
    stop("No hay genes significativos con los cortes indicados.")
  }
  
  significant = significant %>%
    slice_head(n = min(top_n, nrow(significant)))
  
  significant_genes = significant$Geneid
  significant_genes = significant_genes[
    significant_genes %in% rownames(expr_mat)
  ]
  
  if (length(significant_genes) == 0) {
    stop("Los genes significativos no están presentes en expr_mat.")
  }
  
  significant = significant %>%
    filter(Geneid %in% significant_genes)
  
  expr_sig = expr_mat[significant$Geneid, , drop = FALSE]
  
  zscore_significant = t(scale(t(expr_sig)))
  
  keep = complete.cases(zscore_significant)
  zscore_significant = zscore_significant[keep, , drop = FALSE]
  significant = significant[keep, , drop = FALSE]
  
  row_labels = significant$gene_label
  names(row_labels) = significant$Geneid
  
  metadata = metadata[colnames(zscore_significant), , drop = FALSE]
  
  annotation_cols = list(
    condition = c(
      "NS" = "#4C72B0",
      "PS" = "#DD8452"
    ),
    Instrument = c(
      "Illumina HiSeq 2000" = "#404040",
      "Illumina HiSeq 4000" = "#999999"
    )
  )
  
  annotation_df = metadata %>%
    select(any_of(c("condition", "Instrument")))
  
  top_annotation = HeatmapAnnotation(
    df = annotation_df,
    col = annotation_cols[names(annotation_df)],
    annotation_name_side = "right"
  )
  
  heatmap_colors = colorRamp2(
    c(-2, 0, 2),
    c("#377EB8", "white", "#E41A1C")
  )
  
  p = Heatmap(
    zscore_significant,
    col = heatmap_colors,
    cluster_rows = cluster_rows,
    cluster_columns = cluster_columns,
    show_row_names = show_row_names,
    show_column_names = show_column_names,
    row_labels = row_labels,
    row_split = significant$direction,
    name = heatmap_name,
    column_title = heatmap_title,
    top_annotation = top_annotation,
    row_names_gp = grid::gpar(fontsize = 7),
    column_names_gp = grid::gpar(fontsize = 8),
    heatmap_legend_param = list(
      title = heatmap_name
    )
  )
  
  return(p)
}

save_heatmap_png = function(
  heatmap_obj,
  filename_png,
  width = 9,
  height = 8,
  dpi = 400
) {
  
  png(
    filename = filename_png,
    width = width,
    height = height,
    units = "in",
    res = dpi
  )
  
  draw(heatmap_obj)
  dev.off()
}

# ===================== DIRECTORIOS =====================

results_dir = "results/diffExp"
images_dir = "docs/images"
metadata_path = "data/metadata/Runs_metadata.csv"

dir.create(images_dir, recursive = TRUE, showWarnings = FALSE)

# ===================== PARÁMETROS =====================

FDR_cutoff = 1e-5
LFC_cutoff = 1
top_n = 60

# ===================== CARGA DE RESULTADOS =====================

res_all = read.csv(
  file.path(results_dir, "all_hiseq_PS_vs_NS_all_results.csv"),
  row.names = 1,
  check.names = FALSE
)

res_hiseq2000 = read.csv(
  file.path(results_dir, "hiseq2000_PS_vs_NS_all_results.csv"),
  row.names = 1,
  check.names = FALSE
)

# ===================== CARGA DE MATRICES VST =====================

expr_all = read.csv(
  file.path(results_dir, "inmoose-vst-matrix_all_samples.csv"),
  row.names = 1,
  check.names = FALSE
)

expr_hiseq2000 = read.csv(
  file.path(results_dir, "inmoose-vst-matrix_hiseq2000_only.csv"),
  row.names = 1,
  check.names = FALSE
)

expr_all = as.matrix(expr_all)
expr_hiseq2000 = as.matrix(expr_hiseq2000)

# ===================== METADATA =====================

metadata_raw = read.csv(
  metadata_path,
  check.names = FALSE
)

metadata = metadata_raw %>%
  mutate(
    condition = case_when(
      source_name == "Control skin" ~ "NS",
      source_name == "Conventional psoriatic skin" ~ "PS",
      source_name == "conventional psoriatic skin" ~ "PS",
      TRUE ~ NA_character_
    )
  ) %>%
  select(SRR_file_name, condition, Instrument) %>%
  column_to_rownames("SRR_file_name")

metadata_all = metadata[colnames(expr_all), , drop = FALSE]
metadata_hiseq2000 = metadata[colnames(expr_hiseq2000), , drop = FALSE]

# ===================== HEATMAPS =====================

hm_all = make_deg_heatmap(
  res = res_all,
  expr_mat = expr_all,
  metadata = metadata_all,
  lfc_col = "log2FoldChange",
  fdr_col = "adj_pvalue",
  gene_col = "gene_name",
  LFC_cutoff = LFC_cutoff,
  FDR_cutoff = FDR_cutoff,
  top_n = top_n,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  heatmap_title = "Top DEGs: PS vs NS ajustado por instrumento",
  heatmap_name = "Z score"
)

hm_hiseq2000 = make_deg_heatmap(
  res = res_hiseq2000,
  expr_mat = expr_hiseq2000,
  metadata = metadata_hiseq2000,
  lfc_col = "log2FoldChange",
  fdr_col = "adj_pvalue",
  gene_col = "gene_name",
  LFC_cutoff = LFC_cutoff,
  FDR_cutoff = FDR_cutoff,
  top_n = top_n,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  heatmap_title = "Top DEGs: PS vs NS solo Illumina HiSeq 2000",
  heatmap_name = "Z score"
)

# ===================== GUARDADO =====================

save_heatmap_png(
  heatmap_obj = hm_all,
  filename_png = file.path(images_dir, "heatmap_all_hiseq_adjusted.png"),
  width = 9,
  height = 8,
  dpi = 400
)

save_heatmap_png(
  heatmap_obj = hm_hiseq2000,
  filename_png = file.path(images_dir, "heatmap_hiseq2000_only.png"),
  width = 9,
  height = 8,
  dpi = 400
)

rm(list = ls())
invisible(gc())