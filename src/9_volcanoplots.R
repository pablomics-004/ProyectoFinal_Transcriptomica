#!/usr/bin/env Rscript

# ===================== LIBRERÍAS =====================

packages = c("ggplot2", "ggrepel", "tidyverse")

for (p in packages) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p)
  }
  library(p, character.only = TRUE)
}

# ===================== FUNCIONES =====================

make_volcano_plot = function(
  res,
  lfc_col = "log2FoldChange",
  fdr_col = "adj_pvalue",
  gene_col = "gene_name",
  LFC_cutoff = 0.5,
  FDR_cutoff = 1e-4,
  title = "Volcano plot",
  volcano_LFC_limit = 12,
  volcano_FDR_limit = 20,
  n_labels = 4
) {
  
  res_plot = res %>%
    as.data.frame() %>%
    rownames_to_column("Geneid") %>%
    mutate(
      gene_label = ifelse(
        is.na(.data[[gene_col]]) | .data[[gene_col]] == "",
        Geneid,
        .data[[gene_col]]
      ),
      lfc_value = as.numeric(.data[[lfc_col]]),
      fdr_value = as.numeric(.data[[fdr_col]]),
      fdr_plot = case_when(
        is.na(fdr_value) ~ NA_real_,
        fdr_value == 0 ~ .Machine$double.xmin,
        TRUE ~ fdr_value
      ),
      neg_log10_FDR = -log10(fdr_plot),
      DE = case_when(
        lfc_value > LFC_cutoff & fdr_value < FDR_cutoff ~ "UP",
        lfc_value < -LFC_cutoff & fdr_value < FDR_cutoff ~ "DOWN",
        TRUE ~ "NO"
      )
    )
  
  top_up = res_plot %>%
    filter(DE == "UP") %>%
    arrange(fdr_value, desc(abs(lfc_value))) %>%
    slice_head(n = n_labels)
  
  top_down = res_plot %>%
    filter(DE == "DOWN") %>%
    arrange(fdr_value, desc(abs(lfc_value))) %>%
    slice_head(n = n_labels)
  
  genes_to_label = bind_rows(top_up, top_down)
  
  ggplot(res_plot, aes(x = lfc_value, y = neg_log10_FDR)) +
    geom_point(
      data = res_plot %>% filter(DE == "NO"),
      color = "gray75",
      size = 0.7,
      alpha = 0.08
    ) +
    geom_point(
      data = res_plot %>% filter(DE == "DOWN"),
      color = "#377EB8",
      size = 1.8,
      alpha = 0.9
    ) +
    geom_point(
      data = res_plot %>% filter(DE == "UP"),
      color = "#E41A1C",
      size = 1.8,
      alpha = 0.9
    ) +
    geom_label_repel(
      data = genes_to_label,
      aes(label = gene_label, fill = DE),
      color = "white",
      size = 3.3,
      fontface = "bold",
      box.padding = 0.35,
      point.padding = 0.2,
      label.padding = unit(0.18, "lines"),
      label.r = unit(0.12, "lines"),
      label.size = 0,
      max.overlaps = Inf,
      min.segment.length = 0,
      segment.color = "gray35",
      segment.linewidth = 0.35,
      show.legend = FALSE
    ) +
    scale_fill_manual(
      values = c(
        "UP" = "#E41A1C",
        "DOWN" = "#377EB8"
      )
    ) +
    geom_vline(
      xintercept = c(-LFC_cutoff, LFC_cutoff),
      color = "black",
      linetype = "longdash"
    ) +
    geom_hline(
      yintercept = -log10(FDR_cutoff),
      color = "black",
      linetype = "longdash"
    ) +
    coord_cartesian(
      xlim = c(-volcano_LFC_limit, volcano_LFC_limit),
      ylim = c(0, volcano_FDR_limit),
      clip = "off"
    ) +
    labs(
      title = title,
      x = "log2 Fold Change",
      y = "-log10(FDR)"
    ) +
    theme_classic(base_size = 13, base_line_size = 1) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.margin = margin(10, 25, 10, 25)
    )
}

# ===================== DIRECTORIOS =====================

results_dir = "results/diffExp"
images_dir = "docs/images"

dir.create(images_dir, recursive = TRUE, showWarnings = FALSE)

# ===================== PARÁMETROS =====================

FDR_cutoff = 1e-4
LFC_cutoff = 1

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

# ===================== VOLCANO PLOTS =====================

p1 = make_volcano_plot(
  res = res_all,
  lfc_col = "log2FoldChange",
  fdr_col = "adj_pvalue",
  gene_col = "gene_name",
  LFC_cutoff = LFC_cutoff,
  FDR_cutoff = FDR_cutoff,
  title = "PS vs NS ajustado por instrumento",
  volcano_LFC_limit = 30,
  volcano_FDR_limit = 70,
  n_labels = 4
)

p2 = make_volcano_plot(
  res = res_hiseq2000,
  lfc_col = "log2FoldChange",
  fdr_col = "adj_pvalue",
  gene_col = "gene_name",
  LFC_cutoff = LFC_cutoff,
  FDR_cutoff = FDR_cutoff,
  title = "PS vs NS solo Illumina HiSeq 2000",
  volcano_LFC_limit = 20,
  volcano_FDR_limit = 70,
  n_labels = 4
)

# ===================== GUARDADO =====================

ggsave(
  filename = file.path(images_dir, "volcano_all_hiseq_adjusted.png"),
  plot = p1,
  width = 8,
  height = 6,
  dpi = 400
)

ggsave(
  filename = file.path(images_dir, "volcano_hiseq2000_only.png"),
  plot = p2,
  width = 8,
  height = 6,
  dpi = 400
)

rm(list = ls())
invisible(gc())