library(qs)
library(ggplot2)
library(Seurat)
library(SCP)

sce <- qread("/home/data/t030632/AD/AD/sub_all.qs")
df  <- sce@meta.data
colors <- list(
  
  clusters = c(
    Ast   = '#feb662',
    ExN   = "#70a4c8",
    InN   = '#9cd2ed',
    Micro = '#be7fbf',
    Oligo = '#86c7b4',
    OPC   = '#94c58f',
    Endo  = '#eb5b5e'
  ),
  
  slice = c(
    "#7986CB","#FFAB91","#64B5F6","#AED581",
    "#DCE775","#F48FB1","#9575CD"
  ),
  
  genotype = c(
    WT = '#c4daec',
    AD = '#e6bac5'
  ),
  
  month = c(
    '#efd2c9','#e0cfda','#ffee72','#e6e2a3','#cbdaa9',
    '#64a776','#2d3462','#e6b884','#bc9a7f','#c49abc',
    '#927c9a','#3674a2','#9f8d89','#63a3b8','#c4daec',
    '#61bada','#e29eaf','#4490c4','#de8b36','#c4612f',
    '#76a2be','#408444','#c6adb0','#9d3b62','#2d3462'
  )
)

theme_umap <- function() {
  theme_classic() +
    theme(
      plot.background = element_blank(),
      panel.grid = element_blank(),
      axis.title = element_blank(),
      axis.text  = element_blank(),
      axis.ticks = element_blank(),
      axis.line = element_line(colour = "black", size = 0.3,
                               arrow = arrow(length = unit(0.1,"cm"))),
      legend.position = "none",
      aspect.ratio = 1
    )
}
plot_umap <- function(df, group, cols) {
  ggplot(df, aes(UMAP_1, UMAP_2, color = .data[[group]])) +
    geom_point(size = 0.1, stroke = 0) +
    scale_color_manual(values = cols) +
    scale_x_continuous(limits = range(df$UMAP_1)) +
    scale_y_continuous(limits = range(df$UMAP_2)) +
    theme_umap()
}
plot_legend <- function(df, group, cols, ncol = 2) {
  ggplot(df, aes(UMAP_1, UMAP_2, color = .data[[group]])) +
    geom_point(size = 0) +
    theme_void() +
    scale_color_manual(values = cols, name = '') +
    guides(color = guide_legend(ncol = ncol,
                                override.aes = list(size = 2))) +
    theme(
      legend.text = element_text(size = 5),
      legend.position = 'bottom',
      legend.box.spacing = unit(0, 'cm')
    )
}
p <- plot_umap(df, "clusters", colors$clusters)
ggsave("AD/Fig1/Fig1_CellDimPlot.pdf", p, width = 5, height = 4)

p_leg <- plot_legend(df, "clusters", colors$clusters)
ggsave("AD/Fig1/Fig1_CellDimPlot_legend.pdf", p_leg, width = 3.5, height = 6)

p <- plot_umap(df, "slice", colors$slice)
ggsave("AD/Fig1/Fig1_CellDimPlot_slice.pdf", p, width = 5, height = 4)

p_leg <- plot_legend(df, "slice", colors$slice)
ggsave("AD/Fig1/Fig1_CellDimPlot_slice_legend.pdf", p_leg, width = 3.5, height = 6)

p <- plot_umap(df, "month", colors$month)
ggsave("AD/Fig1/Fig1_CellDimPlot_month.pdf", p, width = 5, height = 4)

p_leg <- plot_legend(df, "month", colors$month, ncol = 3)
ggsave("AD/Fig1/Fig1_CellDimPlot_month_legend.pdf", p_leg, width = 3.5, height = 6)

p <- plot_umap(df, "genotype", colors$genotype)
ggsave("AD/Fig1/Fig1_CellDimPlot_genotype.pdf", p, width = 5, height = 4)

p_leg <- plot_legend(df, "genotype", colors$genotype)
ggsave("AD/Fig1/Fig1_CellDimPlot_genotype_legend.pdf", p_leg, width = 3.5, height = 6)

#####Fig1C#####
sce_ref <- qread("sub_all.qs")
gene_list <- c(
  "Camk2a","Snap25","Slc17a7","Slc30a3",
  "Gad1","Gad2","Sst","Npy","Dpp6",
  "Aqp4","Slc1a2","Gfap","Gja1",
  "Plp1","Mog","Mbp","Cldn11",
  "Olig1","Olig2","Vcan","Pdgfra",
  "Hexb","P2ry12","C1qa","Cx3cr1",
  "Cldn5","Rgs5","Ebf1","Flt1"
)
cluster_cols <- c(
  'Ast'   = '#feb662',
  'ExN'   = "#70a4c8",
  'InN'   = '#9cd2ed',
  'Micro' = '#be7fbf',
  'OL'    = '#86c7b4',
  'OPC'   = '#94c58f',
  'Str'   = '#eb5b5e'
)
mean_gene_exp <- AverageExpression(
  sce_ref,
  features = unique(gene_list),
  group.by = "clusters",
  assays   = "RNA",
  slot     = "data"
)$RNA
cluster_order <- c("ExN","InN","Ast","OL","OPC","Micro","Str")

mean_gene_exp <- mean_gene_exp[, cluster_order]
htdf <- t(scale(t(mean_gene_exp)))
library(ComplexHeatmap)
library(circlize)

col_fun <- colorRamp2(c(-2, 0, 2), c("#0099CC", "white", "#CC0033"))
column_ha <- HeatmapAnnotation(
  Cluster = cluster_order,
  col = list(Cluster = cluster_cols)
)
pdf("Fig1/Fig1C_top5marker_Heatmap.pdf", width = 6, height = 8)

Heatmap(
  htdf,
  name = "Z-score",
  col  = col_fun,
  
  cluster_rows    = FALSE,
  cluster_columns = FALSE,
  
  column_order = cluster_order,
  show_column_names = FALSE,
  
  row_title    = "Cluster marker genes",
  column_title = "Clusters",
  
  row_names_gp   = gpar(fontface = "italic", fontsize = 10),
  row_names_side = "left",
  
  rect_gp = gpar(col = "white", lwd = 0.5),
  
  top_annotation = column_ha,
  
  width  = ncol(htdf) * unit(0.25, "inch"),
  height = nrow(htdf) * unit(0.16, "inch")
)

dev.off()

####Fig1F-G####
library(SCP)
library(dplyr)
library(ggplot2)
library(tidyr)

celltype_cols <- c(
  Ast   = '#feb662',
  ExN   = "#70a4c8",
  InN   = '#9cd2ed',
  Micro = '#be7fbf',
  OL    = '#86c7b4',
  OPC   = '#94c58f',
  Str   = '#eb5b5e'
)


sample_table <- data.frame(
  sample = c(
    "HZ12M01","HZ12M02","HZ12M03","HZ12M04","HZ12M05","HZ12M06",
    "HZ2M01","HZ2M02","HZ2M03","HZ2M04","HZ2M05","HZ2M06",
    "HZ3M01","HZ3M02","HZ3M03","HZ3M04","HZ3M05","HZ3M06",
    "HZ4M03","HZ4M05","HZ4MP101","HZ4MP102","HZ4MP202",
    "HZ5M01","HZ5M02","HZ5M03","HZ5M04","HZ5M05","HZ5M06",
    "HZ6M02","HZ6M03","HZ6M04","HZ6M05","HZ6M06",
    "HZ6MP101","HZ6MP102","HZ6MP202",
    "HZ8M01","HZ8M02","HZ8M03"
  ),
  plaque_id = c(
    "HZ12M01","HZ12M02","HZ12M03","HZ12M04","HZ12M05","HZ12M06",
    "HZ2M01","HZ2M02","HZ2M03","HZ2M04","HZ2M05","HZ2M06",
    "HZ3M01","HZ3M02","HZ3M03","HZ3M04","HZ3M05","HZ3M06",
    "HZ4M03","HZ4M05","HZ4MP101","HZ4MP102","HZ4MP202",
    "HZ5M01","HZ5M02","HZ5M03","HZ5M04","HZ5M05","HZ5M06",
    "HZ6M02","HZ6M03","HZ6M04","HZ6M05","HZ6M06",
    "HZ6MP101","HZ6MP102","HZ6MP202",
    "HZ8M01","HZ8M02","HZ8M03"
  )
)


for (i in seq_len(nrow(sample_table))) {
  
  sample    <- sample_table$sample[i]
  plaque_id <- sample_table$plaque_id[i]
  
  message("Processing ", sample)
  
  sce <- readRDS(paste0("/home/data/t030632/AD/result/ST_result/ST/", sample, "_ST_umap.rds"))
  df.res <- read.csv(paste0("/home/data/t030632/AD/ST/all/", sample, "_obs_data.csv"), row.names = 1)
  
  sce$celltype <- factor(df.res$celltype,
                         levels = c("Ast","ExN","InN","Micro","OL","OPC","Str"))
  
  ## -------- Fig1F : Celltype spatial --------
  p1 <- CellDimPlot(
    sce, group.by = "celltype", reduction = "spatial",
    raster = TRUE, pt.size = 0.12, pt.alpha = 0.6,
    palcolor = celltype_cols, theme_use = "theme_blank"
  )
  
  ggsave(paste0("AD/Fig1/Fig1F_", sample, "_celltype.pdf"), p1, width = 6, height = 5)
  
  
  
  ## -------- Plaque data --------
  plaque <- read.csv(
    paste0("/home/data/t030632/AD/ST/plaque/", plaque_id,
           "_PlaqueCircle/", plaque_id, "_PlaqueSum.csv"),
    row.names = 1
  )
  
  plaque$plaque <- "plaque"
  
  
  ## -------- merge plaque with ST --------
  meta <- sce@meta.data %>%
    tibble::rownames_to_column("cellid") %>%
    left_join(plaque, by = c("row", "col")) %>%
    mutate(plaque = ifelse(is.na(plaque), "NO", "YES"))
  
  sce <- AddMetaData(sce, meta)
  
  
  ## -------- Fig1G : plaque overlay --------
  p2 <- CellDimPlot(
    sce, group.by = "plaque", reduction = "spatial",
    raster = TRUE, pt.size = 0.12, pt.alpha = 0.4,
    palcolor = c(NO = "gray85", YES = "#E31A1C"),
    theme_use = "theme_blank"
  ) +
    stat_density_2d(
      data = plaque,
      aes(x = row, y = col, fill = ..level..),
      geom = "polygon", alpha = 0.3
    ) +
    scale_fill_gradientn(colors = c("white","pink","#E31A1C")) +
    geom_point(
      data = plaque,
      aes(x = row, y = col),
      size = 0.4, color = "#E31A1C", alpha = 0.4
    )
  
  ggsave(paste0("AD/Fig1/Fig1G_", sample, "_plaque.pdf"), p2, width = 6, height = 5)
  
}
