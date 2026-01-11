library(tidyverse)
library(patchwork)
library(ComplexHeatmap)
library(circlize)
library(vegan)
library(openxlsx)
library(ggplot2)
library(qs)
library(scCustomize)

sub_Olis <- qread(file = "Olisjw.qs")
####Fig3A-B####

CellDimPlot(
  srt = sub_Olis, group.by = c("genotype"),raster = F,pt.size = 0.01,pt.alpha = 0.1,
  palcolor = c( "#242467","#f20d0d"),
  reduction = "HarmonyUMAP2D", theme_use = "theme_blank"
)
ggsave("result/sub_OL_genetype.pdf",width = 6,height = 5)
CellDimPlot(
  srt = sub_Olis, group.by = c("month"),raster = F,
  palcolor = c("#fab37f","#e98741","#967568","#794976","#82c785","#edeaa4","#cdaa9f"),
  reduction = "HarmonyUMAP2D", theme_use = "theme_blank"
)
ggsave("result/sub_OL_month.pdf",width = 6,height = 5)


CellDimPlot(
  srt = sub_Olis, group.by = c("slice"),raster = F,
  palcolor = c("#dc8e97","#e3d1db","#74a893","#ac9141","#5ac6e9","#ebce8e","#e5c06e","#7587b1"),
  reduction = "HarmonyUMAP2D", theme_use = "theme_blank"
)
ggsave("result/sub_OL_slice.pdf",width = 6,height = 5)

sub_Olis$subclusters <- factor(sub_Olis$subclusters,)
CellDimPlot(
  srt = sub_Olis, group.by = c("subclusters"),raster = F,
  label = F,
  reduction = "HarmonyUMAP2D", theme_use = "theme_blank"
)
ggsave("result/sub_OL_subclusters.pdf",width = 6,height = 5)

maker <- c(
  "Il33","Ptgds","Qdpr","Apoe","Cryab",
  "Serpina3n","H2-K1","C4b","H2-D1","Cd74",
  "Npas3","Frmd5","Pex5l","Pakap","Pde4b",
  "S100b","Selenop","Klk6","Enpp6","Anln",
  "Ppp1r1b","Drd1","Meis2","Gpr88","Penk",
  "Ctnna3","Robo1","Prr5l","Ctnna2","Klk6",
  "Man1a","Synpr","Prickle1","Ablim2","Plekha1",
  "Nfasc","Man1a","Ctps","Prom1","Fam214a"
  
)

mean_gene_exp<-AverageExpression(sub_Olis,
                                 features=unique(maker),
                                 group.by='subclusters',assays = 'RNA',
                                 slot='data')%>%
  data.frame()%>%
  as.matrix()

#addcolnames
colnames(mean_gene_exp)<-names(B_color)
names(B_color)
#Z-score
htdf<-t(scale(t(mean_gene_exp),scale=T,center=T))

#color
col_fun=colorRamp2(c(-2,0,2),c("#0099CC","white","#CC0033"))

#topannotation
column_ha=HeatmapAnnotation(cluster=colnames(htdf),
                            col=list(cluster=B_color))
pdf("AD/result/Fig2B_top5maker_Heatmap.pdf",width = 7,height = 8)
Heatmap(htdf,
        name="Z-score",
        cluster_columns=F,cluster_rows=F,show_column_names = F,
        row_title="Clustertop5Markergenes",
        column_title="Clusters",
        row_names_gp=gpar(fontface='italic',fontsize=10),
        row_names_side='left',
        border=T,
        rect_gp=gpar(col="white",lwd=1),
        column_names_side='top',
        column_names_rot=45,
        top_annotation=column_ha,
        # column_names_gp = grid::gpar(fontsize = 7),
        # row_names_gp = grid::gpar(fontsize = 7),
        width = ncol(htdf) * unit(0.2, "inch"),
        height = nrow(htdf) * unit(0.12, "inch"),
        #column_split=paste('clsuter',0:8,sep=''),
        col=col_fun)
dev.off()

####Fig3C_miloR####
library(stringr)
library(Seurat)
library(miloR)
library(scater)
library(scran)
library(dplyr)
library(patchwork)
library(qs)
library(BiocParallel)
register(MulticoreParam(workers=10,progressbar=TRUE))

scRNA <- qread("sub_Olisjw_analysis_pro.qs")
scRNA_pre<-as.SingleCellExperiment(scRNA)

scRNA_pre<-Milo(scRNA_pre)
scRNA_pre

reducedDims(scRNA_pre)

scRNA_pre<-buildGraph(scRNA_pre,k=30,d=15,reduced.dim="HARMONYPCA")

scRNA_pre<-makeNhoods(scRNA_pre,prop=0.05,k=30,d=15,
                      refined=TRUE,reduced_dims="HARMONYPCA")
plotNhoodSizeHist(scRNA_pre)

# colData(scRNA_pre)
scRNA_pre<-countCells(scRNA_pre,
                      meta.data=as.data.frame(colData(scRNA_pre)),
                      sample="orig.ident")

head(nhoodCounts(scRNA_pre))

scRNA_design<-data.frame(colData(scRNA_pre))[,c("orig.ident","genotype")]


scRNA_design<-distinct(scRNA_design)
rownames(scRNA_design)<-scRNA_design$orig.ident

scRNA_design


scRNA_pre<-calcNhoodDistance(scRNA_pre,d=15,reduced.dim="HARMONYPCA")


da_results<-testNhoods(scRNA_pre,design=~genotype,reduced.dim = 'HARMONYPCA',design.df=scRNA_design)

head(da_results)

da_results%>%
  arrange(SpatialFDR)%>%
  head()

ggplot(da_results,aes(PValue))+geom_histogram(bins=50)

ggplot(da_results,aes(logFC,-log10(SpatialFDR)))+
  geom_point()+
  geom_hline(yintercept=1)##Marksignificancethreshold(10%FDR)


scRNA_pre<-buildNhoodGraph(scRNA_pre)

##Plotsingle-cellUMAP
umap_pl<-plotReducedDim(scRNA_pre,dimred="HARMONYUMAP2D",
                        colour_by="genotype",text_by="subclusters",
                        text_size=3,point_size=0.5)+guides(fill="none")
umap_pl

##Plotneighbourhoodgraph
nh_graph_pl<-plotNhoodGraphDA(scRNA_pre,da_results,
                              layout="HARMONYUMAP2D",alpha=0.9)

umap_pl+nh_graph_pl+
  plot_layout(guides="collect")

da_results<-annotateNhoods(scRNA_pre,
                           da_results,
                           coldata_col="subclusters")
head(da_results)

ggplot(da_results,aes(subclusters_fraction))+geom_histogram(bins=50)

str(da_results)
table(da_results$subclusters)
range(da_results$SpatialFDR)

da_results$subclusters <- factor(da_results$subclusters)
ggsave('result/Fig3C_miloR_Group.pdf',width = 5,height = 4)  

####Fig3E####
suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(ggplot2)
  library(cowplot)
  library(clusterProfiler)
})

# =========================
# 0) Parameters
# =========================
geneset_xlsx <- "/home/data/t030632/result/18PCDgenes.xlsx"
markers_rds  <- "/home/data/t030632/AD/AD/result/OL_function/Markers.rds"

out_dir_csv  <- "/home/data/t030632/AD/AD/18PCDgenes_res_csv"
out_dir_rds  <- "/home/data/t030632/AD/AD/18PCDgenes_res_rds"
out_fig_pdf  <- "AD/8_PCD_allclusters.pdf"

dir.create(out_dir_csv, showWarnings = FALSE, recursive = TRUE)
dir.create(out_dir_rds, showWarnings = FALSE, recursive = TRUE)

# 是否强制基因名格式（默认 FALSE：更安全，避免把 HLA-A 变成 Hla-a）
fix_gene_case <- FALSE

# 你最终想展示的死亡通路
focus_terms <- c(
  "Anoikis","Apoptosis","Autophagy","Ferroptosis",
  "Lysosome-dependentcelldeath","Necroptosis"
)

# Cluster颜色
B_color <- c(
  "C1_OL_Il33"="#A6CEE3","C2_OL_MHC" = "#1F78B4","C3_OL_Pakap" = "#B2DF8A",
  "C4_OL_Enpp6"="#33A02C","C5_OL_Penk"="#FDBF6F",
  "C6_OL_Prr5l"="#FF7F00","C7_OL_Synpr"="#FB9A99","C8_OL_Ctps"="#E31A1C"
)

# =========================
# 1) Read gene sets (TERM2GENE)
# =========================
geneset <- readxl::read_xlsx(geneset_xlsx)

term2gene <- as.data.frame(geneset) %>%
  pivot_longer(
    cols = everything(),
    names_to = "term",
    values_to = "gene"
  ) %>%
  drop_na() %>%
  mutate(
    term = as.character(term),
    gene = str_trim(as.character(gene))
  )

if (fix_gene_case) {
  # 注意：不建议对人基因这么做（HLA-A会被破坏），所以默认关掉
  term2gene$gene <- str_to_title(str_to_lower(term2gene$gene))
}

# 去重（TERM2GENE 最好不要重复）
term2gene <- distinct(term2gene, term, gene)

# =========================
# 2) Read markers and prepare cluster list
# =========================
Markers <- readRDS(markers_rds)
stopifnot(all(c("cluster","gene","avg_log2FC") %in% colnames(Markers)))

clusters <- sort(unique(Markers$cluster))

# GSEA 需要 named numeric vector
build_geneList <- function(df_cluster) {
  gl <- df_cluster$avg_log2FC
  names(gl) <- df_cluster$gene
  gl <- sort(gl, decreasing = TRUE)
  
  # 去掉 NA / 重复名字（重复会影响 GSEA）
  gl <- gl[!is.na(gl)]
  gl <- gl[!duplicated(names(gl))]
  gl
}

# =========================
# 3) Run GSEA per cluster and save results
# =========================
set.seed(111)

run_one_cluster <- function(cl) {
  message("Running cluster: ", cl)
  
  df_cl <- Markers %>% filter(cluster == cl)
  geneList <- build_geneList(df_cl)
  
  if (length(geneList) < 30) {
    warning("Cluster ", cl, " has too few genes for stable GSEA: ", length(geneList))
  }
  
  gsea_res <- GSEA(
    geneList,
    TERM2GENE = term2gene,
    pvalueCutoff = 1
  )
  
  res_df <- as.data.frame(gsea_res@result)
  
  # 输出
  write.csv(res_df, file = file.path(out_dir_csv, paste0("18PCDgenes_res_", cl, ".csv")), row.names = FALSE)
  saveRDS(gsea_res, file = file.path(out_dir_rds, paste0("18PCDgenes_res_", cl, ".rds")))
  
  res_df$cluster <- cl
  res_df
}

all_res <- map_dfr(clusters, run_one_cluster)

# =========================
# 4) Harmonize term column name and filter pathways
# =========================
# clusterProfiler 的 GSEA result 通常有 ID / Description 两列
term_col <- if ("ID" %in% colnames(all_res)) "ID" else if ("Description" %in% colnames(all_res)) "Description" else NA
stopifnot(!is.na(term_col))

plot_df <- all_res %>%
  mutate(term = .data[[term_col]]) %>%
  filter(term %in% focus_terms)

# 保证 cluster 顺序、颜色可用
plot_df$cluster <- factor(plot_df$cluster, levels = intersect(names(B_color), unique(plot_df$cluster)))
plot_df$term <- factor(plot_df$term, levels = focus_terms)

# =========================
# 5) Plot: per term, barplot of NES across clusters
# =========================
make_one_term_plot <- function(one_term) {
  df <- plot_df %>% filter(term == one_term)
  
  ggplot(df, aes(x = cluster, y = NES, fill = cluster)) +
    geom_col(width = 0.6, linewidth = 0.3) +
    geom_hline(yintercept = 0, linewidth = 0.8) +
    scale_fill_manual(values = B_color, drop = FALSE) +
    scale_y_continuous(
      breaks = seq(-2.5, 2.5, 0.5),
      limits = c(-2.8, 2.8),
      expand = c(0, 0)
    ) +
    labs(x = NULL, y = "NES", title = one_term) +
    theme_classic(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.8),
      axis.text.x = element_blank(),
      axis.text.y = element_text(size = 12, colour = "black"),
      axis.title.y = element_text(size = 16),
      axis.ticks.length = unit(0.2, "cm"),
      axis.ticks = element_line(color = "black", linewidth = 1),
      legend.position = "none",
      plot.title = element_text(size = 14, face = "bold")
    )
}

plots <- lapply(levels(plot_df$term), make_one_term_plot)

final_plot <- plot_grid(plotlist = plots, ncol = 3, align = "h")
ggsave(out_fig_pdf, plot = final_plot, width = 16, height = 8)

####Fig3I-K####
suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(qs)
  library(dplyr)
  library(readr)
  library(tibble)
  library(SingleCellExperiment)
  library(scater)
  library(scran)
  library(RColorBrewer)
  library(ComplexHeatmap)
  library(grid)
})

# =========================
# 0) Parameters / files
# =========================
human_qs   <- "AD/rawdata/human_sub_Oli.qs"
mouse_qs   <- "AD/sub_Olisjw_analysis.qs"
orth_csv   <- "AD/rawdata/mouse_human_marmoset_macaque_orthologs_20231113.csv"

out_figK   <- "FigK.pdf"
out_figMN  <- "Fig3J_MetaNeighborUS.pdf"

assay_use  <- "RNA"
hvg_n      <- 3000

# human/mouse subcluster labels used for heatmap
human_levels <- c("human|HOL1","human|HOL2","human|HOL3","human|HOL4")
mouse_levels <- c(
  "mouse|C1_OL_Il33","mouse|C2_OL_MHC","mouse|C3_OL_Pakap","mouse|C4_OL_Enpp6",
  "mouse|C5_OL_Penk","mouse|C6_OL_Prr5l","mouse|C7_OL_Synpr","mouse|C8_OL_Ctps"
)

# =========================
# 1) Human: module scores + plotting
# =========================
human_OL <- qread(human_qs)

genes_list <- list(
  MHC_I  = c("HLA-A","HLA-B","HLA-C","B2M"),
  MHC_II = c("CD74","HLA-DMA","HLA-DPA1","HLA-DPB1","HLA-DRA")
)

human_OL <- AddModuleScore(human_OL, features = genes_list, name = names(genes_list))
# AddModuleScore 会生成 MHC_I1 / MHC_II1（一般如此），这里做个更直观的别名
score_cols <- grep("^MHC_I|^MHC_II", colnames(human_OL@meta.data), value = TRUE)
# 保险：取每个 name 的第一个 score 列
mhcI_col  <- grep("^MHC_I", score_cols, value = TRUE)[1]
mhcII_col <- grep("^MHC_II", score_cols, value = TRUE)[1]

# 你原来画 Cluster2（不确定是哪个列），这里示例画 MHC_I score
pK <- FeaturePlot(
  human_OL,
  features = mhcI_col,
  pt.size = 0.8,
  raster = FALSE
) + theme_void()

ggsave(out_figK, plot = pK, width = 5, height = 4)

# =========================
# 2) Mouse: map mouse genes -> human symbols (counts-level)
# =========================
sce_mouse <- qread(mouse_qs)
DefaultAssay(sce_mouse) <- assay_use

orth <- read.csv(orth_csv, stringsAsFactors = FALSE)
stopifnot(all(c("mouse_Symbol","human_Symbol") %in% colnames(orth)))

counts <- GetAssayData(sce_mouse, assay = assay_use, slot = "counts")
mouse_genes <- rownames(counts)

my_map <- orth %>%
  transmute(mouse_Symbol = as.character(mouse_Symbol),
            human_Symbol = as.character(human_Symbol)) %>%
  filter(mouse_Symbol %in% mouse_genes) %>%
  filter(!is.na(mouse_Symbol), mouse_Symbol != "") %>%
  filter(!is.na(human_Symbol), human_Symbol != "") %>%
  distinct(mouse_Symbol, .keep_all = TRUE)

idx <- match(mouse_genes, my_map$mouse_Symbol)
mapped <- my_map$human_Symbol[idx]

new_names <- mapped
miss <- is.na(new_names) | new_names == ""
new_names[miss] <- mouse_genes[miss]  # unmapped: keep original mouse symbol

# 合并重复 human symbol（稀疏矩阵求和）
grp <- factor(new_names, levels = unique(new_names))
G <- sparse.model.matrix(~ grp - 1)
colnames(G) <- levels(grp)
rownames(G) <- mouse_genes

counts_humanized <- t(G) %*% counts
stopifnot(identical(rownames(counts_humanized), levels(grp)))

# 替换 assay（只替换 counts；后续会重新 Normalize/Scale）
sce_mouse[[assay_use]] <- CreateAssayObject(counts = counts_humanized)

rm(counts, counts_humanized, G)

sce_mouse$batch <- "mouse"
human_OL$batch  <- "human"

# =========================
# 3) Merge + preprocessing (must assign!)
# =========================
# 确保 subclusters 字段存在
stopifnot("subclusters" %in% colnames(sce_mouse@meta.data))
stopifnot("subclusters" %in% colnames(human_OL@meta.data))

cca.results <- merge(sce_mouse, human_OL)

DefaultAssay(cca.results) <- assay_use
cca.results <- cca.results %>%
  NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(verbose = FALSE) %>%
  ScaleData(verbose = FALSE)

# =========================
# 4) MetaNeighborUS
# =========================
cca.sce <- as.SingleCellExperiment(cca.results)

# HVG (MetaNeighborUS 推荐使用跨 study 的 HVGs)
global_hvgs1 <- variableGenes(dat = cca.sce, exp_labels = cca.sce$batch)
global_hvgs  <- global_hvgs1[seq_len(min(hvg_n, length(global_hvgs1)))]

Aurocs_matrix <- MetaNeighborUS(
  var_genes   = global_hvgs,
  dat         = cca.sce,
  study_id    = cca.sce$batch,
  cell_type   = cca.sce$subclusters,
  fast_version = TRUE
)

# =========================
# 5) Heatmap export
# =========================
# 取子矩阵（如果某些 level 不存在会报错，这里先做交集）
rows_use <- intersect(human_levels, rownames(Aurocs_matrix))
cols_use <- intersect(mouse_levels, colnames(Aurocs_matrix))

mat_plot <- t(Aurocs_matrix[rows_use, cols_use, drop = FALSE])

pdf(out_figMN, width = 4, height = 3.5)
Heatmap(
  mat_plot,
  name = "Similarity",
  col = colorRampPalette(brewer.pal(9, "Blues"))(100),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_gp = gpar(fontsize = 10, fontface = "bold", col = "black"),
  column_names_gp = gpar(fontsize = 10, fontface = "bold", col = "black"),
  heatmap_legend_param = list(title = "Similarity", at = seq(0, 0.6, by = 0.1)),
  border = TRUE,
  column_title = ""
)
dev.off()



