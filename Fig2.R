library(Seurat)
library(dplyr)
library(purrr)
sub_Olis <- qread("/home/data/caixing/AD/sub_Olisjw_analysis.qs")
sub_Olis$group <- paste(sub_Olis$slice, sub_Olis$genotype, sep = "_")
Idents(sub_Olis) <- "group"

slices <- c("TH","STR","ENT","HPF","RS","BF","PFC")

deg_list <- map(slices, function(sl) {
  FindMarkers(
    subset(sub_Olis, slice == sl),
    ident.1 = paste0(sl, "_5xFAD"),
    ident.2 = paste0(sl, "_WT"),
    test.use = "MAST",
    latent.vars = c("sample","month"),
    logfc.threshold = 0,
    min.pct = 0
  ) %>%
    tibble::rownames_to_column("gene") %>%
    mutate(slice = sl)
})

slice_deg <- bind_rows(deg_list)
write.csv(slice_deg, "result/slice_deg/all_slice_DEG.csv", row.names = FALSE)

deg_filt <- slice_deg %>% 
  filter(p_val_adj < 0.05)

up_genes   <- deg_filt %>% filter(avg_log2FC >  0.1)
down_genes <- deg_filt %>% filter(avg_log2FC < -0.1)

up_list   <- split(up_genes$gene,   up_genes$slice)
down_list <- split(down_genes$gene, down_genes$slice)

library(ggVennDiagram)

p_up <- ggVennDiagram(up_list) +
  scale_fill_gradient(low="white", high="#E31A1C") +
  theme_void()

p_down <- ggVennDiagram(down_list) +
  scale_fill_gradient(low="white", high="#1F78B4") +
  theme_void()

ggsave("result/slice_deg/fig2A_venn_up_down.pdf",
       ggpubr::ggarrange(p_up, p_down), width=8, height=4)

library(ComplexHeatmap)
library(circlize)
library(tidyr)

make_fc_mat <- function(df){
  df %>%
    select(gene, slice, avg_log2FC) %>%
    pivot_wider(names_from = slice, values_from = avg_log2FC, values_fill = 0) %>%
    column_to_rownames("gene") %>%
    as.matrix()
}

up_mat   <- make_fc_mat(up_genes)
down_mat <- make_fc_mat(down_genes)
fc_col <- colorRamp2(c(-0.5,0,0.5), c("blue","white","red"))

pdf("result/slice_deg/Fig2B_heatmap_up.pdf",5,8)
Heatmap(up_mat, name="logFC", col=fc_col, cluster_rows=FALSE, cluster_columns=FALSE)
dev.off()

pdf("result/slice_deg/Fig2B_heatmap_down.pdf",5,8)
Heatmap(down_mat, name="logFC", col=fc_col, cluster_rows=FALSE, cluster_columns=FALSE)
dev.off()


dat <- read_xlsx("AD/result/slice_deg/enrichGO_TH_filter.xlsx", sheet=2)

dat <- dat %>%
  mutate(logP = -log10(pvalue),
         group = Cluster)

top_go <- dat %>%
  group_by(group) %>%
  slice_max(logP, n=5)

ggplot(top_go, aes(logP, reorder(Description, logP), fill=group)) +
  geom_col() +
  scale_fill_manual(values=c(ADup="#E31A1C", WTup="#1F78B4")) +
  theme_classic() +
  xlab("-log10(P)") + ylab("")

ggsave("AD/result/slice_deg/Fig2C_TH_GO_barplot.pdf", width=5, height=4)

#####Fig2E####
knitr::opts_chunk$set(echo = TRUE, message = FALSE, warning = FALSE)

suppressPackageStartupMessages({
  library(qs)
  library(Seurat)
  library(data.table)
  library(tidyverse)
  library(ComplexHeatmap)
  library(circlize)
  library(SNFtool)
})

# =========================
# 0) Paths / Parameters
# =========================
base_dir <- "/home/data/caixing/AD"
code_dir <- file.path(base_dir, "code")
res_dir  <- file.path(base_dir, "result", "Traj_Oli")
csv_dir  <- file.path(base_dir, "results", "Traj_Oli")  # 你原代码用 results（注意区别）
dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(csv_dir, showWarnings = FALSE, recursive = TRUE)

themes_file <- file.path(code_dir, "figure_themes.R")
if (file.exists(themes_file)) source(themes_file)

months_order <- c("2M","3M","4M","5M","6M","8M","12M")
region_levels <- c("TH", "STR", "ENT", "HPF", "RS", "BF", "PFC")

deg_cutoff_padj <- 0.05
deg_cutoff_fc   <- 0.16  # 你后面用的 >0.16（约 log2(1.12)）

# =========================
# 1) Load Seurat object
# =========================
sub_Olis <- qread(file.path(base_dir, "sub_Olisjw_analysis.qs"))
DefaultAssay(sub_Olis) <- "RNA"

message("Genotype table:")
print(table(sub_Olis$genotype))
message("Month table:")
print(table(sub_Olis$month))

# 若你想分 WT/AD（你原脚本里 sub_Olis_WT/sub_Olis_AD 直接用了但没定义）
sub_Olis_WT <- subset(sub_Olis, subset = genotype == "WT")
sub_Olis_AD <- subset(sub_Olis, subset = genotype != "WT")  # 这里按你数据情况自行调整

# =========================
# 2) Helper: Sample-level average expression (sample_slice as identity)
# =========================
calc_sample_avg_expr <- function(obj, out_rdata, out_meta_rds, region_levels) {
  Idents(obj) <- paste(obj$sample, obj$slice, sep = "_")
  ave <- AverageExpression(obj, slot = "data", assays = "RNA")$RNA
  
  ave_dt <- as.data.frame(ave) %>%
    rownames_to_column("gene") %>%
    as.data.table()
  
  ave_melt <- melt(ave_dt, id.vars = "gene",
                   variable.name = "donor_region",
                   value.name = "ave.exp")
  
  # donor_region = sample_slice -> 拆分出 sample 和 slice
  ave_melt[, ids := str_split(donor_region, "_")]
  ave_melt[, Donor.ID := map_chr(ids, ~ .x[1])]
  ave_melt[, Unified_region := map_chr(ids, ~ .x[2])]
  ave_melt[, ids := NULL]
  
  ave_melt[, Region := factor(Unified_region, levels = region_levels)]
  
  save(ave_melt, file = out_rdata)
  saveRDS(obj@meta.data, file = out_meta_rds)
  
  invisible(ave_melt)
}

# WT sample-level avg
wt_rdata <- file.path(res_dir, "WT_OL_sample.level.average.expression_filtered.Rdata")
wt_meta  <- file.path(res_dir, "WT_metadata.rds")
wt_ave   <- calc_sample_avg_expr(sub_Olis_WT, wt_rdata, wt_meta, region_levels)

# AD sample-level avg（可选）
ad_rdata <- file.path(res_dir, "AD_OL_sample.level.average.expression_filtered.Rdata")
ad_meta  <- file.path(res_dir, "AD_metadata.rds")
ad_ave   <- calc_sample_avg_expr(sub_Olis_AD, ad_rdata, ad_meta, region_levels)

# =========================
# 3) Helper: Run adjacent-month DE (MAST, latent sample)
# =========================
run_adjacent_deg <- function(obj, months_order, out_dir, latent = "sample") {
  DefaultAssay(obj) <- "RNA"
  Idents(obj) <- "month"
  
  pairs <- tibble(
    ident.1 = months_order[-1],
    ident.2 = months_order[-length(months_order)]
  )
  
  for (i in seq_len(nrow(pairs))) {
    m1 <- pairs$ident.1[i]
    m0 <- pairs$ident.2[i]
    tag <- paste0("p", gsub("M","",m1), "vsp", gsub("M","",m0))
    
    message("Running FindMarkers: ", m1, " vs ", m0, " -> ", tag)
    
    deg <- FindMarkers(
      obj, ident.1 = m1, ident.2 = m0,
      min.pct = 0, logfc.threshold = 0, min.cells.feature = 0,
      test.use = "MAST", return.thresh = 1,
      latent.vars = latent
    )
    
    deg_dt <- as.data.table(deg, keep.rownames = "rn")
    fwrite(deg_dt, file.path(out_dir, paste0(tag, ".csv")))
  }
}

# 对全 sub_Olis 做月龄相邻 DE（你原逻辑）
run_adjacent_deg(sub_Olis, months_order, csv_dir, latent = "sample")

# =========================
# 4) Helper: compute path_group zscore per gene
# =========================
attach_path_group_and_zscore <- function(ave_dt, meta_rds, months_order) {
  ind_meta <- readRDS(meta_rds)
  ind_meta$Donor.ID     <- as.character(ind_meta$sample)
  ind_meta$Path..Group. <- ind_meta$month
  
  ind_meta <- ind_meta[, c("Donor.ID", "Path..Group.")] %>%
    unique()
  
  # data.table join: add path_group
  ave_dt <- copy(ave_dt)
  setDT(ind_meta)
  ave_dt[ind_meta, on = .(Donor.ID), path_group := i.Path..Group.]
  
  ave_dt[, zscore_gene := as.numeric(scale(ave.exp)), by = .(gene)]
  ave_path <- ave_dt[, .(ave_zscore_path = mean(zscore_gene, na.rm = TRUE)),
                     by = .(gene, path_group)]
  ave_path[, path_group := factor(path_group, levels = months_order)]
  
  ave_path
}

# =========================
# 5) Helper: Select significant genes from multiple adjacent DE files
# =========================
read_deg_and_sig_genes <- function(csv_dir, months_order, padj = 0.05, fc = 0.16) {
  # 生成文件名 p3vsp2 ... p12vsp8
  tags <- paste0(
    "p", gsub("M","",months_order[-1]),
    "vsp",
    gsub("M","",months_order[-length(months_order)])
  )
  
  files <- file.path(csv_dir, paste0(tags, ".csv"))
  names(files) <- tags
  
  get_sig <- function(dt) dt[p_val_adj < padj & avg_log2FC > fc, rn]
  
  sig <- Reduce(union, lapply(files, function(f) {
    dt <- fread(f)
    get_sig(dt)
  }))
  
  # 过滤掉噪音基因（和你原来一致）
  bad <- grep("^Rps|^Rpl|^mt-|^H2-|^Gm|^Hbb-|^Hba-|-", sig, ignore.case = TRUE)
  if (length(bad) > 0) sig <- sig[-bad]
  
  unique(sig)
}

# =========================
# 6) Helper: TT clustering (SNF spectral clustering)
# =========================
gt_TT_clustering <- function(dt, k = 6, seed = 9) {
  mtx <- dt[, .(gene, ave_zscore_path, path_group)] %>% unique()
  mtx <- dcast(mtx, gene ~ path_group, value.var = "ave_zscore_path")
  mtx <- mtx %>%
    as.data.frame() %>%
    column_to_rownames("gene") %>%
    as.matrix()
  
  set.seed(seed)
  diss <- dist(mtx)
  sim  <- 1 - as.matrix(diss) / max(diss)
  
  clust <- SNFtool::spectralClustering(sim, K = k)
  clustLab <- as.factor(clust)
  
  annot <- data.table(gene = rownames(mtx), cluster = clustLab)
  dt <- copy(dt)
  dt[annot, on = .(gene), cluster := i.cluster]
  
  ht <- ComplexHeatmap::Heatmap(
    mtx,
    show_row_names = FALSE,
    cluster_columns = FALSE,
    row_split = clustLab
  )
  
  list(dt = dt, clustLab = clustLab, mtx = mtx, p = ht)
}

# =========================
# 7) AD: build TT gene sets and plot trajectories
# =========================
# Load AD average expression dt
load(ad_rdata)  # loads ave_melt as object name used inside calc_sample_avg_expr; better: explicitly re-load:
ad_ave_dt <- qread(NULL) # placeholder (下面直接用 attach_path_group_and_zscore 输出)
ad_path <- attach_path_group_and_zscore(ad_ave, ad_meta, months_order)

# Get significant genes across adjacent stages
path_sig_genes <- read_deg_and_sig_genes(csv_dir, months_order, deg_cutoff_padj, deg_cutoff_fc)
message("n sig genes: ", length(path_sig_genes))

TT_cluster <- gt_TT_clustering(ad_path[gene %in% path_sig_genes, ], k = 6)
TT_cluster$p

TT_gene_clusters <- TT_cluster$dt[, .(gene, cluster)] %>% unique()
write_csv(TT_gene_clusters, file.path(res_dir, "cluster_geneset.csv"))
saveRDS(TT_gene_clusters, file.path(res_dir, "TT_cluster_gene_clusters.rds"))

TT_genetraj <- copy(TT_cluster$dt)
TT_genetraj <- TT_genetraj[!grepl("mt-", gene) & !grepl("^Rpl-", gene)]
TT_genetraj <- TT_genetraj[, .(gene, cluster, path_group, ave_zscore_path)] %>% unique()
TT_genetraj[, ave_zscore_cluster := mean(ave_zscore_path), by = .(cluster, path_group)]

# 名称（带每个cluster在2M时的基因数）
TT_N <- TT_genetraj[path_group == "2M", .N, by = .(cluster)][order(cluster)]$N
rename_dt <- data.table(
  cluster = factor(1:6),
  new_name = paste0("gene set #", 1:6, " (n = ", TT_N, ")")
)
TT_genetraj[rename_dt, on = .(cluster), new_name := i.new_name]
TT_genetraj[, path_group := factor(path_group, levels = months_order)]

# 画线图（需要 my_border_theme()，如果没有就用 theme_classic()）
base_theme <- if (exists("my_border_theme")) my_border_theme() else theme_classic()

p_ad <- ggplot(TT_genetraj, aes(x = path_group, y = ave_zscore_path)) +
  geom_line(color = "gray90", aes(group = gene)) +
  facet_wrap(new_name ~ ., ncol = 1) +
  geom_line(
    data = TT_genetraj[, .(path_group, cluster, ave_zscore_cluster, new_name)] %>% unique(),
    aes(y = ave_zscore_cluster, group = new_name, color = new_name),
    linewidth = 1
  ) +
  base_theme +
  labs(x = "Increasing Pathology (Month)", y = "Standardized gene expression") +
  theme(legend.position = "none", strip.text = element_text(size = 17)) +
  scale_color_manual(values = c("#e9a3c9", "#91bfdb", "#FBA949", "#8BD448", "#FAE442", "#9C4F96"))

ggsave(file.path(res_dir, "Fig2E_AD_traj.pdf"), p_ad, width = 3.5, height = 15)
qsave(TT_genetraj, file.path(res_dir, "TT_genetraj_month_AD.qs"))

# =========================
# 8) WT: apply AD clusters and plot trajectories
# =========================
wt_path <- attach_path_group_and_zscore(wt_ave, wt_meta, months_order)
dt_wt <- wt_path[gene %in% path_sig_genes, ]

# 用 AD 的 cluster 注释（复用 gene set）
annot <- TT_gene_clusters
setDT(annot)
dt_wt[annot, on = .(gene), cluster := i.cluster]
dt_wt <- dt_wt[, .(gene, cluster, path_group, ave_zscore_path)] %>% unique()
dt_wt[, ave_zscore_cluster := mean(ave_zscore_path), by = .(cluster, path_group)]

TT_N_wt <- dt_wt[path_group == "2M", .N, by = .(cluster)][order(cluster)]$N
rename_dt_wt <- data.table(
  cluster = factor(1:6),
  new_name = paste0("gene set #", 1:6, " (n = ", TT_N_wt[1:6], ")")
)
dt_wt[rename_dt_wt, on = .(cluster), new_name := i.new_name]
dt_wt[, path_group := factor(path_group, levels = months_order)]

p_wt <- ggplot(dt_wt, aes(x = path_group, y = ave_zscore_path)) +
  geom_line(color = "gray90", aes(group = gene)) +
  facet_wrap(new_name ~ ., ncol = 1) +
  geom_line(
    data = dt_wt[, .(path_group, cluster, ave_zscore_cluster, new_name)] %>% unique(),
    aes(y = ave_zscore_cluster, group = new_name, color = new_name),
    linewidth = 1
  ) +
  base_theme +
  labs(x = "Increasing Pathology (Month)", y = "Standardized gene expression") +
  theme(legend.position = "none", strip.text = element_text(size = 17)) +
  scale_color_manual(values = c("#e9a3c9", "#91bfdb", "#FBA949", "#8BD448", "#FAE442", "#9C4F96"))

ggsave(file.path(res_dir, "Fig2E_WT_traj.pdf"), p_wt, width = 3.5, height = 15)
qsave(dt_wt, file.path(res_dir, "TT_genetraj_month_WT.qs"))

# =========================
# 9) Wide format export (optional)
# =========================
TT_wide <- reshape(as.data.frame(TT_genetraj), idvar = "gene", timevar = "path_group", direction = "wide")
write_csv(TT_wide, file.path(res_dir, "TT_genetraj_wide_AD.csv"))


#####
celltype_color = c(
  "geneset1"="#66C2A5","geneset2"="#FC8D62","geneset3"="#8DA0CB",
  "geneset4"="#E78AC3",
  "geneset5"="#A6D854","geneset6"="#FFD92F"
)

# dat$logP <- -log10(dat$pvalue)

for (i in c('geneset1','geneset2','geneset3','geneset4','geneset5','geneset6')) {
  dat <- read_xlsx(paste0('AD/result/Traj_Oli/Geneset.xlsx'),
                   sheet = i)
  dat$group_type <- i
  dat$log10FDR <- as.numeric(dat$log10FDR)
  dat2 <- dat %>%
    dplyr::group_by(group_type) %>%
    dplyr::arrange((log10FDR))
  
  cmap <- c("viridis", "magma", "inferno", "plasma", "cividis", "rocket", "mako", "turbo")
  # dat2$group <- 1
  
  # dat2 <- dplyr::filter(dat2,dat2$group_type = i)
  # dat2 <- dat2[order(dat2$Gene_Number, decreasing = F), ]
  dat2
  dat2$Description<-factor(dat2$Description,levels = unique(dat2$Description),ordered = T)
  
  ggplot(data = dat2, aes(x = log10FDR, y = Description,fill=group_type)) +
    geom_bar(width = 0.7,stat = 'identity') +
    theme_classic() + 
    scale_x_continuous(expand = c(0,0)) +
    scale_fill_manual(values = alpha(celltype_color, 0.66)) +
    theme(axis.text.y = element_blank()) + 
    geom_text(data = dat2,
              aes(x = 0.1, y = Description, label = Description),
              size = 4.8,
              hjust = 0)+
    geom_vline(xintercept = -log10(0.05))
  
  ggsave(paste0('AD/Fig2E',i,'_enrichr_filter.pdf'),width = 5,height = 2.5)
  
}


