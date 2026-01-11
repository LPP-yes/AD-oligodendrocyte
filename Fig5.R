library(SCP)
sub_Olis1 <- qread("/home/data/caixing/AD/sub_Olisjw_analysis.qs")

# sub_Olis1 <- NULL
sub_Olis1 <- RunSlingshot(srt = sub_Olis1, group.by = "subclusters", 
                          start = c('C8_OL_Ctps','C7_OL_Synpr'),#end = c('C2_OL_Serpina3n','C3_OL_Pakap','C4_OL_Plin3','C5_OL_Penk'
                          # ),
                          reduction = "Harmonypca")
colnames(sub_Olis1@meta.data)
table(sub_Olis1$subclusters)
Traj1 <- subset(sub_Olis2,Lineage == 'Lineage1')

qsave(Traj1,'result/Slingshot/Traj1.qs')

FeatureDimPlot(Traj1,features = 'Monocle3_Pseudotime', 
               reduction = "umap", raster = F,
               pt.size = 0.01,pt.alpha =  0.5,
               theme_use = "theme_blank")
ggsave('result/Fig5B_Traj1_monocle3_pseudotime.pdf',width = 5,height = 6)

CellDimPlot(Traj1,reduction = 'umap',group.by = "subclusters",label = F,
            pt.size = 0.0000000000000001,pt.alpha = 0.5,
            theme_use = "theme_blank")
ggsave('result/Fig5B_Traj4_subclusters.pdf',width = 5,height = 6)


#####pyscenic####
library(Seurat)
library(qs)
library(tidyverse)
library(SCENIC)
library(ClusterGVis)

sce <- qread("/home/data/caixing/AD/result/TrajectoryAnalysis/sub_Olis_TrajectoryAnalysis.qs")
loom <- open_loom("/home/data/caixing/AD/result/pyscenic/sample_SCENIC.loom")

regulon_incidence <- get_regulons(loom, column.attr.name = "Regulons")
regulon_list      <- regulonsToGeneLists(regulon_incidence)

regulonAUC <- get_regulons_AUC(loom, column.attr.name = "RegulonsAUC")
regulon_thresholds <- get_regulon_thresholds(loom)
embeddings <- get_embeddings(loom)

close_loom(loom)
sub_regulonAUC <- regulonAUC[, match(colnames(sce), colnames(regulonAUC))]

stopifnot(identical(colnames(sub_regulonAUC), colnames(sce)))

cell_meta <- data.frame(
  celltype = sce$subclusters,
  row.names = colnames(sce)
)

save(sub_regulonAUC, cell_meta, sce,
     file = "/home/data/caixing/AD/result/TrajectoryAnalysis/for_rss_and_visual.Rdata")

sub_regulonAUC <- sub_regulonAUC[
  onlyNonDuplicatedExtended(rownames(sub_regulonAUC)), 
]

auc_matrix <- getAUC(sub_regulonAUC)
qsave(auc_matrix, "result/pyscenic_auc.qs")
auc_mat <- qread("result/pyscenic_auc.qs")
traj <- qread("result/Traj4.qs")

meta <- traj@meta.data %>%
  select(Monocle3_Pseudotime, Lineage1, subclusters) %>%
  arrange(Monocle3_Pseudotime)

common_cells <- intersect(rownames(meta), colnames(auc_mat))

auc_mat <- auc_mat[, common_cells]
meta    <- meta[common_cells, ]
meta$Pseudotime_group <- cut(
  meta$Monocle3_Pseudotime,
  breaks = quantile(meta$Monocle3_Pseudotime, probs = seq(0,1,0.01)),
  include.lowest = TRUE,
  labels = paste0("Pse_", 1:100)
)
auc_long <- as.data.frame(t(auc_mat))
auc_long$Pseudotime_group <- meta$Pseudotime_group

tf_long <- auc_long %>%
  pivot_longer(-Pseudotime_group,
               names_to="TF",
               values_to="AUC")

tf_pt_mat <- tf_long %>%
  group_by(TF, Pseudotime_group) %>%
  summarise(meanAUC = mean(AUC), .groups="drop") %>%
  pivot_wider(names_from = Pseudotime_group,
              values_from = meanAUC)

tf_mat <- as.data.frame(tf_pt_mat)
rownames(tf_mat) <- tf_mat$TF
tf_mat <- tf_mat[,-1]

qsave(tf_mat, "result/T4_TF_score.qs")
zscore_clip <- function(mat, clip=2){
  z <- t(scale(t(mat)))
  z[z > clip] <- clip
  z[z < -clip] <- -clip
  return(z)
}

tf_z <- zscore_clip(tf_mat, 2)
ck_all <- clusterData(tf_z, cluster.method="kmeans", cluster.num=5)

pdf('result/Fig5D_T4_pysenic_fiterall.pdf',height=5,width=12,onefile=F)
visCluster(
  ck_all,
  plot.type="heatmap",
  show_row_dend=FALSE,
  ht.col = list(col_range=c(-2,0,2),
                col_color=c("#66C2A5","#EEEEEE","#B2182B"))
)
dev.off()
#####scFEA#####
sce <- qread("/home/data/caixing/AD/result/scMetabolism/sub_Olis_countexp.Seurat.qs")
metabolism_mat <- sce@assays[["METABOLISM"]][["score"]]
up_pathways <- c(
  'Citrate cycle (TCA cycle)','Glycolysis / Gluconeogenesis',
  'Cysteine and methionine metabolism','Glutathione metabolism',
  'Oxidative phosphorylation','Pentose phosphate pathway',
  'Phenylalanine, tyrosine and tryptophan biosynthesis',
  'Alanine, aspartate and glutamate metabolism',
  'Tryptophan metabolism','Tyrosine metabolism','Pyruvate metabolism'
)

down_pathways <- c(
  'Fatty acid biosynthesis','Fatty acid elongation',
  'Lysine degradation','Sphingolipid metabolism',
  'Synthesis and degradation of ketone bodies',
  'Terpenoid backbone biosynthesis',
  'N-Glycan biosynthesis'
)

dol_pathways <- c(
  'Amino sugar and nucleotide sugar metabolism',
  'Glycosphingolipid biosynthesis - globo and isoglobo series',
  'Taurine and hypotaurine metabolism',
  'Alanine, aspartate and glutamate metabolism',
  'Arginine and proline metabolism',
  'Arginine biosynthesis',
  'Biosynthesis of unsaturated fatty acids',
  'Folate biosynthesis',
  'Glycine, serine and threonine metabolism'
)
analyze_metabolism_trajectory <- function(traj_file, tag){
  
  traj <- qread(traj_file)
  
  meta <- traj@meta.data %>%
    select(Monocle3_Pseudotime, Lineage1, subclusters) %>%
    arrange(Monocle3_Pseudotime)
  
  meta$Pseudotime_group <- cut(
    meta$Monocle3_Pseudotime,
    breaks = quantile(meta$Monocle3_Pseudotime, probs=seq(0,1,0.01)),
    include.lowest = TRUE,
    labels = paste0("Pse_",1:100)
  )
  
  common_cells <- intersect(rownames(meta), colnames(metabolism_mat))
  meta <- meta[common_cells,]
  metab <- metabolism_mat[,common_cells]
  
  ## pseudotime bin → pathway mean
  df <- as.data.frame(t(metab))
  df$Pseudotime_group <- meta$Pseudotime_group
  
  long <- df %>%
    pivot_longer(-Pseudotime_group, names_to="Pathway", values_to="Score")
  
  pt_mat <- long %>%
    group_by(Pathway, Pseudotime_group) %>%
    summarise(meanScore=mean(Score), .groups="drop") %>%
    pivot_wider(names_from=Pseudotime_group, values_from=meanScore)
  
  mat <- as.data.frame(pt_mat)
  rownames(mat) <- mat$Pathway
  mat <- mat[,-1]
  
  ## Z-score
  zmat <- t(scale(t(mat)))
  zmat[zmat > 2] <- 2
  zmat[zmat < -2] <- -2
  
  ## 全通路动态
  ck <- clusterData(zmat, cluster.method="kmeans", cluster.num=5)
  
  pdf(paste0("result//",tag,"_metabolism_all.pdf"),4,10)
  visCluster(ck, plot.type="heatmap", show_row_dend=FALSE,
             ht.col=list(col_range=c(-2,0,2),
                         col_color=c("#ABDDA4","#EEEEEE","#F46D43")))
  dev.off()
  
  ## 重点代谢子集
  key_mat <- zmat[c(up_pathways, down_pathways, dol_pathways),]
  
  ck2 <- clusterData(key_mat, cluster.method="kmeans", cluster.num=1)
  
  pdf(paste0("result//",tag,"_metabolism_key.pdf"),4,10)
  visCluster(ck2, plot.type="heatmap", show_row_dend=FALSE,
             ht.col=list(col_range=c(-2,0,2),
                         col_color=c("#4ea64a","#EEEEEE","#8e4c99")),
             markGenes=rownames(key_mat))
  dev.off()
  
  qsave(mat, paste0("result//",tag,"_M_score.qs"))
}
analyze_metabolism_trajectory("result//Traj1.qs","T1")
analyze_metabolism_trajectory("result//Traj2.qs","T2")
analyze_metabolism_trajectory("result//Traj3.qs","T3")
analyze_metabolism_trajectory("result//Traj4.qs","T4")

