suppressPackageStartupMessages({
  library(qs)
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
  library(ggplot2)
  library(cowplot)
  library(RColorBrewer)
  library(nichenetr)
  library(magrittr)
})
sub_Olis <- qread(file = "sub_Olisjw_analysis.qs")

Olig <- sub_Olis@tools$DEtest_subclusters$AllMarkers_wilcox
ribo_genes=rownames(sub_Olis)[grep("^Rps|^Rpl|^mt-|^H2-|^Gm|^Hbb-|^Hba-|-",rownames(sub_Olis),ignore.case=T)]
print(ribo_genes)

Olig <- dplyr::filter(Olig,!Olig$gene %in%ribo_genes)
unique(Olig$gene)

table(Olig$group1)

Olig <- dplyr::filter(Olig,Olig$p_val_adj < 0.05)
unique(Olig$gene)
DEGs <- Olig[with(Olig, avg_log2FC > 0.5 & p_val_adj < 0.05), ]

table(DEGs$group1)

DEGs <- DEGs[,c('gene','group1')]

# df <- list(DEGs)
result <- split(DEGs$gene, DEGs$group1)
result[['background']] <- unique(Olig$gene)
saveRDS(result,'rawdata/OL_markers.rds')
# Heatmap of regulatory potential for top-ranked ligands
lr_network <- readRDS(("/home/data/caixing/nichenetr/nichenet-mouse-network/lr_network_mouse_21122021.rds"))
ligand_target_matrix <- readRDS(("/home/data/caixing/nichenetr/nichenet-mouse-network/ligand_target_matrix_nsga2r_final_mouse.rds"))
weighted_networks <- readRDS(("/home/data/caixing/nichenetr/nichenet-mouse-network/weighted_networks_nsga2r_final_mouse.rds"))

# load("./data/fib_markers.RData")
DElist.ls <- readRDS("rawdata/OL_markers.rds")
figdir <- 'result/nichenet'


getNichenet <- function(gene.oi, gene.bg, out.prefix) {
  str(ligand_target_matrix)
  lr_network_expressed <- lr_network %>% filter(to %in% gene.bg)
  potential_ligands <- lr_network_expressed %>%
    pull(from) %>%
    unique()
  ligand_activities <- predict_ligand_activities(
    geneset = gene.oi,
    background_expressed_genes = gene.bg,
    ligand_target_matrix = ligand_target_matrix,
    potential_ligands = potential_ligands
  )
  write_csv(ligand_activities, sprintf("%s.ligand_activities.csv", out.prefix))
  return(ligand_activities)
}

plotNichenet <- function(gene.oi, gene.bg, out.prefix, ligand_activities,
                         width.target = 7, width.pearson = 1.3, width.comb = 8, height = 7,
                         rel_widths = c(0.19, 0.81), comb.rel.height = c(10, 2),
                         n.top = 20, pearson.max = 0.2) {
  # ligand target
  best_upstream_ligands <- ligand_activities %>%
    top_n(n.top, pearson) %>%
    arrange(-pearson) %>%
    pull(test_ligand)
  active_ligand_target_links_df <- best_upstream_ligands %>%
    lapply(get_weighted_ligand_target_links,
           geneset = gene.oi,
           ligand_target_matrix = ligand_target_matrix, n = 250
    ) %>%
    bind_rows()
  active_ligand_target_links_df <- active_ligand_target_links_df %>% filter(!is.na(weight))
  active_ligand_target_links <- prepare_ligand_target_visualization(
    ligand_target_df = active_ligand_target_links_df,
    ligand_target_matrix = ligand_target_matrix,
    cutoff = 0.25
  )
  order_ligands <- intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
  order_targets <- active_ligand_target_links_df$target %>% unique()
  order_targets <- intersect(order_targets, rownames(active_ligand_target_links))
  active_ligand_target_links.debug <<- active_ligand_target_links
  order_targets.debug <<- order_targets
  order_ligands.debug <<- order_ligands
  vis_ligand_target <- active_ligand_target_links[order_targets, order_ligands] %>% t()
  colnames(vis_ligand_target) <- make.names(colnames(vis_ligand_target))
  gene.unexplained <- setdiff(gene.oi, colnames(active_ligand_target_links))
  vis_ligand_target.debug <<- vis_ligand_target
  vis_ligand_target[vis_ligand_target > 0.01] <- 0.01
  
  # ligand pearson
  ligand_pearson_matrix <- ligand_activities %>%
    select(pearson) %>%
    as.matrix() %>%
    magrittr::set_rownames(ligand_activities$test_ligand)
  
  vis_ligand_pearson <- ligand_pearson_matrix[order_ligands, ] %>%
    as.matrix(ncol = 1) %>%
    magrittr::set_colnames("Pearson")
  vis_ligand_pearson.debug <<- vis_ligand_pearson
  vis_ligand_pearson[vis_ligand_pearson > pearson.max] <- pearson.max
  
  # ligand target
  lr_network_top <- lr_network %>%
    filter(from %in% best_upstream_ligands & to %in% gene.bg) %>%
    distinct(from, to)
  best_upstream_receptors <- lr_network_top %>%
    pull(to) %>%
    unique()
  lr_network_top_df <- weighted_networks$lr_sig %>%
    filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors)
  lr_network_top_df <- lr_network_top_df %>% spread("from", "weight", fill = 0)
  lr_network_top_matrix <- lr_network_top_df %>%
    select(-to) %>%
    as.matrix() %>%
    magrittr::set_rownames(lr_network_top_df$to)
  
  dist_receptors <- dist(lr_network_top_matrix, method = "binary")
  hclust_receptors <- hclust(dist_receptors, method = "ward.D2")
  order_receptors <- hclust_receptors$labels[hclust_receptors$order]
  
  dist_ligands <- dist(lr_network_top_matrix %>% t(), method = "binary")
  hclust_ligands <- hclust(dist_ligands, method = "ward.D2")
  order_ligands_receptor <- hclust_ligands$labels[hclust_ligands$order]
  
  vis_ligand_receptor_network <- lr_network_top_matrix[order_receptors, order_ligands_receptor]
  
  p_ligand_pearson <- vis_ligand_pearson %>% make_heatmap_ggplot("Prioritized ligands", "Ligand activity",
                                                                 color = brewer.pal(9, "Oranges")[6], legend_position = "top",
                                                                 x_axis_position = "top",
                                                                 legend_title = "Pearson correlation coefficient\ntarget gene prediction ability)"
  )
  p_ligand_pearson <- p_ligand_pearson + theme(legend.key.width = unit(2, "cm"), axis.text.y = element_text(face = "italic", size = 12, color = "black"))
  ggsave(file = sprintf("%s.ligand_pearson.pdf", out.prefix), width = width.pearson, height = height)
  
  # ligand target
  p_ligand_target_network <- vis_ligand_target %>%
    make_heatmap_ggplot("Prioritized ligands", "genes in receiver cells",
                        color = "#602ca3", legend_position = "top",
                        x_axis_position = "top",
                        legend_title = "Regulatory potential"
    ) +
    theme(axis.text.x = element_text(face = "italic", size = 12), legend.key.width = unit(2, "cm"))
  ggsave(file = sprintf("%s.ligand_target.pdf", out.prefix), width = width.target, height = height)
  
  # ligand combine
  figures_without_legend <- plot_grid(
    p_ligand_pearson +
      theme(legend.position = "none", axis.ticks = element_blank()) +
      theme(axis.title.x = element_text()),
    p_ligand_target_network + theme(
      legend.position = "none",
      axis.text.y = element_blank(),
      axis.title.y = element_blank(),
      axis.ticks = element_blank()
    ),
    align = "h", nrow = 1,
    rel_widths = rel_widths,
    rel_heights = c(
      nrow(vis_ligand_pearson),
      nrow(vis_ligand_target) + 3
    )
  )
  legends <- plot_grid(as_ggplot(get_legend(p_ligand_pearson)),
                       as_ggplot(get_legend(p_ligand_target_network)),
                       nrow = 2, align = "h"
  )
  pp <- plot_grid(figures_without_legend, legends, rel_heights = comb.rel.height, nrow = 2, align = "hv")
  ggsave(file = sprintf("%s.comb.pdf", out.prefix), pp, width = width.comb, height = height)
}

myoF.clusters <- unique(Olig$group1)
for (i in 1:length(myoF.clusters)) {
  gene.oi <- DElist.ls[[myoF.clusters[i]]]
  gene.bg <- DElist.ls[["background"]]
  
  ligand_activities <- getNichenet(
    gene.oi = gene.oi,
    gene.bg = gene.bg,
    out.prefix = sprintf("%s/%s.nichenet", figdir, myoF.clusters[i])
  )
  plotNichenet(gene.oi, gene.bg,
               out.prefix = sprintf("%s/%s.nichenet", figdir, myoF.clusters[i]),
               ligand_activities = ligand_activities, n.top = 20,
               width.target = 13.1,
               width.pearson = 1.3,
               width.comb = 19, height = 6
  )
}

