#####Ro/e#####
sub_Olis <- qread(file = "sub_Olisjw_analysis_pro.qs")

meta = sub_Olis@meta.data
table(meta$slice,meta$subclusters)

#----- 03. Ro/e -----#
if(TRUE){
  divMatrix <- function(m1, m2){
    ## Divide each element in turn in two same dimension matrixes
    ##
    ## Args:
    #' @m1: the first matrix
    #' @m2: the second matrix
    ##
    ## Returns:
    ## a matrix with the same dimension, row names and column names as m1. 
    ## result[i,j] = m1[i,j] / m2[i,j]
    dim_m1 <- dim(m1)
    dim_m2 <- dim(m2)
    if( sum(dim_m1 == dim_m2) == 2 ){
      div.result <- matrix( rep(0,dim_m1[1] * dim_m1[2]) , nrow = dim_m1[1] )
      row.names(div.result) <- row.names(m1)
      colnames(div.result) <- colnames(m1)
      for(i in 1:dim_m1[1]){
        for(j in 1:dim_m1[2]){
          div.result[i,j] <- m1[i,j] / m2[i,j]
        }
      }   
      return(div.result)
    }
    else{
      warning("The dimensions of m1 and m2 are different")
    }
  }
  
  ROIE <- function(crosstab){
    ## Calculate the Ro/e value from the given crosstab
    ##
    ## Args:
    #' @crosstab: the contingency table of given distribution
    ##
    ## Return:
    ## The Ro/e matrix 
    rowsum.matrix <- matrix(0, nrow = nrow(crosstab), ncol = ncol(crosstab))
    rowsum.matrix[,1] <- rowSums(crosstab)
    colsum.matrix <- matrix(0, nrow = ncol(crosstab), ncol = ncol(crosstab))
    colsum.matrix[1,] <- colSums(crosstab)
    allsum <- sum(crosstab)
    roie <- divMatrix(crosstab, rowsum.matrix %*% colsum.matrix / allsum)
    row.names(roie) <- row.names(crosstab)
    colnames(roie) <- colnames(crosstab)
    return(roie)
  }
}

meta$month
colnames(meta)
##--- Figure 4A; Tissue preference of each B cell subset evaluated by the Ro/e index
plot_df = meta #%>% filter(Tissue %in% c("Blood","Adjacent non-tumor","Tumor"))

plot_df$Annotation <- plot_df$subclusters
plot_df$Annotation = as.character(plot_df$Annotation)
plot_df$Tissue = as.character(plot_df$slice)
table(plot_df$slice)

summary <- table(plot_df[,c('Annotation','Tissue')])
roe <- as.data.frame(ROIE(summary))
# roe <- roe[,c("WT","5xFAD")]

require(ComplexHeatmap)
require(circlize)
col_fun <- colorRamp2(c(0, 1, 1.5, 2),
                      c("#fffde7", "#ffe0b2","#ff9800", "#e65100"))

Heatmap(t(roe),
        col = col_fun,
        cluster_rows = F,
        cluster_columns = T,
        clustering_distance_rows = "euclidean",
        clustering_method_rows = "ward.D",
        column_names_gp = grid::gpar(fontsize = 6),
        row_names_gp = grid::gpar(fontsize = 6))

pdf(paste0("result/Fig4A_roe_slice.pdf"),width = 4,height = 3.5)
Heatmap(
  # t(roe[,c('ENT','PFC','RS','STR','HPF','BF','TH')]),
  roe,
  col = col_fun,
  cluster_rows = F,
  cluster_columns = T,
  clustering_distance_rows = "euclidean",
  clustering_method_rows = "ward.D",
  column_names_gp = grid::gpar(fontsize = 6),
  row_names_gp = grid::gpar(fontsize = 6),
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(sprintf("%.1f", roe[i, j]), x, y, gp = gpar(fontsize = 6))
  },
  width = ncol(roe) * unit(0.3, "inch"),
  height = nrow(roe) * unit(0.15, "inch"),
  name = "Ro/e"
)
dev.off()

#####