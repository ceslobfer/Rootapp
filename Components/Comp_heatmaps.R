library("pheatmap")


heatMapRoots <- function(table,clusterCut){
  table[table==0] <- Inf
  logMinValue <- log(min(table))
  table <- log10(table)
  table <- table + abs(logMinValue)
  table[table==Inf] <- 0
  if (length(row.names(table))<=55) {
    if (clusterCut == 0) {
      pheatmap(table,clustering_method = "ward.D2", cluster_cols = F, cluster_rows = F)
    }else{
      pheatmap(table, cutree_rows = clusterCut,clustering_method = "ward.D2", cluster_cols = F)
    }
    
  }else{
    if (clusterCut == 0) {
      pheatmap(table,clustering_method = "ward.D2", cluster_cols = F,cluster_rows = F,show_rownames = F)
    }else{
      pheatmap(table, cutree_rows = clusterCut,clustering_method = "ward.D2", cluster_cols = F,show_rownames = F)
    }
    
  }
}



