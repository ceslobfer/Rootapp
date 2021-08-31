
tablaSeqGoR <- read.csv("www/seq_gos_raiz.txt", sep = "\t",col.names = c("seqName","GoID"),stringsAsFactors = F)
tablaSeqDesRoots <- read.csv("www/seq_description_raiz.txt", sep = "\t",col.names = c("seq","description"), row.names = 1,stringsAsFactors = F)
RM_data <- read.csv("www/CPM_raizMicrodiseccion_RM_mean.txt", sep = "\t", row.names = 1, stringsAsFactors = F )



goTableDesR <- function(seq){
  
  goTable <- data.frame(stringsAsFactors = F)
  tablaSearch <- data.frame(stringsAsFactors = F)
  if (seq %in% tablaSeqGoR$seqName) {
    tableSearch <- tablaSeqGoR[which(tablaSeqGoR$seqName == seq),]
    for (go in tableSearch$GoID) {
      goTable <- rbind(goTable,c(go,tablaGoDes[go,]),stringsAsFactors = F)
    }
    colnames(goTable) <- c("GO id", "GO descriptor")
  }
  return(goTable)
}




searchGoSeqR <- function(go){
  res <- NULL
  tablaSearch <- data.frame(stringsAsFactors = F)
  if (go %in% tablaSeqGoR$GoID) {
    tablaSearch <- tablaSeqGoR[which(tablaSeqGoR$GoID == go),]
    
    for (seq in tablaSearch$seqName) {
      if (seq %in% row.names(RM_data)) {
        res <- rbind.data.frame(res,c(seq,tablaSeqDesRoots[seq,]),stringsAsFactors = F)
      }
      
    }
    colnames(res) <- c("Sequence ID", "NCBI description blast")
  }
  
  return(res)
}

