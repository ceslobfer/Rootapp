RC_data <- read.csv("www/CPM_raizMicrodiseccion_RC_mean.txt", sep = "\t", row.names = 1, stringsAsFactors = F )
RDV_data <- read.csv("www/CPM_raizMicrodiseccion_RDV_mean.txt", sep = "\t", row.names = 1, stringsAsFactors = F )
RDC_data <- read.csv("www/CPM_raizMicrodiseccion_RDC_mean.txt", sep = "\t", row.names = 1, stringsAsFactors = F )
RM_data <- read.csv("www/CPM_raizMicrodiseccion_RM_mean.txt", sep = "\t", row.names = 1, stringsAsFactors = F )
tablaSeqDesRoots <- read.csv("www/seq_description_raiz.txt", sep = "\t",col.names = c("seq","description"), row.names = 1,stringsAsFactors = F)
totalRootExpression <- read.csv("www/CPM_totalRaiz.txt", sep = " ",row.names = 1, stringsAsFactors = F)
microExpression <- read.csv("www/CPM_microdiseccion.txt", sep = " ",row.names = 1, stringsAsFactors = F)


seqCorRoots <- function(seq){
  tablaCor <- data.frame(stringsAsFactors = F)
  seqQuery <- as.numeric(microExpression[seq,])
  for (seq2 in row.names(microExpression)) {
    if (seq2!=seq) {
      correlation <- cor.test(seqQuery,as.numeric(microExpression[seq2,]))
      if (correlation$p.value < 0.05){
        tablaCor <- rbind.data.frame(tablaCor,c(seq2, tablaSeqDesRoots[seq2,],correlation$estimate,correlation$p.value),stringsAsFactors = F)
      }
    }
  }
  colnames(tablaCor) <- c("Sequence.ID","Blast.description","Correlation.value","Correlation.pvalue")
  tablaCor <- transform(tablaCor, Correlation.value = as.numeric(Correlation.value))
  tablaCor <- transform(tablaCor, Correlation.pvalue = as.numeric(Correlation.pvalue))
  tablaCor <- tablaCor[order(abs(tablaCor$Correlation.value), decreasing = TRUE),]
  return(tablaCor)
}