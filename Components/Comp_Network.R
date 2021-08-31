library(shiny)
library(xtable)
module_dataRoots <- read.csv("www/roots_module.txt", sep = "\t",col.names = c("seqName","module"), stringsAsFactors = F )
tablaDiffRoots <-read.csv("www/microdiseccion_diffExpression.txt", sep = "\t", col.names = c("seqName","Tissue","Condition","Log Fold Change"), stringsAsFactors = F )
tablaSeqDesRoots <- read.csv("www/seq_description_raiz.txt", sep = "\t",col.names = c("seq","description"), row.names = 1,stringsAsFactors = F)


# listGenes <- list("ppt_57269" = list("Name" = "ppt_57269","expression" = datosEmbriones["ppt_57269",]),
#                  "ppt_251470" = list("Name" = "ppt_251470","expression" = datosEmbriones["ppt_251470",]),
#                   "ppt_57268" = list("Name" = "ppt_57267","expression" = datosEmbriones["ppt_57267",]))
#print(listGenes)

moduleList <- function(){
  resList <- list()
  resList[["Module selection"]] <- "Module"
  for (module in unique(module_dataRoots$module)) {
    resList[module]<-module
  }
  return(resList)
}

corGenerationRoots <- function(listGenes, pvalueFilter, corFilter){
  NameVector <- c()
  sourceV <- c()
  targetV <- c()
  corValuesV <- c()
  signCor <- c()
  colorGene <- c()
  labelGene <- c()
  toolTipGene <- c()
  groupGene <- c()
  for (seq1 in names(listGenes)) {
    #Node configuration
    if (!(seq1 %in% module_dataRoots$seqName)) {
      colorGene <-c(colorGene,"gray")
      groupGene <- c(groupGene,"gray")
    }else{
      colorGene <-c(colorGene,module_dataRoots[which(module_dataRoots$seqName==seq1),]$module)
      groupGene <- c(groupGene,module_dataRoots[which(module_dataRoots$seqName==seq1),]$module)
    }
    print(paste0("<b>", listGenes[[seq1]][["Name"]],"</b>"))
    labelGene <- c(labelGene,listGenes[[seq1]][["Name"]])
    diffTable <- tablaDiffRoots[which(tablaDiffRoots$seqName ==  seq1),c(2,3,4)]
    htmlTable <- print(xtable(diffTable, align="llll"), include.rownames=FALSE,
          type="html")
    toolTipGene <- c(toolTipGene,paste0("<p><b>", seq1,"</b><br>",tablaSeqDesRoots[seq1,],"<br>",
                                        htmlTable,"</p>"))
    for (seq2 in names(listGenes)) {
      pasteSeq <- ""
      if (seq1 > seq2) {
        pasteSeq <- paste(seq1,seq2,sep = "")
      }else{
        if (seq2 > seq1) {
          pasteSeq <- paste(seq2,seq1,sep = "")
          
        }
      }
      
      if (seq1 != seq2 && !(pasteSeq %in% NameVector)) {
        NameVector <- c(NameVector,pasteSeq)
        corTest <- cor.test(as.numeric(listGenes[[seq1]]$expression),as.numeric(listGenes[[seq2]]$expression))
        corValue <- corTest$estimate
        pvalue <- corTest$p.value
        if (pvalue < pvalueFilter && abs(as.numeric(corValue))>corFilter) {
          sourceV <- c(sourceV, seq1 )
          targetV <- c(targetV, seq2 )
          corValuesV <- c(corValuesV, as.character(round(as.numeric(corValue), digits = 3)))
          if (as.numeric(corValue > 0)) {
            signCor <- c(signCor, "red")
          }else{
            signCor <- c(signCor, "blue")
          }
          
        }
        

      }
      
    }
  }
  
  nodes <- data.frame(id = names(listGenes), 
                      color = list(background = colorGene, border = rep("balck",length(names(listGenes)))), 
                      label= labelGene,
                      title = toolTipGene)
  edges <- data.frame(from = sourceV, to = targetV,
                      label = corValuesV,
                      color = signCor)
  return(list(nodes, edges))
}

