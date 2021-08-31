library(ggplot2)
totalRootExpression <- read.csv("www/CPM_totalRaiz.txt", sep = " ",row.names = 1, stringsAsFactors = F)


RC_data <- read.csv("www/CPM_raizMicrodiseccion_RC_mean.txt", sep = "\t", row.names = 1, stringsAsFactors = F )
RDV_data <- read.csv("www/CPM_raizMicrodiseccion_RDV_mean.txt", sep = "\t", row.names = 1, stringsAsFactors = F )
RDC_data <- read.csv("www/CPM_raizMicrodiseccion_RDC_mean.txt", sep = "\t", row.names = 1, stringsAsFactors = F )
RM_data <- read.csv("www/CPM_raizMicrodiseccion_RM_mean.txt", sep = "\t", row.names = 1, stringsAsFactors = F )
barplotRootsFunc <- function(seq){
  RC_C <- RC_data[seq,"RC_C"]
  RC_N <- RC_data[seq,"RC_N"]
  RDC_C <- RDC_data[seq,"RDC_C"]
  RDC_N <- RDC_data[seq,"RDC_N"]
  RDV_C <- RDV_data[seq,"RDV_C"]
  RDV_N <- RDV_data[seq,"RDV_N"]
  RM_C <- RM_data[seq,"RM_C"]
  RM_N <- RM_data[seq,"RM_N"]
  expression <- c(RC_C,RDC_C,RDV_C,RM_C,RC_N,RM_N,RDC_N,RDV_N)
  barData <- data.frame(stringsAsFactors = F)
  barData <- rbind.data.frame(c(RC_N,RM_N,RDC_N,RDV_N), stringsAsFactors = F)
  colnames(barData) <- c("RC_N-RC_C","RM_N-RM_C","RDC_N-RDC_C","RDV_N-RDV_C")
  barData <- rbind.data.frame(barData,c(RC_C,RM_C,RDC_C,RDV_C), stringsAsFactors = F)
  row.names(barData)<-c("Ammonium (3mM)","Control")
  rangeExpression <- 25*((max(expression)-min(expression))/100)
  barplot(as.matrix(barData), main="CPM expression comparative",font.lab = 2,cex.lab=1.2,cex.names = 0.6, font = 2,
          xlab="Comparatives",ylab = "CPM expression",ylim = c(0,(max(expression) + (rangeExpression))), density=c(30,30),angle =c(45) ,col=c("darkblue","red"),
          legend = rownames(barData),args.legend = list(text.font = 2),beside=TRUE)
}

barplotRootsTotalFunc <- function(seq){
  C <- mean(as.numeric(totalRootExpression[seq,c(1,3,5)]))
  N <- mean(as.numeric(totalRootExpression[seq,c(2,4,6)]))
  expression <- c(C,N)
  barData <- data.frame(stringsAsFactors = F)
  barData <- rbind.data.frame(expression, stringsAsFactors = F)
  colnames(barData) <- c("Control","Ammonium (3mM)")
  row.names(barData)<-c("Total Root")
  rangeExpression <- 25*((max(expression)-min(expression))/100)
  barplot(as.matrix(barData), main="CPM expression comparative",font.lab = 2,cex.lab=1.2,cex.names = 0.6, font = 2,
          xlab="Comparatives",ylab = "CPM expression",ylim = c(0,(max(expression) + (rangeExpression))), density=c(30,30),angle =c(45) ,col=c("darkblue","red"),
          legend = NULL,args.legend = list(text.font = 2),beside=TRUE)
}
