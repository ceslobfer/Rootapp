library(shiny)
library(shinyjs)
library(visNetwork)
library(shinythemes)
library(shinycssloaders)
library(DT)
library(ggplot2)
library(car)
library(nortest)
library(tseries)
library(RcmdrMisc)
library(lmtest)

RM_data <- read.csv("www/CPM_raizMicrodiseccion_RM_mean.txt", sep = "\t", row.names = 1, stringsAsFactors = F )
microMeanExpression <- read.csv("www/CPM_raizMicrodiseccion_mean_total.txt", sep = "\t", row.names = 1, stringsAsFactors = F )
tablaSeqDesRoots <- read.csv("www/seq_description_raiz.txt", sep = "\t",col.names = c("seq","description"), row.names = 1,stringsAsFactors = F)
tablaDiffRoots <-read.csv("www/microdiseccion_diffExpression.txt", sep = "\t", col.names = c("seqName","Tissue","Condition","Log Fold Change"), stringsAsFactors = F )
tablaDiffTotalRoots <-read.csv("www/diff_totalRaiz_file.txt", sep = "\t", col.names = c("seqName","Tissue","Condition","Log Fold Change"), stringsAsFactors = F )
tablaSeqGoR <- read.csv("www/seq_gos_raiz.txt", sep = "\t",col.names = c("seqName","GoID"),stringsAsFactors = F)
keggTable <- read.csv("www/kegg_total_transcriptome.txt", col.names = c("seqName","KeggID"),sep = "\t", stringsAsFactors = F )
totalRootExpression <- read.csv("www/CPM_totalRaiz.txt", sep = " ",row.names = 1, stringsAsFactors = F)
hubGenes_roots <- read.csv("www/hubGenes_microdisec.txt", sep = "\t", col.names = c("seqName","module"))
module_dataRoots <- read.csv("www/roots_module.txt", sep = "\t",col.names = c("seqName","module"), stringsAsFactors = F )
seqJoinGoTable <- read.csv("www/seq_gos_raiz_join.txt", sep = "\t",col.names = c("seqName","GoID"), row.names = 1,stringsAsFactors = F)
microExpression <- read.csv("www/CPM_microdiseccion.txt", sep = " ",row.names = 1, stringsAsFactors = F)


source("Components/Comp_Tables.R")
source("Components/Comp_barplot.R")
source("Components/Comp_corTable.R")
source("Components/Comp_roots.R")
source("Components/Comp_heatmaps.R")
source("Components/Comp_Network.R")
shinyServer(function(input, output) {
    
    
#Network gene Roots
    values <- reactiveValues()
    values$list <- list()
    edgeNodes <- reactiveValues()
    edgeNodes$list <- list()
    observeEvent(input$buttonAddNetwork, {
        if (input$seqID_R_N %in% row.names(RM_data)) {
            if (input$geneName_R_N != "") {
                values$list[[input$seqID_R_N]] <- list("Name" = input$geneName_R_N,"expression" = microExpression[input$seqID_R_N,])
            }else{
                values$list[[input$seqID_R_N]] <- list("Name" = input$seqID_R_N,"expression" = microExpression[input$seqID_R_N,])
            }
        
            edgeNodes$list <- corGenerationRoots(values$list,0.05,0)
        }
    })
    
    observeEvent(input$buttonDeleteNetwork, {
        if (input$seqID_R_N %in% names(values$list)) {
            values$list[[input$seqID_R_N]] <- NULL
            edgeNodes$list <- corGenerationRoots(values$list,0.05,0)
        }
    })
    
    observeEvent(input$fileNetRoots,{
        inFile <- input$fileNetRoots
        listSeq <- read.csv(inFile$datapath, header = F,stringsAsFactors = F)
        if (ncol(listSeq) == 1) {
            for (i in c(1:nrow(listSeq))) {
                seqName <- listSeq$V1[i]
                if (seqName %in% row.names(RM_data)) {
                    values$list[[seqName]] <- list("Name" = seqName,"expression" = microExpression[seqName,])
                }
            }
        }
        if (ncol(listSeq) == 2) {
            for (i in c(1:nrow(listSeq))) {
                seqName <- listSeq$V1[i]
                nameGene <- listSeq$V2[i]
                if (seqName %in% row.names(RM_data)) {
                    values$list[[seqName]] <- list("Name" = nameGene,"expression" = microExpression[seqName,])
                }
            }
        }
        print(values$list)
        edgeNodes$list <- corGenerationRoots(values$list,0.05,0)
    })
    
    observeEvent(input$buttonEmptyNetwork, {
        for (seqID in names(values$list)) {
            values$list[[seqID]] <- NULL
        }
        edgeNodes$list <- corGenerationRoots(values$list,0.05,0)
    })
    observeEvent(input$buttonFilterNetwork,{
        edgeNodes$list <- corGenerationRoots(values$list,input$CortePvalueNetR,input$CorteCorNetR)
    })

    output$networkRoots <- renderVisNetwork({
        resNet <- edgeNodes$list
        if (length(edgeNodes$list) != 0) {
            nodes <- resNet[[1]]
            edges <- resNet[[2]]
            visNetwork(nodes, edges)  %>% 
                visOptions(highlightNearest = TRUE, height ="650px",nodesIdSelection = TRUE) %>% 
                visNodes(borderWidth = 2,physics = FALSE, font = list(size = 17, face="arial black"),
                         shadow = list(enabled = TRUE, size = 8), shape = "ellipse")%>% 
                visInteraction(navigationButtons = TRUE) %>% 
                visPhysics(enabled = FALSE,solver = "repulsion", repulsion =list(centralGravity=0,nodeDistance=400) ) %>%
                visIgraphLayout(physics = FALSE, randomSeed = 12) %>% 
                visEdges(shadow = TRUE,width = 4,length = 400, smooth = FALSE,physics = FALSE)
                
        }
        
    })

#General table Roots
    generalTableR <- eventReactive(input$buttonGeneralR, {
        seqList <- c()
        #Differential Condition
        if (!is.null(input$checkGroupR) && input$checkGroupR != "" ) {
            for (condition in input$checkGroupR) {
                condition <- strsplit(condition,split = ",")[[1]]
                searchCondition <- tablaDiffRoots[which(tablaDiffRoots$Tissue == condition[1] & tablaDiffRoots$Condition == condition[2]),1]
                if (length(seqList) == 0) {
                    seqList <- searchCondition
                }else{
                    seqList <- intersect(seqList, searchCondition)
                }
            }
        }
        
        #Go Condition
        if (!is.null(input$goRTotal)) {
            goSearch <- strsplit(input$goRTotal,split = ",")[[1]]
            for (go in goSearch) {
                seqsGo <- tablaSeqGoR[which(tablaSeqGoR$GoID == go),1]
                if (length(seqList) == 0) {
                    seqList <- seqsGo
                }else{
                    seqList <- intersect(seqList, seqsGo)
                }
            }
        }
        if (input$moduleRTotal != "Module") {
            seqsModule <- module_dataRoots[which(module_dataRoots$module==input$moduleRTotal),]$seqName
            if (length(seqList) == 0) {
                seqList <- seqsModule
            }else{
                seqList <- intersect(seqList, seqsModule)
            }
        }
        tablaSearch <- data.frame(stringsAsFactors = F)
        if (!is.null(seqList)) {
            print(seqList)
            for (seq in seqList) {
                if (seq %in% hubGenes_roots$seqName) {
                    tablaSearch <- rbind.data.frame(tablaSearch,c(seq, tablaSeqDesRoots[seq,],seqJoinGoTable[seq,],"HUB"),stringsAsFactors = F)
                }else{
                    tablaSearch <- rbind.data.frame(tablaSearch,c(seq, tablaSeqDesRoots[seq,],seqJoinGoTable[seq,],"NO-H"),stringsAsFactors = F)
                }
                
            }
        }
        if (input$goRTotal == "" && input$moduleRTotal == "Module" && is.null(input$checkGroupR)) {
            seqList <- row.names(RM_data)
            for (seq in seqList) {
                
                if (seq %in% hubGenes_roots$seqName) {
                    tablaSearch <- rbind.data.frame(tablaSearch,c(seq, tablaSeqDesRoots[seq,],seqJoinGoTable[seq,],"HUB"),stringsAsFactors = F)
                }else{
                    tablaSearch <- rbind.data.frame(tablaSearch,c(seq, tablaSeqDesRoots[seq,],seqJoinGoTable[seq,],"NO-H"),stringsAsFactors = F)
                }
                
            }
        }
        
        
        tablaSearch <- tablaSearch[which(tablaSearch[,1] %in% row.names(RM_data)),]
        colnames(tablaSearch) <- c("Seq ID","Seq Description","GO ID", "HUBGENE")
        tablaSearch
        
        
        
    })
    
    output$tableGeneralRoots <- DT::renderDataTable({
        table <- generalTableR()
        if (length(rownames(table)) > 0) {
            DT::datatable(table,extensions = 'Buttons',options = list(initComplete = JS(
                "function(settings, json) {",
                "$(this.api().table().header()).css({'background-color': '#1c1b1b', 'color': '#ffffff'});",
                "}"),dom = 'Bfrtip',pageLength = length(row.names(table)),buttons = c('copy', 'csv', 'excel', 'pdf', 'print')),
                filter = "top",
                selection = 'multiple',
                style = 'bootstrap',
                class = 'cell-border stripe',
                rownames = FALSE,
            )
        }
    })


    heatmapValuesR <- reactiveValues()
    observeEvent(input$fileHeatmapRoots,{
        inFile <- input$fileHeatmapRoots
        listSeq <- read.csv(inFile$datapath, header = T,stringsAsFactors = F)
        table <- microMeanExpression[listSeq[,1],]
        heatmapValuesR$table <- table
    })
    #Descarga heatmap Roots
    output$downloadHeatRoots <- downloadHandler(
        filename = function(){paste("heatMapRoots",'.png',sep='')},
        content = function(file){
            ggsave(filename = file, plot = heatMapRoots(heatmapValuesR$table,input$clusterHeatR) )
        }
    )
    
    newCorTableR <- eventReactive(input$buttonCorR, {
        tablaCorResR <- NULL
        if (input$seqID_R %in% row.names(RM_data)){
            tablaCorResR<- seqCorRoots(input$seqID_R)
            
        }
        tablaCorResR
    })
    
    output$tableCorR <- DT::renderDataTable({
        if (input$seqID_R %in% row.names(RM_data)){
            table <- newCorTableR()
            DT::datatable(table,extensions = 'Buttons',options = list(initComplete = JS(
                "function(settings, json) {",
                "$(this.api().table().header()).css({'background-color': '#1c1b1b', 'color': '#ffffff'});",
                "}"),dom = 'Bfrtip',pageLength = length(row.names(table)),buttons = c('copy', 'csv', 'excel', 'pdf', 'print')),
                filter = "top",
                selection = 'multiple',
                style = 'bootstrap',
                class = 'cell-border stripe',
                rownames = FALSE,
            )
        }
    })
    
    
    output$searchGO_R <- DT::renderDataTable({
        table<-searchGoSeqR(input$goR)
        DT::datatable(table,extensions = 'Buttons',options = list(initComplete = JS(
            "function(settings, json) {",
            "$(this.api().table().header()).css({'background-color': '#1c1b1b', 'color': '#ffffff'});",
            "}"),dom = 'Bfrtip',pageLength = length(row.names(table)),buttons = c('copy', 'csv', 'excel', 'pdf', 'print')),
            filter = "top",
            selection = 'multiple',
            style = 'bootstrap',
            class = 'cell-border stripe',
            rownames = FALSE,
        )
    })
    
    output$tableKeggSeq <- renderTable({
        if (input$keggID %in% keggTable$KeggID) {
            keggTable[which(keggTable$KeggID == input$keggID & keggTable$seqName %in% row.names(datosEmbriones)),]
        }else{
            if (input$seqID %in% row.names(datosEmbriones) & input$seqID %in% keggTable$seqName) {
                keggTable[which(keggTable$seqName == input$seqID),]
            }
        }
        
    })

    output$tableGOR <- renderTable({
        if (input$seqID_R %in% row.names(RM_data) &&
            input$seqID_R %in% tablaSeqGoR$seqName) {
            goTableDesR(input$seqID_R)
        }
        
    })

    output$descriptionSeq_R <- renderUI({
        if (input$seqID_R %in% row.names(RM_data)) {
            column(
                h5(p(tablaSeqDesRoots[input$seqID_R,],style="color:white;text-align:center")),
                width=12,style="background-color:black;border-radius: 8px")
        }
    })


#Barplot Raiz
    hide("barplotRoots")
    output$barplotRoots <- renderPlot({
        res <- plot.new()
        if (input$seqID_R %in% row.names(RM_data)) {
            res <- barplotRootsFunc(input$seqID_R)
        }
        res
    })
    
    observeEvent(input$seqID_R,{
        if (input$seqID_R %in% row.names(RM_data)) {
            show("barplotRoots")
        }else{
            hide("barplotRoots")
        }
        
    })

#Heatmap plot
hide("heatmapPlot")
output$heatmapPlot <- renderPlot({
    if (unique(rownames(heatmapValuesR$table) %in% row.names(RM_data))) {
        heatMapRoots(heatmapValuesR$table,input$clusterHeatR)
    }
})

observeEvent(input$fileHeatmapRoots,{
    if (unique(rownames(heatmapValuesR$table) %in% row.names(RM_data))) {
        show("heatmapPlot")
    }else{
        hide("heatmapPlot")
    }
})


#Legend plot Roots
    hide("legendtHeatR")
    output$legendtHeatR <- renderImage({
        if (input$seqID_R %in% row.names(RM_data)) {
            legendColorsR(input$seqID_R)
        }
        list(src="www/legendR.svg",width = "310px",
             height = "180px")
    },deleteFile = TRUE)
    
    observeEvent(input$seqID_R,{
        if (input$seqID_R %in% row.names(RM_data)) {
            show("legendtHeatR")
        }else{
            hide("legendtHeatR")
        }
        
    })

    #Legend plot Roots
    hide("legendtHeatE")
    output$legendtHeatE <- renderImage({
        if (input$seqID %in% row.names(datosEmbriones)) {
            legendColorsE(input$seqID)
        }
        list(src="www/legendE.svg",width = "310px",
             height = "180px")
    },deleteFile = TRUE)
    
    observeEvent(input$seqID,{
        if (input$seqID %in% row.names(datosEmbriones)) {
            show("legendtHeatE")
        }else{
            hide("legendtHeatE")
        }
        
    })    
    
    
    output$diffRootTable <- renderTable({
        if (input$seqID_R %in% tablaDiffRoots$seqName) {
            tableDiffExpressionRoots(input$seqID_R)
        }
        
    })
   

    #Roots Heatmaps
    hide("heatRoots_N")
    output$heatRoots_N <- renderImage({
        if (input$seqID_R %in% row.names(RM_data)) {
            heatMapRootsN(input$seqID_R)
        }
        list(src="www/raiz_heat_modN.svg",width = "280px",
             height = "660px",style="margin-right:-40px")
    },deleteFile = TRUE)
    
    observeEvent(input$seqID_R,{
        if (input$seqID_R %in% row.names(RM_data)) {
            show("heatRoots_N")
        }else{
            hide("heatRoots_N")
        }
        
    })
    
    hide("heatRoots_C")
    output$heatRoots_C <- renderImage({
        if (input$seqID_R %in% row.names(RM_data)) {
            heatMapRootsC(input$seqID_R)
        }
        list(src="www/raiz_heat_modC.svg",width = "280px",
             height = "660px",style="margin-right:-40px")
    },deleteFile = TRUE)
    
    observeEvent(input$seqID_R,{
        if (input$seqID_R %in% row.names(RM_data)) {
            show("heatRoots_C")
        }else{
            hide("heatRoots_C")
        }
        
    })

})