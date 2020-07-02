#Plots for Mouse Embryo Trajectory Analysis

library(SummarizedExperiment)
library(SingleCellExperiment)
library(tradeSeq)
library(ggplot2)
library(fishpond)
library(readr)
library(dplyr)
library(tikzDevice)
library(gtools)
base_dir <- "/Users/Scott/Documents/Dissertation Data/SingleCellProject/MouseEmbryoData/"
source('~/Documents/Dissertation/SingleCellProject/code/Simulation Code/SingleCellProjectFunctions.R')

PlotsSaveDir <- paste0(base_dir, "EmVsNoEMPlots/")
if(!dir.exists(PlotsSaveDir)){dir.create(PlotsSaveDir)}

tikzLockFile <- paste0(PlotsSaveDir, "/TikzTemp___LOCK")
if(file.exists(tikzLockFile)){system(paste0("rm ", tikzLockFile))}

options("tikzDocumentDeclaration" = "\\documentclass[12pt,letterpaper,twoside]{report}\n")
options("tikzMetricsDictionary" = paste(PlotsSaveDir, "/TikzTemp", sep = ""))
options("tikzLatexPackages" = c(
  "\\usepackage{tikz}",
  "\\usepackage[active,tightpage,psfixbb]{preview}",
  "\\PreviewEnvironment{pgfpicture}",
  "\\setlength\\PreviewBorder{0pt}",
  "\\usepackage{bm}",
  "\\usepackage{newtxtext}")
)


pdflatex_loc <- '/Library/TeX/texbin/pdflatex'

stdal <- TRUE

gene_mappingT <- data.frame(read_tsv(file.path(base_dir, "GeneNamesMapping/gnames.txt"),
                                     col_names = F))
colnames(gene_mappingT) <- c("gene_id", "MGI")
gene_mapping <- distinct(gene_mappingT)

#This is a function written by Mike that extracts hack to get time, yhat, lineage from tradeSeq for one gene to be able to construct expression
  #profiles across pseudotime
  #see https://gist.github.com/mikelove/c91d293312b1cd2f9d71657819dc41a9
getCurvesMike <- function(models, counts, gene, nPoints = 100, lwd = 2, size = 2/3, alpha = 2/3) {
  if (is(gene, "character")) {
    id <- which(names(models) %in% gene)
  }
  else id <- gene
  if(length(id)==0){
    stop("No fitted values for this gene were found")
  }
  dm <- colData(models)$tradeSeq$dm
  y <- unname(counts[id, ])
  X <- colData(models)$tradeSeq$X
  slingshotColData <- colData(models)$slingshot
  pseudotime <- slingshotColData[, grep(x = colnames(slingshotColData),
                                        pattern = "pseudotime")]
  if (is.null(dim(pseudotime)))
    pseudotime <- matrix(pseudotime, ncol = 1)
  nCurves <- length(grep(x = colnames(dm), pattern = "t[1-9]"))
  betaMat <- rowData(models)$tradeSeq$beta[[1]]
  beta <- betaMat[id, ]
  lcol <- timeAll <- rep(0, nrow(dm))
  for (jj in seq_len(nCurves)) {
    for (ii in seq_len(nrow(dm))) {
      if (dm[ii, paste0("l", jj)] == 1) {
        timeAll[ii] <- dm[ii, paste0("t", jj)]
        lcol[ii] <- jj
      }
      else {
        next
      }
    }
  }
  df <- data.frame(time = timeAll, gene_count = y,
                   lineage = as.character(lcol))
  rows <- sample(seq_len(nrow(df)), nrow(df), replace = FALSE)
  df <- df[rows, ]
  out <- list()
  for (jj in seq_len(nCurves)) {
    df <- tradeSeq:::.getPredictRangeDf(dm, jj, nPoints = nPoints)
    Xdf <- tradeSeq:::predictGAM(lpmatrix = X, df = df, pseudotime = pseudotime)
    yhat <- c(exp(t(Xdf %*% t(beta)) + df$offset))
    out[[jj]] <- data.frame(time=df[,paste0("t",jj)], yhat, lineage=as.character(jj))
  }
  return(do.call(rbind, out))
}


SameLineages <- TRUE

if(SameLineages==TRUE){
    load("/Users/Scott/Documents/Dissertation Data/SingleCellProject/MouseEmbryoData/tradeSeqMouseDataResultsEMandNoEMSameLineages.RData")
    AdjPvals <- loadRData("/Users/Scott/Documents/Dissertation Data/SingleCellProject/MouseEmbryoData/AllAdjPvalsAssociationSameLineages20InfReps.RData")
    SamLin <- "SameLineages"
}else if(SameLineages==FALSE){
  load("/Users/Scott/Documents/Dissertation Data/SingleCellProject/MouseEmbryoData/tradeSeqMouseDataResultsEMandNoEMSeparateLineagesForNoEMSameCellClusters.RData")
  #load("/Users/Scott/Documents/Dissertation Data/SingleCellProject/MouseEmbryoData/tradeSeqMouseDataResultsEMandNoEMSeparateLineagesForNoEM.RData")
    AdjPvals <- loadRData("/Users/Scott/Documents/Dissertation Data/SingleCellProject/MouseEmbryoData/AllAdjPvalsAssociationSeparateLineages20InfReps.RData")
    SamLin <- "SeparateLineages"
}





print(AdjPvals[AdjPvals$MGI=="Nme2",])
print(AdjPvals[AdjPvals$MGI=="Rpl36a",])

#Load the data for each time point as a tximeta object
load("/Users/Scott/Documents/Dissertation Data/SingleCellProject/MouseEmbryoData/QuantObj8.5UsingEMTximetaObj.RData")
load("/Users/Scott/Documents/Dissertation Data/SingleCellProject/MouseEmbryoData/QuantObj8.0UsingEMTximetaObj.RData")
load("/Users/Scott/Documents/Dissertation Data/SingleCellProject/MouseEmbryoData/QuantObj8.25UsingEMTximetaObj.RData")

#Modify cell name to correspond with the names used in the tradeSeq analysis
colnames(QuantObj8.0UsingEM) <- paste0(colnames(QuantObj8.0UsingEM), "_1")
colnames(QuantObj8.25UsingEM) <- paste0(colnames(QuantObj8.25UsingEM), "_2")
colnames(QuantObj8.5UsingEM) <- paste0(colnames(QuantObj8.5UsingEM), "_3")


Counts8.25EMFull <- assays(QuantObj8.25UsingEM)$counts

#Get the 500 cells used in the trajectory analysis
#Get pseudotime results for EM and NoEM results
pt1 <- data.frame(sce_res[[1]]$slingshot)
pt1$cell_id <- rownames(pt1)

pt1_NoEM <- data.frame(sce_res[[2]]$slingshot)
pt1_NoEM$cell_id <- rownames(pt1_NoEM)

#Use the same 500 cells at each time point that were used in the trajectory analysis
QuantObj8.0ToUse2 <- QuantObj8.0UsingEM[,colnames(QuantObj8.0UsingEM)%in%pt1$cell_id]
QuantObj8.25ToUse2 <- QuantObj8.25UsingEM[,colnames(QuantObj8.25UsingEM)%in%pt1$cell_id]
QuantObj8.5ToUse2 <- QuantObj8.5UsingEM[,colnames(QuantObj8.5UsingEM)%in%pt1$cell_id]

QuantObj <- cbind(QuantObj8.0ToUse2, QuantObj8.25ToUse2, QuantObj8.5ToUse2)

Counts8.0EM <- as.matrix(assays(QuantObj8.0ToUse2)$counts)
Counts8.25EM <- as.matrix(assays(QuantObj8.25ToUse2)$counts)
Counts8.5EM <- as.matrix(assays(QuantObj8.5ToUse2)$counts)
# 
# #Modify cell name to correspond with the names used in the tradeSeq analysis
# colnames(Counts8.0EM) <- paste0(colnames(Counts8.0EM), "_1")
# colnames(Counts8.25EM) <- paste0(colnames(Counts8.25EM), "_2")
# colnames(Counts8.5EM) <- paste0(colnames(Counts8.5EM), "_3")
# rm(QuantObj8.0EM, QuantObj8.25EM, QuantObj8.5EM)
# gc()

#QuantObj8.0ToUse <- QuantObj8.0UsingEM
#QuantObj8.25ToUse <- QuantObj8.25UsingEM
#QuantObj8.5ToUse <- QuantObj8.5UsingEM

rm(QuantObj8.0UsingEM, QuantObj8.25UsingEM, QuantObj8.5UsingEM)
gc()

#Now, load in the tximport objects for the NoEM counts and convert them to be SummarizedExperiment objects for use with plotInfReps
load("/Users/Scott/Documents/Dissertation Data/SingleCellProject/MouseEmbryoData/QuantObj8.25NoEM.RData")
load("/Users/Scott/Documents/Dissertation Data/SingleCellProject/MouseEmbryoData/QuantObj8.5NoEM.RData")
load("/Users/Scott/Documents/Dissertation Data/SingleCellProject/MouseEmbryoData/QuantObj8.0NoEM.RData")

Counts8.0NoEM <- QuantObj8.0NoEM$counts
Counts8.25NoEM <- QuantObj8.25NoEM$counts
Counts8.5NoEM <- QuantObj8.5NoEM$counts

#Modify cell name to correspond with the names used in the tradeSeq analysis
colnames(Counts8.0NoEM) <- paste0(colnames(Counts8.0NoEM), "_1")
colnames(Counts8.25NoEM) <- paste0(colnames(Counts8.25NoEM), "_2")
colnames(Counts8.5NoEM) <- paste0(colnames(Counts8.5NoEM), "_3")
rm(QuantObj8.0NoEM, QuantObj8.25NoEM, QuantObj8.5NoEM)
gc()


Counts8.25NoEMFull <- Counts8.25NoEM

#Calculate InfRV for the genes based on the cells used in the analysis
  #This is a modification of the computeInfRV function in fishpond to support use of the mean and variance components (which fishpond did not
  #as of writing this)
computeInfRVSVB <- function(y, pc = 5, shift = 0.01){
  infVar <- assays(y)$variance
  mu <- assays(y)$mean
  InfRV <- pmax(infVar - mu, 0)/(mu + pc) + shift
  mcols(y)$meanInfRV <- rowMeans(InfRV)
  y
}

QuantObj <- computeInfRVSVB(QuantObj)

#Add InfRV to the full list of Pvalues loaded above
InfRV <- data.frame(rowData(QuantObj)$meanInfRV)
colnames(InfRV) <- "InfRV"
InfRV$gene_id <- rownames(InfRV)

AdjPvals2 <- merge(AdjPvals, InfRV, by = "gene_id")

#Make sure columns for the NoEM counts are in the same order as for the EM counts just to be safe
Counts8.0NoEM2 <- as.matrix(Counts8.0NoEM[,colnames(QuantObj8.0ToUse2)])
Counts8.25NoEM2 <- as.matrix(Counts8.25NoEM[,colnames(QuantObj8.25ToUse2)])
Counts8.5NoEM2 <- as.matrix(Counts8.5NoEM[,colnames(QuantObj8.5ToUse2)])


rm(Counts8.0NoEM, Counts8.25NoEM, Counts8.5NoEM)
gc()

#Create SummarizedExperiment object and ensure the cells are in the same order as from the object for the EM counts
QuantObjNoEMT <- SummarizedExperiment(assays = SimpleList(counts = cbind(Counts8.0NoEM2, Counts8.25NoEM2, Counts8.5NoEM2)))
QuantObjNoEM <- QuantObjNoEMT[,colnames(QuantObj)]

extractPT <- function(x, Lin = FALSE){
  #print(x)
  w1 <- as.numeric(x["cellWeights.curve1"])
  w2 <- as.numeric(x["cellWeights.curve2"])
  if(w1 > w2){
    y <- as.numeric(x["pseudotime.curve1"])
    lin <- "1"
  }else if(w2 > w1){
    y <- as.numeric(x["pseudotime.curve2"])
    lin <- "2"
  }else{
    pt1 <- as.numeric(x["pseudotime.curve1"])
    pt2 <- as.numeric(x["pseudotime.curve2"])
    #Set the seed like this to ensure cells get the same lineage assignment between EM and NoEM results if the lineages are fixed the ones
      #fit from the EM results always
    set.seed(round(mean(c(pt1,pt2), na.rm=T)*10e6))
    lin <- sample(c("1", "2"), 1)
    if(lin=="1"){
      y <- pt1
    }else if(lin=="2"){
      y <- pt2
    }
    
  }
  
  if(Lin==TRUE){
    return(lin)
  }else{
    return(y)
  }
  
}


pt2 <- subset(pt1, pt1$cell_id %in% colnames(QuantObj))
pt2_NoEM <- subset(pt1_NoEM, pt1_NoEM$cell_id %in% colnames(QuantObjNoEM))

pt2$PT <- apply(pt2, 1, extractPT)
pt2$Lin <- apply(pt2, 1, extractPT, Lin = TRUE)

pt2_NoEM$PT <- apply(pt2_NoEM, 1, extractPT)
pt2_NoEM$Lin <- apply(pt2_NoEM, 1, extractPT, Lin = TRUE)



#Make sure cells are in the same order as in the SummarizedExperiment object prior to adding the PT information
pt3 <- pt2[colnames(QuantObj),]
pt3_NoEM <- pt2_NoEM[colnames(QuantObjNoEM),]
colData(QuantObj)$PT <- pt3$PT
colData(QuantObj)$Lineage <- pt3$Lin

colData(QuantObjNoEM)$PT <- pt3_NoEM$PT
colData(QuantObjNoEM)$Lineage <- pt3_NoEM$Lin
"ENSMUSG00000059040.5"
"ENSMUSG00000002731.6"
plotRes <- function(sce_res, counts_res, gene_mapping, QuantObj, QuantObjNoEM, gene = NULL, MGI = NULL, PlotsSaveDir, 
                    PlotForPaper = FALSE, SamLin){
  if(!is.null(gene)&!is.null(MGI)){
    stop("Only need to specify one of gene and MGI since the other will be determined automatically by MGI")
  }
  if(is.null(gene)){
    gene <- gene_mapping$gene_id[gene_mapping$MGI==MGI]
  }else{
    MGI <- gene_mapping$MGI[gene_mapping$gene_id==gene]
  }
  
  if(length(gene)==0){
    print(paste0("The input MGI ", MGI, " could not be mapped to a full gene_id"))
    return(NULL)
  }
  
  if(!(gene %in% rownames(sce_res[[1]]))){
    print(paste0("EM tradeSeq results were not found for the current gene/MGI", gene, " (", MGI, ")"))
    return(NULL)
  }
  if(!(gene %in% rownames(sce_res[[2]]))){
    return(NULL)
  }
  


  fil_piece <- paste0(MGI, "Plot1", SamLin)
  fil_name <- paste0(PlotsSaveDir, fil_piece, ".tex")
  tikz(file = fil_name, height=5.5, width=8.5, standAlone = stdal, sanitize = F)
  par(cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5, mar=c(5,6,4,1)+.03)
  #fil_name <- paste0(PlotsSaveDir, fil_piece, ".pdf")
  #pdf(file = fil_name, width = 8.5, height = 5.5)
  
  #ylimUpper <- max(as.numeric(assays(QuantObj)[["mean"]][gene,]))
  ylimUpper <- max(as.numeric(assays(QuantObj)[["mean"]][gene,]) + qnorm(.975)*sqrt(as.numeric(assays(QuantObj)[["variance"]][gene,])))
  
  t1 <- getCurvesMike(sce_res[[1]], counts_res[[1]], gene = gene, nPoints = 1000)
  t2 <- getCurvesMike(sce_res[[2]], counts_res[[2]], gene = gene, nPoints = 1000)
  #main = paste0(gene, " (", MGI, ") for Counts Obtained using EM")
  plotInfReps(QuantObj, gene, x = "PT", cov = "Lineage", xlab = "pseudotime", legend = T, legendPos = "topright", legendTitle = TRUE,
              main = paste0(MGI, " for Counts Obtained using EM"), ylim = c(0, ylimUpper))
  
  lines(t1$time[1:1000], t1$yhat[1:1000], col = "dodgerblue", lwd = 2)
  lines(t1$time[1001:2000], t1$yhat[1001:2000], col = "goldenrod4", lwd = 2)
  dev.off()
  tryCatch(dev.off(), error = function(x){})
  SVBCompiletikzPlot(PlotsSaveDir, fil_piece, pdflatex_loc = pdflatex_loc)
  rm(fil_name)
  rm(fil_piece)

  print(paste0("ylim upper is ", ylimUpper))
  fil_piece <- paste0(MGI, "Plot2", SamLin)
  fil_name <- paste0(PlotsSaveDir, fil_piece, ".tex")
  tikz(file = fil_name, height=5.5, width=8.5, standAlone = stdal, sanitize = F)
  par(cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5, mar=c(5,6,4,1)+.03)
  
  #fil_name <- paste0(PlotsSaveDir, fil_piece, ".pdf")
  #pdf(file = fil_name, width = 8.5, height = 5.5)
  #main =  paste0(gene, " (", MGI, ") for Counts Obtained Without EM")
  plotInfReps(QuantObjNoEM, gene, x = "PT", cov = "Lineage", xlab = "pseudotime", legend = T, legendPos = "topright",legendTitle = TRUE,
              main =  paste0(MGI, " for Counts Obtained Without EM"), useMean = FALSE, ylim = c(0, ylimUpper))
  
  lines(t2$time[1:1000], t2$yhat[1:1000], col = "dodgerblue", lwd = 2)
  lines(t2$time[1001:2000], t2$yhat[1001:2000], col = "goldenrod4", lwd = 2)
  dev.off()
  tryCatch(dev.off(), error = function(x){})
  SVBCompiletikzPlot(PlotsSaveDir, fil_piece, pdflatex_loc = pdflatex_loc)
  rm(fil_name)
  rm(fil_piece)
  
  
  t1$type <- "EM"
  t2$type <- "NoEM"
  t3 <- rbind(t1,t2)
  t3$Lineage <- paste0(t3$type, "Lin", t3$lineage)
  #+ylim(0,10)
  #ggtitle(paste0("Predicted Counts by Lineage and Pseudotime for EM and NoEM\nAlevin Counts for ", gene, " (", MGI, ")"))
  ggp <- ggplot(data = t3, aes(x=time, y=yhat, group=Lineage, fill=type, color=Lineage))+
    geom_line(size=1.5)+scale_colour_manual(values = c("red", "black", "blue", "green4"))+
    ylab("Predicted Count")+xlab("Pseudotime")+ggtitle(MGI)+
    theme(axis.text.y=element_text(size=rel(1.5), face = "bold"),
          axis.text.x=element_text(size=rel(1.5), face = "bold"), plot.title=element_text(size=rel(1.75), hjust=.5,face = "bold"), axis.title=element_text(size=rel(1.5), face = "bold"),
          legend.text = element_text(size=rel(1.25), face = "bold"),legend.title = element_text(size=rel(1.5), face = "bold"))+
    scale_y_continuous(limits = c(-0.25, NA))
  
  
  fil_piece <- paste0(MGI, "Plot3", SamLin)
  fil_name <- paste0(PlotsSaveDir, fil_piece, ".tex")
  tikz(file = fil_name, height=4.5, width=8.5, standAlone = stdal, sanitize = F)
  #fil_name <- paste0(PlotsSaveDir, fil_piece, ".pdf")
  #pdf(file = paste0(PlotsSaveDir, "/", MGI, "Plot3.pdf"), width = 8.5, height = 5.5)
  print(ggp)
  dev.off()
  tryCatch(dev.off(), error = function(x){})
  SVBCompiletikzPlot(PlotsSaveDir, fil_piece, pdflatex_loc = pdflatex_loc)
  rm(fil_name)
  rm(fil_piece)
  
}

plotRes(sce_res = sce_res, counts_res = counts_res, gene_mapping = gene_mapping, 
        QuantObj = QuantObj, QuantObjNoEM = QuantObjNoEM, gene = "ENSMUSG00000059040.5", PlotsSaveDir = PlotsSaveDir, SamLin = SamLin)

#Search for genes where EM vs NoEM seems to especially make a difference
#First row are genes Mike identified as potentially useful
#Second row are genes chosen by SVB with high InfRV and very different pvalues between EM and NoEM
#Third row are genes from Supp Figure 2d on Page 35 of the paper, "A single-cell molecular map of mouse gastrulation and early organogenesis"
#Genes on the next line are the ones highlighted in the volcano plot in the figure mentioned on the previous line
  #but did not show a big difference between EM and NoEM counts because the counts were very similar
#"Rhox5", "Trap1a",
# MGIV <- c("Eno1b", "Eno1", "Hmgb1", "Nme2", "Cox7c", "Gfpt2", "Tnnt1", "Has2", "Taf9", "Cited1", "Cdx1", 
#           "Rpl36a", "Gm21596",
#           "B4galt6","Kitl", "Sct", "Hoxa10", "Cxcl12", "2610528A11Rik", "Tfpi", "Trf", "Xlr3a", "Trap1a")

#Temporarily change save directory to differentiate from other save plot
#PlotsSaveDir <- paste0(base_dir, "EmVsNoEMPlotsTesting/")
#if(!dir.exists(PlotsSaveDir)){dir.create(PlotsSaveDir)}
# MGIV <- c("Kpna2", "Fyb", "Gm49405", "Usp50", "Mrpl53", "Gm49322", "Mrpl23", "Gng10","Isca1",
#           "Gm9803","Igf2"     ,     "Taf9"      ,    "Gm10131"    ,   "mt-Co1"   ,     "Gm20634"   ,    "Gm10282"   ,    "1810009A15Rik" ,"Tmem254a",
#           "Grcc10"   ,     "Amd1"    ,      "Gm2000"    ,    "Ube2v1"    ,    "Sec61g"    ,    "Znrd2"    ,     "Pin4"     ,     "Hist1h2ao"  ,   "Hist1h2ap",
#           "Lbhd1"     ,    "Gm43041"    ,   "Gpx4-ps2"   ,   "Gm20075"    ,   "Gm50388"    ,   "Gm50387")
MGIV <- c("Eno1b", "Eno1", "Hmgb1", "Nme1", "Nme2", "Cox7c", "Rpl36a")

for(i in 1:length(MGIV)){
  MGI <- MGIV[i]
  plotRes(sce_res = sce_res, counts_res = counts_res, gene_mapping = gene_mapping, 
          QuantObj = QuantObj, QuantObjNoEM = QuantObjNoEM, MGI = MGI, PlotsSaveDir = PlotsSaveDir, SamLin = SamLin)
}

#Now, create the 4 panel figure for the paper
AdjPvalsForTable <- AdjPvals2[AdjPvals2$MGI%in%MGIV,]


  fil_piece <- paste0("NmeComparisonPlot", SamLin)
  fil_name <- paste0(PlotsSaveDir, fil_piece, ".tex")
  tikz(file = fil_name, height=1.25*9, width=1.25*8.5*1.545455, standAlone = stdal, sanitize = F)
  par(mfrow=c(2,2), cex.lab=3, cex.axis=2.1, cex.main=3, cex.sub=2.25, mar=c(5,6,4,1)+1)
  ylimUpper <- max(as.numeric(assays(QuantObj)[["mean"]]["ENSMUSG00000037601.11",]) + qnorm(.975)*sqrt(as.numeric(assays(QuantObj)[["variance"]]["ENSMUSG00000037601.11",])))
  
  #Nme1 results
  t1 <- getCurvesMike(sce_res[[1]], counts_res[[1]], gene = "ENSMUSG00000037601.11", nPoints = 1000)
  t2 <- getCurvesMike(sce_res[[2]], counts_res[[2]], gene = "ENSMUSG00000037601.11", nPoints = 1000)
  par(cex.lab=3, cex.axis=2.75, cex.main=3, cex.sub=3, mar=(c(5,6,4,1)+.5),xaxt='n', yaxt='n')
  plotInfReps(QuantObj, "ENSMUSG00000037601.11", x = "PT", cov = "Lineage", xlab = "", main = "", legend = F, ylim = c(0, ylimUpper))
  par(xaxt='s', yaxt='s')
  axis(side = 1, at = seq(0,100,20), las=1, font=2, labels = c("","","","","",""))
  mtext(c("0", "20", "40", "60", "80", "100"), side = 1, line = 2, cex = 2.5, at = c(0,20,40,60,80,100))
  mtext("pseudotime", side = 1, line = 4, cex = 2.5)
  axis(side = 2, at = seq(0,100,50))
  title("A\\hspace*{1.10in}Nme1", adj = 0)
  
  lines(t1$time[1:1000], t1$yhat[1:1000], col = "dodgerblue", lwd = 2)
  lines(t1$time[1001:2000], t1$yhat[1001:2000], col = "goldenrod4", lwd = 2)
  par(cex.lab=3, cex.axis=2.75, cex.main=3, cex.sub=3, mar=(c(5,6,4,1)+.5),xaxt='n', yaxt='n')
  plotInfReps(QuantObjNoEM, "ENSMUSG00000037601.11", x = "PT", cov = "Lineage", xlab = "",
              main =  "", useMean = FALSE, ylim = c(0, ylimUpper))
  par(xaxt='s', yaxt='s')
  axis(side = 1, at = seq(0,100,20), las=1, font=2, labels = c("","","","","",""))
  mtext(c("0", "20", "40", "60", "80", "100"), side = 1, line = 2, cex = 2.5, at = c(0,20,40,60,80,100))
  mtext("pseudotime", side = 1, line = 4, cex = 2.5)
  axis(side = 2, at = seq(0,100,50))
  title("B\\hspace*{1.10in}Nme1", adj = 0)
  
  lines(t2$time[1:1000], t2$yhat[1:1000], col = "dodgerblue", lwd = 2)
  lines(t2$time[1001:2000], t2$yhat[1001:2000], col = "goldenrod4", lwd = 2)
  
  #Nme2 results
  ylimUpper <- max(as.numeric(assays(QuantObj)[["mean"]]["ENSMUSG00000020857.11",]) + qnorm(.975)*sqrt(as.numeric(assays(QuantObj)[["variance"]]["ENSMUSG00000020857.11",])))
  t1 <- getCurvesMike(sce_res[[1]], counts_res[[1]], gene = "ENSMUSG00000020857.11", nPoints = 1000)
  t2 <- getCurvesMike(sce_res[[2]], counts_res[[2]], gene = "ENSMUSG00000020857.11", nPoints = 1000)
  par(cex.lab=3, cex.axis=2.75, cex.main=3, cex.sub=3, mar=(c(5,6,4,1)+.50),xaxt='n', yaxt='n')
  plotInfReps(QuantObj, "ENSMUSG00000020857.11", x = "PT", cov = "Lineage", xlab = "", main = "", legend = F, ylim = c(0, ylimUpper))
  par(xaxt='s', yaxt='s')
  axis(side = 1, at = seq(0,100,20), las=1, font=2, labels = c("","","","","",""))
  mtext(c("0", "20", "40", "60", "80", "100"), side = 1, line = 2, cex = 2.5, at = c(0,20,40,60,80,100))
  mtext("pseudotime", side = 1, line = 4, cex = 2.5)
  axis(side = 2, at = seq(0,300,100))
  title("C\\hspace*{1.10in}Nme2", adj = 0)
  
  lines(t1$time[1:1000], t1$yhat[1:1000], col = "dodgerblue", lwd = 2)
  lines(t1$time[1001:2000], t1$yhat[1001:2000], col = "goldenrod4", lwd = 2)
  par(cex.lab=3, cex.axis=2.75, cex.main=3, cex.sub=3, mar=(c(5,6,4,1)+.50),xaxt='n', yaxt='n')
  plotInfReps(QuantObjNoEM, "ENSMUSG00000020857.11", x = "PT", cov = "Lineage", xlab = "", main = "", useMean = FALSE, ylim = c(0, ylimUpper))
  par(xaxt='s', yaxt='s')
  axis(side = 1, at = seq(0,100,20), las=1, font=2, labels = c("","","","","",""))
  mtext(c("0", "20", "40", "60", "80", "100"), side = 1, line = 2, cex = 2.5, at = c(0,20,40,60,80,100))
  mtext("pseudotime", side = 1, line = 4, cex = 2.5)
  axis(side = 2, at = seq(0,300,100))
  title("D\\hspace*{1.10in}Nme2", adj = 0)
  
  lines(t2$time[1:1000], t2$yhat[1:1000], col = "dodgerblue", lwd = 2)
  lines(t2$time[1001:2000], t2$yhat[1001:2000], col = "goldenrod4", lwd = 2)
  
  dev.off()
  tryCatch(dev.off(), error = function(x){})
  SVBCompiletikzPlot(PlotsSaveDir, fil_piece, pdflatex_loc = pdflatex_loc)
  rm(fil_name)
  rm(fil_piece)

  
  
  
#Now line plot that shows the underestimation of the counts when not using EM
    #Use the cells from the 8.25 time
CountsEMT <- Counts8.25EM 
CountsNoEMT <- Counts8.25NoEM2

CountsEMT <- Counts8.25EMFull
CountsNoEMT <- Counts8.25NoEMFull[,colnames(CountsEMT)]
  
  filt_func <- function(x, FixedNGenes = FALSE){
    ncells_high_exp <- sum(x >= 3)
    return(ncells_high_exp)
  }

vals1 <- apply(CountsEMT, 1, filt_func)

#Genes need to have a count of at least 3 in at least 5 cells to pass filtering
filtered_genes <- rownames(CountsEMT)[vals1 > 10]
  
if(sum(rownames(CountsEMT)!=rownames(CountsNoEMT))!=0){
  stop("The EM and NoEM counts genes are not in the same order, fix this to be safe")
}

if(sum(colnames(CountsEMT)!=colnames(CountsNoEMT))!=0){
  stop("The EM and NoEM counts cells are not in the same order, fix this to be safe")
}

RSumEM <- data.frame(rowSums(CountsEMT))
colnames(RSumEM) <- "rowSumEM"
RSumEM$gene_id <- rownames(RSumEM)

RSumNoEM <- data.frame(rowSums(CountsNoEMT))
colnames(RSumNoEM) <- "rowSumNoEM"
RSumNoEM$gene_id <- rownames(RSumNoEM)

RSumAll <- merge(RSumEM, RSumNoEM, by = "gene_id")
RSumAll$Diff <- RSumAll$rowSumEM - RSumAll$rowSumNoEM

#Add .01 to each rowSum to be able to accurately compute a log ratio
RSumAll$rowSumEM <- RSumAll$rowSumEM + .01
RSumAll$rowSumNoEM <- RSumAll$rowSumNoEM + .01

RSumAll$Log2Ratio <- log2(RSumAll$rowSumEM/RSumAll$rowSumNoEM)

#Restrict to genes with a different enough log2 ratio for the plot and that have a count of at least 3 in at least 10 cells
RSumPlotT <- subset(RSumAll, abs(RSumAll$Log2Ratio) > log2(1.5) & RSumAll$gene_id%in%filtered_genes)
RSumPlotTEM <- RSumPlotT[,c("gene_id", "rowSumEM")]
colnames(RSumPlotTEM) <- c("gene_id", "RawCount")
RSumPlotTEM$type <- "EM"

RSumPlotTNoEM <- RSumPlotT[,c("gene_id", "rowSumNoEM")]
colnames(RSumPlotTNoEM) <- c("gene_id", "RawCount")
RSumPlotTNoEM$type <- "NoEM"

RSumPlot <- rbind(RSumPlotTNoEM, RSumPlotTEM)
RSumPlot$Log2TotalCount <- log2(RSumPlot$RawCount)
RSumPlot$Log10TotalCount <- log10(RSumPlot$RawCount)

RSumPlot$type <- factor(RSumPlot$type, levels = c("NoEM", "EM"))
#RSumPlot$RawCount[RSumPlot$RawCount==0] <- .01
#RSumPlot$type[RSumPlot$type=="NoEM"] <- 0
#RSumPlot$type[RSumPlot$type=="EM"] <- 1

violinp <- ggplot(RSumPlot, aes(x=type, y=RawCount))+geom_violin()+ geom_jitter(height=0,width=.05,alpha = 0.3,size=rel(2.5))+geom_line(aes(group = RSumPlot$gene_id), alpha = 0.25)+
  scale_y_log10()+ylab("Total Gene Count")+xlab("Type")+ggtitle("")+
  theme(axis.text.y=element_text(size=rel(2), face = "bold"),
        axis.text.x=element_text(size=rel(2), face = "bold"), plot.title=element_text(size=rel(1.75), hjust=.5,face = "bold"), axis.title=element_text(size=rel(2), face = "bold"),
        legend.text = element_text(size=rel(1.25), face = "bold"),legend.title = element_text(size=rel(1.5), face = "bold"))

fil_piece <- "EMNoEMViolinPlot"
fil_name <- paste0(PlotsSaveDir, fil_piece, ".tex")
tikz(file = fil_name, height=5.5, width=8.5, standAlone = stdal, sanitize = F)
print(violinp)
dev.off()
tryCatch(dev.off(), error = function(x){})
SVBCompiletikzPlot(PlotsSaveDir, fil_piece, pdflatex_loc = pdflatex_loc)
rm(fil_name)
rm(fil_piece)



stop()



CountsEM <- as.matrix(CountsEMT)
CountsNoEM <- as.matrix(CountsNoEMT)
  
Counts8.0EM <- assays(QuantObj8.0ToUse2)$counts
Counts8.25EM <- assays(QuantObj8.25ToUse2)$counts
Counts8.5EM <- assays(QuantObj8.5ToUse2)$counts




Counts8.0EM2 <- as.matrix(Counts8.0EM[,colnames(Counts8.0NoEM2)])
Counts8.25EM2 <- as.matrix(Counts8.25EM[,colnames(Counts8.25NoEM2)])
Counts8.5EM2 <- as.matrix(Counts8.5EM[,colnames(Counts8.5NoEM2)])

RMeans8.0EMT <- rowMeans(Counts8.0EM2)

RMeans8.0NoEMT <- rowMeans(Counts8.0NoEM2)
  
RMeans8.0EM <- data.frame(RMeans8.0EMT)
colnames(RMeans8.0EM) <- "MeanCount"
RMeans8.0EM$type <- "EM"
RMeans8.0EM$gene_id <- rownames(RMeans8.0EM)

RMeans8.0NoEM <- data.frame(RMeans8.0NoEMT)
colnames(RMeans8.0NoEM) <- "MeanCount"
RMeans8.0NoEM$type <- "NoEM"
RMeans8.0NoEM$gene_id <- rownames(RMeans8.0NoEM)



#Make sure gene names are in the proper order
ordered_gene_names <- gtools::mixedsort(RMeans8.0EM$gene_id)
RMeans8.0EM2 <- RMeans8.0EM[ordered_gene_names,]
RMeans8.0NoEM2 <- RMeans8.0NoEM[ordered_gene_names,]

Log2Ratio <- data.frame(log2(RMeans8.0EM2$MeanCount/RMeans8.0NoEM2$MeanCount))
colnames(Log2Ratio) <- "Log2Ratio"
Log2Ratio$gene_id <- rownames(Log2Ratio)
RowsToKeep <- Log2Ratio$gene_id[Log2Ratio$Log2Ratio]
stop()

lines(t1$time[1:1000], t1$yhat[1:1000], col = "dodgerblue", lwd = 2)
lines(t1$time[1001:2000], t1$yhat[1001:2000], col = "goldenrod4", lwd = 2)
lines(t2$time[1:1000], t2$yhat[1:1000], col = "black", lwd = 2)
lines(t2$time[1001:2000], t2$yhat[1001:2000], col = "red", lwd = 2)

plotInfReps(QuantObj, "ENSMUSG00000059040.5", x = "PT", cov = "Lin", xlab = "pseudotime", legend = T, main = paste0("ENSMUSG00000059040.5", " (Eno1b)"))
t1 <- getCurvesMike(sce_res[[2]], counts_res[[2]], gene = "ENSMUSG00000059040.5", nPoints = 1000)
lines(t1$time[1:1000], t1$yhat[1:1000], col = "dodgerblue", lwd = 2)
lines(t1$time[1001:2000], t1$yhat[1001:2000], col = "goldenrod4", lwd = 2)

plotInfReps(QuantObj, "ENSMUSG00000109350.1", x = "PT", cov = "Lin", xlab = "pseudotime", legend = T, ylim = c(0,150))
t1 <- getCurvesMike(sce_res[[1]], counts_res[[1]], gene = "ENSMUSG00000109350.1", nPoints = 1000)
lines(t1$time[1:1000], t1$yhat[1:1000], col = "dodgerblue", lwd = 2)
lines(t1$time[1001:2000], t1$yhat[1001:2000], col = "goldenrod4", lwd = 2)


#Eno1b (ENSMUSG00000059040.5) is a good example
#Also #Gm44805 (ENSMUSG00000109350.1)
t1 <- getCurvesMike(sce_res[[1]], counts_res[[1]], gene = "ENSMUSG00000059040.5", nPoints = 1000)
t2 <- getCurvesMike(sce_res[[2]], counts_res[[2]], gene = "ENSMUSG00000059040.5", nPoints = 1000)
head(t1)

plot(t1$time, t1$yhat, col = t1$lineage, xlab = "Pseudotime", ylab = "YHat")
plot(t2$time, t2$yhat, col = t2$lineage, xlab = "Pseudotime", ylab = "YHat")







t1 <- getCurvesMike(sce_res[[1]], counts_res[[1]], gene = "ENSMUSG00000059040.5", nPoints = 500)

# 
# #Some plots of the trajectories for some of the genes that may be of interest with the biggest differences
# 
# #Genes Mike had identified
# #Eno1b (ENSMUSG00000059040.5)
# print(plot1 <- plotSmoothers(sce_res[[1]], counts = assays(sce_res[[1]])$counts, gene = "ENSMUSG00000059040.5"))
# print(plot2 <- plotSmoothers(sce_res[[2]], counts = assays(sce_res[[2]])$counts, gene = "ENSMUSG00000059040.5"))
# 
# #Taf9 (ENSMUSG00000052293.14)
# print(plot1 <- plotSmoothers(sce_res[[1]], counts = assays(sce_res[[1]])$counts, gene = "ENSMUSG00000052293.14"))
# print(plot2 <- plotSmoothers(sce_res[[2]], counts = assays(sce_res[[2]])$counts, gene = "ENSMUSG00000052293.14"))
# 
# #Some genes that have an exceptionally big difference between EM and NoEM results
# #Gm44805 (ENSMUSG00000109350.1) (This is an especially good one for illustrating a big difference between EM and noEM)
# print(plot1 <- plotSmoothers(sce_res[[1]], counts = assays(sce_res[[1]])$counts, gene = "ENSMUSG00000109350.1"))
# print(plot2 <- plotSmoothers(sce_res[[2]], counts = assays(sce_res[[2]])$counts, gene = "ENSMUSG00000109350.1"))
# 
# #Pebp1 (ENSMUSG00000032959.12)
# #This is one is an especially good example of a general association trend, as expression dips way down and comes
# #back up somewhat for the EM counts but not at all for the non EM counts
# print(plot1 <- plotSmoothers(sce_res[[1]], counts = assays(sce_res[[1]])$counts, gene = "ENSMUSG00000032959.12"))
# print(plot2 <- plotSmoothers(sce_res[[2]], counts = assays(sce_res[[2]])$counts, gene = "ENSMUSG00000032959.12"))
# 
# #Eno1b
# 
# 
# #Taf9
# 
# 
# plotSmoothers(sce_res[[1]], counts = assays(sce_res[[1]])$counts, gene = "ENSMUSG00000029838.11")
# plotSmoothers(sce_res[[2]], counts = assays(sce_res[[2]])$counts, gene = "ENSMUSG00000029838.11")
# 
# plotSmoothers(sce_res[[1]], counts = assays(sce_res[[1]])$counts, gene = "ENSMUSG00000063524.13")
# plotSmoothers(sce_res[[2]], counts = assays(sce_res[[2]])$counts, gene = "ENSMUSG00000063524.13")
# 
# 
# 
# 
# 
# 
# print(plot1 + scale_colour_manual(values = c("red", "blue")))




