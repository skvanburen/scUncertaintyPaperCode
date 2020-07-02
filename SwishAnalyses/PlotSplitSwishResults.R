#Plot split swish results
onCluster <- FALSE

if(onCluster==TRUE){
  #Specify libPaths to ensure the correct one is used
  .libPaths("/nas/longleaf/home/skvanbur/lib64/R/library")
  source("~/SingleCellProject/SingleCellProjectFunctions.R")
  def_wd2 <- "/pine/scr/s/k/skvanbur/SingleCellProject/SimulationResults/"
}else{
  source('/Users/Scott/Documents/Dissertation/SingleCellProject/code/Simulation Code/SingleCellProjectFunctions.R')
  def_wd2 <- "/Users/Scott/Documents/Dissertation Data/SingleCellProject/SimulationResults/"
}


library(fishpond)
library(tximeta)
library(SummarizedExperiment) # load explicitly to allow mcols function to work
library(readr)
library(gtools)
library(iCOBRA)
library(dplyr)
library(ggplot2)
library(tikzDevice)
sim_number <- 11001

PlotsSaveDir <- "/Users/Scott/Documents/Dissertation/SingleCellProject/SwishPlots/"
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

save_dir <- paste0(def_wd2, "Sim", sim_number, "/")

SwishOutDir <- paste0(save_dir, "SwishFiles/")
if(!dir.exists(SwishOutDir)){dir.create(SwishOutDir)}

SplatterSaveDir <- paste0(save_dir, "SplatterFiles/")
#Load condition (Group) information and call it Group
load(file = paste0(SplatterSaveDir, "simObject.RData"))

SwishObject <- loadRData(file.path(SwishOutDir, "SwishObjectSplit.RData"))

SwishObject2 <- addStatsFromCSV(SwishObject, file.path(SwishOutDir, "summary.csv"))

pvalsT <- data.frame(rowData(SwishObject2)[,"pvalue", drop = FALSE])
colnames(pvalsT) <- "splitSwish"


adjpvalsT <- data.frame(apply(pvalsT, 2, function(x){p.adjust(x, method = "fdr")}))

pvalsT$gene_id <- rownames(pvalsT)
adjpvalsT$gene_id <- rownames(adjpvalsT)

#Now, add results from running swish without splitting into parts

SwishFullResults <- loadRData(file.path(SwishOutDir, "SwishFullResults.RData"))
pvalsTFull <- data.frame(rowData(SwishFullResults)[,"pvalue", drop = FALSE])
colnames(pvalsTFull) <- "swish"

adjpvalsTFull <- data.frame(apply(pvalsTFull, 2, function(x){p.adjust(x, method = "fdr")}))

pvalsTFull$gene_id <- rownames(pvalsTFull)
adjpvalsTFull$gene_id <- rownames(adjpvalsTFull)

pvals <- merge(pvalsT, pvalsTFull, by = "gene_id")
adjpvals <- merge(adjpvalsT, adjpvalsTFull, by = "gene_id")
rownames(pvals) <- pvals$gene_id
rownames(adjpvals) <- pvals$gene_id
pvals$gene_id <- NULL
adjpvals$gene_id <- NULL

TruthT <- data.frame(rowData(SplatterSimObject)[,c("de"), drop = F])
TruthT$status <- as.numeric(TruthT$de)
TruthT$de <- NULL

Truth <- TruthT[gtools::mixedorder(rownames(TruthT)), , drop = FALSE]
Truth$gene_id <- rownames(Truth)

#Now, calculate and add InfRV information
CoverageRes <- loadRData(file.path(save_dir, "CoverageResDropMinnowSkippedGenes.RData"))
InfRVVals <- CoverageRes$gene_level_coverage[,"MeanInfRV", drop = F]
InfRVVals2 <- subset(InfRVVals, rownames(InfRVVals)%in%rownames(pvals))
InfRVQuan <- quantile(InfRVVals2$MeanInfRV, probs = c((1/3), 2/3))
InfRVVals2$InfRV <- cut(InfRVVals2$MeanInfRV,
                       breaks=c(0.01, InfRVQuan[1], InfRVQuan[2], 1),
                       include.lowest=TRUE)
InfRVVals2$gene_id <- rownames(InfRVVals2)
Truth2 <- merge(Truth, InfRVVals2, by = "gene_id")
rownames(Truth2) <- Truth2$gene_id
Truth2$gene_id <- NULL

cobraDat <- COBRAData(pval = pvals, padj = adjpvals, truth = Truth2)

#cobraPerfOverall <- calculate_performance(cobraDat, binary_truth = "status")
cobraPerfInfRV <- calculate_performance(cobraDat, binary_truth = "status", splv = "InfRV", maxsplit = Inf)
colorss <- c("black", "blue")
cobraResForPlotting <- prepare_data_for_plot(cobraPerfInfRV, colorscheme = colorss)
#ggtitle("Comparison of Swish and SplitSwish Results by InfRV")
SwishResPlot <- plot_fdrtprcurve(cobraResForPlotting, xaxisrange = c(0,.15), plottype = "points")+ggtitle("")+
  theme(legend.text = element_text(size=rel(1), face = "bold"),legend.title = element_text(size=rel(1), face = "bold"),axis.text.y=element_text(size=rel(1.25), face = "bold"),
        axis.text.x=element_text(size=rel(1.25), face = "bold"), plot.title=element_text(size=rel(1.25), hjust=.5,face = "bold"), axis.title=element_text(size=rel(1.25), face = "bold"),
        strip.text = element_text(size = rel(1.10)))+facet_wrap(~splitval, nrow = 2)

fil_piece <- "SwishVsSplitSwishICOBRAPlot"
fil_name <- paste0(PlotsSaveDir, fil_piece, ".tex")

tikz(file = fil_name, height=5.5, width=8.5, standAlone = TRUE, sanitize = F)
print(SwishResPlot)
dev.off()
pdflatex_loc <- '/Library/TeX/texbin/pdflatex'
SVBCompiletikzPlot(PlotsSaveDir, fil_piece, pdflatex_loc)
