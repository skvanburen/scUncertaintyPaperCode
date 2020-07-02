#File that contains code to plot tradeSeq results

#Specify libPaths to ensure the correct one is used
.libPaths("/nas/longleaf/home/skvanbur/lib64/R/library")

useKnownClusters <- TRUE
if(useKnownClusters==TRUE){
  dir_mod4 <- "UsingKnownClusters"
}else{
  dir_mod4 <- ""
}
onCluster <- T
if(onCluster==TRUE){
  source("~/SingleCellProject/SingleCellProjectFunctions.R")
  base_dir <- "/pine/scr/s/k/skvanbur/SingleCellProject/SimulationResults/"
  save_dir <- "/pine/scr/s/k/skvanbur/SingleCellProject/tradeSeqiCobraRes/"
  PlotsSaveDir <- paste0("/pine/scr/s/k/skvanbur/SingleCellProject/tradeSeqPowerPlots", dir_mod4, "/")
}

array_val <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
if(!dir.exists(PlotsSaveDir)){dir.create(PlotsSaveDir, recursive = T)}
library(iCOBRA)
library(ggplot2)
library(tikzDevice)
library(splatter)
setwd(PlotsSaveDir)


if(onCluster==TRUE){
  pdflatex_loc <- '/nas/longleaf/home/skvanbur/texlive/2020/bin/x86_64-linux/pdflatex'
  options(tikzLatex = pdflatex_loc)
}else{
  pdflatex_loc <- 'pdflatex'
}

options("tikzDocumentDeclaration" = "\\documentclass[12pt,letterpaper,twoside]{report}\n")
options("tikzMetricsDictionary" = paste(PlotsSaveDir, "/TikzTemp", array_val, sep = ""))
options("tikzLatexPackages" = c(
  "\\usepackage{tikz}",
  "\\usepackage[active,tightpage,psfixbb]{preview}",
  "\\PreviewEnvironment{pgfpicture}",
  "\\setlength\\PreviewBorder{0pt}",
  "\\usepackage{bm}",
  "\\usepackage{newtxtext}")
)
tikzLockFile <- paste0(PlotsSaveDir, "/TikzTemp", array_val, "___LOCK")
if(file.exists(tikzLockFile)){system(paste0("rm ", tikzLockFile))}

TotalGenes <- 60179

testTypes <- c("Assoc", "Pattern", "DiffEnd", "StartEnd", "EarlyDE")


  
  SimNumberToUse <- array_val
  print(paste0("Current Simulation Number is ", SimNumberToUse))

  for(j in 1:length(testTypes)){
    curr_type <- testTypes[j]
    GenerateiCobraPlots(SimulatedInfReps = TRUE, OutliersRemoved = FALSE, SimNumberToUse = SimNumberToUse, TotalGenes = TotalGenes, StatType = curr_type, useKnownClusters = useKnownClusters, pdflatex_loc = pdflatex_loc)
    print(gc())
    GenerateiCobraPlots(SimulatedInfReps = FALSE, OutliersRemoved = FALSE, SimNumberToUse = SimNumberToUse, TotalGenes = TotalGenes, StatType = curr_type, useKnownClusters = useKnownClusters, pdflatex_loc = pdflatex_loc)
    print(gc())
  }



#c(15,16,19:27,28)
#c(29,30,32)
#c(1:4,15:32,36:39)
# for(ll in 36:39){
#   
#   SimNumberToUse <- ll
#   
#   GenerateiCobraPlots(TRUE, FALSE, SimNumberToUse, TotalGenes = TotalGenes)
#   GenerateiCobraPlots(TRUE, TRUE, SimNumberToUse, TotalGenes = TotalGenes)
# }

# TableAlevinAndSplatterCountstradeSeqResults(SimNumberToUse = SimNumberToUse, base_dir = base_dir, CompareSplatterAndMeanWaldStatResults = TRUE,
#                                             CompareSplatterAndMeanWaldStatResultsCounts = ComparisonOfRes2, pvalue_touse = "PvalMeanWaldStat")
#I had tested out the CompMI methods, but I'm not sure I expect it to work and Mike didn't seem to like these approaches as much
# CombinePvals <- CombineMIPVals(AllPvals = full_tradeSeqres_boot, pvalToUse = "pvalue")
# CombinePvals$pAdjFDRMI <- p.adjust(lol$pvalt, method = "fdr")
# 
# CombinePvals2 <- lol[,c("gene_id", "pvalt")]
# FulltradeSeqResFinalSim9_2 <- merge(FulltradeSeqResFinalSim9, CombinePvals2, by = "gene_id")
# 
# 
# 
# rA <- rowSums(actual_counts)
# rE <- rowSums(expected_counts)
# cA <- colSums(actual_counts)
# cE <- colSums(expected_counts)
# 
# load("/Users/Scott/Documents/Dissertation Data/SingleCellProject/Simulation Results/Sim9/tradeSeqResults/tradeSeqResultsSplatterExpectedCounts.RData")
# filtered_genes <- rownames(association_test_res)
# actual_counts_filtered <- actual_counts[filtered_genes,]
# expected_counts_filtered <- expected_counts[filtered_genes,]
# 
# rAF <- rowSums(actual_counts_filtered)
# rEF <- rowSums(expected_counts_filtered)
# cAF <- colSums(actual_counts_filtered)
# cEF <- colSums(expected_counts_filtered)