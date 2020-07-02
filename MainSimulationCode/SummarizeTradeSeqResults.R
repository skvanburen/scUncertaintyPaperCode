onCluster <- TRUE
if(onCluster==TRUE){
   
   #Specify libPaths to ensure the correct one is used
   .libPaths("/nas/longleaf/home/skvanbur/lib64/R/library")
   
   source('~/SingleCellProject/SingleCellProjectFunctions.R')
   base_dir <- "/pine/scr/s/k/skvanbur/SingleCellProject/SimulationResults/"
   save_dir <- "/pine/scr/s/k/skvanbur/SingleCellProject/tradeSeqiCobraRes/"
   SimNumberToUse <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
}

#Summarize tradeSeq Results

library(dplyr)
library(iCOBRA)
library(compositions)
library(psych)

useKnownClusters <- TRUE

#Load coverage results to extract information about gene uniqueness and InfRV (averaged across cells)
CoverageResT <- loadRData(paste0(base_dir, "Sim", SimNumberToUse, "/", "CoverageRes.RData"))
CoverageResT2 <- CoverageResT$gene_level_coverage
CoverageResToUse <- CoverageResT2[,c("gene_id", "MeanInfRV", "uniqueness_ratio", "UniqFactor")]

TierRes <- loadRData(paste0(base_dir, "Sim", SimNumberToUse, "/", "AlevinTierInfo.RData"))

#Read in and summarize the significance results for the tradeSeq results
if(SimNumberToUse%in% c(10117,10118,11001)){
   nboot <- 20
}else{
   nboot <- 100
}
 
   SimObject <- loadRData(paste0(base_dir, "Sim", SimNumberToUse, "/SplatterFiles/simObject.RData"))
   GeneLevelDEParamsT <- SimObject$tde_overall
   
   #Now, merge in the genenames to be used that are not just G1, G2, etc
   ENSGNames <- read.table(file = paste0(base_dir, "Sim", SimNumberToUse, "/", "SplatterFiles/", "quants_mat_rows.txt"))
   GeneLevelDEParamsT$gene_id <- ENSGNames$V1
   GeneLevelDEParamsT$feature_id <- NULL
   
   GeneLevelDEParams <- data.frame(GeneLevelDEParamsT, stringsAsFactors = F)
   FilterNoGenes <- TRUE
   
   #Now, get a list of "filtered_genes" that are genes with sufficient count to keep in the analysis
   
   Counts <- read.csv(file = paste0(base_dir, "Sim", SimNumberToUse, "/", "SplatterFiles/", "quants_mat.csv"), header = FALSE)
   Counts2 <- data.frame(Counts)
   rownames(Counts2) <- ENSGNames$V1
   filt_func <- function(x){
      ncells_high_exp <- sum(x >= 10)
      return(ncells_high_exp)
   }
   rows_to_keep <- apply(Counts2, 1, filt_func)
   genes_to_keep <- rows_to_keep > 10
   genenames_to_keep <- rownames(Counts2)[genes_to_keep]
   
   filtered_genes <- genenames_to_keep


iCobraSaveMainDir <- save_dir

   DataTypeNum <- 1
   returnPvalsForiCobra(SimNumberToUse = SimNumberToUse, FilteringtypeNum = FilteringtypeNum, nboot = nboot, CoverageResToUse = CoverageResToUse, 
                                         UseTop3000Genes = UseTop3000Genes, base_dir = base_dir, FilterNoGenes = FilterNoGenes, TierRes = TierRes, GeneLevelDEParams = GeneLevelDEParams,
                                         iCobraSaveMainDir = iCobraSaveMainDir, useKnownClusters = useKnownClusters, filtered_genes = filtered_genes)
   print(gc())
   DataTypeNum <- 2
   returnPvalsForiCobra(SimNumberToUse = SimNumberToUse, FilteringtypeNum = FilteringtypeNum, nboot = nboot, CoverageResToUse = CoverageResToUse, 
                        UseTop3000Genes = UseTop3000Genes, base_dir = base_dir, FilterNoGenes = FilterNoGenes, TierRes = TierRes, GeneLevelDEParams = GeneLevelDEParams,
                        iCobraSaveMainDir = iCobraSaveMainDir, useKnownClusters = useKnownClusters, filtered_genes = filtered_genes)

print(gc())



