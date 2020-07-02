#Code to run the tradeSeq Results
onCluster <- TRUE

#Specify libPaths to ensure the correct one is used
.libPaths("/nas/longleaf/home/skvanbur/lib64/R/library")


useKnownClusters <- TRUE
if(onCluster==TRUE){
  if(version$nickname=="Planting of a Tree"){
    .libPaths("/nas/longleaf/home/skvanbur/bin/R3.6.0")
  }
  
  
  array_val <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
  
  Sys.sleep(sample(1:200, 1))
  source("~/SingleCellProject/SingleCellProjectFunctions.R")
  def_wd2 <- "/pine/scr/s/k/skvanbur/SingleCellProject/SimulationResults/"
}

#library(mclust)
library(slingshot)

#Load version 1.1.21 to avoid a bug in the tradeSeq code
#devtools::install_github(repo = "statOmics/tradeSeq", ref = "1aabc6ed558c865c2b322fb51a0213893e731303", lib = "~/InstalledtradeSeq1.1.21/")
library(tradeSeq, lib = "~/InstalledtradeSeq1.1.21/")
library(SingleCellExperiment)
library(mgcv)

if(packageVersion("tradeSeq", lib = "~/InstalledtradeSeq1.1.21/") < "1.1.20"){
  stop("tradeSeq version must be at least 1.1.20 because important fixes have been introduced.")
}

if(packageVersion("tradeSeq", lib = "~/InstalledtradeSeq1.1.21/") > "1.1.21"){
  stop("A bug seems to be introduced by changes in 1.1.23 that try to capture and deal with errors.  Avoid this version of the package for fitting the fitGAM statement but use it when calculating the test results because it does seem to fix other potential issues")
}
sim_number <- floor((array_val - 1)/1000)
Part <- array_val%%1000
if(Part==0){
  Part <- 1000
}

runResRemovingOutliers <- FALSE

#Set sim_number_to_use to ensure correct files are loaded

sim_number_to_use <- sim_number + 10000 + 100

#Randomly choose a seet to set before running the fitGAM step
current_seed <- (3*array_val) - 12


save_dir <- paste0(def_wd2, "Sim", sim_number_to_use, "/")
print(paste0("Current simulation directory results are being loaded from is ", save_dir))


SplatterSaveDir <- paste0(save_dir, "SplatterFiles/")

MinnowSaveDir <- paste0(save_dir, "MinnowOutput/")

tradeSeqSaveDir <- paste0(save_dir, "tradeSeqResults/")
if(!dir.exists(tradeSeqSaveDir)){dir.create(tradeSeqSaveDir, recursive = TRUE)}

SimObject <- loadRData(paste0(SplatterSaveDir, "simObject.RData"))
  
  Cell_groupsT <- SimObject$prior_information$groups_id
  Cell_groupsT2 <- Cell_groupsT[gtools::mixedorder(Cell_groupsT$cell_id),]
  Cell_groupsT3 <- Cell_groupsT2$group_id
  names(Cell_groupsT3) <- paste0("Cell", 1:length(Cell_groupsT3))
  
  DyntoyStartMilestones <- SimObject$prior_information$start_milestones
  DyntoyEndMilestones <- SimObject$prior_information$end_milestones
  
  
  Minnow_Orig_cell_names <- read.table(paste0(MinnowSaveDir, "alevin/true_cell_names.txt"))
  Minnow_Assigned_cell_names <- read.table(paste0(MinnowSaveDir, "alevin/quants_mat_rows.txt"))
  
  if(nrow(Minnow_Orig_cell_names)!=nrow(Minnow_Assigned_cell_names)){
    stop("The number of rows of the original cell names from Splatter and from the assigned names from minnow don't match even though they should, figure out why")
  }
  
  Minnow_full_cell_info <- cbind(Minnow_Orig_cell_names, Minnow_Assigned_cell_names, stringsAsFactors = FALSE)
  colnames(Minnow_full_cell_info) <- c("Cell", "CellBarcode")
  
  if(sum(names(Cell_groupsT3)!=Minnow_full_cell_info$Cell)==0){
    KnownClusters <- Cell_groupsT3
    names(KnownClusters) <- Minnow_full_cell_info$CellBarcode
  }else{
    stop()
  }


Splatter_countsT <- loadRData(paste0(save_dir, "SplatterCountsScaledToAlevinLibSize.RData"))
Splatter_countsUnscaledT <- loadRData(paste0(save_dir, "SplatterCounts.RData"))
Alevin_countsT <- loadRData(paste0(save_dir, "AlevinCounts.RData"))


#Keep the list of genes to keep the same between different datasets and choose based on the splatter counts (unscaled to Alevin library size)
genenames_to_keep <- rownames(Splatter_countsUnscaledT)



Splatter_counts <- Splatter_countsT[genenames_to_keep,]
Splatter_countsUnscaled <- Splatter_countsUnscaledT[genenames_to_keep,]

Alevin_counts <- Alevin_countsT[genenames_to_keep,]

Mean_Boot_CountsT <- loadRData(paste0(save_dir, "BootSampsMeans.RData"))
Mean_Boot_CountsT2 <- Mean_Boot_CountsT[,colnames(Alevin_counts)]
Mean_Boot_Counts <- Mean_Boot_CountsT2[genenames_to_keep,]


InfRepCountsT <- loadRData(paste0(save_dir, "alevinInfRepData/InfRepLevelFiles/", "infRepDatRep", Part, ".RData"))

#The colnames of the Alevin_counts are in order, put these in order even though it is not critical for now because of the fact that everything comes from the same path
InfRepCountsT2 <- InfRepCountsT[,colnames(Alevin_counts)]

InfRepCounts <- InfRepCountsT2[genenames_to_keep,]


InfRepCountsSimulatedT <- loadRData(paste0(save_dir, "SimulatedInfRepData/", "SimulatedInfRepDataRep", Part, ".RData"))

InfRepCountsSimulatedT2 <- InfRepCountsSimulatedT[,colnames(Alevin_counts)]
InfRepCountsSimulated <- InfRepCountsSimulatedT2[genenames_to_keep,]


#Run the results corresponding to the non-InfRep results under array value 1
  if(Part==1){
    
    FulltradeSeqPipeline(tradeSeqSaveDir = tradeSeqSaveDir, InputCounts =  Splatter_countsUnscaled, Part = NULL, CountsFromSplatter = "Unscaled",
                         current_seed = current_seed, useKnownClusters = useKnownClusters, KnownClusters = KnownClusters,
                         DyntoyStartMilestones = DyntoyStartMilestones, DyntoyEndMilestones = DyntoyEndMilestones)
    
    print(gc())

    
    
    #Now, Alevin Counts
    FulltradeSeqPipeline(tradeSeqSaveDir = tradeSeqSaveDir, InputCounts =  Alevin_counts, Part = NULL, CountsFromSplatter = FALSE,
                         current_seed = current_seed, useKnownClusters = useKnownClusters, KnownClusters = KnownClusters,
                         DyntoyStartMilestones = DyntoyStartMilestones, DyntoyEndMilestones = DyntoyEndMilestones)
    
    print(gc())
    

    
  }
  
  
  #Now run tradeSeq on the infReps, both actual ones drawn from alevin and the simulated ones using Negative Binomial
  #For now, the genes used are not restricted to be the same as using the regular point estimates 


FulltradeSeqPipeline(tradeSeqSaveDir = tradeSeqSaveDir, InputCounts =  InfRepCounts, Part = Part, CountsFromSplatter = FALSE,
                     current_seed = current_seed, SimulatedData = FALSE, useKnownClusters = useKnownClusters, KnownClusters = KnownClusters,
                     DyntoyStartMilestones = DyntoyStartMilestones, DyntoyEndMilestones = DyntoyEndMilestones)
print(gc())
  
  
  #Run and Save results removing potential outliers for the simulated data based on mean/var of infreps
FulltradeSeqPipeline(tradeSeqSaveDir = tradeSeqSaveDir, InputCounts =  InfRepCountsSimulated, Part = Part, CountsFromSplatter = FALSE,
                     current_seed = current_seed, SimulatedData = TRUE, useKnownClusters = useKnownClusters, KnownClusters = KnownClusters,
                     DyntoyStartMilestones = DyntoyStartMilestones, DyntoyEndMilestones = DyntoyEndMilestones)

print(gc())




