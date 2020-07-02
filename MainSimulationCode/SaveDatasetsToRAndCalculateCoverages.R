#Specify libPaths to ensure the correct one is used
.libPaths("/nas/longleaf/home/skvanbur/lib64/R/library")

array_val <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
print(paste0("The current array value is ", array_val))
print(paste0("The current slurm job id is ", as.numeric(Sys.getenv("SLURM_ARRAY_JOB_ID"))))



#Set to FALSE if data has already been saved so no need to waste the time doing it again
RunSaveAlevinDataToRData <- TRUE
RegenerateNewInfReps <- TRUE

CalculateCoverageRes <- TRUE
CalculateMainCoverageRes <- TRUE
#Whether to remove genes that minnow is unable to simulate reads from
DropMinnowSkippedGenesFromCoverageCalc <- TRUE

#Calculate coverage without including zero values
CalculateCoveragesWithoutZeros <- TRUE


library(readr)
library(tximport)
library(fishpond)
library(data.table)
library(gtools)
library(MASS)
#library(BiocParallel)
source("~/SingleCellProject/SingleCellProjectFunctions.R")


if(packageVersion("tximport") < "1.15.13"){
  stop("Please update tximport to be at least version 1.15.13 to properly import the results")
}
sim_number <- array_val

if(sim_number%in%c(10117,10118,11001)){
  nboot <- 20
}else{
  nboot <- 100
}

#Set seed to simulate the pseudo-inferential replicates
GenerateInfRepSeed <- sim_number * 12316

gencode_files_location <- "~/SingleCellProject/gencode_v32_files/"
def_wd2 <- "/pine/scr/s/k/skvanbur/SingleCellProject/SimulationResults/"
save_dir <- paste0(def_wd2, "Sim", sim_number, "/")

if(!dir.exists(save_dir)){dir.create(save_dir, recursive = TRUE)}

SplatterSaveDir <- paste0(save_dir, "SplatterFiles/")
if(!dir.exists(SplatterSaveDir)){dir.create(SplatterSaveDir, recursive = TRUE)}

MinnowSaveDir <- paste0(save_dir, "MinnowOutput/")
if(!dir.exists(MinnowSaveDir)){dir.create(MinnowSaveDir, recursive = TRUE)}

AlevinSaveDir <- paste0(save_dir, "AlevinOutput/")
if(!dir.exists(AlevinSaveDir)){dir.create(AlevinSaveDir, recursive = TRUE)}

if((sim_number >=100 & sim_number < 999) | sim_number%in%c(10117,10118)){
  #This code is for dynverse simulations
  load(file = paste0(SplatterSaveDir, "simObject.RData"))
  numCells <- dim(dataset$counts)[1]
}else{
  #This code is for splatter simulations
  load(file = paste0(SplatterSaveDir, "simObject.RData"))
  numCells <- dim(SplatterSimObject)[2]
}

if(numCells > 1000){
  stop("Check the specification of the number of cells, it appears to be incorrect")
}


##########################################################################################
#Now, use tximport to save the results of the Alevin quantification into R data files

alevin_fil <- paste0(AlevinSaveDir, "alevin/quants_mat.gz")
alevin_resDropInfRepsT1 <- proc.time()
alevin_res <- tximport(alevin_fil, type = "alevin", dropInfReps = TRUE, alevinArgs = list(filterBarcodes = FALSE, tierImport = TRUE))
alevin_resDropInfRepsT2 <- proc.time() - alevin_resDropInfRepsT1

TimeToSaveAlevinResDropInfReps1 <- proc.time()
save(alevin_res, file = paste0(save_dir, "AlevinDataNoInfRepsSim", sim_number, ".RData"))
TimeToSaveAlevinResDropInfReps2 <- proc.time() - TimeToSaveAlevinResDropInfReps1

TimeToSaveAlevinResDropInfRepsToDisk <- TimeToSaveAlevinResDropInfReps2[3]

DropInfRepsAlevinResReadInTime <- alevin_resDropInfRepsT2[3]
DropInfRepsAlevinResMemorySizeBytes <- as.numeric(object.size(alevin_res))

DropInfRepsAlevinResSizeOnDisk <- file.size(paste0(save_dir, "AlevinDataNoInfRepsSim", sim_number, ".RData"))

save(DropInfRepsAlevinResReadInTime, DropInfRepsAlevinResMemorySizeBytes,TimeToSaveAlevinResDropInfRepsToDisk, DropInfRepsAlevinResSizeOnDisk,
     file = paste0(save_dir, "AlevinDataNoInfReps", "SizeAndComputationTimeRes.RData"))

#The function SaveAlevinDataToRData will save the Alevin data object (with and without all bootstrap replicates) and will save separate files 
  #for each cell/replicate for potential use later
if(RunSaveAlevinDataToRData==TRUE){
  SaveAlevinDataToRData(alevin_fil = alevin_fil, save_dir = save_dir, ReSaveAllAlevinFiles = FALSE)
}

##########################################################################################

print(gc())

##########################################################################################
#Save Counts Output by Alevin and counts from Splatter for easy access later
alevin_res <- loadRData(paste0(save_dir, "AlevinDataNoInfRepsSim", sim_number, ".RData"))
Alevin_countsT <- alevin_res$counts

ordered_cell_names <- gtools::mixedsort(colnames(Alevin_countsT))
ordered_gene_names <-  gtools::mixedsort(rownames(Alevin_countsT))


Alevin_countsT2 <- Alevin_countsT[,ordered_cell_names]
Alevin_counts <- as.matrix(Alevin_countsT2[ordered_gene_names,])

save(Alevin_counts, file = paste0(save_dir, "AlevinCounts.RData"))


#Now save the counts output by Splatter, both as output by Splatter and with counts scaled
#such that cells have the same expected lirary size in each case
#The library size for the Alevin counts will always be smaller, so this will serve to scale down the Splatter
#counts to make the counts from Alevin and Splatter more comparable to each other
#Verify the actual counts output by alevin are similar to the expected counts (ie the values simulated by splatter)
#These are the splatter simulations
if(sim_number < 100 | sim_number %in% c(1001)){
  Splatter_countsT <- counts(SplatterSimObject)
  
  if(ncol(Splatter_countsT)!=numCells){
    stop("Check the dimensions of the splatter counts")
  }
}else{
  #These are for the dynverse simulations
  #Load the counts from here now from the actual dataset simulation object because the counts were scaled to have an average library size of
    #1,000,000 and the counts in the object are not scaled to have this library size
  Splatter_countsT <- read.table(file = paste0(SplatterSaveDir, "quants_mat.csv"), quote="", sep=",")
  if(ncol(Splatter_countsT)!=numCells){
    stop("Check the dimensions of the splatter counts")
  }
  row_names <- read.table(paste0(SplatterSaveDir, "quants_mat_rows.txt"))
  rownames(Splatter_countsT) <- row_names$V1
}


col_names <- read.table(paste0(MinnowSaveDir, "alevin/quants_mat_rows.txt"))
col_names2 <- read.table(paste0(MinnowSaveDir, "alevin/true_cell_names.txt"))
col_names3 <- read.table(paste0(SplatterSaveDir, "quants_mat_cols.txt"))
colnames(Splatter_countsT) <- as.character(col_names$V1)

Splatter_countsT2 <- Splatter_countsT[,ordered_cell_names]
Splatter_counts <- Splatter_countsT2[ordered_gene_names,]

save(Splatter_counts, file = paste0(save_dir, "SplatterCounts.RData"))

Lib_Sizes_Splatter <- colSums(Splatter_counts)
Lib_Sizes_Alevin <- colSums(Alevin_counts)
Lib_Sizes_Mult_Factor <- Lib_Sizes_Alevin/Lib_Sizes_Splatter

Splatter_counts_scaledT <- Splatter_counts
for(j in 1:numCells){
  Splatter_counts_scaledT[,j] <- Splatter_counts_scaledT[,j] * Lib_Sizes_Mult_Factor[j]
}
Splatter_counts_scaled <- Splatter_counts_scaledT
#Verify the new "library sizes" are the same
if(!all.equal(colSums(Splatter_counts_scaled), Lib_Sizes_Alevin)){
  stop("Scaling the Splatter Counts to have the same library size as the alevin counts didn't work properly")
}

save(Splatter_counts_scaled, file = paste0(save_dir, "SplatterCountsScaledToAlevinLibSize.RData"))

#Now, save the Tier information from alevin (which is saved per cell/gene)
TierInfoT <- alevin_res$tier

TierInfoT2 <- TierInfoT[,ordered_cell_names]
TierInfo <- as.matrix(TierInfoT2[ordered_gene_names,])

save(TierInfo, file = paste0(save_dir, "AlevinTierInfo.RData"))

##########################################################################################
#Calculate coverage results for the current quantifications, both at the cell level and gene levels
uniqueness_fil_dir <- gencode_files_location
if(CalculateCoverageRes==TRUE){

  if(CalculateMainCoverageRes==TRUE){
    CoverageRes <- SummarizeBootstrapSamples(SplatterSaveDir, MinnowSaveDir, uniqueness_fil_dir, save_dir, useSplatterCountsScaledToAlevinLibrarySize = TRUE, sim_number = sim_number,
                                             DropMinnowSkippedGenesFromCoverageCalc = FALSE)
    
    save(CoverageRes, file = paste0(save_dir, "CoverageRes.RData"))
  }
  
  #Results that drop the genes that minnow has to skip due to having too short of a read length
  if(DropMinnowSkippedGenesFromCoverageCalc==TRUE){
    CoverageRes <- SummarizeBootstrapSamples(SplatterSaveDir, MinnowSaveDir, uniqueness_fil_dir, save_dir, useSplatterCountsScaledToAlevinLibrarySize = TRUE, sim_number = sim_number,
                                             DropMinnowSkippedGenesFromCoverageCalc = TRUE)
    
    save(CoverageRes, file = paste0(save_dir, "CoverageResDropMinnowSkippedGenes.RData"))
  }
  
  #All results with zeros removed from coverage results also drop the minnow skipped genes
  if(CalculateCoveragesWithoutZeros==TRUE){
    CoverageRes <- SummarizeBootstrapSamples(SplatterSaveDir, MinnowSaveDir, uniqueness_fil_dir, save_dir, useSplatterCountsScaledToAlevinLibrarySize = TRUE, sim_number = sim_number,
                                             DropMinnowSkippedGenesFromCoverageCalc = TRUE, CalculateCoveragesWithoutZeros = TRUE)
    
    save(CoverageRes, file = paste0(save_dir, "CoverageResCalculatedWithoutZeros.RData"))
  }
  
  
  
  
}


##########################################################################################
print(gc())

#Now, simulate pseudo-inferential replicates using the mean and variance values
if(RegenerateNewInfReps==TRUE){
  print(gc())
  boot_samps_means <- loadRData(paste0(save_dir, "BootSampsMeans.RData"))
  boot_samps_variances <- loadRData(paste0(save_dir, "BootSampsVariances.RData"))
  
  TimetoGenerateInfReps1 <- proc.time()
  #Ensure bootstrap means and variances are in the same order
  if(sum(colnames(boot_samps_means)!=colnames(boot_samps_variances))!=0){
    stop("The cells are ordered differently for the mean and variance values")
  }
  
  if(sum(rownames(boot_samps_means)!=rownames(boot_samps_variances))!=0){
    stop("The genes are ordered differently for the mean and variance values")
  }
  
  #Generate theta values for use with the MASS package 
  #Theta values here are for var = mu + (mu/theta) so set a maximum value to ensure sufficient extra poisson variation is properly
  #present 
  #Mike Love had used phi = 0.001 in DESeq2 for phi=1/theta, so to correspond to that here set theta = 1/.001 = 1000 such that we specify a maximum
  #This is to ensure variance is always (at least slightly) greater than the mean
  #Theta may also be called phi, and is also equal to the "size" value from the base R functions
  #So, set theta to 1000 if theta is negative or greater than 1000
  cell_names <- colnames(boot_samps_means)
  genenames <- rownames(boot_samps_means)
  ngenes <- nrow(boot_samps_means)
  ordered_cell_names <- gtools::mixedsort(cell_names)
  
  theta_vals <- matrix(NA, nrow = ngenes, ncol = numCells)
  for(z in 1:numCells){
    curr_means <- boot_samps_means[,z]
    curr_variances <- boot_samps_variances[,z]
    theta_vals[,z] <- (curr_means^2)/(curr_variances - curr_means)
  }
  theta_vals[!is.finite(theta_vals)] <- 1000
  theta_vals[theta_vals > 1000] <- 1000
  theta_vals[theta_vals < 0] <- 1000
  rownames(theta_vals) <- genenames
  
  #Be careful here: These columns so far are not ordered, so needs to be cell_names, not ordered_cell_names
  colnames(theta_vals) <- cell_names
  
  
  
  SimulatedInfRepDataDirec <- paste0(save_dir, "SimulatedInfRepData/")
  if(!dir.exists(SimulatedInfRepDataDirec)){dir.create(SimulatedInfRepDataDirec, recursive = TRUE)}
  
  #Ensure the boot_samp_means and theta_vals vectors are in the same order before you generate inf reps
  if(sum(colnames(boot_samps_means)!=colnames(theta_vals))!=0){
    stop("The cells are ordered differently for the mean and theta values")
  }
  
  if(sum(rownames(boot_samps_means)!=rownames(theta_vals))!=0){
    stop("The genes are ordered differently for the mean and theta values")
  }
  
  ninfrep <- nboot
  set.seed(GenerateInfRepSeed)
  for(i in 1:ninfrep){
    print(paste0("Current infrep number is ", i))
    curr_samp_vals <- matrix(NA, nrow = ngenes, ncol = numCells)
    for(j in 1:numCells){
      p1 <- boot_samps_means[,j]
      #curr_samp_vals[,j] <- rpois(ngenes, p1)
      curr_theta <- theta_vals[,j]
      curr_samp_vals[,j] <- MASS::rnegbin(ngenes, mu = p1, theta = curr_theta)
    }
    #Be careful here: These columns so far are not ordered, so needs to be cell_names, not ordered_cell_names
    colnames(curr_samp_vals) <- cell_names
    rownames(curr_samp_vals) <- genenames
    
    curr_samp_vals2 <- curr_samp_vals[,ordered_cell_names]
    
    save(curr_samp_vals2, file = paste0(SimulatedInfRepDataDirec, "SimulatedInfRepDataRep", i, ".RData"))
  }
  
  TimetoGenerateInfReps2 <- proc.time() - TimetoGenerateInfReps1
  TimetoGenerateAllInfReps <- TimetoGenerateInfReps2[3]
  
  InfRepFileSize <- file.size(paste0(SimulatedInfRepDataDirec, "SimulatedInfRepDataRep", 1, ".RData"))
  
  save(TimetoGenerateAllInfReps, InfRepFileSize, file = paste0(save_dir, "TimeandDiskSizeofInfRepsFiles.RData"))
  print(gc())
}