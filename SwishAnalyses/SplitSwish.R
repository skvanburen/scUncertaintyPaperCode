#Specify libPaths to ensure the correct one is used
.libPaths("/nas/longleaf/home/skvanbur/lib64/R/library")

#Run split swish using generated pseudo-infReps
library(fishpond)
library(tximeta)
library(SummarizedExperiment) # load explicitly to allow mcols function to work
library(SingleCellExperiment)
library(DESeq2)
library(scran)
sim_number <- 11001

#source("~/SingleCellProject/SingleCellProjectFunctions.R")
def_wd2 <- "/pine/scr/s/k/skvanbur/SingleCellProject/SimulationResults/"

save_dir <- paste0(def_wd2, "Sim", sim_number, "/")
AlevinSaveDir <- paste0(save_dir, "AlevinOutput/")
if(!dir.exists(AlevinSaveDir)){stop("The directory with specified alevin output does not exist")}

MinnowSaveDir <- paste0(save_dir, "MinnowOutput/")

SwishOutDir <- paste0(save_dir, "SwishFiles/")
if(!dir.exists(SwishOutDir)){dir.create(SwishOutDir)}

SplatterSaveDir <- paste0(save_dir, "SplatterFiles/")
#Load condition (Group) information and call it Group
load(file = paste0(SplatterSaveDir, "simObject.RData"))

#Match cell names from Splatter ("Cell1", "Cell2", etc to the actual cell barcodes names assigned by minnow)
CellBarcodeNamesT <- read.table(paste0(MinnowSaveDir, "alevin/quants_mat_rows.txt"))
SplatterCellNamesT <- read.table(paste0(MinnowSaveDir, "alevin/true_cell_names.txt"))
CellBarcodeNames <- cbind(CellBarcodeNamesT, SplatterCellNamesT)
colnames(CellBarcodeNames) <- c("CellBarcode", "SplatterCellName")

SortedCellNames <- CellBarcodeNames[gtools::mixedorder(CellBarcodeNames$CellBarcode),]

Groups <- SplatterSimObject@colData$Group
Groups2 <- cbind(SortedCellNames, Groups)
condition <- Groups2$Group
alevin_fil <- paste0(AlevinSaveDir, "alevin/quants_mat.gz")

coldata <- data.frame(alevin_fil)
colnames(coldata) <- "files"

#Create a fake temporary name for the file to avoid the error in tximeta
coldata$names <- "1"

#Set location for the tximeta files to avoid needing to rerun this everytime
#Also, make sure the tximeta files exist here before running the results for the computation time so they are not 
#included in the results
setTximetaBFC(def_wd2)

sT <- proc.time()
y <- tximeta(coldata, type = "alevin", dropInfReps = TRUE, skipMeta = TRUE, alevinArgs = list(filterBarcodes = FALSE, tierImport = TRUE))
loadTimeSplit <- (proc.time() - sT)[3]
ObjSizeSplit <- format(object.size(y), units="Mb")

save(loadTimeSplit, ObjSizeSplit, file = file.path(SwishOutDir, "SplitLoadTimeAndMemory.RData"))
#Make sure the object is sorted by CellBarcodeName to get the Group label right
ySorted <- y[,SortedCellNames$CellBarcode]
ySorted$condition <- Groups2$Group

###Code below taken from the "Compress.R" file from Mike
# filter at the top
y2 <- labelKeep(ySorted, minCount=3, minN=5)
table(mcols(y2)$keep)
y3 <- y2[mcols(y2)$keep,]

y3 <- as(y3, "SingleCellExperiment")
y3 <- computeSumFactors(y3)

colData(y3)$sizeFactor <- sizeFactors(y3)

#assays(y3) <- lapply(assays(y3), as.matrix) # make dense matrices
#y4 <- makeInfReps(y3, 100)
#y4 <- scaleInfReps(y4, lengthCorrect=FALSE, sfFun=sfFun)

#Save the current SummarizedExperiment object to be able to load it back into R to summarize the results
save(y3, file = file.path(SwishOutDir, "SwishObjectSplit.RData"))

#Run this command once to get the structure of the Snakefile but you will need to upload the modified version to each folder with swish results
  #This file can be found in the splitSwishSnakeFile folder in the SimulationCode folder
#splitSwish(y3, 8, prefix=file.path(SwishOutDir, "swish"), snakefile=file.path(SwishOutDir, "Snakefile"), overwrite = FALSE)
splitSwish(y3, 8, prefix=file.path(SwishOutDir, "swish"), snakefile=NULL, overwrite = TRUE)


# Now, run snakemake on command line (edit as needed)
