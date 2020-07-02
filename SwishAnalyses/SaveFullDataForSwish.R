#Specify libPaths to ensure the correct one is used
.libPaths("/nas/longleaf/home/skvanbur/lib64/R/library")

#Save Full Data for swish
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

#Create a fake temporary name for the file to avoid a potential error in tximeta
coldata$names <- "1"

#Set location for the tximeta files to avoid needing to rerun this everytime
  #Also, make sure the tximeta files exist here before running the results for the computation time so they are not 
  #included in the results
setTximetaBFC(def_wd2)
sT <- proc.time()
y <- tximeta(coldata, type = "alevin", dropInfReps = FALSE, alevinArgs = list(filterBarcodes = FALSE, tierImport = TRUE))
loadTimeFull <- (proc.time() - sT)[3]

#Make sure the object is sorted by CellBarcodeName to get the Group label right
ySorted <- y[,SortedCellNames$CellBarcode]
ySorted$condition <- Groups2$Group


ObjSizeFull <- format(object.size(ySorted), units="Mb")
save(loadTimeFull, ObjSizeFull, file = file.path(SwishOutDir, "FullLoadTimeAndMemory.RData"))
y2 <- labelKeep(ySorted, minCount=3, minN=5)
table(mcols(y2)$keep)
y2 <- y2[mcols(y2)$keep,]

y2 <- as(y2, "SingleCellExperiment")
y2 <- computeSumFactors(y2)
assays(y2) <- lapply(assays(y2), as.matrix)

y2 <- scaleInfReps(y2, lengthCorrect=FALSE, sfFun=sizeFactors(y2))

save(y2, file = file.path(SwishOutDir, "SwishObjectFull.RData"))