sim_number <- 11001

#Specify libPaths to ensure the correct one is used
.libPaths("/nas/longleaf/home/skvanbur/lib64/R/library")

startTime <- proc.time()
library(fishpond)


source("~/SingleCellProject/SingleCellProjectFunctions.R")
def_wd2 <- "/pine/scr/s/k/skvanbur/SingleCellProject/SimulationResults/"

save_dir <- paste0(def_wd2, "Sim", sim_number, "/")

SwishOutDir <- paste0(save_dir, "SwishFiles/")
if(!dir.exists(SwishOutDir)){dir.create(SwishOutDir)}

SwishDat <- loadRData(file.path(SwishOutDir, "SwishObjectFull.RData"))

set.seed(1)
SwishResults  <- swish(SwishDat, x="condition")
TotalTime <- proc.time() - startTime

RunTimeFullSwish <- TotalTime[3]

save(SwishResults, file = file.path(SwishOutDir, "SwishFullResults.RData"))

#Print maximum memory usage as reported by R
print(gc())
MemoryInfoFullSwish <- gc()

save(RunTimeFullSwish, MemoryInfoFullSwish, file = file.path(SwishOutDir, "SwishFullTimeAndMemoryInfo.RData"))