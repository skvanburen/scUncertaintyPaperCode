#Save PBMC 4K Data Quantification Results

if(version$nickname=="Planting of a Tree"){
  .libPaths("/nas/longleaf/home/skvanbur/bin/R3.6.0")
}

library(tximport)
library(fishpond)
source("~/SingleCellProject/SingleCellProjectFunctions.R")

alevin_fil <- "/pine/scr/s/k/skvanbur/SingleCellProject/PBMC4KAlevinOutput/alevin/quants_mat.gz"
save_dir <- "/pine/scr/s/k/skvanbur/SingleCellProject/PBMC4KAlevinRes/"
if(!dir.exists(save_dir)){dir.create(save_dir, recursive = T)}
SaveAlevinDataToRData(alevin_fil = alevin_fil, save_dir = save_dir)

#PBMC4K_alevin_res <- tximport(alevin_fil, type = "alevin", dropInfReps = F)
#save(PBMC4K_alevin_res, file = paste0(save_dir, "PBMC4KAlevinData.RData"))
gc()
print(gc())
PBMC4K_alevin_res <- loadRData(paste0(save_dir, "AlevinData.RData"))

t1 <- PBMC4K_alevin_res$counts

t3 <- Matrix::rowMeans(t1)
t4 <- data.frame(t3)
colnames(t4) <- "MeanTGE"
t4$gene_id <- rownames(t4)

MeanExpressionByGenePBMC4K <- t4

save(MeanExpressionByGenePBMC4K, file = paste0(save_dir, "MeanExpressionByGenePBMC4K.RData"))

gc()
print(gc())

MeanVarInfRVBootstrapSamps <- CalcMeanVarInfRVBootstrapSamples(save_dir = save_dir)
save(MeanVarInfRVBootstrapSamps, file = paste0(save_dir, "MeanVarInfRVBootstrapSamps.RData"))

gc()
print(gc())
