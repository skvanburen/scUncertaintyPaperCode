#Pipeline code
onCluster <- TRUE

#Specify libPaths to ensure the correct one is used
.libPaths("/nas/longleaf/home/skvanbur/lib64/R/library")

library(splatter)
library(readr)
library(tximport)
library(fishpond)
library(data.table)

source("~/SingleCellProject/SingleCellProjectFunctions.R")
  
  #The "name" of this sim is "11001" to distinguish it from others
  array_val <- 11001

  
  

  tx2gene <- read_tsv("~/SingleCellProject/gencode_v32_files/tx2gene.tsv", col_names = F)
  colnames(tx2gene) <- c("tx_id", "gene_id")
  
  #set.seed(current_seed)
  
  gencode_files_location <- "~/SingleCellProject/gencode_v32_files/"
  gencode_v32_index_location <- "/pine/scr/s/k/skvanbur/SingleCellProject/GENCODEv32SalmonIndex/"
  def_wd2 <- "/pine/scr/s/k/skvanbur/SingleCellProject/SimulationResults/"
  whitelistBarcodesFile <- "~/SingleCellProject/737K-august-2016.txt"
  
  sim_number <- array_val
  print(paste0("sim_number is ", sim_number))
  
  save_dir <- paste0(def_wd2, "Sim", sim_number, "/")
  
  if(!dir.exists(save_dir)){dir.create(save_dir, recursive = TRUE)}
  
  SplatterSaveDir <- paste0(save_dir, "SplatterFiles/")
  if(!dir.exists(SplatterSaveDir)){dir.create(SplatterSaveDir, recursive = TRUE)}
  
  MinnowSaveDir <- paste0(save_dir, "MinnowOutput/")
  if(!dir.exists(MinnowSaveDir)){dir.create(MinnowSaveDir, recursive = TRUE)}
  
  AlevinSaveDir <- paste0(save_dir, "AlevinOutput/")
  if(!dir.exists(AlevinSaveDir)){dir.create(AlevinSaveDir, recursive = TRUE)}
  
  #Number of genes is set to 60179 because the real data quantification of the PBMC 4K data (using GENCODE v32 reference for reference transcripts only)
  #using Alevin results in 60,179 genes being present in the final quantification output
  #and we will use these to be the gene names such that we need the same number of genes
  
  nGenes <- 60179
  
  #Simulate 250 cells and take 100 from each group, just as is done in the swish paper code, to ensure each group has exactly 100 cells
  numCells <- 250

  #Simulation code below is adapted from the swishPaper code
  params <- newSplatParams()
  # LFC centered at 2

  params <- setParam(params, "de.facLoc", log(2) * 3)
  params <- setParam(params, "de.facScale", log(2) * 1)
  params <- setParam(params, "de.prob", 0.1)
  
  SplatterSimObject <- splatSimulate(params,
                       group.prob=c(0.5,0.5),
                       method="groups",
                       nGenes=nGenes,
                       batchCells = numCells,
                       seed=1)
  
  SplatterSimObject <- SplatterSimObject[,order(SplatterSimObject$Group)]
  SplatterSimObject <- SplatterSimObject[,c(1:100,151:250)]
  mcols(SplatterSimObject)$de <- with(mcols(SplatterSimObject), DEFacGroup2/DEFacGroup1 != 1)
  print(table(mcols(SplatterSimObject)$de))
  print(table(SplatterSimObject$Group))
  
  #################################################################################################################################
  #Assign gene names such that the expression ordering from the splatter simulated genes are the same as from the PBMC 4K data
  #ie the most expressed gene from splatter simulation is assigned the name of the most expressed gene from the PBMC 4K data, etc
  load("/pine/scr/s/k/skvanbur/SingleCellProject/MeanExpressionByGenePBMC4K.RData")
  RowMeansSplatter <- data.frame(rowMeans(counts(SplatterSimObject)))
  colnames(RowMeansSplatter) <- "rowMean"
  RowMeansSplatter$Rank <- rank(RowMeansSplatter$rowMean)
  MeanExpressionByGenePBMC4K$Rank <- rank(MeanExpressionByGenePBMC4K$MeanTGE)
  
  RowMeansSplatterOrdered <- RowMeansSplatter[order(RowMeansSplatter$Rank),]
  MeanExpressionByGenePBMC4KOrdered <- MeanExpressionByGenePBMC4K[order(MeanExpressionByGenePBMC4K$Rank),]
  
  RowMeansSplatterOrdered$gene_id <- MeanExpressionByGenePBMC4KOrdered$gene_id
  
  RowMeansSplatterOrdered2 <- RowMeansSplatterOrdered[gtools::mixedorder(rownames(RowMeansSplatterOrdered)),]
  
  if(sum(rownames(RowMeansSplatterOrdered2)!=rownames(SplatterSimObject))!=0){
    stop("The rownames that will be assigned to the splatter simulation object are not in the correct order")
  }
  rownames(SplatterSimObject) <- RowMeansSplatterOrdered2$gene_id
  ###############################################################################################################################
  
  save(SplatterSimObject, file = paste0(SplatterSaveDir, "simObject.RData"))
  
  write.table(rownames(SplatterSimObject), file= file.path(SplatterSaveDir, "quants_mat_rows.txt"), quote=FALSE, col.names=FALSE, row.names=FALSE)
  write.table(colnames(SplatterSimObject), file= file.path(SplatterSaveDir, "quants_mat_cols.txt"), quote=FALSE, col.names=FALSE, row.names=FALSE)
  write.table(counts(SplatterSimObject), file= file.path(SplatterSaveDir, "quants_mat.csv"), quote=FALSE, col.names=FALSE, row.names=FALSE, sep=",") 
  
  
  print(paste0("Maximum library size is ", max(colSums(counts(SplatterSimObject)))))
  
  print(gc()) 
