#GenerateDynToyObjects

onCluster <- TRUE
#Specify libPaths to ensure the correct one is used
.libPaths("/nas/longleaf/home/skvanbur/lib64/R/library")
library(data.table)
library(tidyverse, quietly = TRUE)
library(TSCAN)
library(dyntoy)
library(dynwrap)
source("~/SingleCellProject/SingleCellProjectFunctions.R")

seeds <- c(98,333)

tx2gene <- read_tsv("~/SingleCellProject/gencode_v32_files/tx2gene.tsv", col_names = F)
colnames(tx2gene) <- c("tx_id", "gene_id")


gencode_files_location <- "~/SingleCellProject/gencode_v32_files/"
gencode_v32_index_location <- "/pine/scr/s/k/skvanbur/SingleCellProject/GENCODEv32SalmonIndex/"
def_wd2 <- "/pine/scr/s/k/skvanbur/SingleCellProject/SimulationResults/"
whitelistBarcodesFile <- "~/SingleCellProject/737K-august-2016.txt"

#Number of genes is set to 60179 because the real data quantification of the PBMC 4K data (using GENCODE v32 reference for reference transcripts only)
#using Alevin results in 60,179 genes being present in the final quantification output
#and we will use these to be the gene names such that we need the same number of genes
nGenes <- 60179

#Sim 117 corresponds to the trifurcating simulation with 100 cells and 118 gives the simulation with 250 cells
  #"117" and "118" were used as names for each to ensure distinct indexing from other simulations
  #In the final results, 10117 amd 10118 were used as the main results (that used 20 inferential replicates), and
  #117 and 118 used 100 inferential replicates, though 117 and 10117 used the same simulated data objects (as did 118 and 10118)
for(uu in c(117,118)){
  sim_number <- uu
  current_seed <- seeds[uu - 116]
  
  print(paste0("sim_number is ", sim_number))
  
  save_dir <- paste0(def_wd2, "Sim", sim_number, "/")
  
  if(!dir.exists(save_dir)){dir.create(save_dir, recursive = TRUE)}
  
  #Files are saved in the "SplatterSaveDir" even though these simulations are from dyntoy and not splatter, to match same directory structure as the splatter simulations
  SplatterSaveDir <- paste0(save_dir, "SplatterFiles/")
  if(!dir.exists(SplatterSaveDir)){dir.create(SplatterSaveDir, recursive = TRUE)}
  
  MinnowSaveDir <- paste0(save_dir, "MinnowOutput/")
  if(!dir.exists(MinnowSaveDir)){dir.create(MinnowSaveDir, recursive = TRUE)}
  
  AlevinSaveDir <- paste0(save_dir, "AlevinOutput/")
  if(!dir.exists(AlevinSaveDir)){dir.create(AlevinSaveDir, recursive = TRUE)}
  
  set.seed(current_seed)
  #Set the initial dropout_probability_factor value to be the current default, 100
  dropout_probability_factor <- 100
if(sim_number==117){
    numCells <- 100
    differentially_expressed_rate <- 0.20
    f1 <- function() runif(1, 100, 1000)
    scaleFinalCounts <- TRUE
    model <- model_multifurcating(num_branchpoints = 2, max_degree = 3)
  }else if(sim_number==118){
    numCells <- 250
    differentially_expressed_rate <- 0.20
    f1 <- function() runif(1, 100, 1000)
    scaleFinalCounts <- TRUE
    model <- model_multifurcating(num_branchpoints = 2, max_degree = 3)
  }

  
  dataset <- generate_dataset(
    model = model,
    num_cells = numCells,
    num_features = nGenes,
    differentially_expressed_rate = differentially_expressed_rate,
    sample_mean_count = f1,
    normalise = F,
    dropout_probability_factor = dropout_probability_factor
  )
  
  
  DyntoyStartMilestones <- dataset$prior_information$start_milestones
  DyntoyEndMilestones <- dataset$prior_information$end_milestones
  
  print(paste0("Start Milestones are ", DyntoyStartMilestones))
  print(paste0("End Milestones are ", DyntoyEndMilestones))
  
  print(colSums(table(dataset$prior_information$groups_id)))
  countsT <- t(as.matrix(dataset$counts))
  
  if(scaleFinalCounts==TRUE){
    counts <- apply(countsT, 2, function(x){x * (1000000/sum(x))})
  }else{
    counts <- as.matrix(countsT)
  }
  
  load("/pine/scr/s/k/skvanbur/SingleCellProject/MeanExpressionByGenePBMC4K.RData")
  RowMeansDyntoy <- data.frame(rowMeans(counts))
  colnames(RowMeansDyntoy) <- "rowMean"
  RowMeansDyntoy$Rank <- rank(RowMeansDyntoy$rowMean)
  MeanExpressionByGenePBMC4K$Rank <- rank(MeanExpressionByGenePBMC4K$MeanTGE)
  
  RowMeansDyntoyOrdered <- RowMeansDyntoy[order(RowMeansDyntoy$Rank),]
  MeanExpressionByGenePBMC4KOrdered <- MeanExpressionByGenePBMC4K[order(MeanExpressionByGenePBMC4K$Rank),]
  
  RowMeansDyntoyOrdered$gene_id <- MeanExpressionByGenePBMC4KOrdered$gene_id
  
  RowMeansDyntoyOrdered2 <- RowMeansDyntoyOrdered[gtools::mixedorder(rownames(RowMeansDyntoyOrdered)),]
  
  if(sum(rownames(RowMeansDyntoyOrdered2)!=rownames(counts))!=0){
    stop("The rownames that will be assigned to the simulation object are not in the correct order")
  }
  
  rownames(counts) <- RowMeansDyntoyOrdered2$gene_id
  colnames(counts) <- paste0("Cell", 1:numCells)
  
  save(dataset, file = paste0(SplatterSaveDir, "simObject.RData"))
  
  write.table(rownames(counts), file= file.path(SplatterSaveDir, "quants_mat_rows.txt"), quote=FALSE, col.names=FALSE, row.names=FALSE)
  write.table(colnames(counts), file= file.path(SplatterSaveDir, "quants_mat_cols.txt"), quote=FALSE, col.names=FALSE, row.names=FALSE)
  write.table(counts, file= file.path(SplatterSaveDir, "quants_mat.csv"), quote=FALSE, col.names=FALSE, row.names=FALSE, sep=",") 
  
}


