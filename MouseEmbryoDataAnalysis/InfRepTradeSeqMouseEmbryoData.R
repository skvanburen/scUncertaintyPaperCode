#tradeSeq InfRep results for the Mouse Embryo Data
#Add libpaths to work with RStudio Server (custom installation on cluster)
.libPaths("~/lib64/R/library")
Sys.sleep(sample(1:200,1))
array_val <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
#devtools::install_github("statOmics/tradeSeq", ref = "f18706c7a1f7d1c54bff8439d08b64f785dabb47",lib = "~/installedtradeSeq1.1.24/")
library(tradeSeq, lib.loc = "~/InstalledtradeSeq1.1.24/")
if(packageVersion("tradeSeq")!="1.1.24"){stop("Use tradeSeq version 1.1.24 to avoid potential errors that occurred with prior versions")}
library(Matrix)

current_InfRep_seed <- array_val
current_tradeSeq_seed <- array_val*2

source("~/SingleCellProject/SingleCellProjectFunctions.R")
base_dir <- "/pine/scr/s/k/skvanbur/SingleCellProject/MouseEmbryoData/"

curr_save_dir <- paste0(base_dir, "tradeSeqInfRepResults/")
if(!dir.exists(curr_save_dir)){ 
  dir.create(curr_save_dir, recursive = TRUE)
}


curr_save_fil_name <- paste0(curr_save_dir, "PInfRepNum", array_val, ".RData")
print(paste0("Current save file name is ", curr_save_fil_name))

#Fix the lineages to be the same as in the EM results to ensure results can be easily and accurately combined across replicates
load(file.path(base_dir, "EMCountsLineageObjects.RData"))

#Generate the pseudo-inferential replicates
GeneratePInfReps <- function(AlevinMean, AlevinVar, current_InfRep_seed){
  
  if(sum(colnames(AlevinMean)!=colnames(AlevinVar))!=0){
    stop("AlevinMean and AlevinVar objects must have the same column ordering")
  }
  
  if(sum(rownames(AlevinMean)!=rownames(AlevinVar))!=0){
    stop("AlevinMean and AlevinVar objects must have the same row ordering")
  }
  
  AlevinMean2 <- as.matrix(AlevinMean)
  AlevinVar2 <- as.matrix(AlevinVar)
  
  
  cell_names <- colnames(AlevinMean2)
  genenames <- rownames(AlevinMean2)
  ngenes <- nrow(AlevinMean2)
  numCells <- ncol(AlevinMean2)
  theta_vals <- matrix(NA, nrow = ngenes, ncol = numCells)
  for(z in 1:numCells){
    curr_means <- AlevinMean2[,z]
    curr_variances <- AlevinVar2[,z]
    theta_vals[,z] <- (curr_means^2)/(curr_variances - curr_means)
  }
  theta_vals[!is.finite(theta_vals)] <- 1000
  theta_vals[theta_vals > 1000] <- 1000
  theta_vals[theta_vals < 0] <- 1000
  rownames(theta_vals) <- genenames
  colnames(theta_vals) <- colnames(AlevinMean2)
  
  set.seed(current_InfRep_seed)
  curr_samp_vals <- matrix(NA, nrow = ngenes, ncol = numCells)
  for(j in 1:numCells){
    p1 <- AlevinMean2[,j]
    #curr_samp_vals[,j] <- rpois(ngenes, p1)
    curr_theta <- theta_vals[,j]
    curr_samp_vals[,j] <- MASS::rnegbin(ngenes, mu = p1, theta = curr_theta)
  }
  #Be careful here: These columns so far are not ordered, so needs to be cell_names, not ordered_cell_names
  colnames(curr_samp_vals) <- colnames(AlevinMean2)
  rownames(curr_samp_vals) <- rownames(AlevinMean2)
  
  return(curr_samp_vals)
  
}


dat1 <- loadRData(file.path(base_dir, "QuantObj8.0UsingEM.RData"))
dat2 <- loadRData(file.path(base_dir, "QuantObj8.25UsingEM.RData"))
dat3 <- loadRData(file.path(base_dir, "QuantObj8.5UsingEM.RData"))

Mean1T <- dat1$mean[,CellsUsedDat1]
Mean2T <- dat2$mean[,CellsUsedDat2]
Mean3T <- dat3$mean[,CellsUsedDat3]

Mean1 <- Mean1T[filtered_genes,]
Mean2 <- Mean2T[filtered_genes,]
Mean3 <- Mean3T[filtered_genes,]

Var1T <- dat1$variance[,CellsUsedDat1]
Var2T <- dat2$variance[,CellsUsedDat2]
Var3T <- dat3$variance[,CellsUsedDat3]

Var1 <- Var1T[filtered_genes,]
Var2 <- Var2T[filtered_genes,]
Var3 <- Var3T[filtered_genes,]

PInfReps1 <- GeneratePInfReps(Mean1, Var1, current_InfRep_seed)
PInfReps2 <- GeneratePInfReps(Mean2, Var2, current_InfRep_seed)
PInfReps3 <- GeneratePInfReps(Mean3, Var3, current_InfRep_seed)

colnames(PInfReps1) <- paste0(colnames(PInfReps1), "_1")
colnames(PInfReps2) <- paste0(colnames(PInfReps2), "_2")
colnames(PInfReps3) <- paste0(colnames(PInfReps3), "_3")

countsT <- cbind(PInfReps1, PInfReps2, PInfReps3)
counts <- countsT[,FinalcountsCellNames]

if(ncol(counts)!=ncol(countsT)){
  stop("The number of columns of countsT should be the same as counts but it is not")
}

set.seed(current_tradeSeq_seed)

nKnots <- 6

control <- mgcv::gam.control()
control$maxit <- 1000 #set maximum number of iterations to 1K

print(gc())
print(paste0("tradeSeq results will be run with ", nKnots, " knots"))


sce <- fitGAM(counts = counts, sds = crv, nknots = nKnots, sce = TRUE, control = control)

association_test_res <- associationTest(sce, global = TRUE, lineages = FALSE)
start_vs_end_test_res <- startVsEndTest(sce, global = TRUE, lineages = FALSE)



n_of_fit_lin <- length(crv@lineages)
if(n_of_fit_lin > 1){
  patternTestRes <- tryCatch(patternTest(sce, global = TRUE, pairwise = FALSE), error = function(x){print(x)})
  diffEndTestRes <- tryCatch(diffEndTest(sce, global = TRUE, pairwise = FALSE), error = function(x){print(x)})
  earlyDETestRes <- tryCatch(earlyDETest(sce, knots = c(1, 3), global = TRUE, pairwise = FALSE), error = function(x){print(x)})
}else{
  patternTestRes <- NULL
  diffEndTestRes <- NULL
  earlyDETestRes <- NULL
}


save(crv, sce, counts, n_of_fit_lin, association_test_res, 
     start_vs_end_test_res, diffEndTestRes, patternTestRes, earlyDETestRes, file = curr_save_fil_name)

print(gc())
