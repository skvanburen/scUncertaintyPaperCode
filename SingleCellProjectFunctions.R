#Functions for to reproduce analysis from Compression of quantification uncertainty for scRNA-seq counts

#Neat little function that will load an .RData object and save its contents to whatever as the results of the output
#For example can do d <- loadRData("file.RData") and whatever is loaded is saved as d
#Modified from code from stackexchange user ricardo at this link https://stackoverflow.com/questions/5577221/how-can-i-load-an-object-into-a-variable-name-that-i-specify-from-an-r-data-file
loadRData <- function(fileName, objNameToGet = NULL){
  #loads an RData file, and returns it
  load(fileName)
  #print(ls()[ls() != "fileName"])
  if(is.null(objNameToGet)){
    rm(objNameToGet)
    #print(ls()[ls() != "fileName"])
    return(get(ls()[ls() != "fileName"]))
  }else{
    return(get(objNameToGet))
  }
  
}


SaveAlevinDataToRData <- function(alevin_fil, save_dir, ReSaveAllAlevinFiles = FALSE){
  library(tximport)
  library(fishpond)
  
  direc <- paste0(save_dir, "alevinInfRepData/")
  if(!dir.exists(direc)){dir.create(direc, recursive = TRUE)}
  
  InfRepLevelFilesSaveDir <- paste0(direc, "InfRepLevelFiles/")
  if(!dir.exists(InfRepLevelFilesSaveDir)){dir.create(InfRepLevelFilesSaveDir, recursive = TRUE)}
  
  CellLevelFilesSaveDir <- paste0(direc, "CellLevelFiles/")
  if(!dir.exists(CellLevelFilesSaveDir)){dir.create(CellLevelFilesSaveDir, recursive = TRUE)}
  
  #First, use tximport to read in variance matrices of the inf reps and save that
    #tximport does not load in the mean values for some reason so do those yourself later
  
  
  if(!file.exists(paste0(save_dir, "AlevinData.RData")) | ReSaveAllAlevinFiles==TRUE){
    FullAlevinT1 <- proc.time()
    alevin_res <- tximport(alevin_fil, type = "alevin", dropInfReps = FALSE)
    FullAlevinT2 <- proc.time() - FullAlevinT1
    
    FullAlevinResReadInTime <- FullAlevinT2[3]
    FullAlevinResMemorySizeBytes <- as.numeric(object.size(alevin_res))
    
    TimeToSaveAlevinRes1 <- proc.time()
    save(alevin_res, file = paste0(save_dir, "AlevinData.RData"))
    TimeToSaveAlevinRes2 <- proc.time() - TimeToSaveAlevinRes1
    
    TimeToSaveAlevinResToDisk <- TimeToSaveAlevinRes2[3]
    FullAlevinResSizeOnDiskBytes <- file.size(paste0(save_dir, "AlevinData.RData"))
    
    infR <- alevin_res$infReps
    infRepsListMemorySizeBytes <- as.numeric(object.size(infR))
    rm(infR)
    gc()
    
    save(FullAlevinResMemorySizeBytes, FullAlevinResReadInTime, FullAlevinResSizeOnDiskBytes, TimeToSaveAlevinResToDisk, infRepsListMemorySizeBytes,
         file = paste0(save_dir, "FullAlevinDataSizeAndComputationTimeRes.RData"))
    print("AlevinData.RData file with full infReps saved successfully")
  }else{
    alevin_res <- loadRData(paste0(save_dir, "AlevinData.RData"))
  }
  
  
  
  #Now, save datasets by cell that have gene on the rows and infReps on the columns- done for convenience, should help get results later
  infR <- alevin_res$infReps
  
  rm(alevin_res)
  gc()
  

  

  if(is.null(infR)){
    stop("No inferential replicates found so cannot save those cell by cell or separately for each replicate")
  }else{
    ninfreps <- length(infR)
    
    #First, save variances of bootstrap samples
    if(!file.exists(paste0(save_dir, "BootSampsVariances.RData")) | ReSaveAllAlevinFiles==TRUE){
      
      alevin_res1 <- tximport(alevin_fil, type = "alevin", dropInfReps = TRUE)
      
      #This estimate comes straight from alevin, which uses the biased variance estimator (ie with n in the denom instead of n-1)
      #To match with the var() calls in the function that calculates coverage, I will multiply the variances from alevin by 
      #nboot/(nboot - 1)
      boot_samps_variancesT <- as.matrix(alevin_res1$variance)
      boot_samps_variances <- boot_samps_variancesT * (ninfreps/(ninfreps-1))
      save(boot_samps_variances, file = paste0(save_dir, "BootSampsVariances.RData"))
      print("File with variances of the bootstrap samples for each cell saved successfully")
      
      rm(alevin_res1)
      rm(boot_samps_variances)
      gc()
    }
    
    #Now, save separate files for each inferential replicate to reduce future memory requirements
    AllReplicatesSaved <- TRUE
    TimeToSaveInfRepsDatasets1 <- proc.time()
    for(i in 1:ninfreps){
      if(file.exists(paste0(InfRepLevelFilesSaveDir, "infRepDatRep", i, ".RData")) & ReSaveAllAlevinFiles==FALSE){
        AllReplicatesSaved <- FALSE
        next
      }
      st_rep <- proc.time()
      if(i%%25==0){
        print(paste0("Currently Saving InfReps for Replicate Number ", i, " out of ", ninfreps))
        
      }
      assign(paste0("infRepDatRep", i), as.matrix(infR[[i]]))
      save(list = paste0("infRepDatRep", i), file = paste0(InfRepLevelFilesSaveDir, "infRepDatRep", i, ".RData"))
      rm(list = paste0("infRepDatRep", i))
      gc()
      tt_rep <- proc.time()
      print(paste0("Time to save data for Rep ", i, " is ", (tt_rep[3] - st_rep[3])))
      
    }
    TimeToSaveInfRepsDatasets2 <- proc.time() - TimeToSaveInfRepsDatasets1
    TimeToSaveInfRepLevelDatasets <- TimeToSaveInfRepsDatasets2[3]
    gc()
    rm(i)
    
    if(AllReplicatesSaved==TRUE){
      save(TimeToSaveInfRepLevelDatasets, file = paste0(save_dir, "TimeToSaveInfRepLevelDatasets.RData"))
    }
    #Now, save the inferential replicates separately for each cell to make calculating coverage results easier
    genes <- rownames(as.matrix(infR[[1]]))

    cell_names <- colnames(as.matrix(infR[[1]]))
    ncells <- length(cell_names)
    
    #Additionally save the means of the bootstrap samples across each cell/gene to use later
    boot_samps_means <- matrix(NA, nrow = length(genes), ncol = ncells)
    TimeToSaveCellData <- numeric(length = ncells)
    
    AllCellsSaved <- TRUE
    for(c in 1:ncells){

      rearranged_dat <- matrix(NA, nrow = length(genes), ncol = ninfreps)
      rownames(rearranged_dat) <- genes
      curr_cell <- cell_names[c]
      
      st1 <- proc.time()
      print(paste0("Current Cell Being saved is ", c, " out of ", ncells))
      #rearranged_dat <- lapply(infR, function(x){return(x[,curr_cell])})
      for(i in 1:ninfreps){
        rearranged_dat[,i] <- infR[[i]][,curr_cell]
        #curr_dat <- loadRData(paste0(InfRepLevelFilesSaveDir, "infRepDatRep", i, ".RData"))
        #rearranged_dat[,i] <- curr_dat[,curr_cell]
        #curr_infR <- infR[[i]][,curr_cell]
        # for(g in 1:length(genes)){
        #   rearranged_dat[g,i] <- curr_infR[g]
        # }

        #rearranged_dat[,i] <- curr_infR
      }
      
      boot_samps_means[,c] <- rowMeans(rearranged_dat)
      
      if((file.exists(paste0(CellLevelFilesSaveDir, "infRepDat", curr_cell, ".RData"))  & ReSaveAllAlevinFiles==FALSE)){
        AllCellsSaved <- FALSE
        next
      }
      
      rearranged_dat2 <- data.frame(rearranged_dat)
      colnames(rearranged_dat2) <- as.character(1:ninfreps)
      rownames(rearranged_dat2) <- genes
      assign(paste0("infRepDat", curr_cell), rearranged_dat2)

      save(list = paste0("infRepDat", curr_cell), file = paste0(CellLevelFilesSaveDir, "infRepDat", curr_cell, ".RData"))
      
      rm(list = paste0("infRepDat", curr_cell))
      rm(rearranged_dat)
      rm(rearranged_dat2)
      TimeToSaveCellData[c] <- (proc.time() - st1)[3]
      print(paste0("Time to save data for cell ", c, " is ", TimeToSaveCellData[c]))
      if(c%%1==0){
        print(gc())
      }
      
    }
    rownames(boot_samps_means) <- genes
    colnames(boot_samps_means) <- cell_names
    
    save(boot_samps_means, file = paste0(save_dir, "BootSampsMeans.RData"))
    print("File with means of the bootstrap samples for each cell saved successfully")
    
    if(AllCellsSaved==TRUE){
      MeanTimeToSaveCellLevelFiles <- mean(TimeToSaveCellData)
      CellLevelFilesApproxSizeOnDisk <- file.size(paste0(CellLevelFilesSaveDir, "infRepDat", cell_names[1], ".RData"))
      save(MeanTimeToSaveCellLevelFiles, CellLevelFilesApproxSizeOnDisk, file = paste0(save_dir, "MeanTimeToSaveAndDiskSizeCellLevelFiles.RData"))
    }
  }
}


SummarizeBootstrapSamples <- function(SplatterSaveDir, MinnowSaveDir, uniqueness_fil_dir, save_dir, coverageRes = TRUE, calcMeanSDInfRVBootSamps = TRUE,
                                      useSplatterCountsScaledToAlevinLibrarySize = TRUE, sim_number, DropMinnowSkippedGenesFromCoverageCalc,
                                      CalculateCoveragesWithoutZeros = FALSE){
  alevin_res <- loadRData(paste0(save_dir, "AlevinDataNoInfRepsSim", sim_number, ".RData"))
  
  direc <- paste0(save_dir, "alevinInfRepData/")
  InfRepLevelFilesSaveDir <- paste0(direc, "InfRepLevelFiles/")
  CellLevelFilesSaveDir <- paste0(direc, "CellLevelFiles/")
  
  
  #Load in information from Avi about gene uniqueness to be able to stratify genes by gene uniqueness
  uniqueness_infoT <- read_tsv(paste0(uniqueness_fil_dir, "gene_uniqueness_GENCODE_v32.txt"))
  uniqueness_infoT2 <- subset(uniqueness_infoT, uniqueness_infoT$total!=0)
  
  uniqueness_infoT2$uniqueness_ratio <- uniqueness_infoT2$unique / uniqueness_infoT2$total
  uniqueness_info <- uniqueness_infoT2[,c("gene", "uniqueness_ratio")]
  colnames(uniqueness_info) <- c("gene_id", "uniqueness_ratio")
  
  cell_names <- colnames(alevin_res$counts)
  
  #Read in the "True" Counts (ie the ones splatter generated and minnow is generating reads based on)
  if(useSplatterCountsScaledToAlevinLibrarySize==TRUE){
    print("Results will be calculated using splatter counts scaled to the alevin library size")
    print(paste0("File to be loaded is ", paste0(save_dir, "SplatterCountsScaledToAlevinLibSize.RData")))
    Splatter_counts_scaled <- loadRData(paste0(save_dir, "SplatterCountsScaledToAlevinLibSize.RData"))
    SplatterCountsToUseT <- Splatter_counts_scaled
  }else{
    Splatter_countsT <- read_csv(paste0(SplatterSaveDir, "quants_mat.csv"), col_names = FALSE)
    Splatter_counts <- data.frame(Splatter_countsT)
    tt1 <- read_table(paste0(SplatterSaveDir, "quants_mat_rows.txt"), col_names = FALSE)
    rownames(Splatter_counts) <- as.character(tt1$X1)
    
    #Want the column names to correspond to the cell barcodes, not just Cell1, Cell2, etc so load these from the minnowOutput
    #Confusingly, both the col_names and row_names are called "quants_mats_rows" because Splatter and alevin have flipped what are
    #rows and what are columns
    col_names <- read.table(paste0(MinnowSaveDir, "alevin/quants_mat_rows.txt"))
    colnames(Splatter_counts) <- as.character(col_names$V1)
    
    SplatterCountsToUseT <- Splatter_counts
  }

  

  
  #Minnow will ignore genes that it can't/won't use for its simulation (for example if certain genes are not in the BFG file, etc)
    #such that the gene names cannot be extracted from the minnow output
    #Pull these instead from files saved by splatter directly
  if(DropMinnowSkippedGenesFromCoverageCalc==TRUE){
    MinnowSkippedGenes <- data.frame(read_tsv(paste0(uniqueness_fil_dir, "MinnowSkippedGenes.txt"), col_names = FALSE))
    SplatterCountsToUse <- SplatterCountsToUseT[!(rownames(SplatterCountsToUseT)%in%MinnowSkippedGenes[,1]),]
    #SplatterCountsToUse <- subset(SplatterCountsToUseT, !(rownames(SplatterCountsToUseT)%in%as.character(MinnowSkippedGenes[,1])))
  }else{
    SplatterCountsToUse <- SplatterCountsToUseT
  }
  
  print(paste0("Number of genes used to calculate coverage is ", nrow(SplatterCountsToUse)))
  
  
  ncells <- length(colnames(SplatterCountsToUse))
  cell_level_coverage_res <- vector(mode = "list", length = ncells)
  gene_level_coverage_res <- vector(mode = "list", length = ncells)
  
  for(k in 1:ncells){
    stttt <- proc.time()
    print(paste0("Calculating Coverage for cell ", k, " out of ", ncells))
    curr_datT <- loadRData(paste0(CellLevelFilesSaveDir, "infRepDat", cell_names[k], ".RData"))
    curr_SplatterCountsToUse <- SplatterCountsToUse[,cell_names[k], drop = FALSE]
    
    if(ncol(curr_SplatterCountsToUse)!=1){
      stop("The number of columns of curr_SplatterCountsToUse needs to be 1")
    }
    #This is the cell-level file that has all bootstrap replicates for that cell
    curr_datT2 <- curr_datT[rownames(curr_datT) %in% rownames(curr_SplatterCountsToUse),]
    
    
    #Ensure rows of the current data file and the expected counts are in the same order
    if(sum(rownames(curr_datT2)!=rownames(curr_SplatterCountsToUse))!=0){
      curr_dat <- curr_datT2[rownames(curr_SplatterCountsToUse), , drop = FALSE]
    }else{
      curr_dat <- curr_datT2
    }
    
    #Change counts for this particular cell that has a simulated expression of 0 to NA if this is desired
        #This will result in them not being included in the final coverage calculations
    if(CalculateCoveragesWithoutZeros==TRUE){
      curr_SplatterCountsToUse[curr_SplatterCountsToUse==0] <- NA
    }
    
    
    if(k==1){
      print("Dimensions of the data for the current cell that is being used to calculate coverage is given below")
      print(dim(curr_dat))
      
      print(paste0("The number of elements that are zero are ", sum(curr_SplatterCountsToUse[,1]==0, na.rm=T)))
      
      print(paste0("The number of missing elements are ", sum(is.na(curr_SplatterCountsToUse))))
    }

    
    nboot <- ncol(curr_dat)
    gene_names <- rownames(curr_dat)
    num_genes <- length(gene_names)
    
    #If the simulated count (val) is NA (meaning it was 0 and was changed to NA), return NA for the coverage
    iscontainedin <- function(val, min, max){
      if(is.na(val)){
        return(NA)
      }else{
        covera <- ifelse(min<=val & val<=max & (!is.na(min) & !is.na(max)), 1, 0)
        return(covera)
      }

    }
    
    t_val <- qt(0.975, df = nboot - 1)
    norm_val <- qnorm(.975)
    
    st2 <- proc.time()
    quants <-  t(apply(curr_dat, 1, function(x){quantile(x, probs = c(0.025, 0.05, 0.10, 0.25, 0.40, 0.50, 0.60, 0.75, 0.90, 0.95, 0.975), na.rm = T)}))
    tim2 <- proc.time() - st2
    
    
    
    min_curr <- apply(curr_dat, 1, function(x){min(x, na.rm = TRUE)})
    max_curr <- apply(curr_dat, 1, function(x){max(x, na.rm = TRUE)})
    mean_curr <- rowMeans(curr_dat, na.rm = TRUE)
    median_curr <- apply(curr_dat, 1, function(x){median(x, na.rm = TRUE)})
    sd_curr <- apply(curr_dat, 1, function(x){sd(x, na.rm = TRUE)})
    var_curr <- apply(curr_dat, 1, function(x){var(x, na.rm = TRUE)})
    
    
    
    
    genenames <- rownames(mean_curr)
    num_genes <- length(mean_curr)
    
    #Construct size values for the negative binomial interval
    #Theta values here are for var = mu + (mu/theta) so set a maximum value to ensure sufficient extra poisson variation is properly
    #present 
    #Mike had used phi = 0.001 in DESeq2 for phi=1/theta, so to correspond to that here set theta = 1/.001 = 1000 such that we specify a maximum
    #This is to ensure variance is always (at least slightly) greater than the mean
    #Theta may also be called phi, and is also equal to the "size" value from the base R functions
    #So, set theta to 1000 if theta is negative or greater than 1000
    theta_vals <- (mean_curr^2)/(var_curr - mean_curr)
    theta_vals[!is.finite(theta_vals)] <- 1000
    theta_vals[theta_vals > 1000] <- 1000
    theta_vals[theta_vals < 0] <- 1000
    
    #Code to check the code is working as intended
    #check_data <- rnbinom(10000, size = theta_vals[6], mu = mean_curr[6])
    #check_data <- rnbinom(10000, size = theta_vals, mu = mean_curr)
    #head(mean(check_data))
    #head(var(check_data)) #Var should be mu + (1/theta)*mu^2
    
    QLowerNBinom <- qnbinom(0.025, size = theta_vals, mu = mean_curr)
    QUpperNBinom <- qnbinom(0.975, size = theta_vals, mu = mean_curr)
    
    if(k==1){
      print("The first 5 rows of each objects of interest are given below")
      print("quantile values")
      print(head(quants, n=15))
      
      print("min values")
      print(head(min_curr, n=15))
      
      print("max values")
      print(head(max_curr, n=15))
      
      print("mean values")
      print(head(mean_curr, n=15))
      
      print("median values")
      print(head(median_curr, n=15))
      
      print("sd values")
      print(head(sd_curr, n=15))
      
      print("var values")
      print(head(var_curr, n=15))
      
      print("theta values")
      print(head(theta_vals, n=15))
    }
      
      numerator <- var_curr - mean_curr
      numerator[numerator < 0] <- 0
      denom <- mean_curr + 5
      
      
      InfRV <- (numerator/ denom) + 0.01
      LogInfRV <- log(InfRV)
      
    Norm95ConfInvLower <- mean_curr - (norm_val*sd_curr)
    Norm95ConfInvUpper <- mean_curr + (norm_val*sd_curr)
    
    TDist95ConfInvLower <- mean_curr - (t_val*sd_curr)
    TDist95ConfInvUpper <- mean_curr + (t_val*sd_curr)
    
    
    MinMaxDat <- cbind(curr_SplatterCountsToUse, min_curr, max_curr)
    TwoPointFiveNinetySevenPointFiveDat <- cbind(curr_SplatterCountsToUse, quants[,"2.5%"], quants[,"97.5%"])
    FiveNinetyFiveDat <- cbind(curr_SplatterCountsToUse, quants[,"5%"], quants[,"95%"])
    TenNinetyDat <- cbind(curr_SplatterCountsToUse, quants[,"10%"], quants[,"90%"])
    TwentyFiveSeventyFiveDat <- cbind(curr_SplatterCountsToUse, quants[,"25%"], quants[,"75%"])
    FourtySixtyDat <- cbind(curr_SplatterCountsToUse, quants[,"40%"], quants[,"60%"])
    Norm95ConfInvDat <- cbind(curr_SplatterCountsToUse, Norm95ConfInvLower, Norm95ConfInvUpper)
    TDist95ConfInvDat <- cbind(curr_SplatterCountsToUse, TDist95ConfInvLower, TDist95ConfInvUpper)
    NB95InvDat <- cbind(curr_SplatterCountsToUse, QLowerNBinom, QUpperNBinom)
    
    MinMax <- apply(MinMaxDat, 1, function(x){iscontainedin(x[1], x[2], x[3])})
    TwoPointFiveNinetySevenPointFive <- apply(TwoPointFiveNinetySevenPointFiveDat, 1, function(x){iscontainedin(x[1], x[2], x[3])})
    FiveNinetyFive <- apply(FiveNinetyFiveDat, 1, function(x){iscontainedin(x[1], x[2], x[3])})
    TenNinety <- apply(TenNinetyDat, 1, function(x){iscontainedin(x[1], x[2], x[3])})
    TwentyFiveSeventyFive <- apply(TwentyFiveSeventyFiveDat, 1, function(x){iscontainedin(x[1], x[2], x[3])})
    FourtySixty <- apply(FourtySixtyDat, 1, function(x){iscontainedin(x[1], x[2], x[3])})
    Norm95ConfInv <- apply(Norm95ConfInvDat, 1, function(x){iscontainedin(x[1], x[2], x[3])})
    TDist95ConfInv <- apply(TDist95ConfInvDat, 1, function(x){iscontainedin(x[1], x[2], x[3])})
    NB95Inv <- apply(NB95InvDat, 1, function(x){iscontainedin(x[1], x[2], x[3])})
    
    if(k==1){
      print("MinMaxDat")
      print(head(MinMaxDat, n=15))
      
      print("TwoPointFiveNinetySevenPointFiveDat")
      print(head(TwoPointFiveNinetySevenPointFiveDat, n=15))
      
      print("Norm95ConfInvDat")
      print(head(Norm95ConfInvDat, n=15))
      
      print("NB95InvDat")
      print(head(NB95InvDat, n=15))
    }
    
    
    if(k==1){
      print("MinMax")
      print(head(MinMax, n=15))
      
      print("TwoPointFiveNinetySevenPointFive")
      print(head(TwoPointFiveNinetySevenPointFive, n=15))
      
      print("Norm95ConfInv")
      print(head(Norm95ConfInv, n=15))
      
      print("NB95Inv")
      print(head(NB95Inv, n=15))
    }
    
    


    MinMaxWidth <- MinMaxDat[,3] - MinMaxDat[,2]
    TwoPointFiveNinetySevenPointFiveWidth <- TwoPointFiveNinetySevenPointFiveDat[,3] - TwoPointFiveNinetySevenPointFiveDat[,2]
    FiveNinetyFiveWidth <- FiveNinetyFiveDat[,3] - FiveNinetyFiveDat[,2]
    TenNinetyWidth <- TenNinetyDat[,3] - TenNinetyDat[,2]
    TwentyFiveSeventyFiveWidth <- TwentyFiveSeventyFiveDat[,3] - TwentyFiveSeventyFiveDat[,2]
    FourtySixtyWidth <- FourtySixtyDat[,3] - FourtySixtyDat[,2]
    Norm95ConfInvWidth <- Norm95ConfInvDat[,3] - Norm95ConfInvDat[,2]
    TDist95ConfInvWidth <- TDist95ConfInvDat[,3] - TDist95ConfInvDat[,2]
    NB95InvWidth <- NB95InvDat[,3] - NB95InvDat[,2]
    
    #Need to change width to NA for genes with Zero Counts if CalculateCoveragesWithoutZeros is TRUE
      #Otherwise the width is included in the final width even though the count is treated as NA, 
      #making that count not of interest
    if(CalculateCoveragesWithoutZeros==TRUE){
      MinMaxWidth[is.na(curr_SplatterCountsToUse)] <- NA
      TwoPointFiveNinetySevenPointFiveWidth[is.na(curr_SplatterCountsToUse)] <- NA
      FiveNinetyFiveWidth[is.na(curr_SplatterCountsToUse)] <- NA
      TenNinetyWidth[is.na(curr_SplatterCountsToUse)] <- NA
      TwentyFiveSeventyFiveWidth[is.na(curr_SplatterCountsToUse)] <- NA
      FourtySixtyWidth[is.na(curr_SplatterCountsToUse)] <- NA
      Norm95ConfInvWidth[is.na(curr_SplatterCountsToUse)] <- NA
      TDist95ConfInvWidth[is.na(curr_SplatterCountsToUse)] <- NA
      NB95InvWidth[is.na(curr_SplatterCountsToUse)] <- NA
    }  
    
    
    if(k==1){
      print("MinMaxWidth")
      print(head(MinMaxWidth))
      print("TwoPointFiveNinetySevenPointFiveWidth")
      print(head(TwoPointFiveNinetySevenPointFiveWidth))
      print("Norm95ConfInvWidth")
      print(head(Norm95ConfInvWidth))
      print("NB95InvWidth")
      print(head(NB95InvWidth))
      
      print(paste0("The number of TwoPointFiveNinetySevenPointFiveWidth values that are NA are ", sum(is.na(TwoPointFiveNinetySevenPointFiveWidth))))
      print(paste0("The number of NB95InvWidth values that are NA are ", sum(is.na(NB95InvWidth))))
    }
    
    
    res2 <- data.frame(MinMax = MinMax, TwoPointFiveNinetySevenPointFive = TwoPointFiveNinetySevenPointFive,
                       FiveNinetyFive = FiveNinetyFive, TenNinety = TenNinety, 
                       TwentyFiveSeventyFive = TwentyFiveSeventyFive, FourtySixty = FourtySixty,
                       Norm95ConfInv = Norm95ConfInv, TDist95ConfInv = TDist95ConfInv, NB95Inv = NB95Inv,
                       MinMaxWidth = MinMaxWidth, TwoPointFiveNinetySevenPointFiveWidth = TwoPointFiveNinetySevenPointFiveWidth,
                       FiveNinetyFiveWidth = FiveNinetyFiveWidth, TenNinetyWidth = TenNinetyWidth,
                       TwentyFiveSeventyFiveWidth = TwentyFiveSeventyFiveWidth, FourtySixtyWidth = FourtySixtyWidth,
                       Norm95ConfInvWidth = Norm95ConfInvWidth, TDist95ConfInvWidth = TDist95ConfInvWidth, NB95InvWidth = NB95InvWidth,
                       mean_curr = mean_curr, median_curr = median_curr, sd_curr = sd_curr, var_curr = var_curr,
                       InfRV = InfRV, LogInfRV = LogInfRV)
    
    rownames(res2) <- gene_names
    res2[,"cell_name"] <- cell_names[k]
    res2[,"gene_id"] <- gene_names
    
    gene_level_coverage_res[[k]] <- res2
    
    cell_level_coverage_res[[k]] <- list(MinMaxCoverage = mean(res2$MinMax, na.rm = TRUE), TwoPointFiveNinetySevenPointFiveCoverage = mean(res2$TwoPointFiveNinetySevenPointFive, na.rm = TRUE), FiveNinetyFiveCoverage = mean(res2$FiveNinetyFive, na.rm = TRUE), TenNinetyCoverage = mean(res2$TenNinety, na.rm = TRUE),
                                         TwentyFiveSeventyFiveCoverage = mean(res2$TwentyFiveSeventyFive, na.rm = TRUE), FourtySixtyCoverage = mean(res2$FourtySixty, na.rm = TRUE),
                                         Norm95ConfInvCoverage = mean(res2$Norm95ConfInv, na.rm = TRUE), TDist95ConfInvCoverage = mean(res2$TDist95ConfInv, na.rm = TRUE), NB95InvCoverage = mean(res2$NB95Inv, na.rm = TRUE),
                                         MeanMinMaxWidth = mean(res2$MinMaxWidth, na.rm = TRUE), MeanTwoPointFiveNinetySevenPointFiveWidth = mean(res2$TwoPointFiveNinetySevenPointFiveWidth, na.rm = TRUE), MeanFiveNinetyFiveWidth = mean(res2$FiveNinetyFiveWidth, na.rm = TRUE), MeanTenNinetyWidth = mean(res2$TenNinetyWidth, na.rm = TRUE),
                                         MeanTwentyFiveSeventyFiveWidth = mean(res2$TwentyFiveSeventyFiveWidth, na.rm = TRUE), MeanFourtySixtyWidth = mean(res2$FourtySixtyWidth, na.rm = TRUE),
                                         MeanNorm95ConfInvWidth = mean(res2$Norm95ConfInvWidth, na.rm = TRUE), MeanTDist95ConfInvWidth = mean(res2$TDist95ConfInvWidth, na.rm = TRUE), MeanNB95InvWidth = mean(res2$NB95InvWidth, na.rm = TRUE))
    
    print(paste0("Total time for cell ", k, " is ", (proc.time() - stttt)[3]))
  }
  
  cell_level_coverage_resFinal <- data.frame(data.table::rbindlist(cell_level_coverage_res))
  rownames(cell_level_coverage_resFinal) <- cell_names
  
  gene_level_coverage_res2 <- data.frame(data.table::rbindlist(gene_level_coverage_res))
  
  print(paste0("The dimensions of gene_level_coverage_res2 are given below"))
  print(dim(gene_level_coverage_res2))
  if(nrow(gene_level_coverage_res2)!=(ncells*num_genes)){
    stop("The dimensions of gene_level_coverage_res2 seem to be incorrect")
  }
  gene_names <- gene_level_coverage_res[[1]][,"gene_id"]
  num_genes <- length(gene_names)
  gene_level_coverage_res3 <- vector(mode = "list", length = num_genes)
  for(g in 1:num_genes){
    #1, 1001, 2001, 
    if(g %%10000==0){
      print(paste0("Current gene is ", g, " out of ", num_genes))
    }
    curr_vals <- seq(from = g, to = nrow(gene_level_coverage_res2), by = num_genes)
    dat_curr_gene <- gene_level_coverage_res2[curr_vals,]
    
    gene_level_coverage_res3[[g]] <- list(MinMaxCoverage = mean(dat_curr_gene$MinMax, na.rm = TRUE), TwoPointFiveNinetySevenPointFiveCoverage = mean(dat_curr_gene$TwoPointFiveNinetySevenPointFive, na.rm = TRUE), FiveNinetyFiveCoverage = mean(dat_curr_gene$FiveNinetyFive, na.rm = TRUE), TenNinetyCoverage = mean(dat_curr_gene$TenNinety, na.rm = TRUE),
                                          TwentyFiveSeventyFiveCoverage = mean(dat_curr_gene$TwentyFiveSeventyFive, na.rm = TRUE), FourtySixtyCoverage = mean(dat_curr_gene$FourtySixty, na.rm = TRUE),
                                          Norm95ConfInvCoverage = mean(dat_curr_gene$Norm95ConfInv, na.rm = TRUE), TDist95ConfInvCoverage = mean(dat_curr_gene$TDist95ConfInv, na.rm = TRUE), NB95InvCoverage = mean(dat_curr_gene$NB95Inv, na.rm = TRUE),
                                          MeanMinMaxWidth = mean(dat_curr_gene$MinMaxWidth, na.rm = TRUE), MeanTwoPointFiveNinetySevenPointFiveWidth = mean(dat_curr_gene$TwoPointFiveNinetySevenPointFiveWidth, na.rm = TRUE),
                                          MeanFiveNinetyFiveWidth = mean(dat_curr_gene$FiveNinetyFiveWidth, na.rm = TRUE), MeanTenNinetyWidth = mean(dat_curr_gene$TenNinetyWidth, na.rm = TRUE),
                                          MeanTwentyFiveSeventyFiveWidth = mean(dat_curr_gene$TwentyFiveSeventyFiveWidth, na.rm = TRUE), MeanFourtySixtyWidth = mean(dat_curr_gene$FourtySixtyWidth, na.rm = TRUE),
                                          MeanNorm95ConfInvWidth = mean(dat_curr_gene$Norm95ConfInvWidth, na.rm = TRUE), MeanTDist95ConfInvWidth = mean(dat_curr_gene$TDist95ConfInvWidth, na.rm = TRUE),
                                          MeanNB95InvWidth = mean(dat_curr_gene$NB95InvWidth, na.rm = TRUE),
                                          MeanofMeanBootSamps = mean(dat_curr_gene$mean_curr, na.rm = TRUE), MeanofMedianBootSamps = mean(dat_curr_gene$median_curr, na.rm = TRUE),
                                          MedianofMeanBootSamps = median(dat_curr_gene$mean_curr, na.rm = TRUE), MedianofMedianBootSamps = median(dat_curr_gene$median_curr, na.rm = TRUE),
                                          MeanInfRV = mean(dat_curr_gene$InfRV, na.rm = TRUE), MeanLogInvRV = mean(dat_curr_gene$LogInfRV, na.rm = TRUE),
                                          MedianInfRV = median(dat_curr_gene$InfRV, na.rm = TRUE), MedianLogInvRV = median(dat_curr_gene$LogInfRV, na.rm = TRUE),
                                          MeanSDBootSamps = mean(dat_curr_gene$sd_curr, na.rm = TRUE), MedianSDBootSamps = median(dat_curr_gene$sd_curr, na.rm = TRUE),
                                          MeanVarBootSamps = mean(dat_curr_gene$var_curr, na.rm = TRUE), MedianVarBootSamps = median(dat_curr_gene$var_curr, na.rm = TRUE))
  }
  
  gene_level_coverage_resFinalT <- data.frame(rbindlist(gene_level_coverage_res3))
  rownames(gene_level_coverage_resFinalT) <- gene_names
  gene_level_coverage_resFinalT$gene_id <- gene_names
  
  gene_level_coverage_resFinalT2 <- merge(gene_level_coverage_resFinalT, uniqueness_info, by = "gene_id", all = TRUE)
  rownames(gene_level_coverage_resFinalT2) <- gene_level_coverage_resFinalT2$gene_id
  
  #Only keep genes corresponding to the list of possibly simulated genes
    #ie, remove any extra genes in the uniqueness file
  gene_level_coverage_resFinal <- subset(gene_level_coverage_resFinalT2, gene_level_coverage_resFinalT2$gene_id %in% gene_names)
  
  #Overall Coverage Results will average each column over different cells
  cell_level_averaged_coverage <- data.frame(colMeans(cell_level_coverage_resFinal))
  colnames(cell_level_averaged_coverage) <- "MeanAcrossCells"
  
  #Now, aggregate results based on gene overlap bin
  gene_level_coverage_resFinal[,"UniqFactor"] <- cut(gene_level_coverage_resFinal$uniqueness_ratio,
                                          breaks=c(0,0.20,0.40,0.60,0.80,1),
                                          include.lowest=TRUE,
                                          labels=c("0to0.20", "0.20to0.40", "0.40to0.60", "0.60to0.80", "0.80to1"))
  gene_level_coverage_resFinal$UniqFactor[is.na(gene_level_coverage_resFinal$uniqueness_ratio)] <- NA
  GeneLevelCoverageN <- gene_level_coverage_resFinal[,unlist(lapply(gene_level_coverage_resFinal, is.numeric))]
  
  coverage_res_by_uniqueness_bin <- aggregate(GeneLevelCoverageN, by = list(UniqFactor = gene_level_coverage_resFinal$UniqFactor), FUN = mean)
  coverage_res_by_uniqueness_bin$uniqueness_ratio <- NULL
  
  
  return(list(gene_level_coverage = gene_level_coverage_resFinal, cell_level_coverage = cell_level_coverage_resFinal,
              cell_level_averaged_coverage = cell_level_averaged_coverage,
              coverage_res_by_uniqueness_bin = coverage_res_by_uniqueness_bin))
  
}


CalcMeanVarInfRVBootstrapSamples <- function(save_dir){
  
  load(paste0(save_dir, "AlevinData.RData"))
  regular_counts <- alevin_res$counts
  cell_names <- colnames(alevin_res$counts)
  ncells <- ncol(alevin_res$counts)
  
  InfRepMeanVarRes <- vector(mode = "list", length = ncells)
  for(k in 1:ncells){
    stttt <- proc.time()
    print(paste0("Calculating results for cell ", k, " out of ", ncells))
    curr_datT <- loadRData(paste0(save_dir, "alevinInfRepData/", "infRepDat", cell_names[k], ".RData"))
    curr_dat <- curr_datT
    
    curr_regular_counts <- as.matrix(regular_counts[,k, drop = FALSE])
    mean_reg_count <- mean(curr_regular_counts, na.rm = T)
    
    nboot <- ncol(curr_dat)
    gene_names <- rownames(curr_dat)
    num_genes <- length(gene_names)
    
    mean_curr <- rowMeans(curr_dat, na.rm = TRUE)
    sd_curr <- apply(curr_dat, 1, function(x){sd(x, na.rm = TRUE)})
    var_curr <- apply(curr_dat, 1, function(x){var(x, na.rm = TRUE)})
    coef_of_var_curr_SVB <- (sd_curr/(mean_curr + 5)) + 0.01
    #z_scores <- apply(curr_dat, 1, function(x){scale(x)})
    #iqr_curr <- rowIQRs(as.matrix(curr_dat), na.rm = TRUE)
    #iqr_std_curr 
    
    numerator <- var_curr - mean_curr
    numerator[numerator < 0] <- 0
    denom <- mean_curr + 5
    
    
    InfRV <- (numerator/ denom) + 0.01
    LogInfRV <- log(InfRV)
    
    Ress <- as.matrix(data.frame(infRep_Mean = mean_curr, infRep_SD = sd_curr, infRep_Var = var_curr, infRep_coef_of_var_SVB = coef_of_var_curr_SVB,
                                 InfRV = InfRV, LogInfRV = LogInfRV, mean_reg_count = mean_reg_count))
    rownames(Ress) <- gene_names
    InfRepMeanVarRes[[k]] <- Ress
    print(paste0("Time for Cell ", k, " is ", (proc.time() - stttt)[3]))
  }  
  names(InfRepMeanVarRes) <- cell_names
  
  #Now, average across all cells to get single averages for each gene
  InfRepMeanVarRes2 <- data.frame(Reduce("+", InfRepMeanVarRes)/ncells)
  InfRepMeanVarRes2$gene_id <- rownames(InfRepMeanVarRes2)
  InfRepMeanVarRes3 <- InfRepMeanVarRes2[gtools::mixedorder(InfRepMeanVarRes2$gene_id),]
  return(InfRepMeanVarRes3)
}




FulltradeSeqPipeline <- function(tradeSeqSaveDir, InputCounts, Part = NULL, 
                                 SimulatedData = FALSE, CountsFromSplatter = FALSE, current_seed,
                                 useKnownClusters = FALSE, KnownClusters = NULL,DyntoyStartMilestones = NULL, DyntoyEndMilestones = NULL,
                                 ReRunResults = FALSE){
  
  tryCatch(detach("package:tradeSeq", unload=TRUE), error = function(x){print(x)})
  library(tradeSeq, lib = "~/InstalledtradeSeq1.1.21/")
  
  #Set the seed because tradeSeq results can differ very slightly with different seeds
  set.seed(current_seed)
  ############################################################################################
  ## Code until the next ##### is (mostly) from the vignette from the slingshot package
  CompStart <- proc.time()
  
  #If known clusters are specified and the InputCount coolumns are not in the same order, reorder the input counts to be i the same order just to be safe
  if(useKnownClusters==TRUE){
    if(sum(names(KnownClusters)!=colnames(InputCounts))!=0){
      InputCounts <- InputCounts[,names(KnownClusters)]
    }
  }
  
  
  if(is.null(Part)){
    dir_mod <- ""
  }else{
    dir_mod <- paste0("InfRepNum", Part)
  }
  
  
  if(SimulatedData==TRUE){
    dir_mod2 <- "SimulatedData"
  }else if(SimulatedData==FALSE & CountsFromSplatter=="Unscaled"){
    dir_mod2 <- "CountsFromSplatterUnscaledToAlevinLibSize"
  }else{
    dir_mod2 <- ""
  }
  
  dir_mod3 <- "FilterNoGenes"
  
  
  if(useKnownClusters==TRUE){
    dir_mod4 <- "UsingKnownClusters"
  }else{
    dir_mod4 <- ""
  }
  
  if(!dir.exists(tradeSeqSaveDir)){dir.create(tradeSeqSaveDir, recursive = TRUE)}
  
  print(paste0("TradeSeq Results will be saved in ", tradeSeqSaveDir))
  
  
  curr_save_fil_name <- paste0(tradeSeqSaveDir, "tradeSeqResults", dir_mod2, dir_mod, dir_mod3, dir_mod4,".RData")
  
  
  print(paste0("The file to be saved is ", curr_save_fil_name))
  if(file.exists(curr_save_fil_name) & ReRunResults==FALSE){
    return("This file already exists, so the code will not be run again")
  }
  
  
  print(paste0(nrow(InputCounts), " genes remain after filtering (if this is too low you will want to modify the strictness of the filtering)"))
  
  #FQnorm function taken from the tradeSeq vignette
  FQnorm <- function(counts){
    rk <- apply(counts,2,rank,ties.method='min')
    counts.sort <- apply(counts,2,sort)
    refdist <- apply(counts.sort,1,median)
    norm <- apply(rk,2,function(r){ refdist[r] })
    rownames(norm) <- rownames(counts)
    return(norm)
  }
  norm_countsT <- FQnorm(InputCounts)
  
  pca <- prcomp(t(log1p(norm_countsT)), scale. = FALSE)
  
  numPCsToUse <- 5
  rd1T <- pca$x[,1:numPCsToUse]
  
  
  
  rd <- rd1T
  counts <- as.matrix(InputCounts)
  norm_counts <- norm_countsT
  
  SCEObj <- SingleCellExperiment::SingleCellExperiment(assays = S4Vectors::List(counts = counts, norm_counts = norm_counts))
  
  SingleCellExperiment::reducedDims(SCEObj) <- S4Vectors::SimpleList(PCA = rd)
  
  
  if(useKnownClusters==TRUE){
    SummarizedExperiment::colData(SCEObj)$KnownClusters <- KnownClusters
    GMM <- FALSE
    cl <- KnownClusters
  }
  
  start.clus <- DyntoyStartMilestones
  end.clus <- DyntoyEndMilestones
  if(length(end.clus)==0){
    end.clus <- NULL
  }
  
  if(useKnownClusters==TRUE){
    print("The Known clusters (as input via the argument KnownClusters) are being used and the start and end clusters are being fixed to be the values given by start.clus and end.clus given below")
    slingshot_results <- slingshot(SCEObj, clusterLabels = 'KnownClusters', reducedDim = 'PCA')
    print(paste0("Unique cluster labels are given below:"))
    print(unique(KnownClusters))
    print(paste0("The start cluster(s) are given below:"))
    print(start.clus)
    print(paste0("The end cluster(s) are given below:"))
    print(end.clus)
    lin <- getLineages(SCEObj, clusterLabels = SummarizedExperiment::colData(slingshot_results)$KnownClusters, reducedDim = 'PCA', start.clus = start.clus, end.clus = end.clus)
  }
  
  
  
  
  crv <- SlingshotDataSet(getCurves(lin))
  n_of_fit_lin <- length(crv@lineages)
  
  print(paste0("The number of lineages fit is ", n_of_fit_lin))
  print(paste0("The crv object information is printed below"))
  print(crv)
  
  
  if(is.null(nKnots)){
    nKnots <- 6
  }
  
  
  control <- mgcv::gam.control()
  control$maxit <- 1000 #set maximum number of iterations to 1000
  
  print(gc())
  print(paste0("tradeSeq results will be run with ", nKnots, " knots"))
  st1 <- proc.time()
  
  
  print(paste0("The version of tradeSeq that will be used to fit the GAM models is ", packageVersion("tradeSeq")))
  
  sce <- fitGAM(counts = counts, sds = crv, nknots = nKnots, sce = TRUE, control = control)
  
  GAMFitCompTime <- proc.time() - st1
  
  
  #browser()
  st2 <- proc.time()
  print(gc())
  
  #An important bug fix was introduced in version 1.1.23 of the package that fixes testing but unfortunately another bug that affects the fitting of the model was also introduced
  #So, unload the currently loaded version and load the newest version of the package (which is installed in the usual repo)
  pvers1 <- packageVersion("tradeSeq")
  detach("package:tradeSeq", unload=TRUE)
  library(tradeSeq)
  pvers2 <- packageVersion("tradeSeq")
  print(paste0("The originally loaded version of tradeSeq, ", pvers1, ", has been unloaded and replaced with version ", pvers2, " because the former had a bug that caused errors with calculating test stats and the latter has a bug that can cause the fitGAM statement to fail"))
  
  if(n_of_fit_lin==1){
    #With only 1 lineage the lineage specific test will be the same as the global test
    association_test_res <- associationTest(sce, global = TRUE, lineages = FALSE)
  }else{
    association_test_res <- associationTest(sce, global = TRUE, lineages = FALSE)
  }
  print(gc())
  
  end2 <- proc.time() - st2
  print(paste0("Took ", end2[3], " seconds to calculate associationTest results for all genes"))
  
  st3 <- proc.time()
  
  if(n_of_fit_lin==1){
    start_vs_end_test_res <- startVsEndTest(sce, global = TRUE, lineages = FALSE)
  }else{
    start_vs_end_test_res <- startVsEndTest(sce, global = TRUE, lineages = FALSE)
  }
  print(gc())
  
  
  end3 <- proc.time() - st3
  print(paste0("Took ", end3[3], " seconds to calculate startVsEndTest results for all genes"))
  
  
  
  
  
  if(n_of_fit_lin > 1){
    patternTestRes <- tryCatch(patternTest(sce, global = TRUE, pairwise = FALSE), error = function(x){print(x)})
    print("Calculation of patternTestRes is complete")
    diffEndTestRes <- tryCatch(diffEndTest(sce, global = TRUE, pairwise = FALSE), error = function(x){print(x)})
    print("Calculation of diffEndTestRes is complete")
    earlyDETestRes <- tryCatch(earlyDETest(sce, knots = c(1, 2), global = TRUE, pairwise = FALSE), error = function(x){print(x)})
    print("Calculation of earlyDETestRes is complete")
    
    TotalCompTime <- proc.time() - CompStart
    MemoryInfo <- gc()
    
    save(lin, crv, sce, cl, rd, counts, slingshot_results, n_of_fit_lin, association_test_res, 
         start_vs_end_test_res, diffEndTestRes, patternTestRes, earlyDETestRes, GAMFitCompTime,
         TotalCompTime, MemoryInfo, file = curr_save_fil_name)
  }else{
    
    TotalCompTime <- proc.time() - CompStart
    MemoryInfo <- gc()
    
    save(lin, crv, sce, cl, rd, counts, slingshot_results, n_of_fit_lin, association_test_res, start_vs_end_test_res, GAMFitCompTime,
         TotalCompTime, MemoryInfo, file = curr_save_fil_name)
    
  }
  
  
  
  print(gc())
  return(NULL)
}

extractTrulyDEGenes <- function(base_dir, SimNumberToUse){
  SplatterSimObject <- loadRData(paste0(base_dir, "Sim", SimNumberToUse, "/SplatterFiles/simObject.RData"))
  
  params <- SplatterSimObject@metadata$Params
  TrueDEProb <- params@de.prob
  
  GeneLevelDEParams <- data.frame(SummarizedExperiment::rowData(SplatterSimObject))
  GenesTrulyDEParams <- GeneLevelDEParams[GeneLevelDEParams$DEFacPath1!=1,]
  GenesNotTrulyDEParams <- GeneLevelDEParams[GeneLevelDEParams$DEFacPath1==1,]
  
  GenesTrulyDE <- rownames(GenesTrulyDEParams)
  GenesNotTrulyDE <- rownames(GenesNotTrulyDEParams)
  
  return(list(GenesTrulyDE = GenesTrulyDE, GenesNotTrulyDE = GenesNotTrulyDE))
}




ReturnFulltradeSeqRes <- function(base_dir, SimNumberToUse, type = NULL, infRepNum = NULL, OutliersRemoved, UseTop3000Genes = FALSE,
                                  SimulatedData, CountsFromSplatter = FALSE, RunOnMeanBootSamps = FALSE, useKnownClusters, FilterNoGenes){
  
  if(CountsFromSplatter!=FALSE){
    if(!(CountsFromSplatter %in% c("Scaled", "Unscaled"))){
      stop("Specify CountsFromSplatter as either Scaled or Unscaled")
    }
  }
  
  if(SimulatedData==TRUE & CountsFromSplatter!=FALSE){
    stop("Check function arguments, SimulatedData is true meaning the counts cannot be directly from splatter yet CountsFromSplatter if not FALSE")
  }
  
  if(is.null(type) & is.null(infRepNum)){
    stop("Specify exactly one of type and infRepNum to be non null.  Type should be used for the point estimate Alevin results and the splatter results, and infRepNum for the infRep results")
  }
  
  if(!is.null(type) & !is.null(infRepNum)){
    stop("Specify exactly one of type and infRepNum to be non null.  Type should be used for the point estimate Alevin results and the splatter results, and infRepNum for the infRep results")
  }
  
  if(!is.null(type)){
    if(!(type %in%c("Alevin", "Splatter", "SplatterUnscaled"))){
      stop('Specify type to be either "Alevin", "Splatter", "SplatterUnscaled"')
    }
  }
  
  if(CountsFromSplatter!=FALSE & RunOnMeanBootSamps!=FALSE){
    stop("Specify at most one of CountsFromSplatter and RunOnMeanBootSamps to be TRUE")
  }
  
  
  if(UseTop3000Genes==TRUE & FilterNoGenes==TRUE){
    stop("Specify at most one of UseTop3000Genes and FilterNoGenes to be TRUE")
  }

  
  if(is.null(infRepNum)){
    dir_mod <- ""
  }else{
    dir_mod <- paste0("InfRepNum", infRepNum)
  }
  
  
  if(SimulatedData==TRUE){
    dir_mod2 <- "SimulatedData"
  }else if(SimulatedData==FALSE & CountsFromSplatter=="Scaled"){
    dir_mod2 <- "CountsFromSplatter"
  }else if(SimulatedData==FALSE & CountsFromSplatter=="Unscaled"){
    dir_mod2 <- "CountsFromSplatterUnscaledToAlevinLibSize"
  }else if(SimulatedData==FALSE & RunOnMeanBootSamps==TRUE){
    dir_mod2 <- "MeanBootSamplesCounts"
  }else{
    dir_mod2 <- ""
  }
  
  if(UseTop3000Genes==TRUE){
    dir_mod3 <- "Top3000Genes"
  }else if(FilterNoGenes==TRUE){
    dir_mod3 <- "FilterNoGenes"
  }else{
    dir_mod3 <- ""
  }
  
  if(useKnownClusters==TRUE){
    dir_mod4 <- "UsingKnownClusters"
  }else{
    dir_mod4 <- ""
  }
  
  if(OutliersRemoved==TRUE){
    curr_file <- paste0(base_dir, "Sim", SimNumberToUse, "/tradeSeqResults/", "tradeSeqResults", dir_mod2, "OutliersRemoved", dir_mod, dir_mod3, dir_mod4,".RData")
  }else if(OutliersRemoved==FALSE){
    curr_file <- paste0(base_dir, "Sim", SimNumberToUse, "/tradeSeqResults/", "tradeSeqResults", dir_mod2, dir_mod, dir_mod3, dir_mod4,  ".RData")
  }
  
  if(!file.exists(curr_file)){
    print(paste0("The file ", curr_file,  "does not exist. Returning NULL"))
    return(NULL)
  }
  print(paste0("File to be loaded is ", curr_file))
  load(curr_file)
  assoc_res1 <- association_test_res[,c("waldStat", "df", "pvalue"), drop = FALSE]
  colnames(assoc_res1) <- c("waldStatAssoc", "dfAssoc", "pvalueAssoc")
  if(is.null(infRepNum)){
    assoc_res1$padjFDRAssoc <- p.adjust(assoc_res1$pvalueAssoc, method = "fdr")
  }
  assoc_res1$gene_id <- rownames(assoc_res1)
  
  #Add in the number of kMeans clusters used in the analysis as well
  if(useKnownClusters==TRUE){
    assoc_res1$nClust <- length(unique(slingshot_results@colData$KnownClusters))
  }else if(length(unique(slingshot_results@colData$kMeans))!=0){
    assoc_res1$nClust <- length(unique(slingshot_results@colData$kMeans))
  }else{
    assoc_res1$nClust <- length(unique(slingshot_results@colData$GMM))
  }
  
  
  rnames <- assoc_res1$gene_id
  stats_end_res1 <- start_vs_end_test_res[rnames,c("waldStat", "df", "pvalue"), drop = FALSE]
  colnames(stats_end_res1) <- c("waldStatStartEnd", "dfStartEnd", "pvalueStartEnd")
  if(is.null(infRepNum)){
    stats_end_res1$padjFDRStartEnd <- p.adjust(stats_end_res1$pvalueStartEnd, method = "fdr")
  }
  
  
  if(exists("diffEndTestRes")==TRUE){
    if(is.null(dim(diffEndTestRes))){
      diffEndres1 <- NULL
    }else{
      diffEndres1 <- diffEndTestRes[rnames,c("waldStat", "df", "pvalue"), drop = FALSE]
      colnames(diffEndres1) <- c("waldStatDiffEnd", "dfDiffEnd", "pvalueDiffEnd")
      if(is.null(infRepNum)){
        diffEndres1$padjFDRDiffEnd <- p.adjust(diffEndres1$pvalueDiffEnd, method = "fdr")
      }
    }
  }else{
    diffEndres1 <- NULL
  }
  
  if(exists("earlyDETestRes")==TRUE){
    if(is.null(dim(earlyDETestRes))){
      earlyDEres1 <- NULL
    }else{
      earlyDEres1 <- earlyDETestRes[rnames,c("waldStat", "df", "pvalue"), drop = FALSE]
      colnames(earlyDEres1) <- c("waldStatEarlyDE", "dfEarlyDE", "pvalueEarlyDE")
      if(is.null(infRepNum)){
        earlyDEres1$padjFDREarlyDE <- p.adjust(earlyDEres1$pvalueEarlyDE, method = "fdr")
      }
    }
  }else{
    earlyDEres1 <- NULL
  }
  
  
  if(exists("patternTestRes")==TRUE){
    if(is.null(dim(patternTestRes))){
      patternres1 <- NULL
    }else{
      patternres1 <- patternTestRes[rnames,c("waldStat", "df", "pvalue"), drop = FALSE]
      colnames(patternres1) <- c("waldStatPattern", "dfPattern", "pvaluePattern")
      if(is.null(infRepNum)){
        patternres1$padjFDRPattern <- p.adjust(patternres1$pvaluePattern, method = "fdr")
      }
    }

  }else{
    patternres1 <- NULL
  }
  
  FullResT <- bind_cols(assoc_res1, stats_end_res1, diffEndres1, earlyDEres1, patternres1)
  rownames(FullResT) <- FullResT$gene_id
  FullRes <- FullResT
  FullRes$gene_id <- NULL
  if(!is.null(type)){
    colnames(FullRes) <- paste0(colnames(FullRes), type)
  }else{
    FullRes$infRepNum <- infRepNum
  }
  
  if(is.null(n_of_fit_lin)){
    FullRes$NLin <- NA
  }else{
    FullRes$NLin <- n_of_fit_lin
  }
  
  if(sum(duplicated(colnames(FullRes)))!=0){
    stop("Check the colnames of the full set of results for duplicates, there should be no duplicates")
  }
  
  return(FullRes)
  
}




returnPvalsForiCobra <- function(SimNumberToUse, DataTypeNum, nboot, CoverageResToUse, UseTop3000Genes, base_dir, FilterNoGenes, TierRes, GeneLevelDEParams,
                                 iCobraSaveMainDir, useKnownClusters, filtered_genes){
  j <- DataTypeNum
  
  if(j==1){
    type <- "ActualBoot"
    OutliersRemoved <- FALSE
    OutStr <- ""
    SimulatedData <- FALSE
  }else if(j==2){
    type <- "SimulatedBoot"
    OutliersRemoved <- FALSE
    OutStr <- ""
    SimulatedData <- TRUE
  }
  

    dir_mod <- ""
  
  FullResAlevin <- ReturnFulltradeSeqRes(base_dir = base_dir, SimNumberToUse = SimNumberToUse, type = "Alevin", infRepNum = NULL, OutliersRemoved = OutliersRemoved,
                                         UseTop3000Genes = FALSE, SimulatedData = FALSE, CountsFromSplatter = FALSE, RunOnMeanBootSamps = FALSE, useKnownClusters = useKnownClusters,
                                         FilterNoGenes = FilterNoGenes)
  
  #We had originally additionally run tradeSeq on the splatter counts that had been scaled to have the same library size as the alevin esimates, but this is not needed so set it to NULL

  FullResSplatter <-  NULL

  FullResSplatterUnscaled <-  ReturnFulltradeSeqRes(base_dir = base_dir, SimNumberToUse = SimNumberToUse, type = "SplatterUnscaled", infRepNum = NULL, OutliersRemoved = OutliersRemoved,
                                                    UseTop3000Genes = FALSE, SimulatedData = FALSE, CountsFromSplatter = "Unscaled", RunOnMeanBootSamps = FALSE, useKnownClusters = useKnownClusters,
                                                    FilterNoGenes = FilterNoGenes)
  
  #Load in all pvalues across all inferential replicates (bootstrap or simulated pseudo-inferential reps)
  full_tradeSeqres_bootT <- vector(mode = "list", length = nboot)
  for(i in 1:nboot){
    
    if(i%%25==0){
      print(paste0("i is ", i, " out of ", nboot))
    }
    infRepNum <- i
    curr_res <- ReturnFulltradeSeqRes(base_dir = base_dir, SimNumberToUse = SimNumberToUse, type = NULL, infRepNum = infRepNum, OutliersRemoved = OutliersRemoved,
                                      UseTop3000Genes = FALSE, SimulatedData = SimulatedData, CountsFromSplatter = FALSE, RunOnMeanBootSamps = FALSE, useKnownClusters = useKnownClusters,
                                      FilterNoGenes = FilterNoGenes)
    curr_res$gene_id <- rownames(curr_res)
    full_tradeSeqres_bootT[[i]] <- curr_res
    rm(curr_res)
    rm(infRepNum)
  }
  
  
  full_tradeSeqres_bootTT2 <- data.frame(data.table::rbindlist(full_tradeSeqres_bootT, fill = T))
  full_tradeSeqres_bootT2 <- full_tradeSeqres_bootTT2[order(full_tradeSeqres_bootTT2$gene_id),]
  

  
    print(paste0("These results are for the dyntoy simulations and filter out lowly expressed genes from the power analysis.  Specifically, ", 
                 length(filtered_genes), " genes out of a possible ", nrow(FullResSplatterUnscaled), " are kept and used for the power analysis"))

    
    GenesTrulyDE <- as.character(GeneLevelDEParams$gene_id[GeneLevelDEParams$differentially_expressed==TRUE])
    GenesTrulyNotDE <- as.character(GeneLevelDEParams$gene_id[GeneLevelDEParams$differentially_expressed==FALSE])
    
    GenesTrulyDEAssociation <- GenesTrulyDE
    GenesTrulyNotDEAssociation <- GenesTrulyNotDE
    GenesTrulyDEAssociationFiltered <-  subset(GenesTrulyDEAssociation, GenesTrulyDEAssociation%in%filtered_genes)
    GenesTrulyNotDEAssociationFiltered <-  subset(GenesTrulyNotDEAssociation, GenesTrulyNotDEAssociation%in%filtered_genes)
    
    GenesTrulyDEPattern <- GenesTrulyDE
    GenesTrulyNotDEPattern <- GenesTrulyNotDE
    GenesTrulyDEPatternFiltered <-  subset(GenesTrulyDEPattern, GenesTrulyDEPattern%in%filtered_genes)
    GenesTrulyNotDEPatternFiltered <-  subset(GenesTrulyNotDEPattern, GenesTrulyNotDEPattern%in%filtered_genes)
    
    GenesTrulyDEDiffEnd <- GenesTrulyDE
    GenesTrulyNotDEDiffEnd <- GenesTrulyNotDE
    GenesTrulyDEDiffEndFiltered <-  subset(GenesTrulyDEDiffEnd, GenesTrulyDEDiffEnd%in%filtered_genes)
    GenesTrulyNotDEDiffEndFiltered <-  subset(GenesTrulyNotDEDiffEnd, GenesTrulyNotDEDiffEnd%in%filtered_genes)
    
    GenesTrulyDEEarlyDE<- GenesTrulyDE
    GenesTrulyNotDEEarlyDE <- GenesTrulyNotDE
    GenesTrulyDEEarlyDEFiltered <-  subset(GenesTrulyDEEarlyDE, GenesTrulyDEEarlyDE%in%filtered_genes)
    GenesTrulyNotDEEarlyDEFiltered <-  subset(GenesTrulyNotDEEarlyDE, GenesTrulyNotDEEarlyDE%in%filtered_genes)
    
    GenesTrulyDEStartEnd <- GenesTrulyDE
    GenesTrulyNotDEStartEnd <- GenesTrulyNotDE
    GenesTrulyDEStartEndFiltered <-  subset(GenesTrulyDEStartEnd, GenesTrulyDEStartEnd%in%filtered_genes)
    GenesTrulyNotDEStartEndFiltered <-  subset(GenesTrulyNotDEStartEnd, GenesTrulyNotDEStartEnd%in%filtered_genes)
    
  testingTypes <- c("Assoc", "StartEnd", "DiffEnd", "EarlyDE", "Pattern")


  nTestingTypes <- length(testingTypes)
  for(jjj in 1:nTestingTypes){
    curr_type <- testingTypes[jjj]
    print(paste0("Currently calculating results for testType ", jjj , " (", curr_type, ")", " out of ", nTestingTypes))
    
    if(curr_type=="Assoc"){
      AlevinPvalToUse <- "pvalueAssocAlevin"
      SplatterPvalToUse <- "pvalueAssocSplatter"
      SplatterUnscaledPvalToUse <- "pvalueAssocSplatterUnscaled"
      
      InfRepStat <- "waldStatAssoc"
      InfRepdf <- "dfAssoc"
      InfRepPval <- "pvalueAssoc"
      
      GenesTrulyDE <- GenesTrulyDEAssociationFiltered
      GenesTrulyNotDE <- GenesTrulyNotDEAssociationFiltered
      
      if(useKnownClusters==TRUE){
        iCobraSaveDir <- paste0(iCobraSaveMainDir, "AssocResUsingKnownClusters/")
      }else{
        iCobraSaveDir <- paste0(iCobraSaveMainDir, "AssocRes/")
      }
      
      if(!dir.exists(iCobraSaveDir)){dir.create(iCobraSaveDir, recursive = TRUE)}
      
    }else if(curr_type=="Pattern"){
      AlevinPvalToUse <- "pvaluePatternAlevin"
      SplatterPvalToUse <- "pvaluePatternSplatter"
      SplatterUnscaledPvalToUse <- "pvaluePatternSplatterUnscaled"
      
      InfRepStat <- "waldStatPattern"
      InfRepdf <- "dfPattern"
      InfRepPval <- "pvaluePattern"
      
      GenesTrulyDE <- GenesTrulyDEPatternFiltered
      GenesTrulyNotDE <- GenesTrulyNotDEPatternFiltered
      
      if(useKnownClusters==TRUE){
        iCobraSaveDir <- paste0(iCobraSaveMainDir, "PatternResUsingKnownClusters/")
      }else{
        iCobraSaveDir <- paste0(iCobraSaveMainDir, "PatternRes/")
      }
      
      if(!dir.exists(iCobraSaveDir)){dir.create(iCobraSaveDir, recursive = TRUE)}
    }else if(curr_type=="DiffEnd"){
      AlevinPvalToUse <- "pvalueDiffEndAlevin"
      SplatterPvalToUse <- "pvalueDiffEndSplatter"
      SplatterUnscaledPvalToUse <- "pvalueDiffEndSplatterUnscaled"
      
      InfRepStat <- "waldStatDiffEnd"
      InfRepdf <- "dfDiffEnd"
      InfRepPval <- "pvalueDiffEnd"
      
      GenesTrulyDE <- GenesTrulyDEDiffEndFiltered
      GenesTrulyNotDE <- GenesTrulyNotDEDiffEndFiltered
      
      if(useKnownClusters==TRUE){
        iCobraSaveDir <- paste0(iCobraSaveMainDir, "DiffEndResUsingKnownClusters/")
      }else{
        iCobraSaveDir <- paste0(iCobraSaveMainDir, "DiffEndRes/")
      }
      
      if(!dir.exists(iCobraSaveDir)){dir.create(iCobraSaveDir, recursive = TRUE)}
    }else if(curr_type=="EarlyDE"){
      AlevinPvalToUse <- "pvalueEarlyDEAlevin"
      SplatterPvalToUse <- "pvalueEarlyDESplatter"
      SplatterUnscaledPvalToUse <- "pvalueEarlyDESplatterUnscaled"
      
      InfRepStat <- "waldStatEarlyDE"
      InfRepdf <- "dfEarlyDE"
      InfRepPval <- "pvalueEarlyDE"
      
      GenesTrulyDE <- GenesTrulyDEEarlyDEFiltered
      GenesTrulyNotDE <- GenesTrulyNotDEEarlyDEFiltered
      
      if(useKnownClusters==TRUE){
        iCobraSaveDir <- paste0(iCobraSaveMainDir, "EarlyDEResUsingKnownClusters/")
      }else{
        iCobraSaveDir <- paste0(iCobraSaveMainDir, "EarlyDERes/")
      }
      
      if(!dir.exists(iCobraSaveDir)){dir.create(iCobraSaveDir, recursive = TRUE)}
    }else if(curr_type=="StartEnd"){
      AlevinPvalToUse <- "pvalueStartEndAlevin"
      SplatterPvalToUse <- "pvalueStartEndSplatter"
      SplatterUnscaledPvalToUse <- "pvalueStartEndSplatterUnscaled"
      
      InfRepStat <- "waldStatStartEnd"
      InfRepdf <- "dfStartEnd"
      InfRepPval <- "pvalueStartEnd"
      
      GenesTrulyDE <- GenesTrulyDEStartEndFiltered
      GenesTrulyNotDE <- GenesTrulyNotDEStartEndFiltered
      
      if(useKnownClusters==TRUE){
        iCobraSaveDir <- paste0(iCobraSaveMainDir, "StartEndResUsingKnownClusters/")
      }else{
        iCobraSaveDir <- paste0(iCobraSaveMainDir, "StartEndRes/")
      }
      
      if(!dir.exists(iCobraSaveDir)){dir.create(iCobraSaveDir, recursive = TRUE)}
    }else{
      stop(paste0("Specify a valid tradeSeq testing type"))
    }
    
    print(paste0("Directory the current results will be saved to is ", iCobraSaveDir))
    
    AlevinSplatterTradeSeqResCurrSim <- data.frame(matrix(NA, nrow = nrow(FullResAlevin), ncol = 1))
    rnames <- rownames(FullResAlevin)
    rownames(AlevinSplatterTradeSeqResCurrSim) <- rnames
    
    AlevinSplatterTradeSeqResCurrSim[,"AlevinPval"] <- FullResAlevin[,AlevinPvalToUse]

    AlevinSplatterTradeSeqResCurrSim[,"SplatterUnscaledPval"]<-  FullResSplatterUnscaled[rnames,SplatterUnscaledPvalToUse]
    
    AlevinSplatterTradeSeqResCurrSim[,"AlevinAdjP"] <- p.adjust(AlevinSplatterTradeSeqResCurrSim[,"AlevinPval"], method = "fdr")

    AlevinSplatterTradeSeqResCurrSim[,"SplatterUnscaledAdjP"] <- p.adjust(AlevinSplatterTradeSeqResCurrSim[,"SplatterUnscaledPval"], method = "fdr")
    
    AlevinSplatterTradeSeqResCurrSim[,1] <- NULL
    AlevinSplatterTradeSeqResCurrSim$gene_id <- rownames(AlevinSplatterTradeSeqResCurrSim)
    
    AlevinSplatterTradeSeqResCurrSim2 <- subset(AlevinSplatterTradeSeqResCurrSim, AlevinSplatterTradeSeqResCurrSim$gene_id %in% filtered_genes)
    
    full_tradeSeqres_boot <- subset(full_tradeSeqres_bootT2, (!is.na(full_tradeSeqres_bootT2[,InfRepStat])) & full_tradeSeqres_bootT2$gene_id %in% filtered_genes)
    
    if(SimulatedData==TRUE){
      save_fil_allpvals <- paste0(iCobraSaveDir, "AllInfRepPvaluesSim", SimNumberToUse, "SimulatedInfReps", ".RData")
    }else{
      save_fil_allpvals <- paste0(iCobraSaveDir, "AllInfRepPvaluesSim", SimNumberToUse, ".RData")
    }
    
    AllInfRepPvals <- full_tradeSeqres_boot
    AllInfRepPvals$GeneTrulyDE <- ifelse(AllInfRepPvals$gene_id %in% GenesTrulyDE, 1, 0)
    save(AllInfRepPvals, file = save_fil_allpvals)

    sumNAVals <- sum(is.na(full_tradeSeqres_bootT2[,InfRepStat]))
    if(sumNAVals!=0){
      print(paste0("Removed results from ", sumNAVals, " rows out of ", nrow(full_tradeSeqres_bootT2), " due to NA Wald Statistics"))
    }
    
    
    #Find most common DF value across the statistics for the current test/gene and use this as the DF of the Chi-squared distribution pvalues are back
      #transformed to
    MostCommonDF <- as.numeric(DescTools::Mode(full_tradeSeqres_boot[,InfRepdf], na.rm = T))
    
    InvChiSQTrans <- function(x, MostCommonDF){
      ChiSQStats <- qchisq(1 - x, df = MostCommonDF)
      m_val <- mean(ChiSQStats, na.rm = T)
      FinalP <- 1-pchisq(m_val, df=MostCommonDF)
      return(FinalP)
      
    }


    
    
    #Change any raw pvalues that show up as exactly 0 to a very small number > Machine eps since they technically aren't exactly zero
    full_tradeSeqres_boot[,InfRepPval][full_tradeSeqres_boot[,InfRepPval]==0] <- 2e-16
    WaldPvalAggregate <- do.call(data.frame, aggregate(full_tradeSeqres_boot[,InfRepPval], by = list(gene_id = full_tradeSeqres_boot$gene_id), FUN = function(x){c(InvChiSQ = InvChiSQTrans(x, MostCommonDF), Pval50Perc = quantile(x, probs = 0.50), Pval75Perc = quantile(x, probs = 0.75))}))
    colnames(WaldPvalAggregate) <- c("gene_id", "MeanStatAfterInvChiSq", "Pval50Perc", "Pval75Perc")
    WaldPvalAggregate[,"gene_id"] <- as.character(WaldPvalAggregate[,"gene_id"])
    
    
    
    FullResultsT <- WaldPvalAggregate
    FullResults <- merge(FullResultsT, AlevinSplatterTradeSeqResCurrSim2, by = "gene_id")
    
    FullResults2 <- merge(FullResults, CoverageResToUse, by = "gene_id")
    FullResults2$MeanInfRV <- NULL
    FullResults2$uniqueness_ratio <- NULL
    FullResults2$UniqFactor <- NULL
    
    FullResults2$AlevinAdjP <- NULL
    FullResults2$SplatterAdjP <- NULL
    FullResults2$SplatterUnscaledAdjP <- NULL
    rownames(FullResults2) <- FullResults2$gene_id
    
    FullResults2$gene_id <- NULL
    
    colnames(FullResults2)[colnames(FullResults2)=="AlevinPval"] <- "AlevinPointEst"
    
    colnames(FullResults2)[colnames(FullResults2)=="SplatterPval"] <- "SimulatedCountScaled"
    colnames(FullResults2)[colnames(FullResults2)=="SplatterUnscaledPval"] <- "SimulatedCount"

    
    #Final pvalue matrices are created below to additionally add in in the RATs style "pvalues"
    RawPvaluesT <- data.frame(FullResults2)
    AdjPvaluesT <- data.frame(apply(RawPvaluesT, 2, function(x) {p.adjust(x, method = "fdr")}))
    AdjPvaluesT$gene_id <- rownames(AdjPvaluesT)
    
    
    #Now, compare to RATs style of averaging across infReps, which requires x% of replicates are significant at a particular threshold AFTER FDR correction
      #This approach in general requires keeping the entire set of pvalues across infReps, in contrast to the above approach which only requires keeping
      #the final chosen raw pvalue that is then FDR adjusted
    
    #So first, we need to calculate the FDR adjusted pvalues for each infRep one rep at a time
    full_boot_res_forRATs_approachT <- vector(mode = "list", length = nboot)
    curr_adjP_RATs_approach <- paste0(InfRepPval, "AdjForRATsApproach")
    for(ir in 1:nboot){
      curr_boot_dat <- subset(full_tradeSeqres_boot, full_tradeSeqres_boot$infRepNum==ir)
      curr_boot_dat[,curr_adjP_RATs_approach] <- p.adjust(curr_boot_dat[,InfRepPval], method = "fdr")
      
      full_boot_res_forRATs_approachT[[ir]] <- curr_boot_dat
    }
    
    full_boot_res_forRATs_approach <- data.frame(data.table::rbindlist(full_boot_res_forRATs_approachT))
    
  
    RATsStylePvalAggregate <- do.call(data.frame, aggregate(full_boot_res_forRATs_approach[,curr_adjP_RATs_approach], by = list(gene_id = full_boot_res_forRATs_approach$gene_id), 
                                                            FUN = function(x){c(as.numeric(mean(x < 0.01, na.rm=T) >= 0.50),
                                                                                as.numeric(mean(x < 0.05, na.rm=T) >= 0.50),
                                                                                as.numeric(mean(x < 0.10, na.rm=T) >= 0.50),
                                                                                as.numeric(mean(x < 0.01, na.rm=T) >= 0.75),
                                                                                as.numeric(mean(x < 0.05, na.rm=T) >= 0.75),
                                                                                as.numeric(mean(x < 0.10, na.rm=T) >= 0.75))}))
    
    colnames(RATsStylePvalAggregate) <- c("gene_id", "RATsStyle50PercFDR0.01", "RATsStyle50PercFDR0.05", "RATsStyle50PercFDR0.10",
                                          "RATsStyle75PercFDR0.01", "RATsStyle75PercFDR0.05", "RATsStyle75PercFDR0.10")
                                                                                                                                                                   
    #These are adjusted pvalues already (since the adjustment is done in the loop above) to be able to create iCOBRA comparison plots for the RATs approach
      #This will only work for the "points" plots that consider FDR at the fixed thresholds of 0.01, 0.05, and 0.10
      #To do this, can can assign the FDR adjusted pvalue for a gene to be < 0.01 if the result above is significant (ie = 1) at the FDR 0.01 cutoff above,
      #assign it to be > 0.01 and < 0.05 if is it significant at FDR 0.05 but not 0.01, assign it to be > 0.05 and < 0.10 if it is signicant
      #At FDR cutoffbut not 0.05, and > 0.10 if it is not significant at any of the FDR cutoffs
    RATsStylePvalAggregate$RATsStyle50Perc <- -1
    RATsStylePvalAggregate$RATsStyle50Perc[RATsStylePvalAggregate$RATsStyle50PercFDR0.01==1 & RATsStylePvalAggregate$RATsStyle50PercFDR0.05==1 & RATsStylePvalAggregate$RATsStyle50PercFDR0.10==1] <- 0.005
    RATsStylePvalAggregate$RATsStyle50Perc[RATsStylePvalAggregate$RATsStyle50PercFDR0.01==0 & RATsStylePvalAggregate$RATsStyle50PercFDR0.05==1 & RATsStylePvalAggregate$RATsStyle50PercFDR0.10==1] <- 0.025
    RATsStylePvalAggregate$RATsStyle50Perc[RATsStylePvalAggregate$RATsStyle50PercFDR0.01==0 & RATsStylePvalAggregate$RATsStyle50PercFDR0.05==0 & RATsStylePvalAggregate$RATsStyle50PercFDR0.10==1] <- 0.075
    RATsStylePvalAggregate$RATsStyle50Perc[RATsStylePvalAggregate$RATsStyle50PercFDR0.01==0 & RATsStylePvalAggregate$RATsStyle50PercFDR0.05==0 & RATsStylePvalAggregate$RATsStyle50PercFDR0.10==0] <- 0.25
    
    
    RATsStylePvalAggregate$RATsStyle75Perc <- -1
    RATsStylePvalAggregate$RATsStyle75Perc[RATsStylePvalAggregate$RATsStyle75PercFDR0.01==1 & RATsStylePvalAggregate$RATsStyle75PercFDR0.05==1 & RATsStylePvalAggregate$RATsStyle75PercFDR0.10==1] <- 0.005
    RATsStylePvalAggregate$RATsStyle75Perc[RATsStylePvalAggregate$RATsStyle75PercFDR0.01==0 & RATsStylePvalAggregate$RATsStyle75PercFDR0.05==1 & RATsStylePvalAggregate$RATsStyle75PercFDR0.10==1] <- 0.025
    RATsStylePvalAggregate$RATsStyle75Perc[RATsStylePvalAggregate$RATsStyle75PercFDR0.01==0 & RATsStylePvalAggregate$RATsStyle75PercFDR0.05==0 & RATsStylePvalAggregate$RATsStyle75PercFDR0.10==1] <- 0.075
    RATsStylePvalAggregate$RATsStyle75Perc[RATsStylePvalAggregate$RATsStyle75PercFDR0.01==0 & RATsStylePvalAggregate$RATsStyle75PercFDR0.05==0 & RATsStylePvalAggregate$RATsStyle75PercFDR0.10==0] <- 0.25
    
    if(sum(RATsStylePvalAggregate$RATsStyle50Perc==-1 | RATsStylePvalAggregate$RATsStyle75Perc==-1) > 0){
      stop("Check the assignment of the RATsStyle Adjusted PValues")
    }
    
    RATsStylePvalAggregateAdjPvals <- RATsStylePvalAggregate[,c("gene_id", "RATsStyle50Perc", "RATsStyle75Perc")]
    
    
    AdjPvalues <- merge(AdjPvaluesT, RATsStylePvalAggregateAdjPvals, by = "gene_id")
    rownames(AdjPvalues) <- AdjPvalues$gene_id
    AdjPvalues$gene_id <- NULL
    
    #The "raw Pvalues" for the RATs style approaches won't actually be used (since they are only used in a FDR TPR comparison plot that only uses the Adj Pvalues created above)
      #but need to create non-missing entries to ensure iCOBRA will work properly so assign all raw pvalues to 0.0001 to also ensure
      #all raw pvalues are smaller than the adjusted pvalues created above
    RawPvalues <- RawPvaluesT
    RawPvalues$RATsStyle50Perc <- 0.0001
    RawPvalues$RATsStyle75Perc <- 0.0001
    
    
    CoverageResToUse2 <- subset(CoverageResToUse, CoverageResToUse$gene_id %in% filtered_genes)
    UniquenessQuantiles <- quantile(CoverageResToUse2$uniqueness_ratio, probs = seq(0,1,.01), na.rm = T)
    InfRVQuantiles <- quantile(CoverageResToUse2$MeanInfRV, probs = seq(0,1,.01), na.rm = T)
    InfRVQuantiles2 <- quantile(CoverageResToUse2$MeanInfRV, probs = c((1/3), (2/3)), na.rm = T)
    InfRVQuantiles3 <- quantile(CoverageResToUse2$MeanInfRV, probs = c(0.20,0.40,0.60,0.80), na.rm = T)
    
    minInfRV <- min(CoverageResToUse2$MeanInfRV, na.rm = T)
    maxInfRV <- max(CoverageResToUse2$MeanInfRV, na.rm = T)
    
    CoverageResToUse2$UniquenessOne <- ifelse(CoverageResToUse2$uniqueness_ratio==1, 1, 0)
    CoverageResToUse2$UniquenessOne[is.na(CoverageResToUse2$uniqueness_ratio)] <- NA
    CoverageResToUse2$UniquenessBottomTenPerc <- ifelse(CoverageResToUse2$uniqueness_ratio <= UniquenessQuantiles["10%"], 1, 0)
    CoverageResToUse2$UniquenessBottomTenPerc[is.na(CoverageResToUse2$uniqueness_ratio)] <- NA
    
    CoverageResToUse2$InfRVTopOnePercent <- ifelse(CoverageResToUse2$MeanInfRV>=InfRVQuantiles["99%"], 1, 0)
    CoverageResToUse2$InfRVTopFivePercent <- ifelse(CoverageResToUse2$MeanInfRV>=InfRVQuantiles["95%"], 1, 0)
    CoverageResToUse2$InfRVTopTenPercent <- ifelse(CoverageResToUse2$MeanInfRV>=InfRVQuantiles["90%"], 1, 0)
    
    CoverageResToUse2$InfRV <- cut(CoverageResToUse2$MeanInfRV, c(minInfRV,InfRVQuantiles3[1], InfRVQuantiles3[2], InfRVQuantiles3[3], InfRVQuantiles3[4], maxInfRV), include.lowest = T)
    
    
    TierInfoGeneLevelT <- apply(TierRes, 1, function(x){mean(x[x>0])})
    TierInfoGeneLevel <- data.frame(TierInfoGeneLevelT)
    colnames(TierInfoGeneLevel) <- "GeneLevelTier"
    TierInfoGeneLevel$gene_id <- rownames(TierInfoGeneLevel)
    TierInfoGeneLevelFilteredGenes <- subset(TierInfoGeneLevel, TierInfoGeneLevel$gene_id  %in% filtered_genes)
    TierInfoGeneLevelFilteredGenes$TierGr1Point5 <- ifelse(TierInfoGeneLevelFilteredGenes$GeneLevelTier > 1.5, 1, 0)
    
    CoverageResToUse3 <- merge(CoverageResToUse2, TierInfoGeneLevelFilteredGenes, by = "gene_id")
    
    TruthT <- CoverageResToUse2[,"gene_id", drop = F]
    TruthT$GeneTrulyDE <- ifelse(TruthT$gene_id %in% GenesTrulyDE, 1, 0)
    colnames(TruthT) <- c("gene_id", "status")
    
    Truth <- merge(TruthT, CoverageResToUse3, by = "gene_id")
    Truth$MeanInfRV <- NULL
    Truth$uniqueness_ratio <- NULL
    Truth$GeneLevelTier <- NULL
    rownames(Truth) <- Truth$gene_id
    
    if(SimNumberToUse > 101){
      Truth$gene_id <- NULL
      Truth2 <- Truth
    }else{
      SplatterObj <- loadRData(paste0(base_dir, "Sim", SimNumberToUse, "/SplatterFiles/", "simObject.RData"))
      SplatterRowDat <- as.data.frame(rowData(SplatterObj))
      GeneOutlier <- SplatterRowDat[rownames(Truth),"OutlierFactor", drop = FALSE]
      GeneOutlier$OutlierFactorGr1 <- ifelse(GeneOutlier$OutlierFactor > 1, 1, 0)
      GeneOutlier$gene_id <- rownames(GeneOutlier)
      
      Truth2 <- merge(Truth, GeneOutlier, by = "gene_id")
      Truth2$OutlierFactor <- NULL
      rownames(Truth2) <- Truth2$gene_id
      Truth2$gene_id <- NULL
      
      if(!("GeneOutlier" %in% colnames(Truth2))){
        
        Truth2$GeneOutlier <- Truth2$OutlierFactorGr1
      }
      
    }
    
    
    rm(Truth)
    rm(GenesTrulyDE)
    rm(WaldStatAggregate)
    rm(FullResultsT)
    rm(FullResults)
    rm(FullResults2)
    rm(MostCommonDF)

    cobraDat <- COBRAData(pval = RawPvalues, padj = AdjPvalues, truth = Truth2)
    
    cobraPerfOverall <- calculate_performance(cobraDat, binary_truth = "status")
    
    cobraPerfInfRV <- calculate_performance(cobraDat, binary_truth = "status", splv = "InfRV", maxsplit = Inf)
    cobraPerfUniq <- calculate_performance(cobraDat, binary_truth = "status", splv = "UniquenessBottomTenPerc", maxsplit = Inf)
    cobraPerfTier <- calculate_performance(cobraDat, binary_truth = "status", splv = "TierGr1Point5", maxsplit = Inf)
    
    cobraPerfInfRVTopTenPerc <- calculate_performance(cobraDat, binary_truth = "status", splv = "InfRVTopTenPercent", maxsplit = Inf)
    cobraPerfInfRVTopFivePerc <- calculate_performance(cobraDat, binary_truth = "status", splv = "InfRVTopFivePercent", maxsplit = Inf)
    cobraPerfInfRVTopOnePerc <- calculate_performance(cobraDat, binary_truth = "status", splv = "InfRVTopOnePercent", maxsplit = Inf)
    
    cobraPerfOutlier <- NULL
    
    
    if(SimulatedData==TRUE){
      save_fil <- paste0(iCobraSaveDir, "iCobraResSim", SimNumberToUse, "SimulatedInfReps", ".RData")
    }else{
      save_fil <- paste0(iCobraSaveDir, "iCobraResSim", SimNumberToUse, ".RData")
    }

    print(paste0("File to be saved is ", save_fil))

    iCobraRes <- list(RawPvalues = RawPvalues, AdjPvalues = AdjPvalues, Truth = Truth2, cobraPerfOverall = cobraPerfOverall,
                      cobraPerfInfRV = cobraPerfInfRV, cobraPerfUniq = cobraPerfUniq, cobraPerfOutlier = cobraPerfOutlier, cobraPerfTier = cobraPerfTier,
                      cobraPerfInfRVTopTenPerc = cobraPerfInfRVTopTenPerc, cobraPerfInfRVTopFivePerc = cobraPerfInfRVTopFivePerc, cobraPerfInfRVTopOnePerc = cobraPerfInfRVTopOnePerc)
    save(iCobraRes, file = save_fil)
    rm(cobraDat, cobraPerfOverall, cobraPerfInfRV, cobraPerfUniq, cobraPerfOutlier, cobraPerfTier, cobraPerfInfRVTopTenPerc, cobraPerfInfRVTopFivePerc, cobraPerfInfRVTopOnePerc)
    rm(iCobraRes)
    rm(RawPvalues)
    rm(AdjPvalues)
    rm(Truth2)
    
    gc()
  }
  
}




GenerateiCobraPlots <- function(SimulatedInfReps, OutliersRemoved, SimNumberToUse, TotalGenes, StatType = "", useKnownClusters, pdflatex_loc){

      UseTop3000Genes <- FALSE
      fil_mod <- ""
    
    if(SimulatedInfReps==TRUE){
      dir_m <- "SimulatedInfReps"
    }else{
      dir_m <- ""
    }
    
    if(useKnownClusters==TRUE){
      dir_mod4 <- "UsingKnownClusters"
    }else{
      dir_mod4 <- ""
    }
      
      if(StatType=="Assoc"){
        curr_fil <- paste0(save_dir, "AssocRes", dir_mod4, "/iCobraResSim", SimNumberToUse, dir_m, fil_mod, ".RData")
        print(paste0("Current file loaded is ", curr_fil))
        resToUse <- loadRData(curr_fil)
        
        plotNam <- "Overall Association Test"
      }else if(StatType=="Pattern"){
        curr_fil <- paste0(save_dir, "PatternRes", dir_mod4, "/iCobraResSim", SimNumberToUse, dir_m, fil_mod, ".RData")
        print(paste0("Current file loaded is ", curr_fil))
        resToUse <- loadRData(curr_fil)
        
        plotNam <- "Overall Pattern Test"
      }else if(SimNumberToUse%in% 5:14){
        stop("Specify StatType as Assoc or Pattern for simulations 5 to 14")
      }else if(StatType=="DiffEnd"){
        curr_fil <- paste0(save_dir, "DiffEndRes", dir_mod4, "/iCobraResSim", SimNumberToUse, dir_m, fil_mod, ".RData")
        print(paste0("Current file loaded is ", curr_fil))
        resToUse <- loadRData(curr_fil)
        
        plotNam <- "Overall DiffEnd Test"
      }else if(StatType=="EarlyDE"){
        curr_fil <- paste0(save_dir, "EarlyDERes", dir_mod4, "/iCobraResSim", SimNumberToUse, dir_m, fil_mod, ".RData")
        print(paste0("Current file loaded is ", curr_fil))
        resToUse <- loadRData(curr_fil)
        
        plotNam <- "Overall EarlyDE Test"
      }else if(StatType=="StartEnd"){
        curr_fil <- paste0(save_dir, "StartEndRes", dir_mod4, "/iCobraResSim", SimNumberToUse, dir_m, fil_mod, ".RData")
        print(paste0("Current file loaded is ", curr_fil))
        resToUse <- loadRData(curr_fil)
        
        plotNam <- "Overall StartEnd Test"
      }
      
      lineageType <- "Trifurcating"
      
      if(SimNumberToUse %in%c(117,10117)){
        ncells <- 100
      }else if(SimNumberToUse %in%c(118,10118)){
        ncells <- 250
      }else{
        stop("unable to properly determine number of cells for this simulation")
      }
      
    
    
    
    
    
    
    Truth <- resToUse$Truth
    
    if(mean(Truth$status)==0){
      perc_genes_ch <- 0
    }else{
      perc_genes_ch <- fr(100*mean(Truth$status), 0)
    }
    

      cobraPerf <- resToUse$cobraPerfInfRV
      cobraPerf2 <- resToUse$cobraPerfUniq
      cobraPerf4 <- resToUse$cobraPerfInfRVTopTenPerc
    
      legendpos <- "right"

    if("AlevinPointEst" %in% colnames(resToUse$RawPvalues)){
      Al_nam <- "AlevinPointEst"
    }else{
      Al_nam <- "Alevin"
    }
    
    if("SplatterUnscaled" %in% colnames(resToUse$RawPvalues)){
      Sim_cnt_nam <- "SplatterUnscaled"
    }else if("SimulatedCount" %in% colnames(resToUse$RawPvalues)){
      Sim_cnt_nam <- "SimulatedCount"
    }else{
      Sim_cnt_nam <- "SimulatedCountUnscaled"
    }
    
    print(paste0("Sim_cnt_nam is ", Sim_cnt_nam))
    
    rgb2 <- function(x,y,z) rgb(x,y,z,maxColorValue=255)
    keep_methods <- c(Al_nam, "MeanStatAfterInvChiSq", "Pval50Perc", "Pval75Perc", Sim_cnt_nam)

    colorss <- c("black", "blue", rgb2(86,180,233), "green4", rgb2(213,94,0))
    
    if(!(all(keep_methods%in%colnames(resToUse$RawPvalues)))){
      stop("The specifications of the methods to keep is incorrect")
    }
    cobraResForPlotting <- prepare_data_for_plot(cobraPerf, keepmethods = keep_methods, colorscheme = colorss)
    cobraResForPlotting2 <- prepare_data_for_plot(cobraPerf2, keepmethods = keep_methods, colorscheme = colorss)
    cobraResForPlotting4 <- prepare_data_for_plot(cobraPerf4, keepmethods = keep_methods, colorscheme = colorss)

    
    ngenes1 <- max(cobraResForPlotting@fdrtpr$TOT_CALLED)
    ngenes2 <- max(cobraResForPlotting2@fdrtpr$TOT_CALLED)
    ngenes4 <- max(cobraResForPlotting4@fdrtpr$TOT_CALLED)

    
    maxFDR <- max(cobraResForPlotting@fdrtpr$FDR) + 0.05
    maxFDR2 <- max(cobraResForPlotting2@fdrtpr$FDR) + 0.05
    maxFDR4 <- max(cobraResForPlotting4@fdrtpr$FDR) + 0.05

    
    minTPRT <- min(cobraResForPlotting@fdrtpr$TPR)
    minTPR2T <- min(cobraResForPlotting2@fdrtpr$TPR)
    minTPR4T <- min(cobraResForPlotting4@fdrtpr$TPR)

    
    minTPR <- max(0, minTPRT - 0.025)
    minTPR2 <- max(0, minTPR2T - 0.025)
    minTPR4 <- max(0, minTPR4T - 0.025)

    
    maxTPRT <- max(cobraResForPlotting@fdrtpr$TPR)
    maxTPR2T <- max(cobraResForPlotting2@fdrtpr$TPR)
    maxTPR4T <- max(cobraResForPlotting4@fdrtpr$TPR)

    
    maxTPR <- min(1.005, maxTPRT + 0.025)
    maxTPR2 <- min(1.005, maxTPR2T + 0.025)
    maxTPR4 <- min(1.005, maxTPR4T + 0.025)

    
    
    xlimsInfRVPlot <- c(0, maxFDR)
    ylimsInfRVPlot <- c(minTPR, maxTPR)
    
    xlimsUniqPlot <- c(0, maxFDR2)
    ylimsUniqPlot <- c(minTPR2, maxTPR2)
    

    xlimsInfRVPlot4 <- c(0, maxFDR4)
    ylimsInfRVPlot4 <- c(minTPR4, maxTPR4)
    
    FullPlotTitle <- paste0(plotNam, " For ", lineageType, " Lineage with ", ncells, " Cells")
    
    PlotByInfRVT <- plot_fdrtprcurve(cobraResForPlotting, xaxisrange = xlimsInfRVPlot, yaxisrange = ylimsInfRVPlot, plottype = "points")
    PlotByUniquenessT <- plot_fdrtprcurve(cobraResForPlotting2, xaxisrange = xlimsUniqPlot, yaxisrange = ylimsUniqPlot, plottype = "points")
    
    PlotByInfRVTopTenPercT <- plot_fdrtprcurve(cobraResForPlotting4, xaxisrange = xlimsInfRVPlot4, yaxisrange = ylimsInfRVPlot4, plottype = "points")

    
    PlotByInfRV <- PlotByInfRVT + ggtitle(FullPlotTitle)+
      theme(legend.text = element_text(size=rel(1), face = "bold"),legend.title = element_text(size=rel(1), face = "bold"),axis.text.y=element_text(size=rel(1.25), face = "bold"),
            axis.text.x=element_text(size=rel(1.25), face = "bold"), plot.title=element_text(size=rel(1.25), hjust=.5,face = "bold"), axis.title=element_text(size=rel(1.25), face = "bold"),
            legend.position = legendpos, strip.text = element_text(size = rel(1.10)))+
      facet_wrap(~splitval, nrow = 2)
    
    PlotByUniqueness <- PlotByUniquenessT + ggtitle(FullPlotTitle)+
      theme(legend.text = element_text(size=rel(1), face = "bold"),legend.title = element_text(size=rel(1), face = "bold"),axis.text.y=element_text(size=rel(1.25), face = "bold"),
            axis.text.x=element_text(size=rel(1.25), face = "bold"), plot.title=element_text(size=rel(1.25), hjust=.5,face = "bold"), axis.title=element_text(size=rel(1.25), face = "bold"),
            legend.position = legendpos, strip.text = element_text(size = rel(1.10)))+
      facet_wrap(~splitval, nrow = 2)
    
    PlotByInfRVTopTenPerc <- PlotByInfRVTopTenPercT + ggtitle(FullPlotTitle)+
      theme(legend.text = element_text(size=rel(1), face = "bold"),legend.title = element_text(size=rel(1), face = "bold"),axis.text.y=element_text(size=rel(1.25), face = "bold"),
            axis.text.x=element_text(size=rel(1.25), face = "bold"), plot.title=element_text(size=rel(1.25), hjust=.5,face = "bold"), axis.title=element_text(size=rel(1.25), face = "bold"),
            legend.position = legendpos, strip.text = element_text(size = rel(1.10)))+
      facet_wrap(~splitval, nrow = 2)
    
    
    #Now, plot that compares RATs style results across infReps to the approach we propose in the paper
    methodsRATsPlots <- c("Pval50Perc", "Pval75Perc", "RATsStyle50Perc", "RATsStyle75Perc")
    colorssRATs <- c("black", "red")
    cobraResForPlotting50 <- prepare_data_for_plot(cobraPerf, keepmethods = c("Pval50Perc", "RATsStyle50Perc"), colorscheme = colorssRATs)
    cobraResForPlotting75 <- prepare_data_for_plot(cobraPerf, keepmethods = c("Pval75Perc", "RATsStyle75Perc"), colorscheme = colorssRATs)
    
    PlotByInfRVRATs50T <- plot_fdrtprcurve(cobraResForPlotting50, xaxisrange = xlimsInfRVPlot, yaxisrange = ylimsInfRVPlot, plottype = "points")
    PlotByInfRVRATs75T <- plot_fdrtprcurve(cobraResForPlotting75, xaxisrange = xlimsInfRVPlot, yaxisrange = ylimsInfRVPlot, plottype = "points")
    
    PlotByInfRVRATs50 <- PlotByInfRVRATs50T + ggtitle(FullPlotTitle)+
      theme(legend.text = element_text(size=rel(1), face = "bold"),legend.title = element_text(size=rel(1), face = "bold"),axis.text.y=element_text(size=rel(1.25), face = "bold"),
            axis.text.x=element_text(size=rel(1.25), face = "bold"), plot.title=element_text(size=rel(1.25), hjust=.5,face = "bold"), axis.title=element_text(size=rel(1.25), face = "bold"),
            legend.position = legendpos, strip.text = element_text(size = rel(1.10)))+
      facet_wrap(~splitval, nrow = 2)
    
    PlotByInfRVRATs75 <- PlotByInfRVRATs75T + ggtitle(FullPlotTitle)+
      theme(legend.text = element_text(size=rel(1), face = "bold"),legend.title = element_text(size=rel(1), face = "bold"),axis.text.y=element_text(size=rel(1.25), face = "bold"),
            axis.text.x=element_text(size=rel(1.25), face = "bold"), plot.title=element_text(size=rel(1.25), hjust=.5,face = "bold"), axis.title=element_text(size=rel(1.25), face = "bold"),
            legend.position = legendpos, strip.text = element_text(size = rel(1.10)))+
      facet_wrap(~splitval, nrow = 2)
    

    print(paste0("Current PlotsSaveDir is ", PlotsSaveDir))
    
    curr_plot_dir <- paste0(PlotsSaveDir, "PlotsByInfRV/", "Sim", SimNumberToUse, "/", StatType, "/")
    if(!dir.exists(curr_plot_dir)){dir.create(curr_plot_dir, recursive = T)}
    
    fil_piece <- paste0("Sim", SimNumberToUse, "PlotByInfRV", "SimReps", SimulatedInfReps, "OutliersRemoved", OutliersRemoved, "UseTop3000Genes", UseTop3000Genes)
    fil_name <- paste0(curr_plot_dir, fil_piece, ".tex")
    #fil_name <- paste0(curr_plot_dir, fil_piece, ".pdf")
    tikz(file = fil_name, height=6, width=8.5, standAlone = TRUE, sanitize = F)
    #pdf(file = fil_name, height=5, width=8.5)
    print(PlotByInfRV)
    dev.off()
    # tryCatch(dev.off(), error = function(x){})
    SVBCompiletikzPlot(curr_plot_dir, fil_piece, pdflatex_loc = pdflatex_loc)
    rm(fil_name)
    rm(fil_piece)
    rm(curr_plot_dir)
    
    curr_plot_dir <- paste0(PlotsSaveDir, "PlotsByInfRVRATs/", "Sim", SimNumberToUse, "/", StatType, "/")
    if(!dir.exists(curr_plot_dir)){dir.create(curr_plot_dir, recursive = T)}
    
    fil_piece <- paste0("Sim", SimNumberToUse, "PlotByInfRVRATs50", "SimReps", SimulatedInfReps, "OutliersRemoved", OutliersRemoved, "UseTop3000Genes", UseTop3000Genes)
    fil_name <- paste0(curr_plot_dir, fil_piece, ".tex")
    #fil_name <- paste0(curr_plot_dir, fil_piece, ".pdf")
    tikz(file = fil_name, height=6, width=8.5, standAlone = TRUE, sanitize = F)
    #pdf(file = fil_name, height=5, width=8.5)
    print(PlotByInfRVRATs50)
    dev.off()
    # tryCatch(dev.off(), error = function(x){})
    SVBCompiletikzPlot(curr_plot_dir, fil_piece, pdflatex_loc = pdflatex_loc)
    rm(fil_name)
    rm(fil_piece)
    
    fil_piece <- paste0("Sim", SimNumberToUse, "PlotByInfRVRATs75", "SimReps", SimulatedInfReps, "OutliersRemoved", OutliersRemoved, "UseTop3000Genes", UseTop3000Genes)
    fil_name <- paste0(curr_plot_dir, fil_piece, ".tex")
    #fil_name <- paste0(curr_plot_dir, fil_piece, ".pdf")
    tikz(file = fil_name, height=6, width=8.5, standAlone = TRUE, sanitize = F)
    #pdf(file = fil_name, height=5, width=8.5)
    print(PlotByInfRVRATs75)
    dev.off()
    # tryCatch(dev.off(), error = function(x){})
    SVBCompiletikzPlot(curr_plot_dir, fil_piece, pdflatex_loc = pdflatex_loc)
    rm(fil_name)
    rm(fil_piece)
    rm(curr_plot_dir)
    
    
    curr_plot_dir <- paste0(PlotsSaveDir, "PlotsByGeneUniqueness/", "Sim",SimNumberToUse, "/", StatType, "/")
    if(!dir.exists(curr_plot_dir)){dir.create(curr_plot_dir, recursive = T)}
    
    fil_piece <- paste0("Sim", SimNumberToUse,"PlotsByGeneUniqueness", "SimReps", SimulatedInfReps, "OutliersRemoved", OutliersRemoved, "UseTop3000Genes", UseTop3000Genes)
    fil_name <- paste0(curr_plot_dir, fil_piece, ".tex")
    #fil_name <- paste0(curr_plot_dir, fil_piece, ".pdf")
    tikz(file = fil_name, height=6, width=8.5, standAlone = TRUE, sanitize = F)
    #pdf(file = fil_name, height=5, width=8.5)
    print(PlotByUniqueness)
    dev.off()
    tryCatch(dev.off(), error = function(x){})
    SVBCompiletikzPlot(curr_plot_dir, fil_piece, pdflatex_loc = pdflatex_loc)
    rm(fil_name)
    rm(fil_piece)
    rm(curr_plot_dir)
    
    curr_plot_dir <- paste0(PlotsSaveDir, "PlotsByInfRVTopTenPerc/", "Sim", SimNumberToUse, "/", StatType, "/")
    if(!dir.exists(curr_plot_dir)){dir.create(curr_plot_dir, recursive = T)}
    
    fil_piece <- paste0("Sim", SimNumberToUse, "PlotByInfRVTopTenPerc", "SimReps", SimulatedInfReps, "OutliersRemoved", OutliersRemoved, "UseTop3000Genes", UseTop3000Genes)
    fil_name <- paste0(curr_plot_dir, fil_piece, ".tex")
    #fil_name <- paste0(curr_plot_dir, fil_piece, ".pdf")
    tikz(file = fil_name, height=6, width=8.5, standAlone = TRUE, sanitize = F)
    #pdf(file = fil_name, height=5, width=8.5)
    print(PlotByInfRVTopTenPerc)
    dev.off()
    # tryCatch(dev.off(), error = function(x){})
    SVBCompiletikzPlot(curr_plot_dir, fil_piece, pdflatex_loc = pdflatex_loc)
    rm(fil_name)
    rm(fil_piece)
    rm(curr_plot_dir)
    
    
}



SVBCompiletikzPlot <- function(Dir = NULL, fil_piece, pdflatex_loc = NULL){
  

  Dir2 <- shQuote(Dir)
  fil_name <- paste0(Dir, fil_piece, ".tex")
  
  
  if(!file.exists(fil_name)){
    stop(paste0("The file name specified does not exist in the directory specified.  Current file specified is ", fil_name))
  }
  
  #Use shQuote to ensure any and file paths or file names with spaces are properly accessed
    #This needs to be Dir and not Dir2 so the same file doesn't have multiple quotes in it
  fil_name2 <- shQuote(paste0(Dir, fil_piece, ".tex"))
  
  if(is.null(pdflatex_loc)){
    system(paste0("pdflatex "," -interaction=nonstopmode -output-directory=", Dir2, " ", fil_name2), ignore.stdout = TRUE)
  }else{
    if(!file.exists(pdflatex_loc)){
      stop("The specified pdflatex_loc does not exist")
    }
    system(paste0(pdflatex_loc," -interaction=nonstopmode -output-directory=", Dir2, " ", fil_name2), ignore.stdout = TRUE)
  }

  
  extensions_to_delete <- c(".aux", ".log", ".tex")
  for(zz in 1:length(extensions_to_delete)){
    curr_ext <- extensions_to_delete[zz]
    curr_file <- paste0(Dir, fil_piece, curr_ext)
    curr_file2 <- shQuote(paste0(Dir, fil_piece, curr_ext))
    if(file.exists(curr_file)){
      system(paste0("rm ", curr_file2))
    }
  }
  
}


fr <- function(x, n = getOption("digits"), forLatex = TRUE, NLS = FALSE) {
  y <- rep(0, length(x))
  for (i in 1:length(x)) {
    if (x[i] >=0 & x[i]< 10^-(n) & !is.na(x[i])){
      
      if(NLS==TRUE){
        y[i] <- paste("0", ".", paste(rep("0", n), collapse=""), sep = "")
      }else if(forLatex==TRUE){
        y[i] <- paste("$<$", ".", paste(rep("0", n-1), collapse=""), "1", sep = "")
      }else{
        y[i] <- paste("<", ".", paste(rep("0", n-1), collapse=""), "1", sep = "")
      }
      
    }
    if (is.numeric(x[i]) & !is.na(x[i]) &  (x[i] >= 10^-(n) | x[i] < 0)) {
      y[i] <- format(round(x[i], digits = n), nsmall = n, scientific = FALSE)
    }
    # Only truly missing values should ever show up as 0 in y created above, so change them back to missing
    #Note that this is not the same as having a value of 0 in the input vector (x), as these will correctly 
    #Just be formatted to show as <0.0001, etc
    if(y[i] == "0") {
      y[i] <- " "
    }
  } 
  return(y)
}





GenerateCoveragePlots <- function(CoverageVarName, PlotHeight, pdflatex_loc, PlotsSaveDir, PlotNameModifier, ZerosExcludedFromCoverageModifier, SimNumberToUse){
  stdal <- T
  xlim1 <- c(0.5,5.5)
  ylim1 <- c(0,1.1)
  
  if(PlotNameModifier=="ForPaper"){
    PlotWidth <- 10/2
  }else{
    PlotWidth <- 10 
  }

  
  curr_save_dir <- paste0(PlotsSaveDir, CoverageVarName, "/")
  if(!dir.exists(curr_save_dir)){dir.create(curr_save_dir, recursive = T)}
  setwd(curr_save_dir)
  print(paste0("Current Save directory is ", curr_save_dir))
  
  CoverageVarNamePrint <- strsplit(CoverageVarName, "Coverage")[[1]][1]
  if(CoverageVarNamePrint=="Norm95ConfInv"){
    CoverageVarNamePrint2 <- "95\\% Intervals From Normal Distribution Quantiles (Top)\nand from Empirical Bootstrap Quantiles (Bottom)\n"
  }else if(CoverageVarNamePrint=="TwoPointFiveNinetySevenPointFive"){
    #CoverageVarNamePrint2 <- "95\\% Interval From Bootstrap Empirical Distribution Quantiles"
    CoverageVarNamePrint2 <- "\n\n"
  }else if(CoverageVarNamePrint=="FiveNinetyFive"){
    CoverageVarNamePrint2 <- "90\\% Interval From Bootstrap Empirical Distribution Quantiles\n"
  }else if(CoverageVarNamePrint=="MinMax"){
    CoverageVarNamePrint2 <- "100\\% Interval From Bootstrap Empirical Distribution Quantiles\n"
  }else if(CoverageVarNamePrint=="NB95Inv"){
    CoverageVarNamePrint2 <- "95\\% Intervals From Negative Binomial Distribution Quantiles (Top)\nand from Empirical Bootstrap Quantiles (Bottom)\n"
  }else{
    CoverageVarNamePrint2 <- CoverageVarNamePrint
  }
  
  ####################################################################
  ###Old non ggplots
  # #paste0(CoverageVarNamePrint, " Coverage by Gene Uniqueness Level For All Genes (", nrow(gene_level_coverage3), " Total)")
  # #paste0("Coverage of ", CoverageVarNamePrint2, "\nby InfRV and Gene Expression for Filtered Genes (", nrow(gene_level_coverage3_filtered_genes), " Total)")
  # fil_piece <- "GeneUniqVsCoverageAllGenes"
  # fil_name <- paste0(curr_save_dir, fil_piece, ".tex")
  # tikz(file = fil_name, height=PlotHeight, width=PlotWidth, standAlone = stdal, sanitize = F)
  # boxplot(gene_level_coverage3[,CoverageVarName] ~ gene_level_coverage3$UniqFactor, outline = F, ylab = "Coverage", xlab = "Uniqueness Ratio", cex.main = 1.5, xlim = xlim1, ylim = ylim1, cex.lab = 1.5, cex.axis = 1.5, ann=TRUE,
  #         main = paste0("Coverage of ",  CoverageVarNamePrint2, "\n by Gene Uniqueness Quintile For All Genes (", nrow(gene_level_coverage3), " Total)"))
  # dev.off()
  # #print(paste0("pdflatex ","-output-directory=", curr_save_dir, " ", fil_name))
  # system(paste0("pdflatex ","-output-directory=", curr_save_dir, " ", fil_name), ignore.stdout = TRUE)
  # system(paste0("rm ", paste0(curr_save_dir, fil_piece), ".aux"))
  # system(paste0("rm ", paste0(curr_save_dir, fil_piece), ".log"))
  # system(paste0("rm ", paste0(curr_save_dir, fil_piece), ".tex"))
  # rm(fil_name)
  # rm(fil_piece)
  
  # 
  # 
  # 
  # fil_piece <- "GeneUniqVsCoverageFilteredGenes"
  # fil_name <- paste0(curr_save_dir, fil_piece, ".tex")
  # tikz(file = fil_name, height=PlotHeight, width=PlotWidth, standAlone = stdal, sanitize = F)
  # boxplot(gene_level_coverage3_filtered_genes[,CoverageVarName] ~ gene_level_coverage3_filtered_genes$UniqFactor, outline = F, ylab = "Coverage", xlab = "Uniqueness Ratio", cex.main = 1.5, xlim = xlim1, ylim = ylim1, cex.lab = 1.5, cex.axis = 1.5, ann=TRUE,
  #         main = paste0(CoverageVarNamePrint, " Coverage by 
  #         Gene Uniqueness Level For Filtered Genes (", nrow(gene_level_coverage3_filtered_genes), " Total)"))
  # dev.off()
  # system(paste0("pdflatex ","-output-directory=", curr_save_dir, " ", fil_name), ignore.stdout = TRUE)
  # system(paste0("rm ", paste0(curr_save_dir, fil_piece), ".aux"))
  # system(paste0("rm ", paste0(curr_save_dir, fil_piece), ".log"))
  # system(paste0("rm ", paste0(curr_save_dir, fil_piece), ".tex"))
  # rm(fil_name)
  # rm(fil_piece)
  # 
  # fil_piece <- "GeneInfRVVsCoverageAllGenes"
  # fil_name <- paste0(curr_save_dir, fil_piece, ".tex")
  # tikz(file = fil_name, height=PlotHeight, width=PlotWidth, standAlone = stdal, sanitize = F)
  # xlim1 <- c(0.5,2.5)
  # ylim1 <- c(0,1.1)
  # boxplot(gene_level_coverage3[,CoverageVarName] ~ gene_level_coverage3$HighInfRV, outline = F, ylab = "Coverage", xlab = "Inferential Uncertainty Level", cex.main = 1.5, xlim = xlim1, ylim = ylim1, cex.lab = 1.5, cex.axis = 1.5, ann=TRUE,
  #         main = paste0(CoverageVarNamePrint, " Coverage by InfRV 
  #         Level For All Genes (", nrow(gene_level_coverage3), " Total)"))
  # dev.off()
  # system(paste0("pdflatex ","-output-directory=", curr_save_dir, " ", fil_name), ignore.stdout = TRUE)
  # system(paste0("rm ", paste0(curr_save_dir, fil_piece), ".aux"))
  # system(paste0("rm ", paste0(curr_save_dir, fil_piece), ".log"))
  # system(paste0("rm ", paste0(curr_save_dir, fil_piece), ".tex"))
  # rm(fil_name)
  # rm(fil_piece)
  # 
  # fil_piece <- "GeneInfRVVsCoverageFilteredGenes"
  # fil_name <- paste0(curr_save_dir, fil_piece, ".tex")
  # tikz(file = fil_name, height=PlotHeight, width=PlotWidth, standAlone = stdal, sanitize = F)
  # xlim1 <- c(0.5,2.5)
  # ylim1 <- c(0,1.1)
  # boxplot(gene_level_coverage3_filtered_genes[,CoverageVarName] ~ gene_level_coverage3_filtered_genes$HighInfRV, outline = F, ylab = "Coverage", xlab = "Inferential Uncertainty Level", cex.main = 1.5, xlim = xlim1, ylim = ylim1, cex.lab = 1.5, cex.axis = 1.5, ann=TRUE,
  #         main = paste0(CoverageVarNamePrint, " Coverage by InfRV 
  #         Level For Filtered Genes (", nrow(gene_level_coverage3_filtered_genes), " Total)"))
  # dev.off()
  # system(paste0("pdflatex ","-output-directory=", curr_save_dir, " ", fil_name), ignore.stdout = TRUE)
  # system(paste0("rm ", paste0(curr_save_dir, fil_piece), ".aux"))
  # system(paste0("rm ", paste0(curr_save_dir, fil_piece), ".log"))
  # system(paste0("rm ", paste0(curr_save_dir, fil_piece), ".tex"))
  # rm(fil_name)
  # rm(fil_piece)
  # 
  # fil_piece <- "GeneExprVsCoverageAllGenes"
  # fil_name <- paste0(curr_save_dir, fil_piece, ".tex")
  # tikz(file = fil_name, height=PlotHeight, width=PlotWidth, standAlone = stdal, sanitize = F)
  # #pdf(file = paste0(curr_save_dir, "/GeneExprVsCoverageAllGenes.pdf"), height=PlotHeight, width=PlotWidth)
  # xlim1 <- c(0.5,2.5)
  # ylim1 <- c(0,1.1)
  # boxplot(gene_level_coverage3[,CoverageVarName] ~ gene_level_coverage3$HighExpression, outline = F, ylab = "Coverage", xlab = "Average Gene-Level Expression", cex.main = 1.5, xlim = xlim1, ylim = ylim1, cex.lab = 1.5, cex.axis = 1.5, ann=TRUE,
  #         main = paste0(CoverageVarNamePrint, " Coverage by Expression 
  #         Level For All Genes (", nrow(gene_level_coverage3), " Total)"))
  # dev.off()
  # system(paste0("pdflatex ","-output-directory=", curr_save_dir, " ", fil_name), ignore.stdout = TRUE)
  # system(paste0("rm ", paste0(curr_save_dir, fil_piece), ".aux"))
  # system(paste0("rm ", paste0(curr_save_dir, fil_piece), ".log"))
  # system(paste0("rm ", paste0(curr_save_dir, fil_piece), ".tex"))
  # rm(fil_name)
  # rm(fil_piece)
  # 
  # fil_piece <- "GeneExprVsCoverageFilteredGenes"
  # fil_name <- paste0(curr_save_dir, fil_piece, ".tex")
  # tikz(file = fil_name, height=PlotHeight, width=PlotWidth, standAlone = stdal, sanitize = F)
  # xlim1 <- c(0.5,2.5)
  # ylim1 <- c(0,1.1)
  # boxplot(gene_level_coverage3_filtered_genes[,CoverageVarName] ~ gene_level_coverage3_filtered_genes$HighExpression, outline = F, ylab = "Coverage", xlab = "Average Gene-Level Expression", cex.main = 1.5, xlim = xlim1, ylim = ylim1, cex.lab = 1.5, cex.axis = 1.5, ann=TRUE,
  #         main = paste0(CoverageVarNamePrint, " Coverage by Expression 
  #         Level For Filtered Genes (", nrow(gene_level_coverage3_filtered_genes), " Total)"))
  # dev.off()
  # system(paste0("pdflatex ","-output-directory=", curr_save_dir, " ", fil_name), ignore.stdout = TRUE)
  # system(paste0("rm ", paste0(curr_save_dir, fil_piece), ".aux"))
  # system(paste0("rm ", paste0(curr_save_dir, fil_piece), ".log"))
  # system(paste0("rm ", paste0(curr_save_dir, fil_piece), ".tex"))
  # rm(fil_name)
  # rm(fil_piece)
  
  # PlotHeight <- 5.25
  # PlotWidth <- 9
  # tikz(file = paste0(curr_save_dir, "/CoverageInfRVAndExpressionAllGenes.tex"), height=PlotHeight, width=PlotWidth, standAlone = stdal, sanitize = F)
  # ggplot(data=gene_level_coverage3,
  #        aes(x=HighInfRV,y=gene_level_coverage3_filtered_genes[,CoverageVar],fill=HighExpression))+geom_boxplot()+ggtitle(paste0(CoverageVar, " Coverage by InfRV and Expression for All Genes (", nrow(gene_level_coverage3), " Total)"))+
  #   xlab("InfRV")+ylab("Coverage")+
  #   theme(plot.title = element_text(hjust = 0.5))
  # dev.off()
  
  # tikz(file = paste0(curr_save_dir, "/CoverageInfRVAndExpressionFilteredGenes.tex"), height=PlotHeight, width=PlotWidth, standAlone = stdal, sanitize = F)
  # ggplot(data=gene_level_coverage3_filtered_genes,
  #        aes(x=HighInfRV,y=gene_level_coverage3_filtered_genes[,CoverageVar],fill=HighExpression))+geom_boxplot()+ggtitle(paste0(CoverageVar, " Coverage by InfRV and Expression for Filtered Genes (", nrow(gene_level_coverage3_filtered_genes), " Total)"))+
  #   xlab("InfRV")+ylab("Coverage")+
  #   theme(plot.title = element_text(hjust = 0.5))
  # dev.off()
  
  #Set width for ggplot plots (height is set at the top of the file)
  
  if(GeneratePlotsForPaper==TRUE){
    PlotWidth <- 9/2
    PlotWidthUniqueness <- 9
  }else{
    PlotWidth <- 9
    PlotWidthUniqueness <- 9
  }

  
  

  
  #Extract corresponding elements for all genes
  HighInfRV <- gene_level_coverage3[,"HighInfRV"]
  CoverageVar <- gene_level_coverage3[,CoverageVarName]
  HighExpression <- gene_level_coverage3[,"HighExpression"]
  
  #Do include genes with a NA tier value in the tier value boxplots
  TierValue <- gene_level_coverage3$TierFactor[!is.na(gene_level_coverage3$TierFactor)]
  CoverageVarTier <- gene_level_coverage3[,CoverageVarName][!is.na(gene_level_coverage3$TierFactor)]
  
  #Generate sample sizes corresponding to each of the 4 separate boxplots
  n1 <- sum(HighInfRV=="LowInfRV" & HighExpression=="LowExpr")
  n2 <- sum(HighInfRV=="LowInfRV" & HighExpression=="HighExpr")
  n3 <- sum(HighInfRV=="HighInfRV" & HighExpression=="LowExpr")
  n4 <- sum(HighInfRV=="HighInfRV" & HighExpression=="HighExpr")
  
  
  #Custom labels using the sample sizes calcualated above
  if(PlotNameModifier=="ForPaper"){
    my_xlab <- c(paste0("Low InfRV","\nLow Expr ; High Expr", "\n$\\bm{n=", n1, " ; ", n2, "}$"), paste0("High InfRV","\nLow Expr ; High Expr","\n$\\bm{n=", n3, " ; ", n4, "}$"))
  }else{
    my_xlab <- c(paste0("Low InfRV","\n$\\bm{n=", n1, ";", n2, "}$"), paste0("High InfRV (Top 10\\%)","\n$\\bm{n=", n3, ";", n4, "}$"))
  }

  
  
  #Plot of InfRV by coverage and expression for all genes
  fil_piece <- "CoverageInfRVAndExpressionAllGenes"
  fil_name <- paste0(curr_save_dir, fil_piece, PlotNameModifier, ZerosExcludedFromCoverageModifier,".tex")
  tikz(file = fil_name, height=PlotHeight, width=PlotWidth, standAlone = stdal, sanitize = F)
  #paste0(CoverageVarNamePrint, " Coverage by InfRV and \nGene Expression for All Genes (", nrow(gene_level_coverage3), " Total)")
  curr_title <- paste0("Coverage Results for ", CoverageVarNamePrint2, "by InfRV and Gene Expression for All Genes (", nrow(gene_level_coverage3), " Total)")
  if(CoverageVarNamePrint=="TwoPointFiveNinetySevenPointFive"){
    curr_title2 <- CoverageVarNamePrint2
  }else{
    curr_title2 <- curr_title
  }
  
  if(PlotNameModifier=="ForPaper"){
    curr_title2 <- ""
  }
  pl <-ggplot(data=NULL,
              aes(x=HighInfRV,y=CoverageVar,fill=HighExpression))+geom_boxplot()+ggtitle(curr_title2)+
    xlab("InfRV")+ylab("Coverage")+
    scale_x_discrete(labels=my_xlab)+
    theme(legend.text = element_text(size=rel(1), face = "bold"),legend.title = element_text(size=rel(1), face = "bold"),axis.text.y=element_text(size=rel(1.25), face = "bold"),
          axis.text.x=element_text(size=rel(1.4), face = "bold"), plot.title=element_text(size=rel(1.25), hjust=.5,face = "bold"), axis.title=element_text(size=rel(1.25), face = "bold"))
  #ggsave(file = paste0(curr_save_dir, "/CoverageInfRVAndExpressionAllGenes.pdf"), height=PlotHeight, width=PlotWidth)

  if(PlotNameModifier=="ForPaper"){
    print(pl+ xlab("") + theme(legend.position = "none")+scale_fill_manual(values = c("white", "white")))
  }else{
    print(pl)
  }

  dev.off()
  SVBCompiletikzPlot(curr_save_dir, paste0(fil_piece, PlotNameModifier, ZerosExcludedFromCoverageModifier), pdflatex_loc)
  rm(fil_name)
  rm(fil_piece)
  
  
  #Plot the mean tier from Alevin vs coverage
  my_xlab<-paste0(levels(TierValue),"\n$\\bm{n=",table(TierValue), "}$")
  if(sum(is.na(TierValue)) > 0){
    num_NA <- sum(is.na(TierValue))
    my_xlab <- c(my_xlab, paste0("NA\n$\\bm{n=", num_NA, "}$"))
  }
  fil_piece <- "MeanTierAllGenes"
  fil_name <- paste0(curr_save_dir, fil_piece, PlotNameModifier, ZerosExcludedFromCoverageModifier,".tex") 
  tikz(file = fil_name, height=PlotHeight, width=PlotWidth, standAlone = stdal, sanitize = F)
  curr_title <- paste0("Coverage Results for ", CoverageVarNamePrint2, "by Mean Gene Tier for All Genes (", nrow(gene_level_coverage3), " Total)")
  if(CoverageVarNamePrint=="TwoPointFiveNinetySevenPointFive"){
    curr_title2 <- CoverageVarNamePrint2
  }else{
    curr_title2 <- curr_title
  }
  
  if(PlotNameModifier=="ForPaper"){
    curr_title2 <- ""
  }
  pl <- ggplot(data=NULL,
               aes(x=TierValue,y=CoverageVarTier))+geom_boxplot()+ggtitle(curr_title2)+
    xlab("Mean Gene Tier")+ylab("Coverage")+
    scale_x_discrete(labels=my_xlab)+
    theme(axis.text.y=element_text(size=rel(1.25), face = "bold"),
          axis.text.x=element_text(size=rel(1.2), face = "bold"), plot.title=element_text(size=rel(1.25), hjust=.5,face = "bold"), axis.title=element_text(size=rel(1.25), face = "bold"))
  
  print(pl)
  
  #ggsave(file = paste0(curr_save_dir, "/CoverageInfRVAndExpressionFilteredGenes.pdf"), height=PlotHeight, width=PlotWidth)
  dev.off()
  SVBCompiletikzPlot(curr_save_dir, paste0(fil_piece, PlotNameModifier, ZerosExcludedFromCoverageModifier), pdflatex_loc)
  rm(fil_name)
  rm(fil_piece)
  rm(TierValue)
  
  
  #Updated boxplots of Gene Uniqueness level vs coverage using ggplot and splitting all coverage values less than 1 into equally sized groups for all genes
    #And also results for constant cutoffs to make it easier to compare between different experiments
  for(vv in 1:2){
    if(vv==1){
      GeneUniquenessNew <- gene_level_coverage3[,"UniqFactorNew"]
      
      fil_piece <- "CoverageGeneUniquenessggplot"
      PlotWidthUniqToUse <- PlotWidthUniqueness
    }else if(vv==2){
      GeneUniquenessNew <- gene_level_coverage3[,"UniqFactorNewConsCutoffs"]
      
      fil_piece <- "CoverageGeneUniquenessggplotConsCutoffs"
      PlotWidthUniqToUse <- PlotWidth
    }
    
    my_xlab<-paste0(levels(GeneUniquenessNew),"\n$\\bm{n=",table(GeneUniquenessNew), "}$")
    if(sum(is.na(GeneUniquenessNew)) > 0){
      num_NA <- sum(is.na(GeneUniquenessNew))
      my_xlab <- c(my_xlab, paste0("NA\n$\\bm{n=", num_NA, "}$"))
    }
    
    fil_name <- paste0(curr_save_dir, fil_piece, PlotNameModifier, ZerosExcludedFromCoverageModifier,".tex")  
    tikz(file = fil_name, height=PlotHeight, width=PlotWidthUniqToUse, standAlone = stdal, sanitize = F)
    
    #paste0(CoverageVarNamePrint, " Coverage by InfRV and \nGene Expression for Filtered Genes (", nrow(gene_level_coverage3_filtered_genes), " Total)")
    curr_title <- paste0("Coverage Results for ", CoverageVarNamePrint2, "by Gene Uniqueness for All Genes (", nrow(gene_level_coverage3), " Total)")
    if(CoverageVarNamePrint=="TwoPointFiveNinetySevenPointFive"){
      curr_title2 <- CoverageVarNamePrint2
    }else{
      curr_title2 <- curr_title
    }
    
    if(PlotNameModifier=="ForPaper"){
      curr_title2 <- ""
    }
    pl <- ggplot(data=NULL,
                 aes(x=GeneUniquenessNew,y=CoverageVar))+geom_boxplot()+ggtitle(curr_title2)+
      xlab("Gene Uniqueness Ratio")+ylab("Coverage")+
      scale_x_discrete(labels=my_xlab)+
      theme(axis.text.y=element_text(size=rel(1.25), face = "bold"),
            axis.text.x=element_text(size=rel(1.2), face = "bold"), plot.title=element_text(size=rel(1.25), hjust=.5,face = "bold"), axis.title=element_text(size=rel(1.25), face = "bold"))
    
    print(pl)
    
    #ggsave(file = paste0(curr_save_dir, "/CoverageInfRVAndExpressionFilteredGenes.pdf"), height=PlotHeight, width=PlotWidth)
    dev.off()
    SVBCompiletikzPlot(curr_save_dir, paste0(fil_piece, PlotNameModifier, ZerosExcludedFromCoverageModifier), pdflatex_loc)
    rm(fil_name)
    rm(fil_piece)
    rm(GeneUniquenessNew)
  }
 
  
  
  
  #Extract corresponding elements for filtered genes only
  HighInfRV <- gene_level_coverage3_filtered_genes[,"HighInfRV"]
  CoverageVar <- gene_level_coverage3_filtered_genes[,CoverageVarName]
  HighExpression <- gene_level_coverage3_filtered_genes[,"HighExpression"]
  TierValue <- gene_level_coverage3_filtered_genes$TierFactor[!is.na(gene_level_coverage3_filtered_genes$TierFactor)]
  CoverageVarTier <- gene_level_coverage3_filtered_genes[,CoverageVarName][!is.na(gene_level_coverage3_filtered_genes$TierFactor)]
  
  #Generate sample sizes corresponding to each of the 4 separate boxplots
  n1 <- sum(HighInfRV=="LowInfRV" & HighExpression=="LowExpr")
  n2 <- sum(HighInfRV=="LowInfRV" & HighExpression=="HighExpr")
  n3 <- sum(HighInfRV=="HighInfRV" & HighExpression=="LowExpr")
  n4 <- sum(HighInfRV=="HighInfRV" & HighExpression=="HighExpr")
  
  
  #Custom labels using the sample sizes calculated above
  
  
  if(PlotNameModifier=="ForPaper"){
    my_xlab <- c(paste0("Low InfRV","\nLow Expr ; High Expr", "\n$\\bm{n=", n1, " ; ", n2, "}$"), paste0("High InfRV","\nLow Expr ; High Expr","\n$\\bm{n=", n3, " ; ", n4, "}$"))
  }else{
    my_xlab <- c(paste0("Low InfRV","\n$\\bm{n=", n1, ";", n2, "}$"), paste0("High InfRV (Top 10\\%)","\n$\\bm{n=", n3, ";", n4, "}$"))
  }
  #my_xlab<-paste0(levels(HighInfRV),"\n$\\bm{n=",table(HighInfRV), "}$")
  
  
  
  
  #Boxplots of Coverage by InfRV and Expression being low or high (top 10%) for filtered genes
  fil_piece <- "CoverageInfRVAndExpressionFilteredGenes"
  fil_name <- paste0(curr_save_dir, fil_piece, PlotNameModifier, ZerosExcludedFromCoverageModifier,".tex")  
  tikz(file = fil_name, height=PlotHeight, width=PlotWidth, standAlone = stdal, sanitize = F)
  
  #paste0(CoverageVarNamePrint, " Coverage by InfRV and \nGene Expression for Filtered Genes (", nrow(gene_level_coverage3_filtered_genes), " Total)")
  curr_title <- paste0("Coverage Results for ", CoverageVarNamePrint2, "by InfRV and Gene Expression for Filtered Genes (", nrow(gene_level_coverage3_filtered_genes), " Total)")
  if(CoverageVarNamePrint=="TwoPointFiveNinetySevenPointFive"){
    curr_title2 <- CoverageVarNamePrint2
  }else{
    curr_title2 <- curr_title
  }
  
  if(PlotNameModifier=="ForPaper"){
    curr_title2 <- ""
  }
  pl <- ggplot(data=NULL,
               aes(x=HighInfRV,y=CoverageVar,fill=HighExpression))+geom_boxplot()+ggtitle(curr_title2)+
    xlab("")+ylab("Coverage")+
    scale_x_discrete(labels=my_xlab)+
    scale_fill_discrete(name="Expression Level",
                        labels=c("Low", "High (Top 10\\%)"))+
    theme(legend.text = element_text(size=rel(1), face = "bold"),legend.title = element_text(size=rel(1), face = "bold"),axis.text.y=element_text(size=rel(1.25), face = "bold"),
          axis.text.x=element_text(size=rel(1.4), face = "bold"), plot.title=element_text(size=rel(1.25), hjust=.5,face = "bold"), axis.title=element_text(size=rel(1.25), face = "bold"))

  if(PlotNameModifier=="ForPaper"){
    print(pl+ xlab("") + theme(legend.position = "none")+scale_fill_manual(values = c("white", "white")))
  }else{
    print(pl)
  }
  
  #ggsave(file = paste0(curr_save_dir, "/CoverageInfRVAndExpressionFilteredGenes.pdf"), height=PlotHeight, width=PlotWidth)
  dev.off()
  SVBCompiletikzPlot(curr_save_dir, paste0(fil_piece, PlotNameModifier, ZerosExcludedFromCoverageModifier), pdflatex_loc)
  rm(fil_name)
  rm(fil_piece)
  
  
  
  
  #Now, boxplots of the widths of the intervals stratified by InfRV and total expression being low or high (top 10%)
  CurrWidth <- gene_level_coverage3_filtered_genes[, paste0("Mean", CoverageVarNamePrint, "Width")]

  
  
  fil_piece <- "WidthInfRVAndExpressionFilteredGenes"
  fil_name <- paste0(curr_save_dir, fil_piece, PlotNameModifier, ZerosExcludedFromCoverageModifier, ".tex")  
  tikz(file = fil_name, height=PlotHeight, width=PlotWidth, standAlone = stdal, sanitize = F)
  
  #paste0(CoverageVarNamePrint, " Coverage by InfRV and \nGene Expression for Filtered Genes (", nrow(gene_level_coverage3_filtered_genes), " Total)")
  curr_title <- paste0("Interval Widths of ", CoverageVarNamePrint2, "by InfRV and Gene Expression for Filtered Genes (", nrow(gene_level_coverage3_filtered_genes), " Total)")
  if(CoverageVarNamePrint=="TwoPointFiveNinetySevenPointFive"){
    curr_title2 <- CoverageVarNamePrint2
  }else{
    curr_title2 <- curr_title
  }
  
  if(PlotNameModifier=="ForPaper"){
    curr_title2 <- ""
  }
  
  if(SimNumberToUse%in%c(1001,11001)){
    max_y_lim <- 50.1
  }else{
    max_y_lim <- 100.1
  }
  pl <- ggplot(data=NULL,
               aes(x=HighInfRV,y=CurrWidth,fill=HighExpression))+geom_boxplot()+ggtitle(curr_title2)+
    scale_y_continuous(limits = c(0, max_y_lim))+xlab("")+ylab("Interval Width")+
    scale_x_discrete(labels=my_xlab)+
    scale_fill_discrete(name="Expression Level",
                        labels=c("Low", "High (Top 10\\%)"))+
    theme(legend.text = element_text(size=rel(1), face = "bold"),legend.title = element_text(size=rel(1), face = "bold"),axis.text.y=element_text(size=rel(1.25), face = "bold"),
          axis.text.x=element_text(size=rel(1.4), face = "bold"), plot.title=element_text(size=rel(1.25), hjust=.5,face = "bold"), axis.title=element_text(size=rel(1.25), face = "bold"))

  if(PlotNameModifier=="ForPaper"){
    print(pl+ xlab("") + theme(legend.position = "none")+scale_fill_manual(values = c("white", "white")))
  }else{
    print(pl)
  }
  
  #ggsave(file = paste0(curr_save_dir, "/CoverageInfRVAndExpressionFilteredGenes.pdf"), height=PlotHeight, width=PlotWidth)
  dev.off()
  SVBCompiletikzPlot(curr_save_dir, paste0(fil_piece, PlotNameModifier, ZerosExcludedFromCoverageModifier), pdflatex_loc)
  rm(fil_name)
  rm(fil_piece)
  
  
  

  
  
  
  
  
  #Updated boxplot of Gene Uniqueness Level vs coverage using ggplot and splitting all coverage values less than 1 into equally sized groups for filtered genes (vv=1)
    # as well as for the constant cutoffs of 0.75, 0.95, and 1
  
  for(vv in 1:2){
    if(vv==1){
      GeneUniquenessNew <- gene_level_coverage3_filtered_genes[,"UniqFactorNew"]
      
      fil_piece <- "CoverageGeneUniquenessggplotFilteredGenes"
      fil_piece2 <- "GeneWidthByuniquenessggplotFilteredGenes"
      PlotWidthUniqToUse <- PlotWidthUniqueness
    }else if(vv==2){
      GeneUniquenessNew <- gene_level_coverage3_filtered_genes[,"UniqFactorNewConsCutoffs"]
      
      fil_piece <- "CoverageGeneUniquenessggplotFilteredGenesConsCutoffs"
      fil_piece2 <- "GeneWidthByuniquenessggplotFilteredGenesConsCutoffs"
      PlotWidthUniqToUse <- PlotWidth
    }

    fil_name <- paste0(curr_save_dir, fil_piece, PlotNameModifier, ZerosExcludedFromCoverageModifier, ".tex")  
    tikz(file = fil_name, height=PlotHeight, width=PlotWidthUniqToUse, standAlone = stdal, sanitize = F)
    my_xlab<-paste0(levels(GeneUniquenessNew),"\n$\\bm{n=",table(GeneUniquenessNew), "}$")
    if(sum(is.na(GeneUniquenessNew)) > 0){
      num_NA <- sum(is.na(GeneUniquenessNew))
      my_xlab <- c(my_xlab, paste0("NA\n$\\bm{n=", num_NA, "}$"))
    }
    print(paste0("The my_xlab is given below"))
    print(my_xlab)
    #paste0(CoverageVarNamePrint, " Coverage by InfRV and \nGene Expression for Filtered Genes (", nrow(gene_level_coverage3_filtered_genes), " Total)")
    curr_title <- paste0("Coverage Results for ", CoverageVarNamePrint2, "by Gene Uniqueness for Filtered Genes (", nrow(gene_level_coverage3_filtered_genes), " Total)")
    if(CoverageVarNamePrint=="TwoPointFiveNinetySevenPointFive"){
      curr_title2 <- CoverageVarNamePrint2
    }else{
      curr_title2 <- curr_title
    }
    
    if(PlotNameModifier=="ForPaper"){
      curr_title2 <- ""
    }
    
    pl <- ggplot(data=NULL,
                 aes(x=GeneUniquenessNew,y=CoverageVar))+geom_boxplot()+ggtitle(curr_title2)+
      xlab("Gene Uniqueness Ratio")+ylab("Coverage")+
      scale_x_discrete(labels=my_xlab)+
      theme(axis.text.y=element_text(size=rel(1.4), face = "bold"),
            axis.text.x=element_text(size=rel(1.4), face = "bold"), plot.title=element_text(size=rel(1.25), hjust=.5,face = "bold"), axis.title=element_text(size=rel(1.4), face = "bold"))
    print(pl)
    
    #ggsave(file = paste0(curr_save_dir, "/CoverageInfRVAndExpressionFilteredGenes.pdf"), height=PlotHeight, width=PlotWidth)
    dev.off()
    SVBCompiletikzPlot(curr_save_dir, paste0(fil_piece, PlotNameModifier, ZerosExcludedFromCoverageModifier), pdflatex_loc)
    rm(fil_name)
    rm(fil_piece)
    
    
    
    
    #Now, interval widths
    
    #Now, boxplots the interval width of the current intervals (for filtered genes) by gene uniqueness level
    IntervalWidth <- gene_level_coverage3_filtered_genes[,paste0("Mean", CoverageVarNamePrint, "Width")]
    
    
    
    fil_name <- paste0(curr_save_dir, fil_piece2, PlotNameModifier, ZerosExcludedFromCoverageModifier, ".tex")  
    tikz(file = fil_name, height=PlotHeight, width=PlotWidthUniqToUse, standAlone = stdal, sanitize = F)
    
    #paste0(CoverageVarNamePrint, " Coverage by InfRV and \nGene Expression for Filtered Genes (", nrow(gene_level_coverage3_filtered_genes), " Total)")
    curr_title <- paste0("Interval Widths of ", CoverageVarNamePrint2, "by Gene Uniqueness for Filtered Genes (", nrow(gene_level_coverage3_filtered_genes), " Total)")
    if(CoverageVarNamePrint=="TwoPointFiveNinetySevenPointFive"){
      curr_title2 <- CoverageVarNamePrint2
    }else{
      curr_title2 <- curr_title
    }
    
    if(PlotNameModifier=="ForPaper"){
      curr_title2 <- ""
    }
    
    if(SimNumberToUse%in%c(1001,11001)){
      max_y_lim <- 50.1
    }else{
      max_y_lim <- 100.1
    }
    
    pl <- ggplot(data=NULL,
                 aes(x=GeneUniquenessNew,y=IntervalWidth))+geom_boxplot()+ggtitle(curr_title2)+
      xlab("Gene Uniqueness Ratio")+ylab("Interval Width")+
      scale_y_continuous(limits = c(0, max_y_lim))+scale_x_discrete(labels=my_xlab)+
      theme(axis.text.y=element_text(size=rel(1.4), face = "bold"),
            axis.text.x=element_text(size=rel(1.4), face = "bold"), plot.title=element_text(size=rel(1.25), hjust=.5,face = "bold"), axis.title=element_text(size=rel(1.4), face = "bold"))
    print(pl)
    
    #ggsave(file = paste0(curr_save_dir, "/CoverageInfRVAndExpressionFilteredGenes.pdf"), height=PlotHeight, width=PlotWidth)
    dev.off()
    SVBCompiletikzPlot(curr_save_dir, paste0(fil_piece2, PlotNameModifier, ZerosExcludedFromCoverageModifier), pdflatex_loc)
    rm(fil_name)
    rm(fil_piece2)
  }
  
  
  
  
  #Plot the mean tier from Alevin vs coverage for filtered genes only
  
  my_xlab<-paste0(levels(TierValue),"\n$\\bm{n=",table(TierValue), "}$")
  if(sum(is.na(TierValue)) > 0){
    num_NA <- sum(is.na(TierValue))
    my_xlab <- c(my_xlab, paste0("NA\n$\\bm{n=", num_NA, "}$"))
  }
  fil_piece <- "MeanTierFilteredGenes"
  fil_name <- paste0(curr_save_dir, fil_piece, PlotNameModifier, ZerosExcludedFromCoverageModifier,".tex") 
  tikz(file = fil_name, height=PlotHeight, width=PlotWidth, standAlone = stdal, sanitize = F)
  curr_title <- paste0("Coverage Results for ", CoverageVarNamePrint2, "by Mean Gene Tier for Filtered Genes (", nrow(gene_level_coverage3_filtered_genes), " Total)")
  if(CoverageVarNamePrint=="TwoPointFiveNinetySevenPointFive"){
    curr_title2 <- CoverageVarNamePrint2
  }else{
    curr_title2 <- curr_title
  }
  
  if(PlotNameModifier=="ForPaper"){
    curr_title2 <- ""
  }
  
  pl <- ggplot(data=NULL,
               aes(x=TierValue,y=CoverageVarTier))+geom_boxplot()+ggtitle(curr_title2)+
    xlab("Mean Gene Tier")+ylab("Coverage")+
    scale_x_discrete(labels=my_xlab)+
    theme(axis.text.y=element_text(size=rel(1.5), face = "bold"),
          axis.text.x=element_text(size=rel(1.5), face = "bold"), plot.title=element_text(size=rel(1.25), hjust=.5,face = "bold"), axis.title=element_text(size=rel(1.5), face = "bold"))
  
  print(pl)
  
  #ggsave(file = paste0(curr_save_dir, "/CoverageInfRVAndExpressionFilteredGenes.pdf"), height=PlotHeight, width=PlotWidth)
  dev.off()
  SVBCompiletikzPlot(curr_save_dir, paste0(fil_piece, PlotNameModifier, ZerosExcludedFromCoverageModifier), pdflatex_loc)
  rm(fil_name)
  rm(fil_piece)

  
  
  
  #Now, width plot for the tier of the genes
  CurrWidthTier <- gene_level_coverage3_filtered_genes[, paste0("Mean", CoverageVarNamePrint, "Width")][!is.na(gene_level_coverage3_filtered_genes$TierFactor)]
 
  fil_piece <- "WidthMeanTierFilteredGenes"
  fil_name <- paste0(curr_save_dir, fil_piece, PlotNameModifier, ZerosExcludedFromCoverageModifier,".tex") 
  tikz(file = fil_name, height=PlotHeight, width=PlotWidth, standAlone = stdal, sanitize = F)
  curr_title <- paste0("Interval Widths of ", CoverageVarNamePrint2, "by Mean Gene Tier for Filtered Genes (", nrow(gene_level_coverage3_filtered_genes), " Total)")
  if(CoverageVarNamePrint=="TwoPointFiveNinetySevenPointFive"){
    curr_title2 <- CoverageVarNamePrint2
  }else{
    curr_title2 <- curr_title
  }
  
  if(PlotNameModifier=="ForPaper"){
    curr_title2 <- ""
  }
  
  if(SimNumberToUse%in%c(1001,11001)){
    max_y_lim <- 25.1
  }else{
    max_y_lim <- 100.1
  }
  
  pl <- ggplot(data=NULL,
               aes(x=TierValue,y=CurrWidthTier))+geom_boxplot()+ggtitle(curr_title2)+
    scale_y_continuous(limits = c(0, max_y_lim))+xlab("Mean Gene Tier")+ylab("Interval Width")+
    scale_x_discrete(labels=my_xlab)+
    theme(axis.text.y=element_text(size=rel(1.5), face = "bold"),
          axis.text.x=element_text(size=rel(1.5), face = "bold"), plot.title=element_text(size=rel(1.25), hjust=.5,face = "bold"), axis.title=element_text(size=rel(1.5), face = "bold"))
  
  print(pl)
  
  #ggsave(file = paste0(curr_save_dir, "/CoverageInfRVAndExpressionFilteredGenes.pdf"), height=PlotHeight, width=PlotWidth)
  dev.off()
  SVBCompiletikzPlot(curr_save_dir, paste0(fil_piece, PlotNameModifier, ZerosExcludedFromCoverageModifier), pdflatex_loc)
  rm(fil_name)
  rm(fil_piece)
  
  
  rm(TierValue)
  #Now, plot the difference in coverages between norm95 interval and the current selected iterval
  if(CoverageVarName!="Norm95ConfInvCoverage"){
    CurrCoverage <- gene_level_coverage3_filtered_genes[,CoverageVarName]
    NormCoverage <- gene_level_coverage3_filtered_genes[,"Norm95ConfInvCoverage"]
    
    DiffCoverage <- NormCoverage - CurrCoverage 
    
    
    fil_piece <- "DiffInCoverageByuniquenessggplotFilteredGenes"
    fil_name <- paste0(curr_save_dir, fil_piece, PlotNameModifier, ZerosExcludedFromCoverageModifier, ".tex")  
    tikz(file = fil_name, height=PlotHeight, width=PlotWidth, standAlone = stdal, sanitize = F)
    
    #paste0(CoverageVarNamePrint, " Coverage by InfRV and \nGene Expression for Filtered Genes (", nrow(gene_level_coverage3_filtered_genes), " Total)")
    curr_title <- paste0("Difference in Coverage From the Normal Quantile Interval \nfor ", CoverageVarNamePrint2, "\nby Gene Uniqueness for Filtered Genes (", nrow(gene_level_coverage3_filtered_genes), " Total)")
    pl <- ggplot(data=NULL,
                 aes(x=GeneUniquenessNew,y=DiffCoverage))+geom_boxplot()+ggtitle(curr_title)+
      xlab("Gene Uniqueness Ratio")+ylab(paste0("Coverage of 95\\% Interval from Normal Quantiles - \nCoverage of ", CoverageVarNamePrint2))+
      scale_x_discrete(labels=my_xlab)+
      theme(axis.text.y=element_text(size=rel(1.25), face = "bold"),
            axis.text.x=element_text(size=rel(1.4), face = "bold"), plot.title=element_text(size=rel(1.25), hjust=.5,face = "bold"), axis.title=element_text(size=rel(1.25), face = "bold"))
    print(pl)
    
    #ggsave(file = paste0(curr_save_dir, "/CoverageInfRVAndExpressionFilteredGenes.pdf"), height=PlotHeight, width=PlotWidth)
    dev.off()
    SVBCompiletikzPlot(curr_save_dir, paste0(fil_piece, PlotNameModifier, ZerosExcludedFromCoverageModifier), pdflatex_loc)
    rm(fil_name)
    rm(fil_piece)
  }
  
  #Now, plot the difference in widths between norm95 interval and the current selected iterval
  if(CoverageVarName!="Norm95ConfInvCoverage"){
    CurrWidth <- gene_level_coverage3_filtered_genes[, paste0("Mean", CoverageVarNamePrint, "Width")]
    NormWidth <- gene_level_coverage3_filtered_genes[,paste0("Mean", "Norm95ConfInv", "Width")]
    
    DiffWidth <- NormWidth - CurrWidth
    
    
    fil_piece <- "DiffInWidthByuniquenessggplotFilteredGenes"
    fil_name <- paste0(curr_save_dir, fil_piece, PlotNameModifier, ZerosExcludedFromCoverageModifier, ".tex")  
    tikz(file = fil_name, height=PlotHeight, width=PlotWidth, standAlone = stdal, sanitize = F)
    
    #paste0(CoverageVarNamePrint, " Coverage by InfRV and \nGene Expression for Filtered Genes (", nrow(gene_level_coverage3_filtered_genes), " Total)")
    pl <- ggplot(data=NULL,
                 aes(x=GeneUniquenessNew,y=DiffWidth))+geom_boxplot()+ggtitle(paste0("Difference in Width From the Normal Quantile Interval \nfor ", CoverageVarNamePrint2, "\nby Gene Uniqueness for Filtered Genes (", nrow(gene_level_coverage3_filtered_genes), " Total)"))+
      xlab("Gene Uniqueness Ratio")+ylab(paste0("Width of 95\\% Interval from Normal Quantiles - \nWidth of ", CoverageVarNamePrint2))+
      scale_x_discrete(labels=my_xlab)+
      theme(axis.text.y=element_text(size=rel(1.25), face = "bold"),
            axis.text.x=element_text(size=rel(1), face = "bold"), plot.title=element_text(size=rel(1.25), hjust=.5,face = "bold"), axis.title=element_text(size=rel(1.25), face = "bold"))
    print(pl)
    
    #ggsave(file = paste0(curr_save_dir, "/CoverageInfRVAndExpressionFilteredGenes.pdf"), height=PlotHeight, width=PlotWidth)
    dev.off()
    SVBCompiletikzPlot(curr_save_dir, paste0(fil_piece, PlotNameModifier, ZerosExcludedFromCoverageModifier), pdflatex_loc)
    rm(fil_name)
    rm(fil_piece)
  }
  
  
  
  
  
}
