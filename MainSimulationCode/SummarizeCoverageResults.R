#Summarize Coverage Results
onCluster <- TRUE

if(onCluster==TRUE){
  #SimNumberToUse <- 117
  #Specify libPaths to ensure the correct one is used
  .libPaths("/nas/longleaf/home/skvanbur/lib64/R/library")
  SimNumberToUse <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
  Sys.sleep(2*SimNumberToUse)
  base_dir <- "/pine/scr/s/k/skvanbur/SingleCellProject/SimulationResults/"
  PlotsSaveDir <- paste0("/pine/scr/s/k/skvanbur/SingleCellProject/CoveragePlots/Sim", SimNumberToUse, "/")
  source("~/SingleCellProject/SingleCellProjectFunctions.R")
  
}

library(tikzDevice)
library(ggplot2)

if(!dir.exists(PlotsSaveDir)){dir.create(PlotsSaveDir, recursive = TRUE)}

DropMinnowSkippedGenesFromCoverageCalc <- TRUE


if(SimNumberToUse%in% c(10117,10118,11001)){
  nboot <- 20
}else{
  nboot <- 100
}

#Generate differently formatted plots for different purposes (general evaluation, for presentation slides, and for the paper)
for(lo in 1:3){
  if(lo==1){
    GeneratePlotsForSlides <- TRUE
    GeneratePlotsForPaper <- FALSE
  }else if(lo==2){
    GeneratePlotsForSlides <- FALSE
    GeneratePlotsForPaper <- FALSE
  }else if(lo==3){
    GeneratePlotsForSlides <- FALSE
    GeneratePlotsForPaper <- TRUE
  }

  #Generate Main coverage plots (both including zeros and not including them) for the analysis
  for(zo in 1:2){
    if(zo==1){
      EvaluateCoverageWithoutZeros <- TRUE
    }else if(zo==2){
      EvaluateCoverageWithoutZeros <- FALSE
    }


      if(GeneratePlotsForSlides==TRUE){
        PlotHeight <- 5.25/2
        PlotNameModifier <- "ForSlides"
      }else if(GeneratePlotsForPaper==TRUE){
        PlotHeight <- 5.25/2
        PlotNameModifier <- "ForPaper"
      }else{
        PlotHeight <- 5.25
        PlotNameModifier <- ""
      }
      
      #Load the mean of the bootstrap samples
      load(paste0(base_dir, "Sim", SimNumberToUse, "/BootSampsMeans.RData"))
      
      if(EvaluateCoverageWithoutZeros==TRUE){
        load(paste0(base_dir, "Sim", SimNumberToUse, "/CoverageResCalculatedWithoutZeros.RData"))
        ZerosExcludedFromCoverageModifier <- "NoZeros"
        DropMinnowSkippedGenesFromCoverageCalc <- TRUE
        print(paste0("EvaluateCoverageWithoutZeros is TRUE, meaning DropMinnowSkippedGenesFromCoverageCalc is being set to TRUE automatically"))
      }else{
        if(DropMinnowSkippedGenesFromCoverageCalc==TRUE){
          load(paste0(base_dir, "Sim", SimNumberToUse, "/CoverageResDropMinnowSkippedGenes.RData"))
        }else{
          load(paste0(base_dir, "Sim", SimNumberToUse, "/CoverageRes.RData"))
        }
        ZerosExcludedFromCoverageModifier <- ""
      }
      
      
      
      gene_level_coverage <- CoverageRes$gene_level_coverage
      coverage_res_by_uniqueness_bin <- CoverageRes$coverage_res_by_uniqueness_bin
      
      cell_level_coverage <- CoverageRes$cell_level_coverage
      cell_level_averaged_Coverage <- CoverageRes$cell_level_averaged_coverage
      
      Alevin_counts<- loadRData(paste0(base_dir, "Sim", SimNumberToUse, "/AlevinCounts.RData"))
      
      #Use splatter counts that have been scaled to alevin's library size for coverage comparisons because library size itself is somewhat arbitrary and this makes comparisons more fair
      #Splatter_counts <- loadRData(paste0("/Users/Scott/Documents/Dissertation Data/SingleCellProject/Simulation Results/Sim", SimNumberToUse, "/SplatterCounts.RData"))
      Splatter_counts <- loadRData(paste0(base_dir, "Sim", SimNumberToUse, "/SplatterCountsScaledToAlevinLibSize.RData"))
      
      Alevin_counts2 <- data.frame(Alevin_counts)
      Splatter_counts2 <- data.frame(Splatter_counts)
      
      mean_Alevin_counts <- data.frame(rowMeans(Alevin_counts2))
      colnames(mean_Alevin_counts) <- "mean_Alevin_counts"
      mean_Splatter_counts <- data.frame(rowMeans(Splatter_counts2))
      colnames(mean_Splatter_counts) <- "mean_Splatter_counts"
      
      mean_Alevin_counts$gene_id <- rownames(mean_Alevin_counts)
      mean_Splatter_counts$gene_id <- rownames(mean_Splatter_counts)
      
      
      gene_level_coverage2 <- merge(gene_level_coverage, mean_Alevin_counts, by = "gene_id")
      gene_level_coverage3 <- merge(gene_level_coverage2, mean_Splatter_counts, by = "gene_id")
      
      #Now, add tier information
      TierInfo <- loadRData(paste0(base_dir, "Sim", SimNumberToUse, "/AlevinTierInfo.RData"))
      TierInfoGeneSummary <- data.frame(apply(TierInfo,1, function(x){mean(x[x>0], na.rm=T)}))
      colnames(TierInfoGeneSummary) <- "MeanTier"
      TierInfoGeneSummary$gene_id <- rownames(TierInfoGeneSummary)
      
      gene_level_coverage3 <- merge(gene_level_coverage3, TierInfoGeneSummary, by = "gene_id")
      
      quantilevals_tier_all<- quantile(gene_level_coverage3$MeanTier[gene_level_coverage3$MeanTier!=1], probs = seq(0,1,.01), na.rm=T)
      #labels=c("10\\%", "20\\%", "30\\%", "40\\%", "50\\%","60\\%", "70\\%", "80\\%", "90\\%", "100\\%")
      
      dff <- data.frame(x = gene_level_coverage3$MeanTier,y=gene_level_coverage3$TwoPointFiveNinetySevenPointFiveCoverage)
      ggplot(dff,aes(x=x,y=y)) + geom_point(alpha = 0.3)+xlab("Mean Tier")+ylab("Coverage")
      gene_level_coverage3[,"TierFactor"] <- cut(gene_level_coverage3$MeanTier,
                                                    breaks=c(1,1.5,2,2.5,3),
                                                    include.lowest=TRUE)
      
      filt_func <- function(x){
        ncells_high_exp <- sum(x >= 10)
        return(ncells_high_exp)
      }
      
      #vals1 <- apply(Alevin_counts2, 1, filt_func)
      vals1 <- apply(Splatter_counts2, 1, filt_func)
      
      #genes_filtered <- rownames(Alevin_counts2)[vals1 > 10]
      if(DropMinnowSkippedGenesFromCoverageCalc==TRUE){
        genes_filteredT <- rownames(Splatter_counts2)[vals1 > 10]
        genes_filtered <- genes_filteredT[genes_filteredT%in%gene_level_coverage$gene_id]
      }else{
        genes_filtered <- rownames(Splatter_counts2)[vals1 > 10]
      }
      
      
      
      Alevin_counts_Filtered <- Alevin_counts2[genes_filtered,]
      Splatter_counts_Filtered <- Splatter_counts2[genes_filtered,]
      
      if(onCluster==TRUE){
        pdflatex_loc <- '/nas/longleaf/home/skvanbur/texlive/2020/bin/x86_64-linux/pdflatex'
        options(tikzLatex = pdflatex_loc)
      }else{
        pdflatex_loc <- '/Library/TeX/texbin/pdflatex'
      }
      
      
      tikzLockFile <- paste0(PlotsSaveDir, "/TikzTemp___LOCK")
      if(file.exists(tikzLockFile)){system(paste0("rm ", tikzLockFile))}
      options("tikzDocumentDeclaration" = "\\documentclass[12pt,letterpaper,twoside]{report}\n")
      options("tikzMetricsDictionary" = paste(PlotsSaveDir, "/TikzTemp", sep = ""))
      options("tikzLatexPackages" = c(
        "\\usepackage{tikz}",
        "\\usepackage[active,tightpage,psfixbb]{preview}",
        "\\PreviewEnvironment{pgfpicture}",
        "\\setlength\\PreviewBorder{0pt}",
        "\\usepackage{bm}",
        "\\usepackage{newtxtext}")
      )
      
      stdal <- TRUE
      #Calculate "bias" of alevin counts- that is, how do the results differ from the (scaled) splatter counts
        #Create histograms for both all genes and filtered genes
      Diff_Alevin_Splatter_Counts <- Alevin_counts2 - Splatter_counts2
      BiasAveragedAcrossCells <- rowMeans(Diff_Alevin_Splatter_Counts)
      BiasAveragedAcrossCells2 <- BiasAveragedAcrossCells[abs(BiasAveragedAcrossCells) < 2]
      #pdf(file = paste0(PlotsSaveDir, "HistogramAllGenes.pdf"), height = 5, width = 6.5)
      
      fil_piece <- "HistogramAllGenes"
      fil_name <- paste0(PlotsSaveDir, fil_piece, ".tex")
      
      tikz(file = fil_name, height=5, width=6.5, standAlone = stdal, sanitize = F)
      par(mar=c(5.1, 4.1, 4.1, 2.1))
      hist(BiasAveragedAcrossCells2, xlab = "Alevin Count - (Scaled) Splatter Count\n(Averaged Across Cells)",
           main = "Difference Between Alevin Count and Splatter Count \n(Scaled to Alevin Library Size) for all 60,179 Genes\n(Averaged Across Cells)",
           breaks = 40)
      legend("topleft", c(paste0("Mean: ", fr(mean(BiasAveragedAcrossCells),3, NLS=T)), paste0("Median: ", fr(median(BiasAveragedAcrossCells),3, NLS=T)),
                          paste0("SD: ", fr(sd(BiasAveragedAcrossCells),3, NLS=T)), paste0("Max: ", fr(max(BiasAveragedAcrossCells),3, NLS=T)),
                          paste0("Min: ", fr(min(BiasAveragedAcrossCells),3, NLS=T))))
      dev.off()
      
      SVBCompiletikzPlot(PlotsSaveDir, fil_piece, pdflatex_loc)
      
      
      n_genes_filtered <- length(genes_filtered)
      Diff_Alevin_Splatter_Counts_Filtered <- boot_samps_means[genes_filtered,] - Splatter_counts_Filtered
      BiasAveragedAcrossCellsFiltered <- rowMeans(Diff_Alevin_Splatter_Counts_Filtered)
      BiasAveragedAcrossCellsFiltered2 <- BiasAveragedAcrossCellsFiltered[abs(BiasAveragedAcrossCellsFiltered) < 5]
      #pdf(file = paste0(PlotsSaveDir, "HistogramFilteredGenes.pdf"), height = 5, width = 6.5)
      
      
      fil_piece <- "HistogramFilteredGenes"
      fil_name <- paste0(PlotsSaveDir, fil_piece, ".tex")
      tikz(file = fil_name, height=5, width=6.5, standAlone = stdal, sanitize = F)
      par(mar=c(5.1, 4.1, 5.1, 2.1))
      hist(BiasAveragedAcrossCellsFiltered2, xlab = "Alevin Count - (Scaled) Splatter Count\n(Averaged Across Cells)",
           main = paste0("Difference Between Alevin Count and Splatter Count \n(Scaled to Alevin Library Size) for the ", n_genes_filtered, " Genes That Pass Filtering\n(Count > 10 in at least 10 Cells)\n(Averaged Across Cells)"),
           breaks = 40)
      legend("topleft", c(paste0("Mean: ", fr(mean(BiasAveragedAcrossCellsFiltered),3, NLS=T)), paste0("Median: ", fr(median(BiasAveragedAcrossCellsFiltered),3, NLS=T)),
                          paste0("SD: ", fr(sd(BiasAveragedAcrossCellsFiltered),3, NLS=T)), paste0("Max: ", fr(max(BiasAveragedAcrossCellsFiltered),3, NLS=T)),
                          paste0("Min: ", fr(min(BiasAveragedAcrossCellsFiltered),3, NLS=T))))
      
      dev.off()
      SVBCompiletikzPlot(PlotsSaveDir, fil_piece, pdflatex_loc)
      
      
      
      gene_level_coverage3_filtered_genes <- subset(gene_level_coverage3, gene_level_coverage3$gene_id %in% genes_filtered)
      
      #Create indicator variables for high or low uncertainty, with uncertainty being defined as high if the MeanInfRV is n the top 5% of all genes
        #and low otherwise
      #Could think about doing a cutoff for the MeanInfRV value, but not sure what hte interpretations for that would be
      #and the threshold could change a lot depending on what dataset was being used
      HighCutoffAllGenes <- quantile(gene_level_coverage3$MeanInfRV, probs= 0.90)
      HighCutoffFilteredGenes <- quantile(gene_level_coverage3_filtered_genes$MeanInfRV, probs= 0.90)
      
      gene_level_coverage3$HighInfRV <- factor(ifelse(gene_level_coverage3$MeanInfRV >=HighCutoffAllGenes, 1, 0))
      levels(gene_level_coverage3$HighInfRV) <- c("LowInfRV", "HighInfRV")
      gene_level_coverage3_filtered_genes$HighInfRV <- factor(ifelse(gene_level_coverage3_filtered_genes$MeanInfRV >=HighCutoffFilteredGenes, 1, 0))
      levels(gene_level_coverage3_filtered_genes$HighInfRV) <- c("LowInfRV", "HighInfRV")
      
      
      #Now, split into high and low expression levels, with high expression levels being in the top 5% of values and low otherwise
        #Again, could think about doing a specific cutoff value but the magnitude could change a good bit depending on the dataset
      HighExpressionCutoffAllGenes <-quantile(gene_level_coverage3$MeanofMeanBootSamps, probs = 0.90)
      HighExpressionCutoffFilteredGenes <-quantile(gene_level_coverage3_filtered_genes$MeanofMeanBootSamps, probs = 0.90)
      
      gene_level_coverage3$HighExpression <- factor(ifelse(gene_level_coverage3$MeanofMeanBootSamps >=HighExpressionCutoffAllGenes, 1, 0))
      levels(gene_level_coverage3$HighExpression) <- c("LowExpr", "HighExpr")
      
      gene_level_coverage3_filtered_genes$HighExpression <- factor(ifelse(gene_level_coverage3_filtered_genes$MeanofMeanBootSamps >=HighExpressionCutoffFilteredGenes, 1, 0))
      levels(gene_level_coverage3_filtered_genes$HighExpression) <- c("LowExpr", "HighExpr")
      
      gene_level_coverage3$FourClassValue <- 0
      gene_level_coverage3$FourClassValue[gene_level_coverage3$HighExpression=="LowExpr" & gene_level_coverage3$HighInfRV=="LowInfRV"] <- 1
      gene_level_coverage3$FourClassValue[gene_level_coverage3$HighExpression=="HighExpr" & gene_level_coverage3$HighInfRV=="LowInfRV"] <- 2
      gene_level_coverage3$FourClassValue[gene_level_coverage3$HighExpression=="LowExpr" & gene_level_coverage3$HighInfRV=="HighInfRV"] <- 3
      gene_level_coverage3$FourClassValue[gene_level_coverage3$HighExpression=="HighExpr" & gene_level_coverage3$HighInfRV=="HighInfRV"] <- 4
      
      gene_level_coverage3$FourClassValue <- factor(gene_level_coverage3$FourClassValue)
      levels(gene_level_coverage3$FourClassValue) <- c("LowExprLowInfRV", "HighExprLowInfRV", "LowExprHighInfRV", "HighExprHighInfRV")
      
      #Additionally create Gene Unique factors based on quintiles and deciles of values less than 1
        #Just keep values exactly equal to in one group
      
      # c(0,.25,.5,.75, 1)
      quantilevals_all<- quantile(gene_level_coverage3$uniqueness_ratio[gene_level_coverage3$uniqueness_ratio!=1], probs = seq(0,1.01,.2), na.rm=T)
      quantilevals_filtered <- quantile(gene_level_coverage3_filtered_genes$uniqueness_ratio[gene_level_coverage3_filtered_genes$uniqueness_ratio!=1], probs = seq(0,1.01,.2), na.rm=T)
      #labels=c("10\\%", "20\\%", "30\\%", "40\\%", "50\\%","60\\%", "70\\%", "80\\%", "90\\%", "100\\%")
      
      gene_level_coverage3[,"UniqFactorNew"] <- cut(gene_level_coverage3$uniqueness_ratio,
                                                                   breaks=c(quantilevals_all, 1),
                                                                   include.lowest=TRUE)
      
      gene_level_coverage3_filtered_genes[,"UniqFactorNew"] <- cut(gene_level_coverage3_filtered_genes$uniqueness_ratio,
                                                         breaks=c(quantilevals_filtered, 1),
                                                         include.lowest=TRUE)
      
      
      levels(gene_level_coverage3[,"UniqFactorNew"])[nlevels(gene_level_coverage3[,"UniqFactorNew"])] <- "1"
      levels(gene_level_coverage3_filtered_genes[,"UniqFactorNew"])[nlevels(gene_level_coverage3_filtered_genes[,"UniqFactorNew"])] <- "1"
      
      
      gene_level_coverage3[,"UniqFactorNewConsCutoffs"] <- cut(gene_level_coverage3$uniqueness_ratio,
                                                    breaks=c(0, 0.75, 0.95, 1),
                                                    include.lowest=TRUE)
      
      gene_level_coverage3_filtered_genes[,"UniqFactorNewConsCutoffs"] <- cut(gene_level_coverage3_filtered_genes$uniqueness_ratio,
                                                                   breaks=c(0, 0.75, 0.95, 1),
                                                                   include.lowest=TRUE)
      
      
      #We considered calculating coverage quantiles using a normal distribution instead of NB, but this is not used since the confidence bounds can go negative
      #GenerateCoveragePlots("Norm95ConfInvCoverage", PlotHeight = PlotHeight, pdflatex_loc = pdflatex_loc, PlotsSaveDir = PlotsSaveDir, PlotNameModifier = PlotNameModifier, ZerosExcludedFromCoverageModifier = ZerosExcludedFromCoverageModifier,SimNumberToUse = SimNumberToUse)
      GenerateCoveragePlots("TwoPointFiveNinetySevenPointFiveCoverage", PlotHeight = PlotHeight, pdflatex_loc = pdflatex_loc, PlotsSaveDir = PlotsSaveDir, PlotNameModifier = PlotNameModifier, ZerosExcludedFromCoverageModifier = ZerosExcludedFromCoverageModifier,SimNumberToUse = SimNumberToUse)
      GenerateCoveragePlots("NB95InvCoverage", PlotHeight = PlotHeight, pdflatex_loc = pdflatex_loc, PlotsSaveDir = PlotsSaveDir, PlotNameModifier = PlotNameModifier, ZerosExcludedFromCoverageModifier = ZerosExcludedFromCoverageModifier,SimNumberToUse = SimNumberToUse)
    }


#Now, code below here will generate gene-level plots

cell_level_files_dir <- paste0(base_dir, "Sim", SimNumberToUse, "/", "alevinInfRepData/CellLevelFiles/")
cell_names <- colnames(Alevin_counts)
ncells <- length(cell_names)
#Some example genes
genes_to_use <- c("ENSG00000124208.16", "ENSG00000205542.11", "ENSG00000251562.8", "ENSG00000260537.2", "ENSG00000249437.7",
                  "ENSG00000286834.1", "ENSG00000161939.19", "ENSG00000120913.23", "ENSG00000120071.14", "ENSG00000144445.17",
                  "ENSG00000004059.11", "ENSG00000000938.13", "ENSG00000240972.2", gene_level_coverage3_filtered_genes$gene_id[1:10])

gene_level_datT <- vector(mode = "list", length = ncells)
for(i in 1:length(cell_names)){
  print(i)
  curr_cell <- cell_names[i]
  curr_dat <- loadRData(paste0(cell_level_files_dir, "infRepDat", curr_cell, ".RData"))
  curr_dat2 <- curr_dat[genes_to_use,]
  curr_dat2$Cell <- curr_cell
  curr_dat2$CellNumber <- i
  curr_dat2$gene_id <- rownames(curr_dat2)
  
  #Add in point estimate from alevin and the (scaled) splatter count that is the "true" count we are trying to cover
  curr_dat2$AlevinCount <- Alevin_counts2[genes_to_use, curr_cell]
  #Splatter count loaded above into the Splatter_counts2 object2 is scaled
  curr_dat2$ScaledSplatterCount <- Splatter_counts2[genes_to_use, curr_cell]
  
  gene_level_datT[[i]] <- data.frame(curr_dat2)
}

gene_level_dat <- data.frame(data.table::rbindlist(gene_level_datT))

gene_level_dat$BootMean <- rowMeans(gene_level_dat[,1:nboot])
gene_level_dat$BootSD<- apply(gene_level_dat[,1:nboot], 1, function(x){sd(x)})
gene_level_dat$BootVar<- apply(gene_level_dat[,1:nboot], 1, function(x){var(x)})

q_norm_975 <- qnorm(0.975)

gene_level_dat$NormIntLower <- (gene_level_dat$BootMean) - q_norm_975*(gene_level_dat$BootSD)
gene_level_dat$NormIntUpper <- (gene_level_dat$BootMean) + q_norm_975*(gene_level_dat$BootSD)

#Theta values here are for var = mu + (mu/theta) so set a maximum value to ensure sufficient extra poisson variation is properly
#present
#Mike Love had used phi = 0.001 in DESeq2 for phi=1/theta, so to correspond to that here set theta = 1/.001 = 1000 such that we specify a maximum
#This is to ensure variance is always (at least slightly) greater than the mean
#Theta may also be called phi, and is also equal to the "size" value from the base R functions
#So, set theta to 1000 if theta is negative or greater than 1000
gene_level_dat$theta_val <- (gene_level_dat$BootMean^2)/(gene_level_dat$BootVar - gene_level_dat$BootMean)
gene_level_dat$theta_val[!is.finite(gene_level_dat$theta_val)] <- 1000
gene_level_dat$theta_val[gene_level_dat$theta_val > 1000] <- 1000
gene_level_dat$theta_val[gene_level_dat$theta_val < 0] <- 1000

#Code to check the code is working as intended
#check_data <- rnbinom(10000, size = theta_vals[6], mu = mean_curr[6])
#check_data <- rnbinom(10000, size = theta_vals, mu = mean_curr)
#head(mean(check_data))
#head(var(check_data)) #Var should be mu + (1/theta)*mu^2

gene_level_dat$NBinomIntLower <- qnbinom(0.025, size = gene_level_dat$theta_val, mu = gene_level_dat$BootMean)
gene_level_dat$NBinomIntUpper <- qnbinom(0.975, size = gene_level_dat$theta_val, mu = gene_level_dat$BootMean)


gene_level_dat$TwoPointFiveIntLower <- apply(gene_level_dat[,1:nboot], 1, function(x){quantile(x, probs = 0.025)})
gene_level_dat$TwoPointFiveIntUpper <- apply(gene_level_dat[,1:nboot], 1, function(x){quantile(x, probs = 0.975)})

gene_level_dat$SplatterCountCoveredNormInt <- ifelse((gene_level_dat$ScaledSplatterCount >=gene_level_dat$NormIntLower) & (gene_level_dat$ScaledSplatterCount <=gene_level_dat$NormIntUpper), 1, 0)
gene_level_dat$SplatterCountCoveredPercInt <- ifelse((gene_level_dat$ScaledSplatterCount >=gene_level_dat$TwoPointFiveIntLower) & (gene_level_dat$ScaledSplatterCount <=gene_level_dat$TwoPointFiveIntUpper), 1, 0)
gene_level_dat$SplatterCountCoveredNBInt <- ifelse((gene_level_dat$ScaledSplatterCount >=gene_level_dat$NBinomIntLower) & (gene_level_dat$ScaledSplatterCount <=gene_level_dat$NBinomIntUpper), 1, 0)


PlotWidth <- 9
stdal <- TRUE
curr_save_dir <- paste0(PlotsSaveDir, "GeneLevelCoveragePlots", "/")
if(!dir.exists(curr_save_dir)){dir.create(curr_save_dir, recursive = T)}
setwd(curr_save_dir)


for(g in 1:length(genes_to_use)){
  
  curr_gene <- genes_to_use[g]
  
  curr_gene_uniqueness_val <- gene_level_coverage3[gene_level_coverage3$gene_id==curr_gene, "uniqueness_ratio"]
  
  if(length(curr_gene_uniqueness_val)==0){
    curr_gene_uniqueness_val_Print <- "NA"
  }else if(curr_gene_uniqueness_val==0){
    curr_gene_uniqueness_val_Print <- "0"
  }else if(curr_gene_uniqueness_val==1){
    curr_gene_uniqueness_val_Print <- "1"
  }else{
    curr_gene_uniqueness_val_Print <- fr(curr_gene_uniqueness_val, 3)
  }
    
  
  gene_level_dat_currT <- subset(gene_level_dat, gene_level_dat$gene_id==curr_gene)
  gene_level_dat_curr <- gene_level_dat_currT[order(gene_level_dat_currT$ScaledSplatterCount),]
  gene_level_dat_curr$OrderedCellNumber <- 1:ncells
  
  TwoPointFiveIntLower <- gene_level_dat_curr$TwoPointFiveIntLower
  TwoPointFiveIntUpper <- gene_level_dat_curr$TwoPointFiveIntUpper
  
  NormIntLower <- gene_level_dat_curr$NormIntLower
  NormIntUpper <- gene_level_dat_curr$NormIntUpper
  
  NBIntLower <- gene_level_dat_curr$NBinomIntLower
  NBIntUpper <- gene_level_dat_curr$NBinomIntUpper
  
  SplatterCountCoveredNormInt <- gene_level_dat_curr$SplatterCountCoveredNormInt
  SplatterCountCoveredPercInt <- gene_level_dat_curr$SplatterCountCoveredPercInt
  SplatterCountCoveredNBInt <- gene_level_dat_curr$SplatterCountCoveredNBInt
  
  min_obs_value <- min(c(0, NormIntLower, NBIntLower))
  if(min_obs_value==0){
    overall_min_value <- 0
  }else if(min_obs_value > 0){
    overall_min_value <- 0.90 * min_obs_value
  }else if(min_obs_value < 0){
    overall_min_value <- 1.10 * min_obs_value
  }
  
  #If the normal distribution based approach is used, the minimum value for each interval can go below 0 and distort the y axis limits
    #for the other methods
    #In this case, set the limit to be 0 (such that some of those intervals could potentially get cutoff)
  if(overall_min_value < 0){
    overall_min_value <- 0
  }
  
  overall_max_value <- 1.10*max(c(NormIntUpper, NBIntUpper, TwoPointFiveIntUpper, gene_level_dat_curr$ScaledSplatterCount))
  
  for(ty in 1:3){
    

    if(ty==1){
      IntervalKind <- "NormInt"
      
      Lower <- NormIntLower
      Upper <- NormIntUpper
      
      SplatterCountCovered <- SplatterCountCoveredNormInt
      
      #curr_title <- paste0("Cell-Level Coverage Results for 95\\% Interval from Normal Quantiles\nfor Gene ", curr_gene)
      if(PlotNameModifier=="ForPaper"){
        curr_title <- paste0(curr_gene)
      }else{
        curr_title <- paste0("Cell-Level Coverage Results for 95\\% Interval from Normal Quantiles (Top)\nand from Empirical Bootstrap Quantiles (Bottom)\nfor Gene ", curr_gene, " (Gene Uniqueness Value = ", curr_gene_uniqueness_val_Print, ")")
      }

    }else if(ty==2){
      IntervalKind <- "EmpBootInt"
      
      
      Lower <- TwoPointFiveIntLower
      Upper <- TwoPointFiveIntUpper
      
      SplatterCountCovered <- SplatterCountCoveredPercInt
      
      #curr_title <- paste0("Cell-Level Coverage Results for 95\\% Interval from Empirical Bootstrap Quantiles\nfor Gene ", curr_gene)
      curr_title <- "\n\n"
      
    }else if(ty==3){
      IntervalKind <- "NBInt"
      
      Lower <- NBIntLower
      Upper <- NBIntUpper
      
      SplatterCountCovered <- SplatterCountCoveredNBInt
      if(PlotNameModifier=="ForPaper"){
        curr_title <- paste0(curr_gene)
      }else{
        curr_title <- paste0("Cell-Level Coverage Results for 95\\% Interval from Negative Binomial Quantiles (Top)\nand from Empirical Bootstrap Quantiles (Bottom)\nfor Gene ", curr_gene, " (Gene Uniqueness Value = ", curr_gene_uniqueness_val_Print, ")")
      }

    }
    
    SplatterCountCoveredFactor <- factor(SplatterCountCovered, levels = c("0","1"))
    
    SplatterCountCoveredNo <- mean(SplatterCountCovered==0)*100
    SplatterCountCoveredYes <- mean(SplatterCountCovered==1)*100
    
    #Plot of InfRV by coverage and expression for all genes
    #fil_piece <- paste0("GeneLevelPlot", IntervalKind, curr_gene)
    fil_piece <- paste0("GeneLevelPlot", IntervalKind, g)
    fil_name <- paste0(curr_save_dir, fil_piece, PlotNameModifier, ".tex")
    tikz(file = fil_name, height=PlotHeight, width=PlotWidth, standAlone = stdal, sanitize = F)
    #paste0(CoverageVarNamePrint, " Coverage by InfRV and \nGene Expression for All Genes (", nrow(gene_level_coverage3), " Total)")
    
    #, group=group,color=group
    #,color=group
    #,color=group,shape=group
    #  scale_y_continuous(limits = c(-100, 1000))+
    #coord_cartesian(ylim = c(20, 73))+
    pl <- ggplot(data=gene_level_dat_curr, aes(x=OrderedCellNumber, y=ScaledSplatterCount, color = SplatterCountCoveredFactor)) +
      scale_color_manual(values = c("red", "black"), name="Simulated Count \nCovered", labels=c(paste0("No", " (", SplatterCountCoveredNo, "\\%)"),
                                                                                               paste0("Yes", " (", SplatterCountCoveredYes, "\\%)")), drop = FALSE)+
      geom_errorbar(aes(ymin=Lower, ymax=Upper))+
      scale_x_continuous(labels=NULL, breaks = 1:ncells)+
      scale_y_continuous(limits = c(overall_min_value, overall_max_value))+
      theme(axis.ticks.x = element_line(color = "red", size = ifelse(SplatterCountCovered,0,1)))+
      geom_point()+ylab("(Scaled) Simulated Count")+xlab("Cell")+
      theme(legend.text = element_text(size=rel(1), face = "bold"),legend.title = element_text(size=rel(1), face = "bold"), axis.text.y=element_text(size=rel(1.25), face = "bold"),
            axis.text.x=element_text(size=rel(1), face = "bold"), plot.title=element_text(size=rel(1.25), hjust=.5,face = "bold"), axis.title=element_text(size=rel(1.25), face = "bold"))+
      ggtitle(curr_title)
    print(pl)
    dev.off()
    SVBCompiletikzPlot(curr_save_dir, paste0(fil_piece, PlotNameModifier), pdflatex_loc)
    rm(fil_name)
    rm(fil_piece)
    
  }
  
}


print(gc())
}


