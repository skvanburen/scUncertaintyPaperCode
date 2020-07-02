#Code to Combine the Mouse Embryo Data tradeSeq results
library(data.table)
library(readr)
library(dplyr)
base_dir <- "/pine/scr/s/k/skvanbur/SingleCellProject/MouseEmbryoData/"
curr_save_dir <- paste0(base_dir, "tradeSeqInfRepResults/")

source("~/SingleCellProject/SingleCellProjectFunctions.R")
for(lin_fit_type in 1:2){
  for(test_type in 1:2){
    if(test_type==1){
      test <- "start_vs_end_test_res"
      test2 <- "start_end_res"
      save_nam <- "StartEnd"
    }else if(test_type==2){
      test <- "association_test_res"
      test2 <- "assoc_res"
      save_nam <- "Association"
    }
    
    #Could use only 20 infreps instead of 100, for now just 100 is being used
    #n_infrep_used <- 100
    n_infrep_used <- 20
    
    if(n_infrep_used==20){
      nInfRep_saveMod <- "20InfReps"
    }else{
      nInfRep_saveMod <- ""
    }
    CurrResT <- vector(mode = "list", length = n_infrep_used)
    for(i in 1:n_infrep_used){
      print(paste0("i is ", i, " out of ", n_infrep_used))
      curr_save_fil_name <- paste0(curr_save_dir, "PInfRepNum", i, ".RData")
      curr_res <- loadRData(curr_save_fil_name, objNameToGet = test)
      curr_res$gene_id <- rownames(curr_res)
      curr_res$infRepNum <- i
      CurrResT[[i]] <- curr_res
    }
    
    
    
    InvChiSQTrans <- function(x, MostCommonDF){
      ChiSQStats <- qchisq(1 - x, df = MostCommonDF)
      m_val <- mean(ChiSQStats, na.rm = T)
      FinalP <- 1-pchisq(m_val, df=MostCommonDF)
      return(FinalP)
    }
    
    CurrRes <- data.frame(data.table::rbindlist(CurrResT))
    InfRepPval <- "pvalue"
    MostCommonDF <- as.numeric(DescTools::Mode(full_tradeSeqres_boot[,"df"], na.rm = T))
    WaldPvalAggregate <- do.call(data.frame, aggregate(CurrRes[,InfRepPval], by = list(gene_id = CurrRes$gene_id), FUN = function(x){c(InvChiSQ = InvChiSQTrans(x, MostCommonDF), Pval50Perc = quantile(x, probs = 0.50, na.rm=T), Pval75Perc = quantile(x, probs = 0.75, na.rm=T))}))
    colnames(WaldPvalAggregate) <- c("gene_id", "MeanStatAfterInvChiSq", "Pval50Perc", "Pval75Perc")
    
    #Now, load in the EM and NoEM results
    #These were run on local computer and uploaded to the cluster to enable easier debugging and generation of plots, etc
    if(lin_fit_type==1){
      NonInfRepResT <- loadRData(paste0(base_dir, "tradeSeqMouseDataResultsEMandNoEMSeparateLineagesForNoEMSameCellClusters.RData"), objNameToGet = test2)
      save_mod <- "SeparateLineages"
    }else if(lin_fit_type==2){
      NonInfRepResT <- loadRData(paste0(base_dir, "tradeSeqMouseDataResultsEMandNoEMSameLineages.RData"), objNameToGet =  test2)
      save_mod <- "SameLineages"
    }
    #
    
    colnames(NonInfRepResT[[1]])[colnames(NonInfRepResT[[1]])=="pvalue"] <- "pvalueEM"
    colnames(NonInfRepResT[[2]])[colnames(NonInfRepResT[[2]])=="pvalue"] <- "pvalueNoEM"
    
    NonInfRepResT2EM <- NonInfRepResT[[1]][, "pvalueEM", drop = F]
    NonInfRepResT2NoEM <- NonInfRepResT[[2]][, "pvalueNoEM", drop = F]
    
    NonInfRepResT2EM$gene_id <- rownames(NonInfRepResT2EM)
    NonInfRepResT2NoEM$gene_id <- rownames(NonInfRepResT2NoEM)
    
    
    NonInfRepRes <- merge(NonInfRepResT2EM, NonInfRepResT2NoEM, by = "gene_id")
    
    
    AllRes <- merge(NonInfRepRes, WaldPvalAggregate, by = "gene_id")
    
    
    gene_mappingT <- data.frame(read_tsv(file.path(base_dir, "GeneNamesMapping/gnames.txt"),
                                         col_names = F))
    colnames(gene_mappingT) <- c("gene_id", "MGI")
    gene_mapping <- distinct(gene_mappingT)
    
    UnAdjPvals <- merge(AllRes, gene_mapping, by = "gene_id")
    UnAdjPvals[,2][UnAdjPvals[,2]==0] <- .Machine$double.eps 
    UnAdjPvals[,3][UnAdjPvals[,3]==0] <- .Machine$double.eps 
    UnAdjPvals[,4][UnAdjPvals[,4]==0] <- .Machine$double.eps 
    UnAdjPvals[,5][UnAdjPvals[,5]==0] <- .Machine$double.eps 
    
    
    AdjPvals <- UnAdjPvals
    AdjPvals[,2] <- p.adjust(UnAdjPvals[,2], method = "fdr")
    AdjPvals[,3] <- p.adjust(UnAdjPvals[,3], method = "fdr")
    AdjPvals[,4] <- p.adjust(UnAdjPvals[,4], method = "fdr")
    AdjPvals[,5] <- p.adjust(UnAdjPvals[,5], method = "fdr")
    colnames(AdjPvals)[2:5] <- paste0("Adj", colnames(AdjPvals)[2:5])
    
    assign(paste0("AdjPvals", save_nam, save_mod, nInfRep_saveMod), AdjPvals)
    save(list = paste0("AdjPvals", save_nam, save_mod, nInfRep_saveMod), file = paste0(base_dir, "AllAdjPvals", save_nam, save_mod, nInfRep_saveMod, ".RData"))
    
    #subset(AdjPvals[,-1], AdjPvals$MGI %in% c("Eno1b", "Eno1", "Hmgb1", "Nme2", "Nme2", "Cox7c", "Gfpt2", "Tnnt1", "Has2", "Taf9", "Cited1", "cdx1"))
  }
  
}



#Now, merge in InfRV information  across the three time points
