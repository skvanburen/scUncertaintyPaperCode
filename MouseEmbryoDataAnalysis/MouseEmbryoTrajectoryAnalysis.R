#Code for the trajectory analysis on the mouse embryo data
  #This includes code that fixes the cell clusters and lineages/pseudotimes to be the same across the EM vs noEM counts and code that does not

library(slingshot)
library(tradeSeq)
library(RColorBrewer)
library(tikzDevice)
library(readr)
library(dplyr)
current_tradeSeq_seed <- 2
source("/Users/Scott/Documents/Dissertation/SingleCellProject/code/Simulation Code/SingleCellProjectFunctions.R")
UseDifferentCellClusters <- FALSE

base_dir <- "/Users/Scott/Documents/Dissertation Data/SingleCellProject/MouseEmbryoData/"
PlotsSaveDir <- file.path(base_dir, "EmVsNoEMPlots/")
if(!dir.exists(PlotsSaveDir)){dir.create(PlotsSaveDir)}

gene_mappingT <- data.frame(read_tsv(file.path(base_dir, "GeneNamesMapping/gnames.txt"),
                                     col_names = F))
colnames(gene_mappingT) <- c("gene_id", "MGI")
gene_mapping <- distinct(gene_mappingT)

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

pdflatex_loc <- '/Library/TeX/texbin/pdflatex'

stdal <- TRUE

#Type EM results must be run first for this code to work
for(oute in 1:2){
  if(oute==1){
    FitSeparateLineagesForNoEM <- TRUE
  }else if(oute==2){
    FitSeparateLineagesForNoEM <- FALSE
  }
  
  
  assoc_res <- vector(mode = "list", length = 2)
  start_end_res <- vector(mode = "list", length = 2)
  pattern_res <- vector(mode = "list", length = 2)
  diffEnd_res <- vector(mode = "list", length = 2)
  earlyDE_res <- vector(mode = "list", length = 2)
  
  sce_res <- vector(mode = "list", length = 2)
  crv_res <- vector(mode = "list", length = 2)
  counts_res <- vector(mode = "list", length = 2)
  
  #Type EM results must be run first for this code to work
  for(run in 1:2){
    if(run==1){
      dat1 <- loadRData("/Users/Scott/Documents/Dissertation Data/SingleCellProject/MouseEmbryoData/QuantObj8.0UsingEM.RData")
      dat2 <- loadRData("/Users/Scott/Documents/Dissertation Data/SingleCellProject/MouseEmbryoData/QuantObj8.25UsingEM.RData")
      dat3 <- loadRData("/Users/Scott/Documents/Dissertation Data/SingleCellProject/MouseEmbryoData/QuantObj8.5UsingEM.RData")
      type <- "EM"
    }else if(run==2){ 
      dat1 <- loadRData("/Users/Scott/Documents/Dissertation Data/SingleCellProject/MouseEmbryoData/QuantObj8.0NoEM.RData")
      dat2 <- loadRData("/Users/Scott/Documents/Dissertation Data/SingleCellProject/MouseEmbryoData/QuantObj8.25NoEM.RData")
      dat3 <- loadRData("/Users/Scott/Documents/Dissertation Data/SingleCellProject/MouseEmbryoData/QuantObj8.5NoEM.RData")
      type <- "NoEM"
    }
    
    
    
    if(type=="EM"){
      set.seed(1)
      #sample(1:ncol(QuantObj8.0UsingEM$counts), 100)
      colsDat1 <- sample(1:ncol(dat1$counts), 500)
      colsDat2 <- sample(1:ncol(dat2$counts), 500)
      colsDat3 <- sample(1:ncol(dat3$counts), 500)
      
      t1T <- dat1$counts[,colsDat1]
      t2T <- dat2$counts[,colsDat2]
      t3T <- dat3$counts[,colsDat3]
      
      CellsUsedDat1 <- gtools::mixedsort(colnames(t1T))
      CellsUsedDat2 <- gtools::mixedsort(colnames(t2T))
      CellsUsedDat3 <- gtools::mixedsort(colnames(t3T))
      
      save(CellsUsedDat1, CellsUsedDat2, CellsUsedDat3, file = "/Users/Scott/Documents/Dissertation Data/SingleCellProject/MouseEmbryoData/CellsUsedtradeSeq.RData")
      
      t1T2 <- t1T[,CellsUsedDat1]
      t2T2 <- t2T[,CellsUsedDat2]
      t3T2 <- t3T[,CellsUsedDat3]
    }else{
      t1T2 <- dat1$counts[,CellsUsedDat1]
      t2T2 <- dat2$counts[,CellsUsedDat2]
      t3T2 <- dat3$counts[,CellsUsedDat3]
    }


    

    
    
    ordered_gene_names <- gtools::mixedsort(rownames(t1T2))
    
    t1 <- t1T2[ordered_gene_names,]
    t2 <- t2T2[ordered_gene_names,]
    t3 <- t3T2[ordered_gene_names,]
    
    sum(colnames(t1)%in% colnames(t2))
    sum(colnames(t1)%in% colnames(t3))
    #Rename cells slightly in case the same barcode is present for multiple runs
    colnames(t1) <- paste0(colnames(t1), "_1")
    colnames(t2) <- paste0(colnames(t2), "_2")
    colnames(t3) <- paste0(colnames(t3), "_3")
    rm(dat1)
    rm(dat2)
    rm(dat3)
    gc()
    
    sum(rownames(t1)!=rownames(t2))
    sum(rownames(t1)!=rownames(t3))
    
    
    InputCountsT <- cbind(t1,t2,t3)
    #rm(t1,t2,t3)
    gc()
    
    filt_func <- function(x, FixedNGenes = FALSE){
      ncells_high_exp <- sum(x >= 3)
      return(ncells_high_exp)
    }
    
    vals1 <- apply(InputCountsT, 1, filt_func)
    InputCountsT2 <- InputCountsT[vals1 > 10,]
    #Keep the same genes in the NoEM results as in the EM results
    if(type=="EM" | FitSeparateLineagesForNoEM==TRUE){
      #Filter genes based on needing a count of at least 3 in at least 10 cells
      
      #If type is "NoEM use the same cell order as for the EM results
      if(type=="EM"){
        #ColSIndex <- order(colSums(InputCountsT2))
        #InputCounts <- InputCountsT2[,ColSIndex]
        InputCounts <- as.matrix(InputCountsT2)
        InputCountsEM <- InputCounts
        dimCountsEM <- dim(InputCountsEM)
      }else{
        #Make the cells are in the same order to enable the use of the same clusters
        #cell_n <- colnames(InputCountsEM)
        InputCountsT3 <- as.matrix(InputCountsT2)
        InputCounts <- InputCountsT3[,FinalcountsCellNames]
      }

      
      
      
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
      filtered_genes <- rownames(counts)
      norm_counts <- norm_countsT
      
      #Even if fitting separate lineages for the NoEM results, keep the same cell clusters to ensure the results can be more directly compared
      #And that the start cluster for the NoEM results can be specified based on the start cluster from the EM results
      if(type=="EM" | UseDifferentCellClusters==TRUE){
        nKmeans <- 5
        cl <- kmeans(rd, centers = nKmeans)$cluster
        print(paste0("Dimension reduction is being done using kMeans with ", nKmeans, " clusters"))
        GMM <- FALSE
        #We now would like to use some prior information to determine the start cluster.  This is done using prior information from Figure 2
          #of the paper this data comes from, "A single-cell molecular map of mouse gastrulation and early organogenesis"
          #We specifically have cells from times 8.0, 8.25, and 8.50 in our analysis, and we choose genes to be marker genes that are highly expressed
          #for cells tending to come from time 8.0 (since the cells are expected to develop over real time such that cells
          #collected at time 8.5 should be expected to be further along the developmental trajectory
          #After manual insection of Figure 2b from the paper, we will choose marker genes for "Foregut 2" since these have a higher concentration of 
          #8.0 time point cells than other gut cells
          #Genes shown in the paper to be highly expressed in these "Foregut 2" cells are "Osr1", "Ripply3", "Hoxa1", "Irx2", "Fzd2", "Irx1", and we will choose
          #the cluster of cells with the highest average expression of these genes (in terms of ranks) to be the start cluster
        
        Mark_genes <- c("Osr1", "Ripply3", "Hoxa1", "Irx2", "Fzd2", "Irx1")
        Mark_genes_full_namesT <- gene_mapping[gene_mapping$MGI%in%Mark_genes,]
        Mark_genes_full_names <- Mark_genes_full_namesT$gene_id
        
        #Now, choose the cluster with the highest expression of the Marker genes as the starting cluster
          #Higher rank implies higher expression
        counts_sub <- subset(counts, rownames(counts)%in%Mark_genes_full_names)
        counts_rank <- apply(counts_sub, 1, rank)
        
        #Now, calculate mean rank of each gene across each cluster
        clus_mean_ranks <- matrix(NA, nrow = nKmeans, ncol = ncol(counts_rank))
        for(curr_clus in 1:length(unique(cl))){
          
        curr_cl_vals <- subset(cl, cl==curr_clus)
        counts_rank_sub <- subset(counts_rank, rownames(counts_rank)%in%names(curr_cl_vals))
        clus_mean_ranks[curr_clus,] <- colMeans(counts_rank_sub)
        }
        
        #Now, start cluster will be the one that has the highest average rank of expression across the marker genes
          #that pass filtering and remain present in the final analysis
          #After running, cluster "5" is chosen as the start, and this is the same one slingshot selects if letting it choose by leaving start.clus = NULL
        start.clus <- as.character(which(rowMeans(clus_mean_ranks)==max(rowMeans(clus_mean_ranks))))
        
        #End clusters are not specified because we don't know exactly how many lineages to expect
        end.clus <- NULL
      }else{
        start.clus <- c(crv_res[[1]]@lineages$Lineage1[1])
        end.clus <- c(crv_res[[1]]@lineages$Lineage1[length(crv_res[[1]]@lineages$Lineage1)], crv_res[[1]]@lineages$Lineage2[length(crv_res[[1]]@lineages$Lineage2)])
      }
      
      print(paste0("start.clus is ", start.clus))
      print(paste0("end.clus is ", end.clus))
      
      SCEObj <- SingleCellExperiment::SingleCellExperiment(assays = S4Vectors::List(counts = counts, norm_counts = norm_counts))
      
      SingleCellExperiment::reducedDims(SCEObj) <- S4Vectors::SimpleList(PCA = rd)
      
      
      SummarizedExperiment::colData(SCEObj)$kMeans <- cl
      
      slingshot_results <- slingshot(SCEObj, clusterLabels = 'kMeans', reducedDim = 'PCA', start.clus = start.clus, end.clus = end.clus)
      lin <- getLineages(SCEObj, clusterLabels = SummarizedExperiment::colData(slingshot_results)$kMeans, reducedDim = 'PCA')
      crv <- SlingshotDataSet(getCurves(lin))
      
      n_of_fit_lin <- length(crv@lineages)
      print(paste0("The number of lineages fit is ", n_of_fit_lin))
      print(paste0("The crv object information is printed below"))
      print(crv)
      
      TimeNumber <- as.numeric(as.character(lapply(strsplit(colnames(counts), "_"), function(x){x[2]})))
      names(TimeNumber) <- colnames(counts)
      coloring_var <- cl
      #coloring_var <- TimeNumber
      FullColors <- brewer.pal(9,"Set1")[coloring_var]
      names(FullColors) <- as.character(coloring_var)
      
      leg_col <- vector(mode = "character", length = nKmeans)
      for(jj in 1:length(unique(coloring_var))){
        leg_col[jj] <- unique(subset(FullColors, names(FullColors)==as.character(jj)))
      }
      if(oute==1){
        fil_piece <- paste0(type, "RDPlot")
        fil_name <- paste0(PlotsSaveDir, fil_piece, ".tex")
        tikz(file = fil_name, height=5.5, width=8.5, standAlone = stdal, sanitize = F)
        
        plotTitle <- type
        plot(rd, col = FullColors, pch=16, asp = 1, main = plotTitle)
        #legend(x = 40, y = 30, legend = sort(unique(coloring_var)), col = leg_col, pch = 19)
        legend("topleft", legend = sort(unique(coloring_var)), col = leg_col, pch = 19)
        lines(crv, lwd=2, col="black")
        dev.off()
        tryCatch(dev.off(), error = function(x){})
        SVBCompiletikzPlot(PlotsSaveDir, fil_piece, pdflatex_loc = pdflatex_loc)
        rm(fil_name)
        rm(fil_piece)
      }
     
      
      
      if(type=="EM"){
        FinalcountsCellNames <- colnames(counts)
        save(CellsUsedDat1, CellsUsedDat2, CellsUsedDat3, filtered_genes, FinalcountsCellNames, crv, file = "/Users/Scott/Documents/Dissertation Data/SingleCellProject/MouseEmbryoData/EMCountsLineageObjects.RData")
      }
      
    }else if(type=="NoEM" & FitSeparateLineagesForNoEM==FALSE){
      InputCountsT2 <- InputCountsT[filtered_genes,]
      #InputCounts <- InputCountsT2[,ColSIndex]
      InputCounts <- InputCountsT2
      counts <- as.matrix(InputCounts[,FinalcountsCellNames])
      if(sum(dim(counts)!=dimCountsEM)!=0){
        stop("The dimensions of the counts do not line up even though they should")
      }
    }
    
    
    
    
    set.seed(current_tradeSeq_seed)
    
    
    n_of_fit_lin <- length(crv@lineages)
    nKnots <- 6
    
    control <- mgcv::gam.control()
    control$maxit <- 1000 #set maximum number of iterations to 1K
    
    print(gc())
    print(paste0("tradeSeq results will be run with ", nKnots, " knots"))
    st1 <- proc.time()
    
    
    sce <- fitGAM(counts = counts, sds = crv, nknots = nKnots, sce = TRUE, control = control)
    
    
    association_test_res <- associationTest(sce, global = TRUE, lineages = FALSE)
    start_vs_end_test_res <- startVsEndTest(sce, global = TRUE, lineages = FALSE)
    
    patternTestRes <- tryCatch(patternTest(sce, global = TRUE, pairwise = FALSE), error = function(x){print(x)})
    diffEndTestRes <- tryCatch(diffEndTest(sce, global = TRUE, pairwise = FALSE), error = function(x){print(x)})
    earlyDETestRes <- tryCatch(earlyDETest(sce, knots = c(1, 3), global = TRUE, pairwise = FALSE), error = function(x){print(x)})
    
    
    
    
    assoc_res[[run]] <- association_test_res
    start_end_res[[run]] <- start_vs_end_test_res
    
    pattern_res[[run]] <- patternTestRes
    diffEnd_res[[run]] <- diffEndTestRes
    earlyDE_res[[run]] <- earlyDETestRes
    
    sce_res[[run]] <- sce
    crv_res[[run]] <- crv
    
    counts_res[[run]] <- counts
  }
  
  if(UseDifferentCellClusters==TRUE){
    fil_mod <- ""
  }else{
    fil_mod <- "SameCellClusters"
  }
  
  if(FitSeparateLineagesForNoEM==TRUE){
    save(crv_res, assoc_res, start_end_res, pattern_res, diffEnd_res, earlyDE_res, counts_res, sce_res, CellsUsedDat1, CellsUsedDat2, CellsUsedDat3, file = paste0("/Users/Scott/Documents/Dissertation Data/SingleCellProject/MouseEmbryoData/tradeSeqMouseDataResultsEMandNoEMSeparateLineagesForNoEM", fil_mod, ".RData"))
  }else{
    save(crv_res, assoc_res, start_end_res, pattern_res, diffEnd_res, earlyDE_res, counts_res, sce_res, CellsUsedDat1, CellsUsedDat2, CellsUsedDat3, file = "/Users/Scott/Documents/Dissertation Data/SingleCellProject/MouseEmbryoData/tradeSeqMouseDataResultsEMandNoEMSameLineages.RData")
  }
}
