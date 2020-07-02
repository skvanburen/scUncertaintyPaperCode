#This code provides some sample code for how to calculate the uncertainty-aware pvalues after repeating an analysis on separate 
  #simulated pseudo-inferential replicates

#This code assumes the data frame full_pvalues has results for each gene across all replicates (which each being a separate row) and a column named "gene_id" corresponding to 
  #the gene information

MostCommonDF <- as.numeric(DescTools::Mode(full_pvalues[,"df"], na.rm = T))


InvChiSQTrans <- function(x, MostCommonDF){
  ChiSQStats <- qchisq(1 - x, df = MostCommonDF)
  m_val <- mean(ChiSQStats, na.rm = T)
  FinalP <- 1-pchisq(m_val, df=MostCommonDF)
  return(FinalP)
}


#Change any raw pvalues that show up as exactly 0 to a very small number > Machine eps since they technically aren't exactly zero
full_pvalues[,"pvalue"][full_pvalues[,"pvalue"]==0] <- 2e-16

WaldPvalAggregate <- do.call(data.frame, aggregate(full_pvalues[,"pvalue"], by = list(gene_id = full_pvalues$gene_id), FUN = function(x){c(InvChiSQ = InvChiSQTrans(x, MostCommonDF), Pval50Perc = quantile(x, probs = 0.50), Pval75Perc = quantile(x, probs = 0.75))}))
colnames(WaldPvalAggregate) <- c("gene_id", "MeanStatAfterInvChiSq", "Pval50Perc", "Pval75Perc")
WaldPvalAggregate[,"gene_id"] <- as.character(WaldPvalAggregate[,"gene_id"])


RawPvalues <- data.frame(WaldPvalAggregate)
AdjPvalues <- data.frame(apply(RawPvalues, 2, function(x) {p.adjust(x, method = "fdr")}))
