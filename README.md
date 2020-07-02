This repository contains code to reproduce results from "Compression of quantification uncertainty for scRNA-seq counts."  For simple sample code to calculate the uncertainty aware pvalues, see the file "UncertaintyAwarePvalueSampleCode.R."  For other code to reproduce analyses from the paper, see the rest of the code and the instructions below.


Last updated July 2, 2020


Instructions to reproduce the main analyses from the paper are given below.  Note that much of the code here utilizes functions defined in the "SingleCellProjectFunctions.R" file.

First, code to reproduce the main coverage and trajectory analysis results:

1. First, code in the file `QuantifyPBMC4K.sh' within the subdirectory `QuantifyPBMC4KData' will quantify the PBMC 4K data, and `SavePBMC4KDataToRData.R' will save the results in an R friendly format.  These results are used to assign realistic gene names to simulated data based on the rank of gene expression.

2. Then, generate the simulation objects using the splatter package (for the two-group difference simulations) or the dyntoy package (for the trajectory simulations).  The former is done in the file `GenerateSplatterObjectForSwish.R' within the `SwishAnalyses' subdirectory and the latter is done in the file `GenerateDyntoyObjects.R' within the `MainSimulationCode' subdirectory.

3. Then, run minnow to simulate the reads corresponding to the simulated counts and alevin to quantify the counts and generate the bootstrap replicates and compressed uncertainty estimates.  The file `RunMinnowAndAlevin.bash' runs minnow and alevin with 100 bootstrap replicates, and the file `RunAlevinWith20InfReps.bash' can be run after to repeat alevin with 20 bootstrap replicates.

4. Then, the file `SaveDatasetsToRAndCalculateCoverages.R' will save the alevin quantification data to R, calculate all the coverage results, and generate simulated pseudo inferential replicates.

5. Now, Run `SummarizeCoverageResults.R' to generate all coverage related plots.

6. Next, `RunTradeSeq.R' will run the tradeSeq code.

7. Then, the file `SummarizeTradeSeq.R' will import the full results from tradeSeq and compute the uncertainty aware pvalues.  These are saved in a format to be used for iCOBRA plotting.

8. The file `PlotTradeSeqPowerResults.R' will then generate the iCOBRA plots corresponding to the tradeSeq results.

Now, code for Mouse Embryo trajectory analysis:

9. Run `ImportMouseEmbryoResultsTo.R' to import the quantified data into R.

10. Run the file MouseEmbryoTrajectoryAnalysis.R to run the trajectory analysis, both forcing the EM and NoEM results to have the same cell clusters and lineages/pseudotimes and not forcing them to.

11. Run the file InfRepTradeSeqMouseEmbryoData.R to conduct the trajectory analysis for each pseudo-inferential replicate.

12. Run AnalyzeMouseEmbryoTradeSeqResults.R to combine all trajectory results and calculate the uncertainty-aware pvalues.

13. Run the file PlotsForMouseEmbryoTrajectoryAnalysis.R to generate all plots from the mouse embryo trajectory analyses.

Now, code corresponding to the swish and SplitSwish analyses:

14. For the swish results, first run the file SaveFullDataForSwish.R to save the full data (including bootstrap replicates) for the swish analysis.  Then, run the file `RunSwish.R' to run the full swish results.

15.  For the splitSwish results, first run the file "SplitSwish.R" to save the data (without the full bootstrap replicates) and save the necessary Snakemake file to run the splitSwish method.  The snakeMake file can then be run in command line shell.

16. Lastly, the iCOBRA plot comparing performance of swish to SplitSwish can be plotted using the file `PlotSplitSwishResults.R'.

