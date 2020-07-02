#!/bin/bash

#  
#
#  Created by Scott Van Buren on 1/27/20.
#


#SBATCH --time=23:59:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=16g
#SBATCH --output=/pine/scr/s/k/skvanbur/SingleCellProject/PlotTradeSeqPowerResultsslurmLogs/slurmLog_%a.out
module add gcc/9.1.0
mkdir -p /pine/scr/s/k/skvanbur/SingleCellProject/PlotTradeSeqPowerResultsRLogs/
srun ~/bin/R CMD BATCH --vanilla ~/SingleCellProject/PlotTradeSeqPowerResults.R /pine/scr/s/k/skvanbur/SingleCellProject/PlotTradeSeqPowerResultsRLogs/PlotTradeSeqPowerResults_$SLURM_ARRAY_TASK_ID.Rout
