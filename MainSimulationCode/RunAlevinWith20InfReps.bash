#!/bin/bash

#  RunMinnowAndAlevin.sh
#  
#
#  Created by Scott Van Buren on 2/25/20.
#  

#SBATCH --time=23:59:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --mem=29g
#SBATCH --output=/pine/scr/s/k/skvanbur/SingleCellProject/RunAlevinWith20InfRepsslurmLogs/slurmLog_%a.out

#This code repeats the alevin call while generating only 20 bootstrap replicates (and thus estimating compressed uncertainty based on 20 replicates also)

module load salmon/1.1.0; salmon alevin -l ISR -1 /pine/scr/s/k/skvanbur/SingleCellProject/SimulationResults/Sim$SLURM_ARRAY_TASK_ID/MinnowOutput/hg_100_S1_L001_R1_001.fastq.gz -2 /pine/scr/s/k/skvanbur/SingleCellProject/SimulationResults/Sim$SLURM_ARRAY_TASK_ID/MinnowOutput/hg_100_S1_L001_R2_001.fastq.gz -i /pine/scr/s/k/skvanbur/SingleCellProject/GENCODEv32SalmonIndex/ --chromium -p 2 -o /pine/scr/s/k/skvanbur/SingleCellProject/SimulationResults/Sim$SLURM_ARRAY_TASK_ID/AlevinOutput/ --tgMap ~/SingleCellProject/gencode_v32_files/tx2gene.tsv --whitelist /pine/scr/s/k/skvanbur/SingleCellProject/SimulationResults/Sim$SLURM_ARRAY_TASK_ID/MinnowOutput/alevin/quants_mat_rows.txt --numCellBootstraps 20 --dumpFeatures
