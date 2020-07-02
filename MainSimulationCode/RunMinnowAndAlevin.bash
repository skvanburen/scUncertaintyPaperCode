#!/bin/bash

#  RunMinnowAndAlevin.sh
#  
#
#  Created by Scott Van Buren on 2/25/20.
#  

#SBATCH --time=47:59:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --mem=49g
#SBATCH --output=/pine/scr/s/k/skvanbur/SingleCellProject/RunMinnowAndAlevinslurmLogs/slurmLog_%a.out

let nbootsamps=100

LD_LIBRARY_PATH=~/bin/minnow-latest30Mar2020/lib ~/bin/minnow-latest30Mar2020/bin/minnow simulate --splatter-mode --g2t ~/SingleCellProject/gencode_v32_files/tx2gene.tsv --inputdir /pine/scr/s/k/skvanbur/SingleCellProject/SimulationResults/Sim$SLURM_ARRAY_TASK_ID/SplatterFiles/ --PCR 4 -r ~/SingleCellProject/gencode_v32_files/gencode.v32.transcripts.fa -e 0.01 -p 2 -o /pine/scr/s/k/skvanbur/SingleCellProject/SimulationResults/Sim$SLURM_ARRAY_TASK_ID/MinnowOutput/ --dbg --gfa ~/SingleCellProject/gencode_v32_files/dbg_gencode_v32_SVB_k101.gfa/dbg.gfa -w ~/SingleCellProject/737K-august-2016.txt --bfh ~/SingleCellProject/gencode_v32_files/bfh_gencode_v32.txt --custom --gencode --librarySize 1200000
#### --librarySize above should be set to 1,200,000 for the dynverse simulations (to make sure enough molecules were produced) but only needs to be 200,000 for splatter ones

#Now, validate the minnow output to ensure the resulting library sizes match what is expected
    #First, generate the splatter gene count in .gz format
LD_LIBRARY_PATH=~/bin/minnow-latest30Mar2020/lib ~/bin/minnow-latest30Mar2020/bin/validate dump -f <( gunzip -c /pine/scr/s/k/skvanbur/SingleCellProject/SimulationResults/Sim$SLURM_ARRAY_TASK_ID/MinnowOutput/hg_100_S1_L001_R1_001.fastq.gz ) -t ~/SingleCellProject/gencode_v32_files/tx2gene.tsv -o /pine/scr/s/k/skvanbur/SingleCellProject/SimulationResults/Sim$SLURM_ARRAY_TASK_ID/MinnowOutput/splatter_genecount.gz


#The python script below checks and verifies that the generation of molecules via minnow is correct
#When checking the log for this file, the output should return a check mark and not an x if everything is working correctly
module load python/3.6.6
python ~/SingleCellProject/splatter_simulation_sanity_check.py /pine/scr/s/k/skvanbur/SingleCellProject/SimulationResults/Sim$SLURM_ARRAY_TASK_ID/MinnowOutput/alevin/ /pine/scr/s/k/skvanbur/SingleCellProject/SimulationResults/Sim$SLURM_ARRAY_TASK_ID/MinnowOutput/splatter_genecount.gz


#Now, run Alevin and generate the bootstrap replicates
module load salmon/1.1.0; salmon alevin -l ISR -1 /pine/scr/s/k/skvanbur/SingleCellProject/SimulationResults/Sim$SLURM_ARRAY_TASK_ID/MinnowOutput/hg_100_S1_L001_R1_001.fastq.gz -2 /pine/scr/s/k/skvanbur/SingleCellProject/SimulationResults/Sim$SLURM_ARRAY_TASK_ID/MinnowOutput/hg_100_S1_L001_R2_001.fastq.gz -i /pine/scr/s/k/skvanbur/SingleCellProject/GENCODEv32SalmonIndex/ --chromium -p 2 -o /pine/scr/s/k/skvanbur/SingleCellProject/SimulationResults/Sim$SLURM_ARRAY_TASK_ID/AlevinOutput/ --tgMap ~/SingleCellProject/gencode_v32_files/tx2gene.tsv --whitelist /pine/scr/s/k/skvanbur/SingleCellProject/SimulationResults/Sim$SLURM_ARRAY_TASK_ID/MinnowOutput/alevin/quants_mat_rows.txt --numCellBootstraps $nbootsamps --dumpFeatures
