#!/bin/bash

#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=5
#SBATCH --mem=32g


module load salmon/1.1.0
 salmon alevin -l ISR -1 /pine/scr/s/k/skvanbur/SingleCellProject/PBMC4KData/fastqs/pbmc4k_S1_L001_R1_001.fastq.gz /pine/scr/s/k/skvanbur/SingleCellProject/PBMC4KData/fastqs/pbmc4k_S1_L002_R1_001.fastq.gz  -2 /pine/scr/s/k/skvanbur/SingleCellProject/PBMC4KData/fastqs/pbmc4k_S1_L001_R2_001.fastq.gz /pine/scr/s/k/skvanbur/SingleCellProject/PBMC4KData/fastqs/pbmc4k_S1_L002_R2_001.fastq.gz -i /pine/scr/s/k/skvanbur/SingleCellProject/GENCODEv32SalmonIndex --chromium -p 5 -o /pine/scr/s/k/skvanbur/SingleCellProject/PBMC4KAlevinOutput --tgMap ~/SingleCellProject/gencode_v32_files/tx2gene.tsv  --numCellBootstraps 100 --dumpFeatures
