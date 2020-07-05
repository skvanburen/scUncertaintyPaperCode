#!/bin/bash

#  GenerageDBGFileGENCODEv32.sh
#  
#
#  Created by Scott Van Buren on 1/20/20.
#  

#SBATCH --time=23:59:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=16g

LD_LIBRARY_PATH=~/bin/minnow-latest17Jan2020/lib ~/bin/minnow-latest17Jan2020/bin/minnow index -r ~/SingleCellProject/gencode_v32_files/gencode.v32.transcripts.fa -k 101 -f 20 --tmpdir /pine/scr/s/k/skvanbur/TwoPaCoTemp -p 4 -o ~/SingleCellProject/gencode_v32_files/dbg_gencode_v32_SVB_k101.gfa
