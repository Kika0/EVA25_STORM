#!/bin/bash

#SBATCH -J dc_boot-job
#SBATCH -c 101
#SBATCH -o dc_boot.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=m.speers@lancaster.ac.uk
#SBATCH --mem 300000

srun Rscript /beegfs/client/default/speersm/data_challenge/conditional_extremes_bootstrapping/model_fit.R