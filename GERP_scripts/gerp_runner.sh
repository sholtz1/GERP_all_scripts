#!/bin/bash

#SBATCH --account=beetlegenes
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sholtz1@uwyo.edu
#SBATCH --job-name=gerp_chromosome_NC_007417.3


bash run_gerp.sh NC_007416.3

