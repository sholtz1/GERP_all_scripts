#!/bin/bash

#SBATCH --account=beetlegenes
#SBATCH --time=120:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sholtz1@uwyo.edu
#SBATCH --job-name=Drosophila_melanogaster_genomic.fna
module load kentutils/1.04.0

bash  make_alignments.sh /project/beetlegenes/sholtz1/GERP/data/sequences/query/Asbolus_verrucosus_genomic.fna

