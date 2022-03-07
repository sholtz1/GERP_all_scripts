#!/bin/bash
# This script will create a master mpileup file, from which it will create
#	a master sync file. It will also create a single column .csv file
#	with the column order of the mpilup and sync files


#SBATCH --account=beetlegenes
#SBATCH --time=120:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sholtz1@uwyo.edu
#SBATCH --job-name=make_masters

module load swset/2018.05
module load gcc/7.3.0
module load trimmomatic/0.36
module load bwa/0.7.17
module load samtools/1.9


cd /gscratch/sholtz1




samtools mpileup -B $( ls *.clean.sort.bam ) > /gscratch/sholtz1/AllPops.mpileup


ls *.clean.sort.bam > ../AllPops_ColOrder.csv
