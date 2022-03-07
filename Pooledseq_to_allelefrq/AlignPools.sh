#!/bin/bash


# script to process the DNA sequence files and convert them to .sort files

#SBATCH --account=beetlegenes
#SBATCH --time=120:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sholtz1@uwyo.edu
#SBATCH --job-name=pipeline_test

module load swset/2018.05
module load gcc/7.3.0
module load trimmomatic/0.36
module load bwa/0.7.17
module load samtools/1.9

# Change the directory to the original folder with all the DNA pools. We don't need to change
#	the working directory back after this because the trimmomatic and bwa commands now
#	use full filepaths, but if you wanted you could use the cd command again to change the
#	working directory
cd /project/beetlegenes/2017BeetleData/DNAPools/

# Find all the file names for the forwad and backward pools
Forward=($( find *_R1_*))
Backward=($( find *_R2_*))
NumFiles=${#Forward[@]}

for (( i=0; i<${NumFiles}; i++));
do

# Trim the reads and save the intermediate files in the scratch directory
trimmomatic PE -threads 30 -phred33 /project/beetlegenes/2017BeetleData/DNAPools/${Forward[$i]} /project/beetlegenes/2017BeetleData/DNAPools/${Backward[$i]} /gscratch/sholtz1/${Forward[$i]}_forward_paired.fq.gz /gscratch/sholtz1/${Forward[$i]}_forward_unpaired.fq.gz /gscratch/sholtz1/${Backward[$i]}_reverse_paired.fq.gz /gscratch/sholtz1/${Backward[$i]}_reverse_unpaired.fq.gz ILLUMINACLIP:/project/beetlegenes/sholtz1/beetle_test/NexteraPE-PE.fa:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:15 MINLEN:100

# Align both the forward and backward reads to the reference genome
bwa mem -t 30 /project/beetlegenes/sholtz1/beetle_test/Tribolium_ref_2020.fna  /gscratch/sholtz1/${Forward[$i]}_forward_paired.fq.gz /gscratch/sholtz1/${Backward[$i]}_reverse_paired.fq.gz > /gscratch/sholtz1/${Forward[$i]}.sam
done


# Then you can delete or move intermediate files from your scratch directory depending on if you want to keep them
#mv /gscratch/sholtz1/*.sam /project/beetlegenes/sholtz1/beetle_test/TestFiles 

#rm /gscratch/sholtz1/*.gz


