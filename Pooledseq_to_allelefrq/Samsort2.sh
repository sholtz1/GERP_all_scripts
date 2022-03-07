#!/bin/bash

#SBATCH --account=beetlegenes
#SBATCH --time=120:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sholtz1@uwyo.edu
#SBATCH --job-name=samsort2




module load swset/2018.05
module load gcc/7.3.0
module load bwa/0.7.17
module load samtools/1.9

cd /project/beetlegenes/2017BeetleData/DNAPools/


Forward=($( find *_R1_* ))
#Count number of files in a list
NumFiles=${#Forward[@]}
# For loop do it the number of times in the NumFiles list i++ means to add 1 each loop
for (( i=0; i<${NumFiles}; i++));
do



# Sort it to make it a .sort file
# Sort it to make it a .sort file
samtools view -q 20 -b -o /gscratch/sholtz1/${Forward[$i]}.temp2.bam /gscratch/sholtz1/sort_bam/${Forward[$i]}.sort.bam_path/${Forward[$i]}.sort.bam_realigned.bam
samtools sort -o /gscratch/sholtz1/${Forward[$i]}.clean.sort.bam /gscratch/sholtz1/${Forward[$i]}.temp2.bam

done

#rm /project/beetlegenes/sholtz1/beetle_test/TestFiles/*picard*
#rm -rf /project/beetlegenes/sholtz1/beetle_test/TestFiles/*_path
#rm  /project/beetlegenes/sholtz1/beetle_test/CleanBam/*.sort.bam_realigned.bam
#rm /gscratch/sholtz1/*temp*
