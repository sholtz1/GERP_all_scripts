#!/bin/bash

#SBATCH --account=beetlegenes
#SBATCH --time=120:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sholtz1@uwyo.edu
#SBATCH --job-name=samsort




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
samtools view -q 20 -b -o /gscratch/sholtz1/${Forward[$i]}.temp.bam /gscratch/sholtz1/${Forward[$i]}.sam
samtools sort -o /gscratch/sholtz1/sort_bam/${Forward[$i]}.sort.bam /gscratch/sholtz1/${Forward[$i]}.temp.bam

done


#rm /gscratch/sholtz1/*temp*

#mv  /project/beetlegenes/sholtz1/beetle_test/TestFiles/*.sam /gscratch/sholtz1/

