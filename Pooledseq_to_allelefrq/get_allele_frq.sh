#!/bin/bash



#SBATCH --account=beetlegenes
#SBATCH --time=160:00:00
#SBATCH --mem=200G
#SBATCH --partition=teton-gpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --cpus-per-task=1
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sholtz1@uwyo.edu
#SBATCH --job-name=allele_frequency



module load swset/2018.05
module load gcc/7.3.0
module load jdk/8u172-b11
module load intel/18.0.1
module load perl/5.30.1

##Use the mpileup2sync.jar script to make a sync file from the mpileup file
java -ea -Xmx100g -jar /project/beetlegenes/Paths_software/popoolation2/mpileup2sync.jar --input /gscratch/sholtz1/AllPops.mpileup --output /gscratch/sholtz1/AllPops.sync --fastq-type sanger --min-qual 20 --threads 30

##Use the snp-frequency-diff.pl script from Popoolation2 to get the SNP frequency differences among populations

perl /project/beetlegenes/Paths_software/popoolation2/snp-frequency-diff.pl --input /gscratch/sholtz1/AllPops.sync --output-prefix /gscratch/sholtz1/AllPops_allelefrq --min-count 4 --min-coverage 4 --max-coverage 22
