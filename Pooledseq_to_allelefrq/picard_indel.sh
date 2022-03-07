#!/bin/bash



#SBATCH --account=beetlegenes
#SBATCH --time=120:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sholtz1@uwyo.edu
#SBATCH --job-name=picard_realign2


#bwa mem fastq ref sam

module load swset/2018.05
module load gcc/7.3.0
module load picard/2.20.1
module load samtools/1.9
module load gatk/4.1.8.0-py27
module load jdk/8u172-b11 

cd /gscratch/sholtz1/sort_bam/


InFiles=($( find *.bam ))

picard=/project/beetlegenes/Paths_software/picard.jar
gatk=/project/beetlegenes/Paths_software/
ref=/project/beetlegenes/sholtz1/beetle_test/Tribolium_ref_2020.fna
ref_prefix=`echo $ref | cut -d "." -f1`

#load the nessary tools

NumFiles=${#InFiles[@]}




java -jar $picard CreateSequenceDictionary \
      R=$ref \
      O=${ref_prefix}.dict

samtools faidx $ref

#per bam steps
#loop over list of bam files
for (( i=0; i<${NumFiles}; i++));
do

echo "***************** this is ${InFiles[$i]} *******************"
	#set up bam names and prefixes for subsequent output
	cBam=`echo ${InFiles[$i]} | awk -F "/" '{print $NF}'`
	cPrefix=`echo $cBam | sed 's/.sorted.bam//g'`
	cPath=`echo ${InFiles[$i]}_path`
mkdir $cPath	
echo $cPrefix

echo "***************** path is $cPath *******************"

	## PICARD TOOLS ##
	
	#clean bam
	java -jar $picard  CleanSam VALIDATION_STRINGENCY=LENIENT I=${InFiles[$i]} O=${InFiles[$i]}.picardClean

echo "*********************************"	
echo "************COMMAND 1 ***********"
echo "*********************************"	
	
	#index bam	
	java -jar $picard BuildBamIndex I=${InFiles[$i]}.picardClean 
 	#provide dummy sequence information, RGSM will appear as column name in vcf
echo "*********************************"	
echo "************COMMAND 2 ***********"
echo "*********************************"	
	
	java -jar $picard AddOrReplaceReadGroups VALIDATION_STRINGENCY=LENIENT INPUT=${InFiles[$i]}.picardClean OUTPUT=$cPath/${cPrefix}_headers.sorted.bam SO=coordinate RGLB=NA RGPL=NA RGPU=NA RGSM=$cPrefix
	#index new bam
echo "*********************************"	
echo "************COMMAND 3 ***********"
echo "*********************************"	
	
	java -jar $picard BuildBamIndex VALIDATION_STRINGENCY=LENIENT I=$cPath/${cPrefix}_headers.sorted.bam 

echo "*********************************"	
echo "************COMMAND 4 ***********"
echo "*********************************"	
	
	## GATK ##
	
	#find regions near indels that should be realinged 
	#java -jar $gatk -T RealignerTargetCreator -R $ref -I $cPath/${cPrefix}_duplicates.sorted.bam -o $cPath/${cPrefix}_realignTargets.intervals
	java -jar $gatk/GenomeAnalysisTK.jar -T  RealignerTargetCreator -R $ref -I $cPath/${cPrefix}_headers.sorted.bam -o $cPath/${cPrefix}_realignTargets.intervals 
	#perform realignment	
	java -jar $gatk/GenomeAnalysisTK.jar -T  IndelRealigner -R $ref -I $cPath/${cPrefix}_headers.sorted.bam -targetIntervals $cPath/${cPrefix}_realignTargets.intervals -o $cPath/${cPrefix}_realigned.bam	
	
	#index new bam
	java -jar $gatk/GenomeAnalysisTK.jar -T BuildBamIndex VALIDATION_STRINGENCY=LENIENT I=$cPath/${cPrefix}_realigned.bam 
	
	#cohort approach, build site calls for all possible sites (gVCF table)
	#java -jar $gatk -T HaplotypeCaller -R $ref -I $cPath/${cPrefix}_realigned.bam --emitRefConfidence GVCF --variant_index_type LINEAR --variant_index_parameter 128000 -o $cPath/${cPrefix}.raw.snps.indels.g.vcf
	#mv $cPath/*realigned* /project/beetlegenes/sholtz1/beetle_test/CleanBam	
	## CLEAN UP ##
	#rm $cPath/*duplicates*
	#rm $cPath/*headers*
	#rm *picard*
done

	#example for perfoming joint genotying
#java -jar ../GenomeAnalysisTK.jar -T GenotypeGVCFs -R /home/nkane/Research/Project_Sunflower/Reference_Bronze/HA412.v1.0.bronze.20140814.fasta --variant p3_B7.raw.snps.indels.g.vcf --variant p3_C1.raw.snps.indels.g.vcf --variant p3_F1.raw.snps.indels.g.vcf -o testBams.vcf

