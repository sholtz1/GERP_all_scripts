This pipeline takes raw pooled sequencing reads filters them, preforms quality control, and takes to to a final product with counts of base pairs for each polymorphic location in each pool of sequencing data. 

Required software
trimmomatic/0.36
bwa/0.7.17
picard/2.20.1
samtools/1.9
gatk/4.1.8.0-py27
jdk/8u172-b11
VCFtools
Popoolation2




1.Run AlignPools.sh to trim reads and align the to the reference
2.Run Samsort.sh to create a .sort file\
3.Run picard_indel.sh to to soft-clip alignments beyond the end of reference sequences and to set the mapping quality to 0 for unmapped reads and to target and realign reads around indels.
4. Run Samsort2.sh 
6. Run make_masters.sh  run make_masters.sh to create a master mpileup file and a master sync file
7. Run get_allele_frq.sh to get from the pileup file to a file with base pair counts for each sequenced pool for all polymorphic locations
8. Run R_seperate_alleles.R to seperate the base pairs into colums that make them easier to work with.

Final oupput is called Allele_frequencies.delim

Chromosome Location Ref Landscape Generation Core_Edge S A T C G N del
NC_007416.3 26771 N 22 0 "NA" S31 0 16 0 0 0 0
NC_007416.3 26771 N 22 8 "NA" S38 0 21 0 0 0 0
NC_007416.3 26771 N 31 8 C S21 0 11 0 0 0 0
NC_007416.3 26771 N 31 8 E S25 0 9 0 0 0 0
