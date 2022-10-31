library(caroline)
library(tidyverse)
library(ape)



## Read in fold dataset
sitefolds <- read.table("sitefold.txt", header = TRUE)
##Change column names for combining
colnames(sitefolds) <- c("Chromosome", "Location", "Fold")
##Read in formatted allele frequencies
allele <- read.table("Filtered_FRQs_para.delim", header = TRUE)
## Join fold data with allele frequency dataset

allele_fold <- left_join(allele, sitefolds, by = c("Chromosome", "Location"))


##Next lets attach the full exon information

##Make names match
Full_filtered_data  <- allele_fold

## Read in annotation data for tribolium
Trib_annotation <- read.gff("/project/beetlegenes/sholtz1/GERP/data/sequences/ref/Tribolium_annotation.gff")

## Remove refseq rows since they contain the whole gnome
Trib_annotation <- Trib_annotation %>%
  filter(type != "region")


## We need to find out for all sites  whether they fall inside an exon or not

## Remove duplicates from the annotation
Trib_annotation_unique <- Trib_annotation[!duplicated(Trib_annotation[c('seqid', 'start', 'end')]),]

## Change to chacacter for comparison
Trib_annotation_unique$seqid <- as.character(Trib_annotation_unique$seqid)

## Take out all locations in the annotation that arent exons
Trib_annotation_unique <- Trib_annotation_unique %>%
  filter(type == "exon")
## make a new row to put results into
Full_filtered_data$Exon <- rep(NA, nrow(Full_filtered_data))

for(i in 1:nrow(Full_filtered_data)){
  # Find if a row in our data matches with an exon from the refrence.
  GERP_row <- which(Trib_annotation_unique$seqid == Full_filtered_data$Chromosome[i] &
                      Trib_annotation_unique$start <= Full_filtered_data$Location[i] &
                      Trib_annotation_unique$end >= Full_filtered_data$Location[i])
# If its matches an exon put yes if it doesnt put no into the new column
  if(is_empty(GERP_row)){Yes_No <- "No"} else{Yes_No <- "Yes"}

  Full_filtered_data$Exon[i] <- Yes_No

}

### read in all GERP information so that it can be combined with allele FRQ data
NC_007416.3_rm.rates.bed <-read.delim("/project/beetlegenes/sholtz1/GERP/analyses/gerp/used_inc_tribolium/NC_007416.3_rm.rates.bed", header = FALSE)
NC_007417.3_rm.rates.bed <-read.delim("/project/beetlegenes/sholtz1/GERP/analyses/gerp/used_inc_tribolium/NC_007417.3_rm.rates.bed", header = FALSE)
NC_007418.3_rm.rates.bed <-read.delim("/project/beetlegenes/sholtz1/GERP/analyses/gerp/used_inc_tribolium/NC_007418.3_rm.rates.bed", header = FALSE)
NC_007419.2_rm.rates.bed <-read.delim("/project/beetlegenes/sholtz1/GERP/analyses/gerp/used_inc_tribolium/NC_007419.2_rm.rates.bed", header = FALSE)
NC_007420.3_rm.rates.bed <-read.delim("/project/beetlegenes/sholtz1/GERP/analyses/gerp/used_inc_tribolium/NC_007420.3_rm.rates.bed", header = FALSE)
NC_007421.3_rm.rates.bed <-read.delim("/project/beetlegenes/sholtz1/GERP/analyses/gerp/used_inc_tribolium/NC_007421.3_rm.rates.bed", header = FALSE)
NC_007422.5_rm.rates.bed <-read.delim("/project/beetlegenes/sholtz1/GERP/analyses/gerp/used_inc_tribolium/NC_007422.5_rm.rates.bed", header = FALSE)
NC_007423.3_rm.rates.bed <-read.delim("/project/beetlegenes/sholtz1/GERP/analyses/gerp/used_inc_tribolium/NC_007423.3_rm.rates.bed", header = FALSE)
NC_007424.3_rm.rates.bed <-read.delim("/project/beetlegenes/sholtz1/GERP/analyses/gerp/used_inc_tribolium/NC_007424.3_rm.rates.bed", header = FALSE)
NC_007425.3_rm.rates.bed <-read.delim("/project/beetlegenes/sholtz1/GERP/analyses/gerp/used_inc_tribolium/NC_007425.3_rm.rates.bed", header = FALSE)

#Change column names so they match for joining

colnames(NC_007416.3_rm.rates.bed) <- c("Chromosome","not_important", "Location", "Max_RS", "RS")
colnames(NC_007417.3_rm.rates.bed) <- c("Chromosome","not_important", "Location", "Max_RS", "RS")
colnames(NC_007418.3_rm.rates.bed) <- c("Chromosome","not_important", "Location", "Max_RS", "RS")
colnames(NC_007419.2_rm.rates.bed) <- c("Chromosome","not_important", "Location", "Max_RS", "RS")
colnames(NC_007420.3_rm.rates.bed) <- c("Chromosome","not_important", "Location", "Max_RS", "RS")
colnames(NC_007421.3_rm.rates.bed) <- c("Chromosome","not_important", "Location", "Max_RS", "RS")
colnames(NC_007422.5_rm.rates.bed) <- c("Chromosome","not_important", "Location", "Max_RS", "RS")
colnames(NC_007423.3_rm.rates.bed) <- c("Chromosome","not_important", "Location", "Max_RS", "RS")
colnames(NC_007424.3_rm.rates.bed) <- c("Chromosome","not_important", "Location", "Max_RS", "RS")
colnames(NC_007425.3_rm.rates.bed) <- c("Chromosome","not_important", "Location", "Max_RS", "RS")


#Put all GERP information into one object for joining
rates_bed <- rbind(NC_007416.3_rm.rates.bed, NC_007417.3_rm.rates.bed, NC_007418.3_rm.rates.bed,
                   NC_007419.2_rm.rates.bed, NC_007420.3_rm.rates.bed, NC_007421.3_rm.rates.bed,
                   NC_007422.5_rm.rates.bed, NC_007423.3_rm.rates.bed, NC_007424.3_rm.rates.bed,
                   NC_007425.3_rm.rates.bed)

## Add gerp data for each specific locus in each chromosome
AlleleFreqs <- left_join(Full_filtered_data, rates_bed, by = c("Chromosome", "Location"))

#Remove unnecessary column

AlleleFreqs <- AlleleFreqs %>%
  select(-not_important)

## Write file for next step of analysis
write_delim(AlleleFreqs, "/project/beetlegenes/sholtz1/GERP/Full_data_bigtree.delim")


