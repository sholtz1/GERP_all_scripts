##Attach packages
library(tidyverse)
library(caroline)
library(parallel)
library(ape)


### Most of this code mirrors the Format allele frequency script, but ends differently
## so generate the file necessary for running Ensembl VEP 


##Load data
Allele_frequencies <- read.table("Allele_frequencies.delim", header = TRUE, sep = " ")


## Get the min and max numbers for the base pairs into new columns. 
Base_pairs <-Allele_frequencies %>%
  select(8:13)

Base_names <- names(Allele_frequencies)[8:13]

##This gets the total number or base pairs present at a specific location.
Allele_frequencies$primary <- rowSums(Base_pairs)

##makes a column to tell how many alleles are present for that patch
Allele_frequencies$num_alleles <- apply(Base_pairs, 1, function(x) length(which(x !="0")))


# Replace 0 counts with NA and then identify the frequency and identity of the dominant allele 
# at each location
# NOTE: This is not the true major allele as we calculate that across all populations later
Base_pairs[Base_pairs == 0] <- NA


# Now create two convenience functions to identify the count and identity of the dominant
# allele, then use the apply function to use them across the full dataset
f1 <- function(x){
  max_count <- max(x, na.rm = TRUE)
  return(max_count)
}
f2 <- function(x){
  max_allele <- which.max(x)
  return(Base_names[max_allele[1]])
}
Allele_frequencies$max_count <- apply(X=Base_pairs, MARGIN=1, FUN=f1)
Allele_frequencies$max_allele <- apply(X=Base_pairs, MARGIN=1, FUN=f2)


##Remove unnecessary large object

rm("Base_pairs")

print("Checkpoint 1")
##Remove sites with more than two alleles as we cant deal with polyallelic sites. 
Allele_frequencies <- Allele_frequencies %>%
  filter(num_alleles <= 2)


##Make numeric and

Allele_frequencies$primary <- as.numeric(Allele_frequencies$primary)
Allele_frequencies$alternate <- as.numeric(Allele_frequencies$alternate)


# Now filter to remove any sites that are fixed for the same allele across
# all populations. Importantly, this will still keep any sites that are fixed,
# but fixed for different alleles in different populations.


IncludedSites <- Allele_frequencies %>%
  group_by(Chromosome, Location, Landscape) %>%
  summarize(count_diff = sum(total_reads - max_count), allele_diff = length(unique(max_allele))) %>%
  filter((allele_diff > 1) | (count_diff) > 0)

Allele_frequencies <- Allele_frequencies %>%
  filter(paste(Chromosome, Location, Landscape) %in% 
           paste(IncludedSites$Chromosome,IncludedSites$Location,IncludedSites$Landscape))



print("Checkpoint 2")

# Now we need to identify the major allele for each genomic location by comparing
# across all our populations. This protects against erroneous identification of 
# major alleles due to stochasticity in the pooling and sequencing process.


##create a new dataframe with the total number of each base pair for each landscape/location 
FRQ_test <- Allele_frequencies %>%
  group_by(Chromosome, Location) %>%
  summarize(A = sum(A), `T` = sum(`T`), G = sum(G), 
            C = sum(C), del = sum(del))


# Use a similar approach to what we did before to generate the count and identity
# of the higher frequency allele for genomic location across all landscapes
Base_names <- names(FRQ_test)[3:7]
f3 <- function(x){
  max_count <- max(x, na.rm = TRUE)
  return(max_count)
}
f4 <- function(x){
  max_allele <- which.max(x)
  return(Base_names[max_allele[1]])
}
FRQ_test$Major_count <- apply(FRQ_test[,3:7], 1, FUN = f3)
FRQ_test$Major_allele <- apply(FRQ_test[,3:7], 1, FUN = f4)

# Get the global allele frequency of the major allele.
# Higher global major allele frequency could mean changes away from it are
# potentially more deleterious.
FRQ_test$Global_tot <- rowSums(FRQ_test[,3:7]) 



############ this section is where the code diverges form the format allele frequency code#################


##### filter for annotated areas of the genome and add in strand

trib_gff <- read.gff("Tribolium_annotation.gff")

Trib_annotation <- trib_gff %>%
  filter(type == "gene"| type == "CDS" | type == "exon")



## Remove duplicates from the annotation
Trib_annotation_unique <- Trib_annotation[!duplicated(Trib_annotation[c('seqid', 'start', 'end')]),]

## Change to character for comparison
Trib_annotation_unique$seqid <- as.character(Trib_annotation_unique$seqid)

## make a new row to put results into
FRQ_test$Exon <- rep(NA, nrow(FRQ_test))

for(i in 1:nrow(FRQ_test)){
  # Find if a row in our data matches with an exon from the refrence.
  GERP_row <- which(Trib_annotation_unique$seqid == FRQ_test$Chromosome[i] &
                      Trib_annotation_unique$start <= FRQ_test$Location[i] &
                      Trib_annotation_unique$end >= FRQ_test$Location[i])
  # If its matches an exon put yes if it doesnt put no into the new column
  if(is_empty(GERP_row)){Yes_No <- "No"} else{Yes_No <- "Yes"}
  
  FRQ_test$Exon[i] <- Yes_No
  
}

FRQ_test_coding <- FRQ_test %>%
  filter(Exon =="Yes")

##### add in strand information

FRQ_test_coding$strand <- rep(NA, nrow(FRQ_test_coding))

for(i in 1:nrow(FRQ_test_coding)){
  # Find if a row in our data matches with an exon from the refrence.
  Exon_row <- which(Trib_annotation_unique$seqid == FRQ_test_coding$Chromosome[i] &
                      Trib_annotation_unique$start <= FRQ_test_coding$Location[i] &
                      Trib_annotation_unique$end >= FRQ_test_coding$Location[i])
  # put the stand from the identified column into the df
  
  FRQ_test_coding$strand[i] <- Trib_annotation_unique[Exon_row, 7]
  
}


## use the annotation information to put which strand each coding region is on
FRQ_test_coding <- FRQ_test_coding %>%
  mutate(strand = case_when(strand == "1" ~ "-",
                            strand == "2" ~ "+"
                            
  ))
## format the columns to match the necessary format for ensembl VEP

Ensembl_VEP <- FRQ_test_coding %>%
  select(Chromosome, Location,New, alt, strand)%>%
  mutate(allele = paste(New, "/", alt)) %>%
  mutate(end = Location) %>% 
  select(-New, -alt) %>%
  select(Chromosome, Location, end, allele, strand)

## these are the column names VEP looks for

colnames(Ensembl_VEP) <-c("chromosome", "start", "end", "allele", "strand")


## the chromosome names are different in the reference that ensembl vep uses
Ensembal_VEP <- Ensembal_VEP %>%
  mutate(allele = gsub(" ", "", allele)) %>%
  mutate(chromosome = case_when(chromosome == "NC_007416.3" ~ "LGX",
                                chromosome == "NC_007417.3" ~ "LG2",
                                chromosome == "NC_007418.3" ~ "LG3",
                                chromosome == "NC_007419.2" ~ "LG4",
                                chromosome == "NC_007420.3" ~ "LG5",
                                chromosome == "NC_007421.3" ~ "LG6",
                                chromosome == "NC_007422.5" ~ "LG7",
                                chromosome == "NC_007423.3" ~ "LG8",
                                chromosome == "NC_007424.3" ~ "LG9",
                                chromosome == "NC_007425.3" ~ "LG10"))



## this is the file that can be run unsing the online ensembl VEP for tribolium
write.delim(Ensembl_VEP, "Ensembal_VEP.delim")




