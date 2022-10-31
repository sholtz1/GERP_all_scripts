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

##This gets the total number or base pairs present at a specific location.
Allele_frequencies$primary <- rowSums(Base_pairs)

##makes a column to tell how many alleles are present for that patch
Allele_frequencies$count <- apply(Base_pairs, 1, function(x) length(which(x=="0")))



### This section finds the base pair count of the minor SNP
Base_pairs[Base_pairs == 0] <- NA

f1 <- function(x){
  min(x, na.rm = TRUE)
}

Allele_frequencies$alternate <- apply(X=Base_pairs, MARGIN=1, FUN=f1)


##Remove unnecessary large object

rm("Base_pairs")

print("Checkpoint 1")
##Remove sites with more than two alleles
Allele_frequencies <- Allele_frequencies %>%
  filter(count %in% c(4, 5))


##Make numeric and

Allele_frequencies$primary <- as.numeric(Allele_frequencies$primary)
Allele_frequencies$alternate <- as.numeric(Allele_frequencies$alternate)

##Remove columns that only have 1 allele present.

##This makes a dataframe of every location that is fixed within a landscape
## We want to remove these from the analysis since fixed sites dont tell us anything
One_allele <- Allele_frequencies %>%
  group_by(Chromosome, Location, Landscape) %>%
  summarize(Allele_diff = sum(primary - alternate)) %>%
  filter(Allele_diff == 0)

## Add code to make a DF of locations that are fixed at with two different bases

##Start with only alleles that are fixed somewhere
Fixed_only <- Allele_frequencies %>%
  filter(primary == alternate)

## Next we need to get the BP of each allele into a new column for each landscape
#and generation

## Make new row
Fixed_only$Fixed <- rep(NA, nrow(Fixed_only))

## Create DF so we only look at BP columns
Fixed_only_bases <- Fixed_only %>%
  ungroup() %>%
  select(8:13)




## Code to set up parallel processing for the loop
no_cores <- 48
# Setup cluster
clust <- makeCluster(no_cores) 

print(no_cores)
# export the objects you need for your calculations from your environment to each node's environment
clusterExport(clust, varlist=c("Fixed_only_bases", "Fixed_only"))



## Get the allele BP into a new column for all fixed sites

parellel_allele_BP <- do.call(rbind, clusterApplyLB(clust, 1:nrow(Fixed_only), function(i){
  
  
    temp <- names(Fixed_only_bases)[which(Fixed_only_bases[i,] == Fixed_only$alternate[i])]
    
  return(data.frame(Fixed = temp))  
}))


Fixed_only$Fixed <- parellel_allele_BP$Fixed

rm(parellel_allele_BP)


## This code finds locations and landscapes where two different alleles are fixed.
Fixed_only_sum <- Fixed_only %>%
  group_by(Chromosome, Location, Landscape) %>%
  summarise(diff_BP = length(unique(Fixed))) %>%
  filter(diff_BP > 1)
  
## filter our data removing sites identified above use paste to identify 
#each site uniquely but don't remove sites with diff alleles fixed.

Allele_frequencies <- Allele_frequencies %>%
  filter(!(paste(Chromosome, Location, Landscape) %in% 
             paste(One_allele$Chromosome,One_allele$Location,One_allele$Landscape)) |
           (paste(Chromosome, Location, Landscape) %in% 
              paste(Fixed_only_sum$Chromosome,Fixed_only_sum$Location,Fixed_only_sum$Landscape))
         )

print("Checkpoint 2")

## We need to get the alternate allele nucleotide to keep it constant for each
##Landscape


##create a new dataframe with the total number of each base pair for each landscape/location 
FRQ_test <- Allele_frequencies %>%
  group_by(Chromosome, Location) %>%
  summarize(A = sum(A), `T` = sum(`T`), G = sum(G), 
            C = sum(C), del = sum(del))



##Make sure all base pair columns are numeric to that we can use the min function
FRQ_test$A <- as.numeric(FRQ_test$A)
FRQ_test$`T` <- as.numeric(FRQ_test$`T`)
FRQ_test$C <- as.numeric(FRQ_test$C)
FRQ_test$G <- as.numeric(FRQ_test$G)
FRQ_test$del <- as.numeric(FRQ_test$del)

## Make max function with rm.na
f2 <- function(x){
  max(x, na.rm = TRUE)
}


##Make a new column with the maximum number of alleles for that row
FRQ_test$Major <-apply(FRQ_test[,3:7],1,FUN=f2)
FRQ_test$Major <-as.numeric((FRQ_test$Major))

##Make a new column with the minimum number of alleles for each row
FRQ_test[FRQ_test == 0] <- NA
FRQ_test$minor <-apply(FRQ_test[,3:7],1,FUN=f1)
FRQ_test$minor <-as.numeric((FRQ_test$minor))



##Create a new DF that only has nuclotide information.
Base_only <- FRQ_test %>%
  ungroup() %>%
  select(3:7)

Base_only[is.na(Base_only)] <- 0



## Create a new row to put Minor allele BP letter into
FRQ_test$New <- rep(NA, nrow(FRQ_test))

## Get the major allele BP into a new column

# export the objects you need for your calculations from your environment to each node's environment
clusterExport(clust, varlist=c("FRQ_test", "Base_only"))


parellel_major_BP <- do.call(rbind, clusterApplyLB(clust, 1:nrow(FRQ_test), function(i){
  
  
  bases <- names(Base_only)[which(Base_only[i,] == FRQ_test$Major[i])]
  bases2 <- bases[1]
  return(data.frame(Which_base = bases2)) 
}))


FRQ_test$New <- parellel_major_BP$Which_base



########################################################################################
#########put the latter for the minor allele in a new column


## Create a new row to put Minor allele BP letter into
FRQ_test$alt <- rep(NA, nrow(FRQ_test))

## Get the major allele BP into a new column

# export the objects you need for your calculations from your environment to each node's environment
clusterExport(clust, varlist=c("FRQ_test", "Base_only"))


parellel_minor_BP <- do.call(rbind, clusterApplyLB(clust, 1:nrow(FRQ_test), function(i){
  
  
  bases <- names(Base_only)[which(Base_only[i,] == FRQ_test$minor[i])]
  bases2 <- bases[1]
  return(data.frame(Which_base = bases2))  
}))


FRQ_test$alt <- parellel_minor_BP$Which_base

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




