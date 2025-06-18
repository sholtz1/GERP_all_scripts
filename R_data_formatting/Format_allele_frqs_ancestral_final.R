##Attach/detach packages, ggpubr can mess up tidyverse grouping
library(tidyverse)
library(parallel)

##Load data
Allele_frequencies <- read.table("/project/beetlegenes/sholtz1/beetle_test/Allele_frequencies.delim", header = TRUE, sep = " ")
## ID minor allele in gen 0 and the be sure to keep that allele as the minor allele for the other generations. 


## Get the min and max numbers for the base pairs into new columns. 
Base_pairs <-Allele_frequencies %>%
  select(8:13)


##This gets the total number or base pairs present at a specific location.
Allele_frequencies$total_bases <- rowSums(Base_pairs)

##makes a column to tell how many alleles are absent for that patch
Allele_frequencies$alleles_count <- apply(Base_pairs, 1, function(x) length(which(x=="0")))



## change 0's to NA's so that the min function finds the nucleotide with the least bases present
Base_pairs[Base_pairs == 0] <- NA

## make min function that can ignore NA's
f1 <- function(x){min(x, na.rm = TRUE)}

## Get a vector with just the necueotide names
Allele_names <- names(Base_pairs)


## make a function that can identify which neucleotide has the fewest bases present
f2 <- function(x){
  Allele_index <- which.min(x)
  Allele <- Allele_names[Allele_index]
  return(Allele)
}

## fun function to get the count of the nuceotide with the least bases and also the a column with the nuclotide identified
Allele_frequencies$alternate_count <- apply(X = Base_pairs, MARGIN = 1, FUN = f1)
Allele_frequencies$alternate_allele <- apply(X = Base_pairs, MARGIN = 1, FUN = f2)

## make sure these number are numeric so we can do math with them
Allele_frequencies$total_bases <- as.numeric(Allele_frequencies$total_bases)
Allele_frequencies$alternate_count <- as.numeric(Allele_frequencies$alternate_count)


## Get biallelic sites that arent fixed across our entire population
Allele_frequencies <- Allele_frequencies %>%
  ## this filters for sites that have either 1 or two different alleles
  filter(alleles_count %in% c(4,5)) %>%
  ## should we be gourping by landscape? If we want to find fixation across the whole population fixation within the landscape might not be what we want?
  ## or because gene surfing cant act on sites that are already fixed do we not want to worry about these sites?
  group_by(Chromosome, Location, Landscape) %>%
  # make a new column where the number will be 0 if the site is fixed across an entire landscape
  mutate((Allele_diff = sum(total_bases - alternate_count))) %>%
  ## this finds sites where there are multiple alternate alleles
  mutate(diff_BP = length(unique(alternate_allele))) %>% 
  ungroup() %>%
  ##finally we filter for sites taht arent fixed in a landscape or are fixed for diffierent alleles in diffierent landscapes
  filter(Allele_diff !=0 | diff_BP > 1)




print("Checkpoint 1")


## We need to get the ancestral allele neucleotide to keep it constant for each Landscape
## read in ancestral alleles
Anc_Alleles <- read.delim("/project/beetlegenes/sholtz1/GERP/ancestral_alleles/All_ANC.delim", header = TRUE, sep =  " ")


## Join to add the ancestral allele identifier to each landscape/location
Allele_frequencies <- left_join(Allele_frequencies, Anc_Alleles, by = c("Chromosome", "Location"))

##Filter out sites where we dont have information for ancestral allele
Allele_frequencies <- Allele_frequencies %>% filter(!is.na(Base))


##Remove sites where the ancestral allele doesnt appear at all in our populations
Allele_frequencies <- Allele_frequencies %>%
  group_by(Chromosome, Location) %>%
  mutate(
    summed_value = case_when(
      Base == "A" ~ sum(A),
      Base == "G" ~ sum(G),
      Base == "T" ~ sum(`T`),
      Base == "C" ~ sum(C)
    )
  )  %>% filter(summed_value != 0)


## get the number of ancestral alleles present in each landscape into a new column

clusterExport(clust, varlist=c("Allele_frequencies"))

parellel_Allele_frqs_name <- do.call(rbind, clusterApplyLB(clust, 1:nrow(Allele_frequencies), function(i){
  
  
  # get the number in the column identified by the ancestral allele
  alt <- Allele_frequencies[i, which(names(Allele_frequencies) == Allele_frequencies$Base[i])]
  return(data.frame(Alternate = alt))   # this is the end product for each iteration of the loop
}))

Allele_frequencies$ancestral <- parellel_Allele_frqs_name$Alternate



print("Checkpoint 2")

##By dividing the number of minor neucleotides by the total number present at that location we get the allele frequency of the minor allele.
###mutate new column of frequencies. We do this for both the ancestral and major alleles to keep that information

Allele_frequencies <- Allele_frequencies %>%
  mutate(FRQ_ancestral = ancestral/primary)




write.table(Allele_frequencies, "Derived_FRQs_test.delim")
#########
