# Load packages
library(tidyverse)
library(caroline)
library(parallel)

# Load data
Allele_frequencies <- read.table("Allele_frequencies.delim", header = TRUE, sep = " ")

# Isolate the number of reads per allele for each location in a new data frame
# for subsequent analyses
Base_pairs <-Allele_frequencies %>%
  select(8:13)
Base_names <- names(Allele_frequencies)[8:13]

# Calculate the total number of reads present at each specific location.
Allele_frequencies$total_reads <- as.numeric(rowSums(Base_pairs))

# Calculate how many alleles are present for that location
Allele_frequencies$num_alleles <- apply(Base_pairs, 1, function(x) length(which(x != "0")))

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

# Remove unnecessary large object and print output indicating progress within 
# the script
rm("Base_pairs")
print("Checkpoint 1")

# Remove sites with more than two alleles as our analysis can't account for these 
# polyallelic sites. There are also very few so removing them shouldn't affect the 
# analyses.
Allele_frequencies <- Allele_frequencies %>%
  filter(num_alleles <= 2)

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

# First, create a new dataframe with the total number of each base pair for each
# landscape/location 
FRQ_test <- Allele_frequencies %>%
  group_by(Chromosome, Location) %>%
  filter(Generation == 0) %>%
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

## Filter for sites with 3 golbal alleles
f7 <- function(x){
  min_allele <- length(which(x != 0))
  return(min_allele)
}

FRQ_test$Number_alleles <- apply(FRQ_test[,3:7], 1, FUN = f7)

FRQ_test <- FRQ_test %>%
  filter(Number_alleles == 2) %>%
  select(-Number_alleles)


########### Note from Topher #################
# The way I get the global total above is slightly different from your way because it
# includes counts of deletions (del) whereas you only included counts from A,T,C, and G
# This doesn't lead to large changes for most of the data, but there were about 45 rows
# in my subsetted data where it did make a difference. I'm happy to chat about this more
# but I think we should include the del column here because it's a read that doesn't 
# correspond to the major allele, meaning the global frequency of the major allele should
# be slightly lower because of it (most of the places where our results differed in this
# were because you calculated a global frequency of 1 and I calculated something just below
# one.)
########### Note from Topher #################

FRQ_test <- FRQ_test %>%
  mutate(Global = Major_count/Global_tot) %>%   # Use that to calculate the global (across population) frequency of the major allele
  select(Chromosome, Location, Major_count, Major_allele, Global, Global_tot) # Only keep certain columns to reduce space needs

print("Checkpoint 3")
Allele_frequencies <- Allele_frequencies %>%
  filter(total_reads > 9)

# Join to add the major allele frequency identifier to each landscape/location
Allele_frequencies <- left_join(Allele_frequencies, FRQ_test, by = c("Chromosome", "Location"))

# Calculate the frequency of the identified major allele in each population
Allele_names <- names(Allele_frequencies)[8:13]
# The indices for the x vector in the function below correspond to:
# x[10] --> max_allele
# x[12] --> Major_allele
# x[9] --> max_count
# x[7] --> total_reads
f5 <- function(x){
  if(x[10] == x[12]){
    return(as.numeric(x[9]) / as.numeric(x[7]))
  } else{
    MajorColumn <- which(Allele_names == x[12])
    MajorCount <- as.numeric(x[MajorColumn])
    return(MajorCount / as.numeric(x[7]))
  }
}
Allele_frequencies$FRQ <- apply(X = Allele_frequencies[1:nrow(Allele_frequencies),8:21], MARGIN = 1, FUN = f5)

#Allele_frequencies <- Allele_frequencies %>%
  mutate(FRQ = max_count/total_reads)
# Remove unnecessary columns for analysis
Allele_frequencies <- Allele_frequencies %>%
  select(Chromosome, Location, Ref, Landscape, Generation, Core_Edge, FRQ, Global, Global_tot)

# Write individual info file
write.delim(Allele_frequencies, file = "Filtered_FRQs_individual_para_TWL.delim")

# Pivot wider so that data for each landscape is in a single row and changes
# between generations and between core and edge can be easily calculated
Allele_frequencies <- Allele_frequencies %>%
  pivot_wider(names_from = c(Generation, Core_Edge), values_from = FRQ)

print("Checkpoint 5")

# Calculate the allele frequency changes between generation 0 and gen 8 for all 
# landscapes and for both core and edge populations.
# Looks like there are many cases where there aren't values for each generation.
# i.e there are NA's at certain locations that might not have been sequenced.
Allele_frequencies <- Allele_frequencies %>%
  mutate(Change_stationary = `8_NA` - `0_NA`) %>%
  mutate(Change_core = `8_C` - `0_NA`) %>%
  mutate(Change_edge = `8_E` - `0_NA`)

## Write final file for the next step
write.delim(Allele_frequencies, "Filtered_FRQs_para_TWL.delim")
#########
