##Attachpackages
library(tidyverse)
library(caroline)
library(parallel)

##Load data
Allele_frequencies <- read.table("Allele_frequencies.delim", header = TRUE, sep = " ")
## ID minor allele

## Get the min and max numbers for the base pairs into new columns. 
Base_pairs <-Allele_frequencies %>%
  select(8:13)

##This gets the total number of base pairs present at a specific location.
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

##Remove sites with more than two alleles as our analysis cant account for these polyallelic sites
#there are also very few so removing them shouldnt affect the analysis.
Allele_frequencies <- Allele_frequencies %>%
  filter(count %in% c(4, 5))


##Make numeric to work with

Allele_frequencies$primary <- as.numeric(Allele_frequencies$primary)
Allele_frequencies$alternate <- as.numeric(Allele_frequencies$alternate)

##Remove columns that only have 1 allele present.

##This makes a dataframe of every location that is fixed within a landscape
## We want to remove these from the analysis since fixed sites dont tell us anything
One_allele <- Allele_frequencies %>%
  group_by(Chromosome, Location, Landscape) %>%
  summarize(Allele_diff = sum(primary - alternate)) %>%
  filter(Allele_diff == 0)


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




## Code to set up parallele processing for the loop
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

## We need to get the minor allele nucleotide to keep it constant for each
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

##Create a new DF that only has nucleotide information.
Base_only <- FRQ_test %>%
  ungroup() %>%
  select(3:7)

Base_only[is.na(Base_only)] <- 0



## Create a new row to put Minor allele BP letter into
FRQ_test$New <- rep(NA, nrow(FRQ_test))

## Get the minor allele Base into a new column

# export the objects you need for your calculations from your environment to each node's environment
clusterExport(clust, varlist=c("FRQ_test", "Base_only"))


parellel_major_BP <- do.call(rbind, clusterApplyLB(clust, 1:nrow(FRQ_test), function(i){
  
  
  bases <- names(Base_only)[which(Base_only[i,] == FRQ_test$Major[i])]
  bases2 <- bases[1]
  return(data.frame(Which_base = bases2)) 
}))

## Put the output of the function into a the new column
FRQ_test$New <- parellel_major_BP$Which_base

##Remove large object
rm(parellel_major_BP)

## Get the global allele frequency of the major allele.
## Higher global major allele frequency could mean changes away from it are
## Potentially more deleterious.
FRQ_test_2 <- FRQ_test %>%
  group_by(Chromosome, Location) %>%
  summarise(sum = sum(A,`T`,C,G))

FRQ_test$Global_tot <- FRQ_test_2$sum

FRQ_test <- FRQ_test %>%
  mutate(Global = Major/Global_tot)
## Take out colums that will be duplicates when joining. 
FRQ_test <- FRQ_test %>%
  select(Chromosome, Location, New, Major, Global, Global_tot)


print("Checkpoint 3")



## Join to add the minor allele frequency identifier to each landscape/location
Allele_frequencies <- left_join(Allele_frequencies, FRQ_test, by = c("Chromosome", "Location"))


clusterExport(clust, varlist=c("Allele_frequencies"))

parellel_Allele_frqs_name <- do.call(rbind, clusterApplyLB(clust, 1:nrow(Allele_frequencies), function(i){
  
  
  alt <- Allele_frequencies[i, which(names(Allele_frequencies) == Allele_frequencies$New[i])]
  return(data.frame(Alternate = alt))   # this is the end product for each iteration of the loop
}))

Allele_frequencies$alternate <- parellel_Allele_frqs_name$Alternate



print("Checkpoint 4")

##By dividing the number of minor nucleotides by the total number present at that location we get 
##the allele frequency of the minor allele mutate new column of frequencies.
##this value should be able to be over .5
Allele_frequencies <- Allele_frequencies %>%
  mutate(FRQ = alternate/primary)


##Remove unnecessary columns for analysis
Allele_frequencies <- Allele_frequencies %>%
  select(1,2, 3, 4,5 , 6,  FRQ, Global, Global_tot)

##Write individual info file
write.delim(Allele_frequencies, "Filtered_FRQs_individual_para.delim")



##Pivot wider so that data for each landscape is in a single row

Allele_frequencies <- Allele_frequencies %>%
  pivot_wider(names_from = c(Generation, Core_Edge), values_from = FRQ)

print("Checkpoint 5")

##Calculate the allele frequency changes between generation 0 and gen 8 for all landscapes and for both core and edge populations.
##Looks like there are many cases where there arent values for each generation.
##i.e there are NA's at sertain locations that might not have been sequenced.
Allele_frequencies <- Allele_frequencies %>%
  mutate(Change_stationary = `0_NA` - `8_NA`) %>%
  mutate(Change_core = `0_NA` - `8_C`) %>%
  mutate(Change_edge = `0_NA` - `8_E`)



##Change so that negative changes are a reduction in the major allele
##Positive change represents an increase in the major allele

Allele_frequencies$Change_stationary <- -1*(Allele_frequencies$Change_stationary)
Allele_frequencies$Change_core <- -1*(Allele_frequencies$Change_core)
Allele_frequencies$Change_edge <- -1*(Allele_frequencies$Change_edge)

## Write final file for the next step
write.delim(Allele_frequencies, "Filtered_FRQs_para.delim")
#########
