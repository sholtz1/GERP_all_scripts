##Attach/detach packages, ggpubr can mess up tidyverse grouping
detach("package:ggpubr", unload=TRUE)
library(tidyverse)


##Load data
Allele_frequencies <- read.table("Allele_frq_subset_Topher", header = TRUE, sep = " ")
## right now we arent keeping the minor allele constant between generations. 
## ID minor allele in gen 0 and the be sure to keep that allele as the minor allele for the other generations. 

## Get the min and max numbers for the base pairs into new columns. 
Base_pairs <-Allele_frequencies %>%
  select(8:13)

##This gets the total number or base pairs present at a specific location.
Allele_frequencies$primary <- rowSums(Base_pairs)

##makes a colum to tell how many alleles are present for that patch
Allele_frequencies$count <- apply(Base_pairs, 1, function(x) length(which(x=="0")))



### This section finds the base pair count of the minor SNP
Base_pairs[Base_pairs == 0] <- NA

f1 <- function(x){
  min(x, na.rm = TRUE)
}

Allele_frequencies$alternate <- apply(X=Base_pairs, MARGIN=1, FUN=f1)


##Remove unnecessary large object

rm("Base_pairs")

##Remove sites with more than two alleles
Allele_frequencies <- Allele_frequencies %>%
  filter(count %in% c(4, 5))


##Make numberic and

Allele_frequencies$primary <- as.numeric(Allele_frequencies$primary)
Allele_frequencies$alternate <- as.numeric(Allele_frequencies$alternate)

##Remove columns that only have 1 allele present.
##Are we removing cases when alleles were fixed in one generation but polymorphic in another. 

## This code might need to go after IDing the minor allele in case a new allele 
##is completely fixed

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


## Get the allele BP into a new column for all fixed sites
for (i in 1:nrow(Fixed_only)) {
  Fixed_only$Fixed[i] <- names(Fixed_only_bases)[which(Fixed_only_bases[i,] == Fixed_only$alternate[i])]
}


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

## We need to get the alternate allele neucleotide to keep it constant for each
##Landscape

##Test subset for easy computation
##FRQ_subset <- head(Allele_frequencies, 100000)


##create a new dataframe with the total number of each base pair for each landscape/location 
FRQ_test <- Allele_frequencies %>%
  group_by(Chromosome, Location, Landscape) %>%
  summarize(A = sum(A), `T` = sum(`T`), G = sum(G), 
            C = sum(C), del = sum(del))

##Turn 0 to NA so we can use the minumun function
FRQ_test[FRQ_test == 0] <- NA


##Make sure all base pair columns are numeric to that we can use the min function
FRQ_test$A <- as.numeric(FRQ_test$A)
FRQ_test$`T` <- as.numeric(FRQ_test$`T`)
FRQ_test$C <- as.numeric(FRQ_test$C)
FRQ_test$G <- as.numeric(FRQ_test$G)
FRQ_test$del <- as.numeric(FRQ_test$del)



##Make a new column with the minumium number of alleles for that row
FRQ_test$Minor <-apply(FRQ_test[,4:8],1,FUN=f1)
FRQ_test$Minor <-as.numeric((FRQ_test$Minor))

##Creat a new DF that only has nculotide information.
Base_only <- FRQ_test %>%
  ungroup() %>%
  select(4:8)


## Create a new row to put Minor allele BP letter into
FRQ_test$New <- rep(NA, nrow(FRQ_test))

## Get the minor allele BP into a new column
for (i in 1:nrow(FRQ_test)) {
  FRQ_test$New[i] <- names(Base_only)[which(Base_only[i,] == FRQ_test$Minor[i])]
   }

## Take out colums that will be duplicates when joining. 
FRQ_test <- FRQ_test %>%
  select(1:3, 10)

## Join to add the minor allele frequency identifier to each landscape/location
Allele_frequencies <- left_join(Allele_frequencies, FRQ_test, by = c("Chromosome", "Location", "Landscape"))



for (i in 1:nrow(Allele_frequencies)) {
  Allele_frequencies$alternate[i] <- Allele_frequencies[i, which(names(Allele_frequencies) == Allele_frequencies$New[i])]

  
}


##By dividing the number of minor neucleotides by the total number present at that location we get the allele frequency of the minor allele.
###mutate new column of frequencies.
##this value should be able to be over .5
##we will also probably need to change any vaues of 1 we get from this division to 0. 
##these are going to be locations where a base pair was not found in one of the generations
##but was found in the same landscape in another generation. 
Allele_frequencies <- Allele_frequencies %>%
  mutate(FRQ = alternate/primary)




table(Allele_frequencies$FRQ)




##Remove unnecessary columns for analysis
Allele_frequencies <- Allele_frequencies %>%
  select(1,2, 3, 4,5 , 6,  FRQ)



##Pivot wider so that data for each landscape is in a single row

Allele_frequencies <- Allele_frequencies %>%
  pivot_wider(names_from = c(Generation, Core_Edge), values_from = FRQ)

##Calculate the allele frequency changes between generation 0 and gen 8 for all landscapes and for both core and edge populations.
##Looks like there are many cases where there arent values for each generation.
##i.e there are NA's at sertain locations that might not have been sequenced.
Allele_frequencies <- Allele_frequencies %>%
  mutate(Change_stationary = `0_NA` - `8_NA`) %>%
  mutate(Change_core = `0_NA` - `8_C`) %>%
  mutate(Change_edge = `0_NA` - `8_E`)

##Change to absolute value to avoid the direction of the change affecting the data.

Allele_frequencies$Change_stationary <- abs(Allele_frequencies$Change_stationary)
Allele_frequencies$Change_core <- abs(Allele_frequencies$Change_core)
Allele_frequencies$Change_edge <- abs(Allele_frequencies$Change_edge)

write.delim(Allele_frequencies, "Filtered_FRQs.delim")
#########
