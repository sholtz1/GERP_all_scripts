##Attach/detach packages, ggpubr can mess up tidyverse grouping
detach("package:ggpubr", unload=TRUE)
library(tidyverse)


# TWL: I added several comments throughout, but an overall comment for this script is that
# it still needs to be cleaned up. There are multiple unecessary comments, old parts of code,
# or mismatches (referring to a minor allele in the comment when instead the code is using
# the major allele). The goal for the code review was to have the code close to the version
# that would be published with the paper and this script needs a good amount of revision
# before it's at that point.


##Load data
# TWL: Be sure this is updated with the actual file (not the subset for me) in the
# final version.
Allele_frequencies <- read.table("Allele_frq_subset_Topher", header = TRUE, sep = " ")
## right now we arent keeping the minor allele constant between generations. 
## ID minor allele in gen 0 and the be sure to keep that allele as the minor allele for the other generations. 

## Get the min and max numbers for the base pairs into new columns. 
Base_pairs <-Allele_frequencies %>%
  select(8:13)

##This gets the total number or base pairs present at a specific location.
Allele_frequencies$total <- rowSums(Base_pairs)

##makes a colum to tell how many alleles are present for that patch
# TWL: This comment is misleading because the result isn't "how many alleles are present"
# Instead, it's how many alleles are not present since you're checking for the
# frequencies being 0. I'd suggest changing this to be x != 0 so it matches the comment
# which is more intuitive. Then you'll need to update the filter step below, but it
# will again be more intuitive because it will be filtering for counts to be either 1 or 2.
# Also, I notice the 0 here is a character. I'd suggest changing the data type before
# doing any subsequent calculations to avoid possible issues. That will negate the
# need for the as.numeric steps later on in the code too.
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
# TWL: With my earlier suggestion, this will be count %in% c(1,2) or count <= 2
Allele_frequencies <- Allele_frequencies %>%
  filter(count %in% c(4, 5))


##Make numberic and
# TWL: If you follow my earlier suggestion, you won't need these steps.
Allele_frequencies$total <- as.numeric(Allele_frequencies$total)
Allele_frequencies$alternate <- as.numeric(Allele_frequencies$alternate)

##Remove columns that only have 1 allele present.
##Are we removing cases when alleles were fixed in one generation but polymorphic in another. 

## This code might need to go after IDing the minor allele in case a new allele 
##is completely fixed

##This makes a dataframe of every location that is fixed within a landscape
## We want to remove these from the analysis since fixed sites dont tell us anything
One_allele <- Allele_frequencies %>%
  group_by(Chromosome, Location, Landscape) %>%
  summarize(Allele_diff = sum(total - alternate)) %>%
  filter(Allele_diff == 0)

## Add code to make a DF of locations that are fixed at with two different bases

##Start with only alleles that are fixed somewhere
Fixed_only <- Allele_frequencies %>%
  filter(total == alternate)

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

# TWL: The above process works, but it seems a bit cumbersome. Couldn't it be simplified
# by only creating the Fixed_only data frame and then in the last step filtering for
# diff_BP == 1? Then, you could just filter above for !(paste(Chromosome, Location, Landscape) %in% paste(Fixed_only_sum$Chromosome, ...))
# and it would accomplish the same thing in less steps, right?



## We need to get the alternate allele neucleotide to keep it constant for each
##Landscape


##create a new dataframe with the total number of each base pair for each landscape/location 

##Frequencies to be < .5 but this might no be a problem as the change can still be higher
##Could we use a combination of the two methods below? Of would it help
##to use the major allele frequency that way it can be identified even when the 
#minor is 0 this might work better

##Instead of summing we summarize by the first row for each unique location and landscape
##Check and make sure its always gen 0

###FRQ_test <- Allele_frequencies %>%
#  mutate(combo = paste(Chromosome, Location, Landscape)) %>%
#distinct(combo, .keep_all= TRUE) %>%
#  select(-combo)

## New method summarizes without landscape so we are getting the minor allele
## Frequency based on all imformation for that location across landscapes.

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

# TWL: Why is the below step necessary? the max function should by default
# return numeric values and if you're getting something else, it would be worth
# it to understand why as it could relate to a different issue.
FRQ_test$Major <-as.numeric((FRQ_test$Major))

## Subtract major from total to get the minor allele frequency
#FRQ_test <- FRQ_test %>%
#  mutate(Minor = total - Major)

##Creat a new DF that only has nculotide information.
Base_only <- FRQ_test %>%
  ungroup() %>%
  select(3:7)

Base_only[is.na(Base_only)] <- 0

## Create a new row to put Major allele BP letter into
FRQ_test$New <- rep(NA, nrow(FRQ_test))


## Get the minor allele BP into a new column

for (i in 1:nrow(FRQ_test)) {
  FRQ_test$New[i] <- names(Base_only)[which(Base_only[i,] == FRQ_test$Major[i])]
}


## Get the global allele frequency of the major allele.
## Higher global major allele frqency could mean changes away from it are
## Potentially more deleterious.
FRQ_test_2 <- FRQ_test %>%
  group_by(Chromosome, Location) %>%
  summarise(sum = sum(A,`T`,C,G))

FRQ_test$Global_tot <- FRQ_test_2$sum

FRQ_test <- FRQ_test %>%
  mutate(Global = Major/Global_tot)
## Take out colums that will be duplicates when joining. 
FRQ_test <- FRQ_test %>%
  select(Chromosome, Location, New, Major, Global)





## Join to add the minor allele frequency identifier to each landscape/location
Allele_frequencies <- left_join(Allele_frequencies, FRQ_test, by = c("Chromosome", "Location"))

##Create new column
Allele_frequencies$primary <- rep(NA, nrow(Allele_frequencies))

##Get the primary allele using the New column identifier

for (i in 1:nrow(Allele_frequencies)) {
  Allele_frequencies$primary[i] <- Allele_frequencies[i, which(names(Allele_frequencies) == Allele_frequencies$New[i])]

  
}


##By dividing the number of minor neucleotides by the total number present at that location we get the allele frequency of the minor allele.
###mutate new column of frequencies.
##this value should be able to be over .5
## What does a value of 1 mean?
## Double check to make sure this isnt grouped for the calculation.
Allele_frequencies <- Allele_frequencies %>%
  mutate(FRQ = primary/total)





##Remove unnecessary columns for analysis
Allele_frequencies <- Allele_frequencies %>%
  select(1,2, 3, 4,5 , 6, Major, Global, primary,  FRQ)



## Need to writ a DF to keep individual landscape/generation data before pivoting 
## Since that put the full landscape into one row

##write.delim(Allele_frequencies, "Filtered_FRQs_individual.delim")

##Remove data that prevents pivioting
Allele_frequencies <- Allele_frequencies %>%
  select(1,2, 4,5 , 6,   FRQ, Major, Global)

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

##Think about conserving directionality

##Change so that negative changes are a reduction in the major allele
##Positive change represents an increase in the major allele

# TWL: An easier way to achieve this directionality would be to just reverse the 
# order of the subtraction in the mutate steps above. Then you wouldn't need the code
# below.

Allele_frequencies$Change_stationary <- -1*(Allele_frequencies$Change_stationary)
Allele_frequencies$Change_core <- -1*(Allele_frequencies$Change_core)
Allele_frequencies$Change_edge <- -1*(Allele_frequencies$Change_edge)

##write.delim(Allele_frequencies, "Filtered_FRQs.delim")
#########
