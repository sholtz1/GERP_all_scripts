library(ape)
library(tidyverse)
library(caroline)
library(stringr)
library(grantham)
###### finding all annotated SNP variants available in our data


## load in all data to calculate gramtham load scores for each population
Full_filtered_data <- read.table("Full_data_bigtree.delim", header = TRUE)
Variant_effects <- read.table("GxK1C4chiGUhHJZZ.txt", header = TRUE, sep = "\t")

## filter for rows with an amino acid changes, other information is likely 
Variant_effects <- Variant_effects %>%
  filter(X..7 != "-") %>%
  filter(intron_variant == "missense_variant" | intron_variant == "stop_lost" | intron_variant == "coding_sequence_variant")

## split the amino acid columns into an original and variant amino acids.
Variant_effects[c('initial', 'changed')] <- str_split_fixed(Variant_effects$X..7, '/', 2)

## we only need thee amino acid columns and ID column for further analysis. 

Variant_effects <- Variant_effects %>%
  select(1, initial, changed)


## split the first column into chromosome and location in the chromosome.
Variant_effects[c('chromosome', 'Location', 'allele_change')] <- str_split_fixed(Variant_effects$LG10_54578_A.C, '_', 3)


## change to three letter codes for amino acids

Variant_effects <- Variant_effects %>%
  mutate(initial = case_when(
    initial == "A" ~ "Ala",
    initial == "R" ~ "Arg",
    initial == "N" ~ "Asn",
    initial == "D" ~ "Asp",
    initial == "C" ~ "Cys",
    initial == "E" ~ "Glu",
    initial == "Q" ~ "Gln",
    initial == "G" ~ "Gly",
    initial == "H" ~ "His",
    initial == "I" ~ "Ile",
    initial == "L" ~ "Leu",
    initial == "K" ~ "Lys",
    initial == "M" ~ "Met",
    initial == "F" ~ "Phe",
    initial == "P" ~ "Pro",
    initial == "S" ~ "Ser",
    initial == "T" ~ "Thr",
    initial == "W" ~ "Trp",
    initial == "Y" ~ "Try",
    initial == "V" ~ "Val",
  ))


Variant_effects <- Variant_effects %>%
  mutate(changed = case_when(
    changed == "A" ~ "Ala",
    changed == "R" ~ "Arg",
    changed == "N" ~ "Asn",
    changed == "D" ~ "Asp",
    changed == "C" ~ "Cys",
    changed == "E" ~ "Glu",
    changed == "Q" ~ "Gln",
    changed == "G" ~ "Gly",
    changed == "H" ~ "His",
    changed == "I" ~ "Ile",
    changed == "L" ~ "Leu",
    changed == "K" ~ "Lys",
    changed == "M" ~ "Met",
    changed == "F" ~ "Phe",
    changed == "P" ~ "Pro",
    changed == "S" ~ "Ser",
    changed == "T" ~ "Thr",
    changed == "W" ~ "Trp",
    changed == "Y" ~ "Tyr",
    changed == "V" ~ "Val",
  ))

## use the granthamn function to get the Grantham score between the initial and variatnt amino acids
Variant_effects <- Variant_effects %>%
  mutate(grantham = grantham_distance_exact(initial, changed))

##pull out just the Grantham score
Variant_effects$Grantham <- pull(Variant_effects$grantham[,3])
## when there is a loss of stop change in the code assign it a score of 250 as this should be extremely deleterious
Variant_effects <- Variant_effects %>%
  mutate(Grantham = case_when(is.na(Grantham) ~ 250,
         TRUE ~ Grantham))



## change the chromosome name back so that it can joing with the full allele FRQ data.
Variant_effects <- Variant_effects %>%
  mutate(Chromosome = case_when(chromosome == "LGX"~ "NC_007416.3" ,
                                chromosome ==  "LG2" ~ "NC_007417.3" ,
                                chromosome ==  "LG3" ~ "NC_007418.3" ,
                                chromosome ==  "LG4" ~ "NC_007419.2" ,
                                chromosome ==  "LG5" ~ "NC_007420.3" ,
                                chromosome ==  "LG6" ~ "NC_007421.3" ,
                                chromosome ==  "LG7" ~ "NC_007422.5" ,
                                chromosome ==  "LG8" ~ "NC_007423.3" ,
                                chromosome ==  "LG9" ~ "NC_007424.3" ,
                                chromosome ==  "LG10" ~ "NC_007425.3" ))

## shoose only the colums we nned for calculating a Grantham load

Variant_effect_simple <- Variant_effects %>%
  select(Chromosome, Location, Grantham)
##make numeric for joining

Variant_effect_simple$Location <- as.numeric(Variant_effect_simple$Location)


## addd the grantham score for each variant location we can in the large allele frequency change dataset. 
Grantham_alleles <- left_join(Full_filtered_data, Variant_effect_simple, by = c("Chromosome", "Location"))


## remove all loactions we dont have grantham score data for.
Grantham_alleles <- Grantham_alleles  %>%
  filter(!is.na(Grantham))



# Calculate a load score for all populations

Grantham_alleles <- Grantham_alleles %>%
  ##Take only informative sites
  ## Create a value for the frequency of the minor allele times the Grantham score.
  ## The higher the frequency and Grantham the higher this measure of genetic load
  #This is calculated for each pooled sample
  mutate(Gran_Start = ((1-`X0_NA`) * Grantham*Global), Gran_core = ((1-`X8_C`) * Grantham*Global), 
         Gran_edge = ((1-`X8_E`) * Grantham*Global), Gran_end = ((1-`X8_NA`) * Grantham* Global))


############################## test diffierent filtering############################


##################################################################################

## create a genetic load from allele frequencies combined with grantham scores. 
Grantham_loads <- Grantham_alleles %>%
  group_by(Landscape) %>%
  summarise(Start_load = (sum(Gran_Start, na.rm = TRUE)),
            Change_edge = (sum(Gran_edge, na.rm = TRUE)),
            Change_core = (sum(Gran_core, na.rm = TRUE)),
            `NA` = (sum(Gran_end, na.rm = TRUE))) %>%
  filter(Change_edge != 0) %>%
  select(Landscape, Change_edge, Change_core)




## change the names so it can be combined into the final summary. 

Grantham_loads <- Grantham_loads %>%
  pivot_longer(cols=c('Change_edge', 'Change_core'),
                           names_to='treatment',
                           values_to='Grantham_load') 
  

write.delim(Grantham_loads, "Grantham_loads.delim")




