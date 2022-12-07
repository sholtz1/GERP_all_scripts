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
  filter(intron_variant == "missense_variant" | MODIFIER == "HIGH")

## split the amino acid columns into an original and variant amino acids.
Variant_effects[c('initial', 'changed')] <- str_split_fixed(Variant_effects$X..7, '/', 2)

## we only need thee amino acid columns and ID column for further analysis. 

Variant_effects <- Variant_effects %>%
  select(1, initial, changed, MODIFIER)


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

### ideantify columns iwth high impact

Variant_effect_simple$Grantham <- as.numeric(Variant_effect_simple$Grantham)

Variant_effect_simple$Grantham <-replace_na(Variant_effect_simple$Grantham, "HIGH")



##make numeric for joining

Variant_effect_simple$Location <- as.numeric(Variant_effect_simple$Location)



## addd the grantham score for each variant location we can in the large allele frequency change dataset. 
Grantham_alleles <- left_join(Full_filtered_data, Variant_effect_simple, by = c("Chromosome", "Location"))


## Write the final fille with all data for making summaries
write.delim(Grantham_alleles, "Grantham_GERP_Full.delim")







