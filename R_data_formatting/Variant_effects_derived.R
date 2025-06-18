library(ape)
library(tidyverse)
library(caroline)
library(stringr)
library(grantham)
library(tidyverse)

## load in all data to combine allele frequencies and variant effects for each population
Full_data <- read.table("Full_data_newgerp.delim", header = TRUE)
Variant_effects <- read.table("Variant_effects.txt", header = TRUE, sep = "\t")

##filter for important informant in the allele frq file. We are using ancestral frqs
Filtered_data <- Full_data %>%
  select(Chromosome, Location, Landscape, Core_Edge, FRQ_ancestral,MAX_RS, RS) %>%
  filter(!is.na(Core_Edge))

##change to derived/deleterious allele frequency 
Filtered_data <- Filtered_data %>%
  mutate(FRQ_ancestral = 1-FRQ_ancestral)



## change so core and edge are in the same row
Filtered_data <- Filtered_data %>% pivot_wider(names_from = Core_Edge, values_from = FRQ_ancestral)
colnames(Filtered_data) <- c("Chromosome", "Location", "Landscape", "RS", "MAX_RS", "X8_C", "X8_E")



## filter for rows with important genetic changes
## this is also a step that could be thought about specifically if we want to include any of the
##modifier changes now that we think the more deleterious changes aren't whats causing our fitness differences.
Variant_effects <- Variant_effects %>%
  filter(Consequence == "missense_variant" | IMPACT == "HIGH" | IMPACT == "LOW")



## split the amino acid columns into an original and variant amino acids.
Variant_effects[c('initial', 'changed')] <- str_split_fixed(Variant_effects$Amino_acids, '/', 2)

## we only need thee amino acid columns and ID column for further analysis. 

Variant_effects <- Variant_effects %>%
  select(1, initial, changed, IMPACT)


## split the first column into chromosome and location in the chromosome.
Variant_effects[c('chromosome', 'Location', 'allele_change')] <- str_split_fixed(Variant_effects$Uploaded_variation, '_', 3)


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



## choose only the column we need for calculating a Grantham load

Variant_effect_simple <- Variant_effects %>%
  select(Chromosome, Location, Grantham, IMPACT)



##Make numeric

Variant_effect_simple$Grantham <- as.numeric(Variant_effect_simple$Grantham)


Variant_effect_simple$Location <- as.numeric(Variant_effect_simple$Location)

## use unique to avoid duplicates

Variant_effect_simple <- unique(Variant_effect_simple)



## add the grantham score for each variant location we can in the large allele frequency change dataset. 
Grantham_alleles <- left_join(Filtered_data, Variant_effect_simple, by = c("Chromosome", "Location"))



## Write the final file with all data for making summaries
write.delim(Grantham_alleles, "Grantham_GERP_derived.delim")







