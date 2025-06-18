## Final data creation 
library(tidyverse)
## read in all allele FRQ data with GERP scores and grantham scores
Full_filtered_data <- read.table("Grantham_GERP_derived.delim", header = TRUE)




## Combine with fitness data. Needs to have the same names for combination.

## read in fitness data and change names
fitness <-read.table("Beetle_fitness_filtered.delim", header = TRUE)
colnames(fitness) <- c("Landscape", "treatment", "", "Fitness")
fitness <- fitness %>% 
  select(1,2,4) %>%
  mutate(treatment = case_when(treatment == "core" ~ "Change_core",
                               treatment == "edge" ~ "Change_edge"))

################## GERP genetic load without sites from out other loads ##################################

## create an estimate of genetic load from GERP score. higher GERP scores are
##preserved in the phylogeny and therefore change here should be deleterious.


GERP_load <- Full_filtered_data %>%
  ##Take sites that are predicted to be at least moderately deleterious
  filter(RS > 2) %>%
  ## remove sites that contribute to grantham and high effect load
  filter(IMPACT == "LOW" | is.na(IMPACT)) %>%
  ## Create a value for the frequency of the minor allele multiplies by the RS score.
  ## The higher the frequency and RS the higher this measure of genetic load
  #This is calculated for each pooled sample
  mutate(END_core = ((`X8_C`) * RS), 
         END_edge = ((`X8_E`) * RS))

## Modified  method from Valk et al 2019 making a genetic load by multiplying allele frequency 
## by gerp score to get an estimate of load.
pop_genetic_loads <- GERP_load %>%
  group_by(Landscape) %>%
  ##summing these assumes an additive effect of all of the deleterious change at GERP sites
  summarise(edge = (sum(END_edge, na.rm = TRUE)),
            core = (sum(END_core, na.rm = TRUE)))




## Make is so that all values are in the same column this allows for easier comparison
pop_genetic_loads_long <- pop_genetic_loads %>%
  pivot_longer(cols =c("edge", "core"),
               names_to = "Location" , values_to = "Load")

colnames(pop_genetic_loads_long)[2] <- "treatment"

## change names so that it can be joined into the final summary
pop_genetic_loads_long <- pop_genetic_loads_long %>%
  mutate(treatment = case_when(treatment == "edge" ~ "Change_edge",
                               treatment == "core" ~ "Change_core"))

######################################### GERP load with coding sites ##################################
GERP_load_full <- Full_filtered_data %>%
  ##Take only informative sites
  filter(RS > 2) %>%
  ## Create a value for the frequency of the minor allele times the RS score.
  ## The higher the frequency and RS the higher this measure of genetic load
  #This is calculated for each pooled sample
  mutate(END_core = ((`X8_C`) * RS), 
         END_edge = ((`X8_E`) * RS))

## Modified  method from Valk et al 2019 making a genetic load by multiplying allele frequency 
## by gerp score to get an estimate of load.
pop_genetic_loads_full <- GERP_load_full %>%
  group_by(Landscape) %>%
  ##summing these assumes an additive effect of all of the deleterious change at GERP sites
  summarise(edge = (sum(END_edge, na.rm = TRUE)),
            core = (sum(END_core, na.rm = TRUE)))




## Make is so that all values are in the same column this allows for easier comparison
pop_genetic_loads_long_full <- pop_genetic_loads_full %>%
  pivot_longer(cols =c("edge", "core"),
               names_to = "Location" , values_to = "Load_full")



colnames(pop_genetic_loads_long_full)[2] <- "treatment"

## change names so that it can be joind into the final summary
pop_genetic_loads_long_full <- pop_genetic_loads_long_full %>%
  mutate(treatment = case_when(treatment == "edge" ~ "Change_edge",
                               treatment == "core" ~ "Change_core"))


######################################### Format grantham score data ########################
## remove all locations we dont have grantham score data for.
Grantham_alleles <- Full_filtered_data  %>%
  filter(!is.na(Grantham))

### make numeric

Grantham_alleles$Grantham <- as.numeric(Grantham_alleles$Grantham)

# Calculate a load score for all populations

Grantham_alleles <- Grantham_alleles %>%
  ##Take only informative sites
  ## Create a value for the frequency of the derived allele multiplied by the Grantham score.
  ## The higher the frequency and Grantham the higher this measure of genetic load
  #This is calculated for each pooled sample
  mutate(Gran_core = ((`X8_C`) * Grantham), 
         Gran_edge = ((`X8_E`) * Grantham))




## create a genetic load from allele frequencies combined with grantham scores. 
##Summing assumes additive effect of all the deleterious change
Grantham_loads <- Grantham_alleles %>%
  group_by(Landscape) %>%
  summarise(
            Change_edge = (sum(Gran_edge, na.rm = TRUE)),
            Change_core = (sum(Gran_core, na.rm = TRUE)),
            ) %>%
  filter(Change_edge != 0) %>%
  select(Landscape, Change_edge, Change_core)




## change the names so it can be combined into the final summary. 

Grantham_loads <- Grantham_loads %>%
  pivot_longer(cols=c('Change_edge', 'Change_core'),
               names_to='treatment',
               values_to='Grantham_load') 


########################## High impact sites ###################################
## filter for sites identified that would break the protein and should have a high effect on fitness
High_effect_alleles <- Full_filtered_data  %>%
  filter(IMPACT == "HIGH")

## add up the frequency for the derived alleles with predicted high effects
High_effect_alleles <- High_effect_alleles %>%
  group_by(Landscape) %>%
  summarise(Change_edge = sum(`X8_E`, na.rm = TRUE), Change_core = sum(`X8_C`, na.rm = TRUE))

## change column names for joining
High_effect_alleles <- High_effect_alleles %>%
  pivot_longer(cols=c('Change_edge', 'Change_core'),
               names_to='treatment',
               values_to='High_effect_alleles') 
############# Homyzygous high impact ###################
## filter for sites identified that would break the protein and should have a high effect on fitness
High_effect_homo <- Full_filtered_data  %>%
  filter(IMPACT == "HIGH")

##square the frequency to estimate homozygosity occording to HWE
##add up the frequency for the derived alleles with predicted high effects
High_effect_homo <- High_effect_homo %>%
  group_by(Landscape) %>%
  summarise(Change_edge = sum((`X8_E`)^2, na.rm = TRUE), Change_core = sum((`X8_C`)^2, na.rm = TRUE))

##change names for joining
High_effect_homo <- High_effect_homo %>%
  pivot_longer(cols=c('Change_edge', 'Change_core'),
               names_to='treatment',
               values_to='High_effect_homo') 


######################################### GERP load squared to have only homozygous sites count ##################################
GERP_load_full <- Full_filtered_data %>%
  ##Take only informative sites
  filter(RS > 2) %>%
  ## Create a value for the frequency of the derived allele times squared times the RS score.
  ## The higher the frequency and RS the higher this measure of genetic load
  #This is calculated for each pooled sample
  #again we are estimating homozygosity according to HWE
  mutate(END_core = ((`X8_C`)^2 * RS ), 
         END_edge = ((`X8_E`)^2 * RS ))

## Modified  method from Valk et al 2019 making a genetic load by multiplying allele frequency 
## by gerp score to get an estimate of load.
pop_genetic_loads_homo <- GERP_load_full %>%
  group_by(Landscape) %>%
  ##summing these assumes an additive effect of all of the deleterious change at GERP sites
  summarise(
            edge = (sum(END_edge, na.rm = TRUE)),
            core = (sum(END_core, na.rm = TRUE)))




## Make is so that all values are in the same column this allows for easier comparison
pop_genetic_loads_homo <- pop_genetic_loads_homo %>%
  pivot_longer(cols =c("edge", "core"),
               names_to = "Location" , values_to = "Load_homozygous")
## filter for just core and edge loads we will use in our final analyses. 
pop_genetic_loads_homo <- pop_genetic_loads_homo %>%
  filter(Load_homozygous != 0)

colnames(pop_genetic_loads_homo)[2] <- "treatment"

## change names so that it can be joined into the final summary
pop_genetic_loads_homozygous <- pop_genetic_loads_homo %>%
  mutate(treatment = case_when(treatment == "edge" ~ "Change_edge",
                               treatment == "core" ~ "Change_core"))


######### GERP load squared to hace only homozygous sites count filtered for VEP sites ##################################
GERP_load_full <- Full_filtered_data %>%
  ##Take only informative sites
  filter(RS > 2) %>%
  filter(IMPACT == "LOW" | is.na(IMPACT)) %>%
  ## remove sites identified in our coding regions
  ## Create a value for the frequency of the minor allele times the RS score.
  ## The higher the frequency and RS the higher this measure of genetic load
  #This is calculated for each pooled sample
  ## square to get an estimate for proportion of homozygous sites according to hardy weinberg
  mutate( END_core = ((`X8_C`)^2 * RS), 
         END_edge = ((`X8_E`)^2 * RS ))

## Modified  method from Valk et al 2019 making a genetic load by multiplying allele frequency 
## by gerp score to get an estimate of load.
pop_genetic_loads_homo <- GERP_load_full %>%
  group_by(Landscape) %>%
  ##summing these assumes an additive effect of all of the deleterious change at GERP sites
  summarise(
            edge = (sum(END_edge, na.rm = TRUE)),
            core = (sum(END_core, na.rm = TRUE)))




## Make is so that all values are in the same column this allows for easier comparison
pop_genetic_loads_homo <- pop_genetic_loads_homo %>%
  pivot_longer(cols =c("edge", "core"),
               names_to = "Location" , values_to = "Load_homo_filtered")
## filter for just core and edge loads we will use in our final analyses. 
pop_genetic_loads_homo <- pop_genetic_loads_homo %>%
  filter(Load_homo_filtered != 0)


## change names so that it can be joined into the final summary
colnames(pop_genetic_loads_homo)[2] <- "treatment"


pop_genetic_loads_homozygous_filtered <- pop_genetic_loads_homo %>%
  mutate(treatment = case_when(treatment == "edge" ~ "Change_edge",
                               treatment == "core" ~ "Change_core"))


############################## Format grantham score data for homozygous sites only ########################
## remove all locations we dont have grantham score data for.
Grantham_alleles <- Full_filtered_data  %>%
  filter(!is.na(Grantham))

### make numeric

Grantham_alleles$Grantham <- as.numeric(Grantham_alleles$Grantham)

# Calculate a load score for all populations

Grantham_homo <- Grantham_alleles %>%
  ##Take only informative sites
  ## Create a value for the frequency of the minor allele times the Grantham score.
  ## The higher the frequency and Grantham the higher this measure of genetic load
  #This is calculated for each pooled sample
  ## square to get an estimate for proportion of homozygous sites according to hardy weinberg
  mutate( Gran_core = ((`X8_C`)^2 * Grantham), 
         Gran_edge = ((`X8_E`)^2 * Grantham))




## create a genetic load from allele frequencies combined with grantham scores. 
##Summing assumes additive effect of all the deleterious change
Grantham_homo <- Grantham_homo %>%
  group_by(Landscape) %>%
  summarise(
            Change_edge = (sum(Gran_edge, na.rm = TRUE)),
            Change_core = (sum(Gran_core, na.rm = TRUE))) %>%
  filter(Change_edge != 0) %>%
  select(Landscape, Change_edge, Change_core)




## change the names so it can be combined into the final summary. 

Grantham_homo <- Grantham_homo %>%
  pivot_longer(cols=c('Change_edge', 'Change_core'),
               names_to='treatment',
               values_to='Grantham_homo') 



###add in diversity values pi

Landscape_diversity <- read.table("final_diversity.delim", header = TRUE)

## think about the number of deleterious GERP sites at high frequencies
HF_GERP <- Full_filtered_data %>%
  #get core and edge in a single colum 
  pivot_longer(cols =c("X8_C", "X8_E"),
               names_to = "E_C" , values_to = "FRQ") %>%
  ##get the total number of deleterious sites identified by gerp that have over a .5 frequency 
  ##for the dirived allele these sites likely have individuals homozygous for the derived allele and 
  ##therefore will be expressed if they are recessive
  filter(FRQ > .5) %>%
  filter(RS > 2) %>%
  group_by(Landscape, E_C) %>%
  summarise(tot_fixed = n())

##change names for joining
colnames(HF_GERP)<-c("Landscape", "treatment", "HF_GERP")

HF_GERP  <- HF_GERP  %>%
  mutate(treatment = case_when(treatment == "X8_E" ~ "Change_edge",
                               treatment == "X8_C" ~ "Change_core"))

## think about the number of deleterious VEP sites at high frequencies
HF_VEP <- Full_filtered_data %>%
  ##get the total number of deleterious sites identified by gerp that have over a .5 frequency 
  ##for the dirived allele these sites likely have individuals homozygous for the derived allele and 
  ##therefore will be expressed if they are recessive
  pivot_longer(cols =c("X8_C", "X8_E"),
               names_to = "E_C" , values_to = "FRQ") %>%
  filter(FRQ > .5) %>%
  filter(IMPACT %in% c("MODERATE","HIGH")) %>%
  group_by(Landscape, E_C) %>%
  summarise(tot_fixed = n())

colnames(HF_VEP)<-c("Landscape", "treatment", "HF_VEP")

HF_VEP  <- HF_VEP  %>%
  mutate(treatment = case_when(treatment == "X8_E" ~ "Change_edge",
                               treatment == "X8_C" ~ "Change_core"))




############# Join data #######################

## this joins the summaries generated above one by one into a final summary dataframe with everything in it

Final_summary_2 <- left_join(pop_genetic_loads_long, Grantham_loads, by = c("Landscape", "treatment"))

Final_summary_2 <- left_join(Final_summary_2, fitness, by = c("Landscape", "treatment"))

Final_summary_2 <- left_join(Final_summary_2, pop_genetic_loads_long_full, by = c("Landscape", "treatment"))

Final_summary_2 <- left_join(Final_summary_2, High_effect_alleles, by = c("Landscape", "treatment"))

Final_summary_2 <- left_join(Final_summary_2, High_effect_homo, by = c("Landscape", "treatment"))

Final_summary_2 <- left_join(Final_summary_2, pop_genetic_loads_homozygous, by = c("Landscape", "treatment"))

Final_summary_2 <- left_join(Final_summary_2, Grantham_homo, by = c("Landscape", "treatment"))

Final_summary_2 <- left_join(Final_summary_2, pop_genetic_loads_homozygous_filtered, by = c("Landscape", "treatment"))

Final_summary_2 <- left_join(Final_summary_2, Landscape_diversity, by = c("Landscape", "treatment"))

Final_summary_2 <- left_join(Final_summary_2, HF_GERP, by = c("Landscape", "treatment"))

Final_summary_2 <- left_join(Final_summary_2, HF_VEP, by = c("Landscape", "treatment"))

write_delim(Final_summary_2 , "Final_summary_derived.delim")





