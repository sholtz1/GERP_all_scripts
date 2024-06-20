## Final data creation 
library(tidyverse)
## read in all allele FRQ data with GERP scores and grantham scores
Full_filtered_data <- read.table("Grantham_GERP_Full.delim", header = TRUE)

## create colums with changes in frequency from generation 0 to generation 8 for the edge and core
Full_filtered_data <- Full_filtered_data %>%
  mutate(Change_core = `X8_C`- `X0_NA`) %>%
  mutate(Change_edge = `X8_E`- `X0_NA`)



###################Calculate 0 fold change squared for each landscape and treatment#########################

allele_fold <- Full_filtered_data %>%
  filter(!is.na(Fold)) %>%
  filter(is.na(Change_stationary))

##calculate 0 fold change from global allele frequency as that should be what stabilizing selection
## is keeping the frequency near


allele_fold <- allele_fold %>%
  mutate(Change_core = `X8_C`- Global) %>%
  mutate(Change_edge = `X8_E`- Global)

allele_fold$Fold <- as.factor(allele_fold$Fold)


###make dataframe long so that we can compare core and edge treatments
allele_fold_long <- allele_fold %>%
  pivot_longer(cols = c("Change_stationary", "Change_core", "Change_edge"), 
               names_to = "treatment" , values_to = "frq_change") %>%
  mutate(Fold_treat = paste(Fold, treatment, sep = "_"))


## Combine with fitness data. Needs to have the same names for combination.

## read in fitness data and chnage names
fitness <-read.table("Beetle_fitness_filtered.delim", header = TRUE)
colnames(fitness) <- c("Landscape", "treatment", "", "Fitness")
fitness <- fitness %>% 
  select(1,2,4) %>%
  mutate(treatment = case_when(treatment == "core" ~ "Change_core",
                               treatment == "edge" ~ "Change_edge"))


## put fitness data onto the allele FRQ dataset
allele_fold_long_fitness <- right_join(allele_fold_long, fitness,
                                       by = c("Landscape", "treatment"))

## change names for joining
allele_fold_long_fitness <- allele_fold_long_fitness %>%
  mutate(final_frq = case_when(treatment == "Change_core" ~ `X8_C`,
                               treatment == "Change_edge" ~ `X8_E`))


## square allele frequency assuming stablizing selection so any change should be deleterious
allele_fold_change_squared <- allele_fold_long_fitness %>%
  mutate(change_squared = frq_change^2) %>%
  filter(Fold == 0) %>%
  group_by(treatment, Landscape) %>%
  ## get the average allele frequency change for edge and core for each landscape
  summarize(change_0fold_squared = mean(change_squared , na.rm =TRUE), Fitness = mean(Fitness))




##Look at major allele frq for 0 fold sites
## this assumes that 0 fold sites are under tabilizing selection 
#and higher major allele frq should be less deleterious

allele_fold_major_frq <- allele_fold_long_fitness %>%
  filter(Fold == 0) %>%
  group_by(treatment, Landscape) %>%
  ##
  summarize(major_0fold_frq = mean(final_frq , na.rm =TRUE), Fitness = mean(Fitness))

## Just look at non squared change. assumes purifying selction so negative chagne is deleterious
## and positive change should be beneficial.
allele_fold_change_directional <- allele_fold_long_fitness %>%
  filter(Fold == 0) %>%
  group_by(treatment, Landscape) %>%
  ## using mean frq change give an estimate for how much negative vs beneficial change there could be
  summarize(directional_0fold = mean(frq_change , na.rm =TRUE), Fitness = mean(Fitness))


## get allele frequency change for all sites including neutral sites 

Full_data_long <- Full_filtered_data %>%
  pivot_longer(cols = c("Change_stationary", "Change_core", "Change_edge"), 
               names_to = "treatment" , values_to = "frq_change") %>%
  filter(!is.na(frq_change)) %>%
  filter(treatment != "Change_stationary")

## square the change so that negative and positive change dont just cancel each other out
Full_data_squared <- Full_data_long %>%
  mutate(change_squared = frq_change^2)
  
## summarize by landscape and treatment to get final summary stats
Full_data_squared_sum <- Full_data_squared %>%
  group_by(Landscape, treatment) %>%
  summarise(change_squared = (mean(change_squared, na.rm = TRUE)))



## look absolute allele frequency chagnge away from the global total
Full_data_global <- Full_data_long %>%
  mutate(global_change = case_when(
    treatment == "Change_core" ~ abs(Global - X8_C),
    treatment == "Change_edge" ~ abs(Global - X8_E),
  ))

## summarize by landscape and treatment to get final summary stats
Full_data_global_sum <- Full_data_global %>%
  group_by(Landscape, treatment) %>%
  summarise(change_global = (mean(global_change, na.rm = TRUE)))





##### look at grantham scores assuming stabablizing selection

Grantham_global <- Full_data_long %>%
  filter(!is.na(Grantham)) %>%
  mutate(Gran_global_change = case_when(
    treatment == "Change_core" ~ abs(Global - X8_C)*Grantham,
    treatment == "Change_edge" ~ abs(Global - X8_E)*Grantham,
  ))


## summarize by landscape and treatment to get final summary stats
Grantham_global_sum <- Grantham_global %>%
  group_by(Landscape, treatment) %>%
  summarise(Grantham_global = (mean(Gran_global_change, na.rm = TRUE)))


################## GERP genetic load without coding sites ########################################

## create an estimate of genetic load from GERP score. higher GERP scores are
##preserved in the phylogeny and therefore change here should be deleterious.


GERP_load <- Full_filtered_data %>%
  ##Take only informative sites
  filter(RS > 2) %>%
  ## remove sites identified in our coding regions
  filter(is.na(Grantham)) %>%
  ## Create a value for the frequency of the minor allele times the RS score.
  ## The higher the frequency and RS the higher this measure of genetic load
  #This is calculated for each pooled sample
  mutate(Start = ((1-`X0_NA`) * RS), END_core = ((1-`X8_C`) * RS), 
         END_edge = ((1-`X8_E`) * RS), END_Stat = ((1-`X8_NA`) * RS))

## Modified  method from Valk et al 2019 making a genetic load by multiplying allele frequency 
## by gerp score to get an estimate of load.
pop_genetic_loads <- GERP_load %>%
  group_by(Landscape) %>%
  ##summing these assumes an additive effect of all of the deleterious change at GERP sites
  summarise(Start_load = (sum(Start, na.rm = TRUE)),
            edge = (sum(END_edge, na.rm = TRUE)),
            core = (sum(END_core, na.rm = TRUE)),
            `NA` = (sum(END_Stat, na.rm = TRUE)))




## Make is so that all values are in the same column this allows for easier comparison
pop_genetic_loads_long <- pop_genetic_loads %>%
  pivot_longer(cols =c("Start_load", "edge", "core", "NA"),
               names_to = "Location" , values_to = "Load")
## filter for just core and edge loads we will use in our final analyses. 
pop_genetic_loads_long <- pop_genetic_loads_long %>%
  filter(Load != 0) %>%
  filter(Location != "Start_load" & Location != "NA")

colnames(pop_genetic_loads_long)[2] <- "treatment"

## change names so that it can be joind into the final summary
pop_genetic_loads_long <- pop_genetic_loads_long %>%
  mutate(treatment = case_when(treatment == "edge" ~ "Change_edge",
                               treatment == "core" ~ "Change_core"))

######################################### GERP load with coding sites ##################################
GERP_load_full <- Full_filtered_data %>%
  ##Take only informative sites
  filter(RS > 2) %>%
  ## remove sites identified in our coding regions
  ## Create a value for the frequency of the minor allele times the RS score.
  ## The higher the frequency and RS the higher this measure of genetic load
  #This is calculated for each pooled sample
  mutate(Start = ((1-`X0_NA`) * RS), END_core = ((1-`X8_C`) * RS), 
         END_edge = ((1-`X8_E`) * RS), END_Stat = ((1-`X8_NA`) * RS))

## Modified  method from Valk et al 2019 making a genetic load by multiplying allele frequency 
## by gerp score to get an estimate of load.
pop_genetic_loads_full <- GERP_load_full %>%
  group_by(Landscape) %>%
  ##summing these assumes an additive effect of all of the deleterious change at GERP sites
  summarise(Start_load = (sum(Start, na.rm = TRUE)),
            edge = (sum(END_edge, na.rm = TRUE)),
            core = (sum(END_core, na.rm = TRUE)),
            `NA` = (sum(END_Stat, na.rm = TRUE)))




## Make is so that all values are in the same column this allows for easier comparison
pop_genetic_loads_long_full <- pop_genetic_loads_full %>%
  pivot_longer(cols =c("Start_load", "edge", "core", "NA"),
               names_to = "Location" , values_to = "Load_full")
## filter for just core and edge loads we will use in our final analyses. 
pop_genetic_loads_long_full <- pop_genetic_loads_long_full %>%
  filter(Load_full != 0) %>%
  filter(Location != "Start_load" & Location != "NA")

colnames(pop_genetic_loads_long_full)[2] <- "treatment"

## change names so that it can be joind into the final summary
pop_genetic_loads_long_full <- pop_genetic_loads_long_full %>%
  mutate(treatment = case_when(treatment == "edge" ~ "Change_edge",
                               treatment == "core" ~ "Change_core"))

######################################### GERP load with coding sites filtered for meduim RS ##################################
GERP_load_full_filtered <- Full_filtered_data %>%
  ##Take only informative sites
  filter(RS > 2 & RS <6) %>%
  ## remove sites identified in our coding regions
  ## Create a value for the frequency of the minor allele times the RS score.
  ## The higher the frequency and RS the higher this measure of genetic load
  #This is calculated for each pooled sample
  mutate(Start = ((1-`X0_NA`) * RS), END_core = ((1-`X8_C`) * RS), 
         END_edge = ((1-`X8_E`) * RS), END_Stat = ((1-`X8_NA`) * RS))

## Modified  method from Valk et al 2019 making a genetic load by multiplying allele frequency 
## by gerp score to get an estimate of load.
pop_genetic_loads_full_filtered <- GERP_load_full_filtered %>%
  group_by(Landscape) %>%
  ##summing these assumes an additive effect of all of the deleterious change at GERP sites
  summarise(Start_load = (sum(Start, na.rm = TRUE)),
            edge = (sum(END_edge, na.rm = TRUE)),
            core = (sum(END_core, na.rm = TRUE)),
            `NA` = (sum(END_Stat, na.rm = TRUE)))




## Make is so that all values are in the same column this allows for easier comparison
pop_genetic_loads_long_full_filtered <- pop_genetic_loads_full_filtered %>%
  pivot_longer(cols =c("Start_load", "edge", "core", "NA"),
               names_to = "Location" , values_to = "Load_full_filtered")
## filter for just core and edge loads we will use in our final analyses. 
pop_genetic_loads_long_full_filtered <- pop_genetic_loads_long_full_filtered %>%
  filter(Load_full_filtered != 0) %>%
  filter(Location != "Start_load" & Location != "NA")

colnames(pop_genetic_loads_long_full_filtered)[2] <- "treatment"

## change names so that it can be joind into the final summary
pop_genetic_loads_long_full_filtered <- pop_genetic_loads_long_full_filtered %>%
  mutate(treatment = case_when(treatment == "edge" ~ "Change_edge",
                               treatment == "core" ~ "Change_core"))


######################################### Format grantam score data ########################
## remove all loactions we dont have grantham score data for.
Grantham_alleles <- Full_filtered_data  %>%
  filter(!is.na(Grantham))

### make numeric

Grantham_alleles$Grantham <- as.numeric(Grantham_alleles$Grantham)

# Calculate a load score for all populations

Grantham_alleles <- Grantham_alleles %>%
  ##Take only informative sites
  ## Create a value for the frequency of the minor allele times the Grantham score.
  ## The higher the frequency and Grantham the higher this measure of genetic load
  #This is calculated for each pooled sample
  mutate(Gran_Start = ((1-`X0_NA`) * Grantham), Gran_core = ((1-`X8_C`) * Grantham), 
         Gran_edge = ((1-`X8_E`) * Grantham), Gran_end = ((1-`X8_NA`) * Grantham))


### test diffierent filtering?



## create a genetic load from allele frequencies combined with grantham scores. 
##Summing assumes additive effect of all the deleterious change
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
##################grantham averages #############################################

Grantham_averages <- Grantham_alleles %>%
  group_by(Landscape) %>%
  summarise(Start_load = (mean(Gran_Start, na.rm = TRUE)),
            Change_edge = (mean(Gran_edge, na.rm = TRUE)),
            Change_core = (mean(Gran_core, na.rm = TRUE)),
            `NA` = (mean(Gran_end, na.rm = TRUE))) %>%
  filter(Change_edge != 0) %>%
  select(Landscape, Change_edge, Change_core)




## change the names so it can be combined into the final summary. 

Grantham_averages <- Grantham_averages %>%
  pivot_longer(cols=c('Change_edge', 'Change_core'),
               names_to='treatment',
               values_to='Grantham_average') 



########################## High impact sites ###################################

High_effect_alleles <- Full_filtered_data  %>%
  filter(!is.na(Grantham)) %>%
  filter(IMPACT == "HIGH")

High_effect_alleles <- High_effect_alleles %>%
  group_by(Landscape) %>%
  summarise(Change_edge = sum(`X8_E`, na.rm = TRUE), Change_core = sum(`X8_C`, na.rm = TRUE)) %>%
  filter(Change_edge != 0)

High_effect_alleles <- High_effect_alleles %>%
  pivot_longer(cols=c('Change_edge', 'Change_core'),
               names_to='treatment',
               values_to='High_effect_alleles') 
############# Average high impact, low frq sites may not be imporant to fitness ###################
High_effect_averages <- Full_filtered_data  %>%
  filter(!is.na(Grantham)) %>%
  filter(IMPACT == "HIGH")

High_effect_averages <- High_effect_averages %>%
  group_by(Landscape) %>%
  summarise(Change_edge = mean(`X8_E`, na.rm = TRUE), Change_core = mean(`X8_C`, na.rm = TRUE)) %>%
  filter(Change_edge != 0)

High_effect_averages <- High_effect_averages %>%
  pivot_longer(cols=c('Change_edge', 'Change_core'),
               names_to='treatment',
               values_to='High_effect_averages') 



############# Join data #######################

## this joins the summaries generated above one by one into a final summary dataframe with everything in it

Final_summary_2 <- left_join(pop_genetic_loads_long, Full_data_squared_sum, by = c("Landscape", "treatment"))

Final_summary_2 <- left_join(Final_summary_2, allele_fold_change_squared, by = c("Landscape", "treatment"))

allele_fold_major_frq <- allele_fold_major_frq %>%
  select(-Fitness)

Final_summary_2 <- left_join(Final_summary_2, allele_fold_major_frq, by = c("Landscape", "treatment")) %>%
  filter(!is.na(Fitness))

allele_fold_change_directional <- allele_fold_change_directional %>%
  select(-Fitness)

Final_summary_2 <- left_join(Final_summary_2, allele_fold_change_directional, by = c("Landscape", "treatment")) %>%
  filter(!is.na(Fitness))

Final_summary_2 <- left_join(Final_summary_2, Grantham_loads, by = c("Landscape", "treatment"))

Final_summary_2 <- left_join(Final_summary_2, pop_genetic_loads_long_full, by = c("Landscape", "treatment"))

Final_summary_2 <- left_join(Final_summary_2, High_effect_alleles, by = c("Landscape", "treatment"))

Final_summary_2 <- left_join(Final_summary_2, High_effect_averages, by = c("Landscape", "treatment"))

Final_summary_2 <- left_join(Final_summary_2, Grantham_averages, by = c("Landscape", "treatment"))

Final_summary_2 <- left_join(Final_summary_2, Full_data_global_sum, by = c("Landscape", "treatment"))

Final_summary_2 <- left_join(Final_summary_2, Grantham_global_sum, by = c("Landscape", "treatment"))

Final_summary_2 <- left_join(Final_summary_2, pop_genetic_loads_long_full_filtered, by = c("Landscape", "treatment"))

write_delim(Final_summary_2 , "Final_summary.delim")





