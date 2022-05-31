## Calculate a reletive genetic load for the population based on the GERP scores
library(tidyverse)

## Read in frequency data where we have freqency for the final population
test_single_pop <- read.table("Filtered_FRQs_1chromosom.delim", header = TRUE)


## lets find the relative load for a single chromosome
## load in single chromosome subset
NC_007416.3_rm.rates.bed <-read.delim("NC_007416.3_rm.rates.bed", header = FALSE)
colnames(NC_007416.3_rm.rates.bed) <- c("Chromosome","not_important", "Location", "Max_RS", "RS")
 



## Add the RS values to the allele frequency file
## must be done on the frequency file that stil had the origonal freuencies not
## just the change
test_single_pop <- left_join(test_single_pop, NC_007416.3_rm.rates.bed, by = c("Chromosome", "Location"))

#rm(NC_007416.3_rm.rates.bed)

## Vectorized way to calculate a genetic load for ach population/generation

test_single_pop <- test_single_pop %>%
  ##Take only informative sites
  filter(Max_RS > 0) %>%
  ## Create a value for the frequency of the minor allele times the RS score.
  ## The higher the frequency and RS the higher this measure of genetic load
  #This is calculated for each pooled sample
  mutate(Start = ((1-`X0_NA`) * RS), END_core = ((1-`X8_C`) * RS), 
         END_edge = ((1-`X8_E`) * RS), END_Stat = ((1-`X8_NA`) * RS))

## This is the method from Valk et al 2019 making the GERP score relative to
##The total frequency. What exactly is this doing and why?
pop_genetic_loads <- test_single_pop %>%
  group_by(Landscape) %>%
  summarise(Start_load = (sum(Start, na.rm = TRUE)/ sum(`X0_NA`, na.rm = TRUE)),
            edge = (sum(END_edge, na.rm = TRUE)/ sum(`X8_E`, na.rm = TRUE)),
            core = (sum(END_core, na.rm = TRUE)/ sum(`X8_C`, na.rm = TRUE)),
            `NA` = (sum(END_Stat, na.rm = TRUE)/ sum(`X8_NA`, na.rm = TRUE))
            )




## Make is so that all values are in the same column this allows for easier comparison
pop_genetic_loads_long <- pop_genetic_loads %>%
  pivot_longer(cols =c("Start_load", "edge", "core", "NA"),
                       names_to = "Location" , values_to = "Load")

## Add in fitness to compare fitness for each replicate to the calculated genetic load
Loads_VS_Fitness <- left_join(pop_genetic_loads_long , beetle_fitness, by = c("Landscape", "Location")) %>%
  filter(!is.na(Fitness))

## SImple plot of genetic load vs fitness. Seems to be some correlation?
Loads_VS_Fitness %>%
  ggplot(aes(x = Load, y = Fitness)) +
  geom_point()+
  geom_smooth(method = "lm")


library(ggpubr)
pop_genetic_loads_long %>%
  filter(generation == "END_edge_load" |
           generation == "END_core_load" ) %>%
  ggplot(aes(x=generation, y= Load )) +
  geom_boxplot() +
  theme_classic() +
  stat_compare_means(label.x = 1.5, label.y = 1.3)+
  stat_compare_means(method = "t.test",label.x = 2, label.y = 1.3 )

detach(package:ggpubr,unload=TRUE)
