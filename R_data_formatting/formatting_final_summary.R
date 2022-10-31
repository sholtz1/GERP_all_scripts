## Final data creation 


Full_filtered_data <- read.table("Full_data_bigtree.delim", header = TRUE)


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


###make datafram long so that we can compare core and edge treatments
allele_fold_long <- allele_fold %>%
  pivot_longer(cols = c("Change_stationary", "Change_core", "Change_edge"), 
               names_to = "treatment" , values_to = "frq_change") %>%
  mutate(Fold_treat = paste(Fold, treatment, sep = "_"))


## Combine with fitness data. Needs to have the same names for combanation.
fitness <-read.table("Beetle_fitness_filtered.delim", header = TRUE)
colnames(fitness) <- c("Landscape", "treatment", "", "Fitness")
fitness <- fitness %>% 
  select(1,2,4) %>%
  mutate(treatment = case_when(treatment == "core" ~ "Change_core",
                               treatment == "edge" ~ "Change_edge"))



allele_fold_long_fitness <- right_join(allele_fold_long, fitness,
                                       by = c("Landscape", "treatment"))

allele_fold_long_fitness <- allele_fold_long_fitness %>%
  mutate(final_frq = case_when(treatment == "Change_core" ~ `X8_C`,
                               treatment == "Change_edge" ~ `X8_E`))


allele_fold_change_squared <- allele_fold_long_fitness %>%
  mutate(change_squared = frq_change^2) %>%
  filter(Fold == 0) %>%
  group_by(treatment, Landscape) %>%
  summarize(change_0fold_squared = mean(change_squared , na.rm =TRUE), Fitness = mean(Fitness))


##Look at mjaor allele frq for 0 fold sites

allele_fold_major_frq <- allele_fold_long_fitness %>%
  #mutate(change_squared = frq_change^2) %>%
  filter(Fold == 0) %>%
  group_by(treatment, Landscape) %>%
  summarize(major_0fold_frq = mean(final_frq , na.rm =TRUE), Fitness = mean(Fitness))

## Just look at non squared change
allele_fold_change_directional <- allele_fold_long_fitness %>%
  #mutate(change_squared = frq_change^2) %>%
  filter(Fold == 0) %>%
  group_by(treatment, Landscape) %>%
  summarize(directional_0fold = mean(frq_change , na.rm =TRUE), Fitness = mean(Fitness))


######################## Calculate ############################

Full_data_long <- Full_filtered_data %>%
  pivot_longer(cols = c("Change_stationary", "Change_core", "Change_edge"), 
               names_to = "treatment" , values_to = "frq_change") %>%
  filter(!is.na(frq_change)) %>%
  filter(treatment != "Change_stationary")

Full_data_squared <- Full_data_long %>%
  mutate(change_squared = frq_change^2)
  

Full_data_squared_sum <- Full_data_squared %>%
  group_by(Landscape, treatment) %>%
  summarise(change_squared = (mean(change_squared, na.rm = TRUE)))


################## GERP genetic load ########################################




GERP_load <- Full_filtered_data %>%
  ##Take only informative sites
  filter(Max_RS > 0) %>%
  ## Create a value for the frequency of the minor allele times the RS score.
  ## The higher the frequency and RS the higher this measure of genetic load
  #This is calculated for each pooled sample
  mutate(Start = ((1-`X0_NA`) * RS), END_core = ((1-`X8_C`) * RS), 
         END_edge = ((1-`X8_E`) * RS), END_Stat = ((1-`X8_NA`) * RS))

## This is the method from Valk et al 2019 making the GERP score relative to
##The total frequency. What exactly is this doing and why?
pop_genetic_loads <- GERP_load %>%
  group_by(Landscape) %>%
  summarise(Start_load = (sum(Start, na.rm = TRUE)),
            edge = (sum(END_edge, na.rm = TRUE)),
            core = (sum(END_core, na.rm = TRUE)),
            `NA` = (sum(END_Stat, na.rm = TRUE)))




## Make is so that all values are in the same column this allows for easier comparison
pop_genetic_loads_long <- pop_genetic_loads %>%
  pivot_longer(cols =c("Start_load", "edge", "core", "NA"),
               names_to = "Location" , values_to = "Load")

pop_genetic_loads_long <- pop_genetic_loads_long %>%
  filter(Load != 0) %>%
  filter(Location != "Start_load" & Location != "NA")

colnames(pop_genetic_loads_long)[2] <- "treatment"


pop_genetic_loads_long <- pop_genetic_loads_long %>%
  mutate(treatment = case_when(treatment == "edge" ~ "Change_edge",
                               treatment == "core" ~ "Change_core"))

## load in grantham score sdataframe
#Grantham_loads <- read.table("Grantham_loads.delim", header = TRUE)





############# Join data #######################

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
write_delim(Final_summary_2 , "Final_summary.delim")





