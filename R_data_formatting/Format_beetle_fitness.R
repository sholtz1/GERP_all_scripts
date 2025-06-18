library(caroline)
library(tidyverse)


## Load in fitness data
beetle_growth_census <- read_csv("F2_census.csv")
beetle_growth_census <- beetle_growth_census %>%
  select(3:9)



##pipeline to get the fitness from the weights,get beetle count estimates from the weights
beetle_growth_census <- beetle_growth_census %>%
  mutate(weight_count = total_weight/(weight_50/50))

## Need to use the actual counts when we have them, if not replace with estimates based on weight.
beetle_growth_census$Count <- ifelse((is.na(beetle_growth_census$Count)), beetle_growth_census$weight_count,beetle_growth_census$Count )


### Get fitness for each column
beetle_growth_census <-beetle_growth_census %>%
  mutate(fitness = Count/Density) %>%
  filter(!is.na(fitness))

## get fitness for each landscape and location
beetle_fitness <- beetle_growth_census %>%
  group_by(Landscape, Location, Treatment) %>%
  summarize(Fitness = mean(fitness))

## WE need to combine the means for the shuffled edge and core since these 
##Should be biologically equivalent 

#make row to put mean into
beetle_fitness$Shuff_mean <- rep(NA, length(beetle_fitness$Landscape))

#Calculate the means for every pair of fitness data points in a new column
for (i in 1:length(beetle_fitness$Landscape)) {
  means <- mean(c(beetle_fitness$Fitness[i], beetle_fitness$Fitness[i+1]))
  beetle_fitness$Shuff_mean[i] <- means
  
}

## Make the shuffled mean the mean of the two fitnesses in the same landscape.
for (i in 1:length(beetle_fitness$Landscape)) {
ifelse(
  beetle_fitness[i, "Landscape"] == beetle_fitness[i+1, "Landscape"], 
  beetle_fitness[i+1, "Shuff_mean"] <- beetle_fitness[i, "Shuff_mean"], print("")
)
}

##Use the shuffled mean when the treatment is shuffled.
##Keep the core and edge fitness data fro the structured treatment. 
beetle_fitness <- beetle_fitness %>%
  mutate(Fitness = case_when(Treatment == "S" ~ Shuff_mean,
         TRUE ~ Fitness)) 

## Turn  the location column to NA for shuffle treatment since this doesn't matter
beetle_fitness <- beetle_fitness %>%
  mutate(Location = case_when(Treatment == "S" ~ "NA",
                               TRUE ~ Location))


## Only keep 1 row per landscape for the shuffled treatment since this is only 1 datapoint
beetle_fitness <- distinct(beetle_fitness)

## keep only the necessary columns
beetle_fitness <- beetle_fitness %>%
  select(1:4)


##Write final fitness file to use in final analyses

write.delim(beetle_fitness, "Beetle_fitness_filtered.delim")
