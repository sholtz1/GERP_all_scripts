library(caroline)
library(tidyverse)


## Load in fitness data
beetle_growth_census <- read_csv("F2_census.csv")
beetle_growth_census <- beetle_growth_census %>%
  select(3:9)

## Since all beetles aren't counted final population needs to be calculated by weight. 
# TWL: Actually, this isn't the best way to calculate this. We know there are size differences
# among populations (particularly among edge and core populations), so this will not give you
# accurate information if you try to average across all populations. You need to do this for
# each row of the data frame to calculate the estimated beetle numbers from each patch. Let me
# know if you have questions about this, but you need to redo the fitness calculations here because
# the data you've gotten from this isn't correct. I'm sure it's pretty close, but it's not the
# correct values.
mean_beetle_weight <- mean(beetle_growth_census$weight_50, na.rm = TRUE)/50


##pipeline to get the fitness from the weights, Change to use weights when we actually have them
# TWL: Also, we don't need to round the estimates here. The weight estimates are inherently estimates,
# so we know they don't truly represent exact counts. If we use the estimates directly rather than
# rounding them, it will keep a slightly more accurate estimate of fitness in turn.
beetle_growth_census <- beetle_growth_census %>%
  mutate(weight_count = total_weight/mean_beetle_weight) %>%
  mutate(weight_count = round(weight_count))

## Need to use the actual counts when we have them, if not replace with estimates based on weight.
beetle_growth_census$Count <- ifelse((is.na(beetle_growth_census$Count)), beetle_growth_census$weight_count,beetle_growth_census$Count )

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

#Calculate the means for every pair of fitness datapoints in a new column
# TWL: This isn't going to work properly. But it looks like you already figured that
# out based on the code below. Why is this still in the script? There are much easier
# and more direct ways to do this calculation than how you currently have it structured.
# Try to clean up this section of the code to be cleaner and more intuitive.
for (i in 1:length(beetle_fitness$Landscape)) {
  means <- mean(c(beetle_fitness$Fitness[i], beetle_fitness$Fitness[i+1]))
  beetle_fitness$Shuff_mean[i] <- means
  
}

## Make the shuffeled mean the mean of the two fitnesses in the same landscape.
# TWL: Rather than using an ifelse statement that prints nothing to the console, 
# why not just use a traditional if statement?
for (i in 1:length(beetle_fitness$Landscape)) {
ifelse(
  beetle_fitness[i, "Landscape"] == beetle_fitness[i+1, "Landscape"], 
  beetle_fitness[i+1, "Shuff_mean"] <- beetle_fitness[i, "Shuff_mean"], print("")
)
}

##Use the shuffeled mean wehen the treatment is shuffeled.
##Keep the core and edge fitness data fro the structured treatment. 
beetle_fitness <- beetle_fitness %>%
  mutate(Fitness = case_when(Treatment == "S" ~ Shuff_mean,
         TRUE ~ Fitness)) 

## Turn  the location column to NA for shuffele treanment since this doesnt matter
beetle_fitness <- beetle_fitness %>%
  mutate(Location = case_when(Treatment == "S" ~ "NA",
                               TRUE ~ Location))


## Only keep 1 row per landscape for the shuffeled treatment since this is only 1 datapoint
beetle_fitness <- distinct(beetle_fitness)

## keep only the necessary columns
beetle_fitness <- beetle_fitness %>%
  select(1:4)


##Load library for start compare means after this is loaded the above functions wont
##work due to the fact grouping is somehow masked
library(ggpubr)

##Core vs edge fitness

## Lets visualize
beetle_fitness %>%
  filter(Treatment == "C") %>%
  ggplot(aes(x=Location, y= Fitness )) +
  geom_boxplot() +
  theme_classic() +
  stat_compare_means(label.x = 1.5, label.y = 8)+
  stat_compare_means(method = "t.test",label.x = 2, label.y = 8 )


write.delim(beetle_fitness, "Beetle_fitness_filtered")
