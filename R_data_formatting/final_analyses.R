## Flnal population summary data with analyses

library(RColorBrewer)
library(tidyverse) 
library(AICcmodavg)
library(gridExtra)
library(lme4)
library(bbmle)




### Read in final data
Final_summary <- read.delim("Final_summary.delim", sep = "")

##Scale all data for comparison



Scaled_summary <- Final_summary %>%
  mutate(Load = c(scale(Load)), 
         Load_full = c(scale(Load_full)), 
          High_effect_alleles = c(scale(High_effect_alleles)),
         High_effect_averages = c(scale(High_effect_averages)),
         change_squared = c(scale(change_squared)),
         Grantham_average = c(scale(Grantham_average)),
         major_0fold_frq = c(scale(major_0fold_frq)),
         Grantham_load = c(scale(Grantham_load))) %>% 
  select(treatment, Landscape, Load,Load_full,High_effect_alleles, change_squared, Grantham_load, 
         Grantham_average, High_effect_averages, Fitness, major_0fold_frq)

## List of all potential models to compare with AIC
## Here is an explanation of the column names since they can be a little confusing
## Load = GERP load with coding sites filtered out
## Load_full = GERP load for all GERP sites
## High_effect_alleles = VEP high effect alleles times the allele frequency summed together
##High_effect_averages = VEP high effect alleles times the allele frequency averaged
## Change squared = all allele frequency changes across the genome.
## Grantham_load = grantham score X allele frequency summed across the population
## Grantham_average = grantham score X allele frequency  averaged across the population
## major_0fold_frq = average of the major allele for all 0 fold sites in a population


all_interactions <- lmer(Fitness ~ Load*Grantham_load*High_effect_alleles + (1 | Landscape), data = Scaled_summary)
all_interactions_averages <- lmer(Fitness ~ Load*Grantham_average*High_effect_averages + (1 | Landscape), data = Scaled_summary)
GERP_interactions <- lmer(Fitness ~ Load*Grantham_load + Load*High_effect_alleles + (1 | Landscape), data = Scaled_summary)
no_interactions <- lmer(Fitness ~ Load +Grantham_load + High_effect_alleles + (1 | Landscape), data = Scaled_summary)
high_GERP_interaction <- lmer(Fitness ~ Load*High_effect_alleles + (1 | Landscape), data = Scaled_summary)
grantam_GERP_interaction <- lmer(Fitness ~ Load*Grantham_load + (1 | Landscape), data = Scaled_summary)
grantam_GERP_interaction_averages <- lmer(Fitness ~ Load*Grantham_average + (1 | Landscape), data = Scaled_summary)
GERP_only <- lmer(Fitness ~ Load_full + (1 | Landscape), data = Scaled_summary)
total_change <- lmer(Fitness ~ change_squared + (1 | Landscape), data = Scaled_summary)
high_only <- lmer(Fitness ~ High_effect_alleles + (1 | Landscape), data = Scaled_summary)
grantham_only <- lmer(Fitness ~ Grantham_load + (1 | Landscape), data = Scaled_summary)
high_grantham <- lmer(Fitness ~ Grantham_load + High_effect_alleles + (1 | Landscape), data = Scaled_summary)
GERP_grantham <- lmer(Fitness ~ Grantham_load + Load + (1 | Landscape), data = Scaled_summary)
GERP_high <- lmer(Fitness ~ Load + High_effect_alleles + (1 | Landscape), data = Scaled_summary)
high_only_averages <- lmer(Fitness ~ High_effect_averages + (1 | Landscape), data = Scaled_summary)
GERP_grantham_averages <- lmer(Fitness ~ Grantham_average + Load + (1 | Landscape), data = Scaled_summary)
GERP_0fold_interaction <- lmer(Fitness ~ Load*major_0fold_frq + (1 | Landscape), data = Scaled_summary)

# Look at AIC of models

models <- list(all_interactions, GERP_interactions, no_interactions, high_GERP_interaction, 
               GERP_only, total_change, grantam_GERP_interaction, high_only, grantham_only, high_grantham,
               GERP_grantham, GERP_high, all_interactions_averages, grantam_GERP_interaction_averages,
               high_only_averages, GERP_grantham_averages, GERP_0fold_interaction)
mod.names <- c( "all_interactions", "GERP_interactions", "no_interactions", "high_GERP_interaction", 
                "GERP_only", "total_change", "grantam_GERP_interaction", "high_only", "Grantham_only", 
                "high_grantham", "GERP_grantham", "GERP_high", "all_interactions_averages",
                "grantam_GERP_interaction_averages", "high_only_averages", "GERP_grantham_averages",
                "GERP_0fold_interaction")

AICctab(models, mnames = mod.names)

########### are loads better predictors of fitness at the edge than the core?#############
GERP_only_edge <- lm(Fitness ~ Load_full , data = Scaled_summary %>% filter(treatment == "Change_edge"))
GERP_only_core <- lm(Fitness ~ Load_full , data = Scaled_summary %>% filter(treatment == "Change_core"))
GERP_only_test <- lm(Fitness ~ Load_full , data = Scaled_summary)

models <- list(GERP_only_edge, GERP_only_core)
mod.names <- c("GERP_only_edge", "GERP_only_core")

AICctab(models, mnames = mod.names)


### look at final models summaries and residual plots
summary(GERP_only_edge)
summary(GERP_only_core)
summary(GERP_only_test)

res <- resid(GERP_only)
plot(fitted(GERP_only), res)



