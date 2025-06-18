## Flnal population summary data with analyses

library(RColorBrewer)
library(tidyverse) 
library(AICcmodavg)
library(gridExtra)
library(lme4)
library(bbmle)
library(plotly)
library(MuMIn)
library(boot)
library(lmerTest)
library(patchwork)
library(emmeans)




## Flnal population summary data with analyses



### Theme for GG plots

My_Theme = theme(   axis.title.x = element_text(size = 16),   
                    axis.title.y = element_text(size = 16),
                    legend.text = element_text(size=17), 
                    axis.text= element_text(size=15),
                    legend.title=element_blank())



### Read in final data
Final_summary <- read.delim("Final_summary_derived.delim", sep = "")

##Scale all data for comparison



Scaled_summary <- Final_summary %>%
  mutate(Load = c(scale(Load)), 
         Load_full = c(scale(Load_full)),
         Load_homozygous = c(scale(Load_homozygous)),
         Load_homo_filtered = c(scale(Load_homo_filtered)),
         High_effect_alleles = c(scale(High_effect_alleles)),
         High_effect_homo = c(scale(High_effect_homo)),
         Grantham_load = c(scale(Grantham_load)),
         Grantham_homo = c(scale(Grantham_homo)),
         Diversity = c(scale(Diversity)),
         HF_GERP = c(scale(HF_GERP)),
         HF_VEP = c(scale(HF_VEP))
         )%>%
  select(treatment, Fitness, Load_full, Load_homozygous, High_effect_alleles, 
         High_effect_homo, Diversity,
         Grantham_load, Grantham_homo, Load,Load_homo_filtered, Landscape, HF_GERP, HF_VEP)

## add in the block numbers so we can include them as a random effect.
Scaled_summary <- Scaled_summary %>%
  mutate(Block = case_when(
    Landscape >= 21 & Landscape <= 40 ~ 2,
    Landscape >= 41 & Landscape <= 60 ~ 3,
    Landscape >= 61 & Landscape <= 80 ~ 4
  ))



## List of all potential models to compare with AIC
## Here is an explanation of the column names since they can be a little confusing
## Load = GERP load with coding sites filtered out
## Load_full = GERP load for all GERP sites
## High_effect_alleles = VEP high effect alleles times the allele frequency summed together
## Grantham_load = grantham score X allele frequency summed across the population
## All things suffixed with _homo are using hardy weinberg to estimate the number frequency of homozygous sites
all_interactions <- lmer(Fitness ~ Load*Grantham_load*High_effect_alleles + (1 | Block), data = Scaled_summary)
Load_Grantham_interactions <- lmer(Fitness ~ Load*Grantham_load + High_effect_alleles + (1 | Block) , data = Scaled_summary)
Grantham_High_interactions <- lmer(Fitness ~ Load + Grantham_load*High_effect_alleles + (1 | Block)  , data = Scaled_summary)
Load_high_interactions <- lmer(Fitness ~ Load*High_effect_alleles+ Grantham_load + (1 | Block)  , data = Scaled_summary)
No_interactions <- lmer(Fitness ~ Load+High_effect_alleles+ Grantham_load + (1 | Block)  , data = Scaled_summary)
Only_Load_high_interactions <- lmer(Fitness ~ Load*High_effect_alleles + (1 | Block)  , data = Scaled_summary)
Only_Load_high <- lmer(Fitness ~ Load+High_effect_alleles + (1 | Block)  , data = Scaled_summary)
Only_Load_Grantham_interactions <- lmer(Fitness ~ Load*Grantham_load + (1 | Block)  , data = Scaled_summary)
Only_Load_Grantham <- lmer(Fitness ~ Load+Grantham_load + (1 | Block)  , data = Scaled_summary)
Only_High_Grantham_interactions <- lmer(Fitness ~ High_effect_alleles*Grantham_load + (1 | Block)  , data = Scaled_summary)
Only_High_Grantham <- lmer(Fitness ~ High_effect_alleles+Grantham_load + (1 | Block)  , data = Scaled_summary)
Only_Load <- lmer(Fitness ~ Load_full + (1 | Block)  , data = Scaled_summary)
Only_Grantham <- lmer(Fitness ~ Grantham_load + (1 | Block)  , data = Scaled_summary)
Only_High <- lmer(Fitness ~ High_effect_alleles + (1 | Block)  , data = Scaled_summary)
Only_Diversity <- lmer(Fitness ~ Diversity + (1 | Block)  , data = Scaled_summary)


## look at the subset of models without homozygosity


models_nohomo <- list(all_interactions, Load_Grantham_interactions, Grantham_High_interactions,
               Load_high_interactions, No_interactions, Only_Load_high_interactions,
               Only_Load_high, Only_Load_Grantham_interactions, Only_Load_Grantham,
               Only_High_Grantham_interactions, Only_High_Grantham, Only_Load, Only_Grantham,
               Only_High,Only_Diversity)


mod.names_nohomo <- c(
  "all_interactions", "Load_Grantham_interactions","Grantham_High_interactions","Load_high_interactions",
  "No_interactions","Only_Load_high_interactions","Only_Load_high",
  "Only_Load_Grantham_interactions","Only_Load_Grantham","Only_High_Grantham_interactions",
  "Only_High_Grantham"," Only_Load","Only_Grantham","Only_High", "Only_Diversity")

AICctab(models_nohomo, mnames = mod.names_nohomo)

summary(Only_Load_high)

## make a summary DF that has the different in load for each landscape between core and edge 
##with a positive value being an increase at the edge
Pairwise_summary <- Scaled_summary %>%
  filter(Landscape != 75) %>%
  group_by(Landscape) %>%
  summarise(Grantham_load = Grantham_load[1]-Grantham_load[2], 
            High_effect_alleles = High_effect_alleles[1]- High_effect_alleles[2],
            Load_full = Load_full[1]- Load_full[2],
            Diversity = Diversity[1]- Diversity[2],
            Fitness = Fitness[1]- Fitness[2])




# Calculate P values to include in the figures
Grantham_load_T <- t.test(Pairwise_summary$Grantham_load, mu = 0, alternative = "greater")
Load_full_T <- t.test(Pairwise_summary$Load_full, mu = 0, alternative = "greater")
High_effect_alleles_T <- t.test(Pairwise_summary$High_effect_alleles, mu = 0, alternative = "greater")
Diversity_T <- t.test(Pairwise_summary$Diversity, mu = 0, alternative = "less")
Fitness_T <- t.test(Pairwise_summary$Fitness, mu = 0, alternative = "less")

#Pairwise summary of GERP load
Pairwise_GERP_load <-Pairwise_summary %>%
  ggplot(aes(x = Load_full)) +
  geom_histogram(bins = 10, fill = "black")+
  theme_classic() +My_Theme +
  xlab("GERP load")+
  ylab("")+
  #coord_cartesian(xlim=c(-2,4), ylim = c(0,5))+
  geom_text(aes(x = Inf, y = Inf, label = paste("p =", signif(Load_full_T$p.value, 3))),
            hjust = 1, vjust = 1, size = 5) +
  theme(plot.margin = margin(0,0,1.5,0, "cm"))



#Pairwise summary of High Effect load
Pairwise_Load_High <-Pairwise_summary %>%
  ggplot(aes(x = High_effect_alleles)) +
  geom_histogram(bins = 10 , fill = "black")+
  theme_classic() +My_Theme +
  xlab("High Effect Load")+
  ylab("")+
  #coord_cartesian(xlim=c(-2,4), ylim = c(0,5))+
  geom_text(aes(x = Inf, y = Inf, label = paste("p =", signif(High_effect_alleles_T$p.value, 3))),
            hjust = 1, vjust = 1, size = 5)+
  theme(plot.margin = margin(0,0,1.5,0, "cm"))


#Pairwise summary of Grantham load
Pairwise_Grantham_Load <-Pairwise_summary %>%
  ggplot(aes(x = Grantham_load)) +
  geom_histogram(bins = 15 , fill = "black")+
  theme_classic() +My_Theme +
  xlab("Grantham Load")+
  ylab("")+
  #coord_cartesian(xlim=c(-2,4), ylim = c(0,5))+
  geom_text(aes(x = Inf, y = Inf, label = paste("p =", signif(Grantham_load_T$p.value, 3))),
            hjust = 1, vjust = 1, size = 5) +
  theme(plot.margin = margin(0,0,1.5,0, "cm"))

#Pairwise summary of Diversity
Pairwise_Diversity <-Pairwise_summary %>%
  ggplot(aes(x = Diversity)) +
  geom_histogram(bins = 15 , fill = "springgreen3")+
  theme_classic() +My_Theme +
  xlab("Diversity")+
  ylab("")+
  #coord_cartesian(xlim=c(-2,4), ylim = c(0,5))+
  geom_text(aes(x = Inf, y = Inf, label = paste("p =", signif(Diversity_T$p.value, 3))),
            hjust = 1, vjust = 1, size = 5) +
  theme(plot.margin = margin(0,0,1.5,0, "cm"))


#Pairwise summary of Fitness
Pairwise_Fitness <-Pairwise_summary %>%
  ggplot(aes(x = Fitness)) +
  geom_histogram(bins = 15 , fill = "springgreen3")+
  theme_classic() +My_Theme +
  xlab("Fitness")+
  ylab("")+
  #coord_cartesian(xlim=c(-2,4), ylim = c(0,5))+
  geom_text(aes(x = Inf, y = Inf, label = paste("p =", signif(Fitness_T$p.value, 3))),
            hjust = 1, vjust = 1, size = 5) +
  theme(plot.margin = margin(0,0,1.5,0, "cm"))



## combine pairwise plots into a single figure
Pairwise_plots_combined <-( Pairwise_GERP_load + Pairwise_Load_High + Pairwise_Grantham_Load) / (Pairwise_Diversity+ Pairwise_Fitness) + plot_annotation(tag_levels = 'A') & 
  theme(plot.tag.position  = c(.25, .96), plot.tag = element_text(size = 25))

Pairwise_plots_combined

pdf(file = "Pairwise_plots_combined.pdf", width = 12, height = 6, onefile = FALSE, paper = "special")
Pairwise_plots_combined
dev.off()


###### Linear comparisons between fitness and loads/diversity



# calculate R squared to use add to figures
R2_Only_Load <- r.squaredGLMM(Only_Load)
R2_Only_High <- r.squaredGLMM(Only_High)
R2_Only_Grantham <- r.squaredGLMM(Only_Grantham)
R2_Diversity <- r.squaredGLMM(Only_Diversity)


#GERP load vs fitness
GERP_linear_plot <- Scaled_summary %>%
  ggplot(aes(x= Load_full, y= Fitness)) +
  geom_point()+
  geom_smooth(method = "lm") +
  xlab("GERP Load") + ylim(3.5, 8.5)+
  theme_classic() +My_Theme +
  geom_text(aes(x = Inf, y = Inf, label = paste("R²m =", signif(digits = 4,R2_Only_Load[1,1]))),
            hjust = 1, vjust = 1.5, size = 6) +
  geom_text(aes(x = Inf, y = Inf, label = paste("R²c =", signif(digits = 4, R2_Only_Load[1,2]))),
            hjust = 1, vjust = 3, size = 6)

GERP_linear_plot


#Grantham load vs fitness
Grantham_linear_plot <- Scaled_summary %>%
  ggplot(aes(x= Grantham_load, y= Fitness)) +
  geom_point()+
  geom_smooth(method = "lm") +
  xlab("Grantham Load") + ylab("Fitness")+
  theme_classic() +My_Theme + ylim(3.5, 8.5)+
  geom_text(aes(x = Inf, y = Inf, label = paste("R²m =",  signif(digits = 4, R2_Only_Grantham[1,1]))),
            hjust = 1, vjust = 1.5, size = 6) +
  geom_text(aes(x = Inf, y = Inf, label = paste("R²c =", signif(digits = 4,R2_Only_Grantham[1,2]))),
            hjust = 1, vjust = 3, size = 6)

Grantham_linear_plot

#High Effect load vs fitness
HighE_linear_plot <- Scaled_summary %>%
  ggplot(aes(x= High_effect_alleles, y= Fitness)) +
  geom_point()+
  geom_smooth(method = "lm") +
  xlab("High Effect Load") + ylab("")+
  theme_classic() +My_Theme + ylim(3.5, 8.5)+
  geom_text(aes(x = Inf, y = Inf, label = paste("R²m =", signif(digits = 4,R2_Only_High[1,1]))),
            hjust = 1, vjust = 1.5, size = 6) +
  geom_text(aes(x = Inf, y = Inf, label = paste("R²c =", signif(digits = 4,R2_Only_High[1,2]))),
            hjust = 1, vjust = 3, size = 6) 

#diversity vs fitness
Diversity_linear_plot <- Scaled_summary %>%
  ggplot(aes(x= Diversity, y= Fitness)) +
  geom_point()+
  geom_smooth(method = "lm", fill = "springgreen3") +
  xlab("Diversity") + ylab("")+
  theme_classic() +My_Theme + ylim(3.5, 8.5)+
  geom_text(aes(x = Inf, y = Inf, label = paste("R²m =", signif(digits = 4,R2_Diversity[1,1]))),
            hjust = 1, vjust = 1.5, size = 6) +
  geom_text(aes(x = Inf, y = Inf, label = paste("R²c =", signif(digits = 4,R2_Diversity[1,2]))),
            hjust = 1, vjust = 3, size = 6) 



## combine all plots together to make the final figure

Combined_linear_plot <-  GERP_linear_plot + HighE_linear_plot + Grantham_linear_plot + Diversity_linear_plot&
  plot_annotation(tag_levels = 'A') & 
  theme(plot.tag.position  = c(.25, .96), plot.tag = element_text(size = 25))

Combined_linear_plot

pdf(file = "Combined_linear_plot.pdf", width = 12, height = 6, onefile = FALSE, paper = "special")
Combined_linear_plot
dev.off()






################# Split data into bins for analysis  ################################################


### Read in full data
Full_filtered_data <- read.table("Grantham_GERP_derived.delim", header = TRUE)


### Make GERP bins####


##  Filter for sites with a positive GERP score 
GERP_bins <- Full_filtered_data %>%
  filter(RS > 0) %>%
  filter(!is.na(X8_C)) %>%
  filter(!is.na(X8_E)) %>%
  select(Chromosome, Location, Landscape ,X8_C, X8_E, RS) %>%
  pivot_longer(cols = c(X8_C, X8_E), names_to = "Treatment", values_to = "frequency") 

## get a new column that identifies data by gerp score and by edge vs core
GERP_bins <- GERP_bins %>%
  mutate( GERP_bin = case_when(
    RS > 0 & RS <=2 & Treatment == "X8_E" ~ "0-2 Edge" ,
    RS > 0 & RS <=2 & Treatment == "X8_C" ~ "0-2 Core",
    RS > 2 & RS <=4 & Treatment == "X8_E" ~ "2-4 Edge",
    RS > 2 & RS <=4 & Treatment == "X8_C" ~ "2-4 Core",
    RS > 4 & Treatment == "X8_E" ~ "4+ Edge",
    RS > 4 & Treatment == "X8_C" ~ "4+ Core",
  ))

## make sure there arent any duplicate data points
GERP_bins <- GERP_bins %>%
  unique()

GERP_bins$GERP_bin <- as.character(GERP_bins$GERP_bin)


############################# Fixed sites GERP models

## get the total number of sites fixed for the derived allele in each landscape and bin

GERP_bins_fixed <- GERP_bins %>%
  group_by(Landscape, GERP_bin) %>%
  dplyr::summarise(fixed = sum(frequency == 1)) %>%
  mutate(Block = case_when(
    Landscape >= 21 & Landscape <= 40 ~ 2,
    Landscape >= 41 & Landscape <= 60 ~ 3,
    Landscape >= 61 & Landscape <= 80 ~ 4
  ))

## split the column so that the Bin identifier and the core vs edge identifier are in separate columns
GERP_bins_fixed[c('Bin', 'Location')] <- str_split_fixed(GERP_bins_fixed$GERP_bin, " ", n = 2)

# Fit a generalized linear model with a Poisson distribution to the data
GERP_GLMM_fixed <- glmer(fixed ~ Bin*Location + (1 | Block), data = GERP_bins_fixed,
                   family = poisson(link = "log"))
summary(GERP_GLMM_fixed)

## use emmeans to get pariwise P values
em_gerp <- emmeans(GERP_GLMM_fixed, ~ Bin *Location)
pairs(em_gerp, by = "Bin")

# Use the boot library to bootstrap confidence intervals for the model. First, make
# a dataframe for the predictions containing all the combinations of the fixed
# effects. This approach allows us to calculate the true confidence interval for
# just the fixed effects. Because of the random effects, we need to use a bootstrap
# approach to approximate the intervals as they can't be computed analytically.
newdat <- expand.grid(Bin = c("0-2", "2-4", "4+"), Location = c("Core", "Edge"))

# Make a function to generate predictions based on the fitted model and the
#    newdat data frame
preds <- function(fitobj){
  preds <- predict(fitobj, newdata=newdat, re.form=NA)
  return(preds)
}

# Now run the bootstrap analysis
bootpreds <- bootMer(GERP_GLMM_fixed, preds, nsim=50)

# Calculate the CIs
conf_ints <- vector(mode='list', length=6)  # 6 combinations of fixed effects
for(i in 1:6){
  conf_ints[[i]] <- boot.ci(bootpreds, type='perc', index=i)
}

# Store everything together for plotting
GERP_results <- NULL
for(i in 1:6){
  GERP_results <- rbind(GERP_results, c(conf_ints[[i]]$t0, conf_ints[[i]]$percent[1,4:5]))
}
colnames(GERP_results) <- c('mean', 'lower', 'upper')
GERP_results_fixed <- cbind(newdat, GERP_results)

# Exponentiate the values to undue the log transformation used in the link function
GERP_results_fixed[,3:5] <- exp(GERP_results_fixed[,3:5])


## graph the data so that each point represents a landscape and bin 
##and use the mean and 95% confidence intervals from the model
GERP_fixed_plot <- ggplot(GERP_results_fixed, aes(x=Bin, y=mean, color = Location)) + 
  geom_point(shape = 19, size = 4, position = position_dodge(0.9)) + 
  geom_point(data = GERP_bins_fixed, aes(y = fixed), shape = 21, size = 2, 
             position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.9),
             stroke = 0.5) +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=0.5, position=position_dodge(0.9)) +
  theme_classic() + My_Theme +
  xlab("GERP Score")+  theme(legend.position = "none") +
  ylab("Number of Fixed sites") + scale_color_manual(values=c("orangered4", "#56B4E9"))

GERP_fixed_plot

# Now save the plot as a pdf of the desired dimensions for a publication quality figure
#pdf(file = "GERP_fixed_plot.pdf", width = 5, height = 3, onefile = FALSE, paper = "special")
#GERP_fixed_plot
#dev.off()


############################### Fixed sites VEP model ####################################
## Get a new column to identify the VEP effect and core vs edge for each data point
#square to get the proportion of homozygous sites

VEP_bins <- Full_filtered_data %>%
  filter(!is.na(IMPACT)) %>%
  filter(!is.na(X8_C)) %>%
  filter(!is.na(X8_E)) %>%
  select(Chromosome, Location, Landscape ,X8_C, X8_E, IMPACT) %>%
  pivot_longer(cols = c(X8_C, X8_E), names_to = "Treatment", values_to = "frequency")

##change names for analysis
VEP_bins <- VEP_bins %>%
  mutate(VEP_bin = case_when(
    IMPACT == "MODIFIER" ~ "Modifier" ,
    IMPACT == "LOW" ~ "Low",
    IMPACT == "MODERATE" ~ "Moderate",
    IMPACT == "HIGH" ~ "High"
  ))

VEP_bins <- VEP_bins %>%
  mutate(Location = case_when(
    Treatment == "X8_E" ~ "Edge" ,
    Treatment == "X8_C" ~ "Core",
  ))


## remove modifier sites so that we have only Low high and moderate effect sites
VEP_bins <- VEP_bins%>%
  filter(VEP_bin != "Modifier")


## get the number of sites fixed for the derived allele in each bin
VEP_bins_fixed <- VEP_bins %>%
  group_by(Location, Landscape, VEP_bin) %>%
  dplyr::summarise(fixed = sum(frequency == 1)) %>%
  mutate(Block = case_when(
    Landscape >= 21 & Landscape <= 40 ~ 2,
    Landscape >= 41 & Landscape <= 60 ~ 3,
    Landscape >= 61 & Landscape <= 80 ~ 4
  ))



# Fit a generalized linear model with a Poisson distribution to the data
VEP_GLMM_fixed <- glmer(fixed ~ VEP_bin*Location + (1 | Block), data = VEP_bins_fixed,
                  family = poisson(link = "log"))

summary(VEP_GLMM_fixed)

## use em means to get pairwise P values for more specific comparisons
em_vep <- emmeans(VEP_GLMM_fixed, ~ VEP_bin *Location)
pairs(em_vep, by = "VEP_bin")

# Now use the same process to generate confidence intervals
newdat <- expand.grid(VEP_bin = c("Low", "Moderate", "High"), Location = c("Core", "Edge"))
bootpreds <- bootMer(VEP_GLMM_fixed, preds, nsim=50)
conf_ints <- vector(mode='list', length=6)  # 6 combinations of fixed effects
for(i in 1:6){
  conf_ints[[i]] <- boot.ci(bootpreds, type='perc', index=i)
}

# Store everything together for plotting
VEP_results <- NULL
for(i in 1:6){
  VEP_results <- rbind(VEP_results, c(conf_ints[[i]]$t0, conf_ints[[i]]$percent[1,4:5]))
}
colnames(VEP_results) <- c('mean', 'lower', 'upper')
VEP_results_fixed <- cbind(newdat, VEP_results)

# Exponentiate the values to undue the log transformation used in the link function
VEP_results_fixed[,3:5] <- exp(VEP_results_fixed[,3:5])

VEP_fixed_plot <- ggplot(VEP_results_fixed, aes(x=VEP_bin, y=mean, color = Location)) + 
  geom_point(shape = 19, size = 4, position = position_dodge(0.9)) + 
  geom_point(data = VEP_bins_fixed, aes(y = fixed), shape = 21, size = 2, 
             position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.9),
             stroke = 0.5) +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=0.5, position=position_dodge(0.9)) +
  theme_classic() + My_Theme + 
  xlab("VEP Category")+
  ylab("") + scale_color_manual(values=c("orangered4", "#56B4E9"))

VEP_fixed_plot 
##### Combine fixed plots into 1


Combined_fixed_plots <- GERP_fixed_plot + VEP_fixed_plot + plot_annotation(tag_levels = 'A') & 
  theme(plot.tag.position  = c(.25, .96), plot.tag = element_text(size = 25))
Combined_fixed_plots 



pdf(file = "Combined_fixed_plots.pdf", width = 12, height = 6, onefile = FALSE, paper = "special")
Combined_fixed_plots
dev.off()




####################### Linear models on GERP bins ######################################


## Get the average allele frequency for each bin/location
Average_GERP_bins <- GERP_bins %>% 
  group_by(GERP_bin, Landscape) %>%
  #filter(frequency !=0) %>%
  summarise(frequency = mean(frequency))

## get the block from the landscape numbers
Average_GERP_bins <- Average_GERP_bins %>%
  mutate(Block = case_when(
    Landscape >= 21 & Landscape <= 40 ~ 2,
    Landscape >= 41 & Landscape <= 60 ~ 3,
    Landscape >= 61 & Landscape <= 80 ~ 4
  ))


## split the idenftifieres so that we have the bin and the location in seperate columns 
Average_GERP_bins[c('Bin', 'Location')] <- str_split_fixed(Average_GERP_bins$GERP_bin, " ", n = 2)



GERP_bins_model<- lmer(frequency ~ Bin*Location + (1 | Block), data = Average_GERP_bins)

summary(GERP_bins_model)
## use em to look at comparisions between core and edge for each bin
em_gerp_average <- emmeans(GERP_bins_model, ~ Bin *Location)
pairs(em_gerp_average, by = "Bin")

# Use the boot library to bootstrap confidence intervals for the model. First, make
# a dataframe for the predictions containing all the combinations of the fixed
# effects. This approach allows us to calculate the true confidence interval for
# just the fixed effects. Because of the random effects, we need to use a bootstrap
# approach to approximate the intervals as they can't be computed analytically.
newdat <- expand.grid(Bin = c("0-2", "2-4", "4+"), Location = c("Core", "Edge"))

# Make a function to generate predictions based on the fitted model and the
#    newdat data frame
preds <- function(fitobj){
  preds <- predict(fitobj, newdata=newdat, re.form=NA)
  return(preds)
}

# Now run the bootstrap analysis
bootpreds <- bootMer(GERP_bins_model, preds, nsim=50)

# Calculate the CIs
conf_ints <- vector(mode='list', length=6)  # 6 combinations of fixed effects
for(i in 1:6){
  conf_ints[[i]] <- boot.ci(bootpreds, type='perc', index=i)
}

# Store everything together for plotting
GERP_results <- NULL
for(i in 1:6){
  GERP_results <- rbind(GERP_results, c(conf_ints[[i]]$t0, conf_ints[[i]]$percent[1,4:5]))
}
colnames(GERP_results) <- c('mean', 'lower', 'upper')
GERP_results <- cbind(newdat, GERP_results)

####### plot using the confidence intervals from the bootstrap

GERP_bins_plot <- ggplot(GERP_results, aes(x=Bin, y=mean, color = Location)) + 
  geom_point(shape = 19, size = 4, position = position_dodge(0.9)) + 
  geom_point(data = Average_GERP_bins, aes(y = frequency), shape = 21, size = 2, 
             position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.9),
             stroke = 0.5) +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=0.5, position=position_dodge(0.9)) +
  theme_classic() + My_Theme + theme(legend.position = "none") +
  xlab("GERP Score")+  scale_color_manual(values=c("orangered4", "#56B4E9"))+
  ylab("Avegrage Allele Frequency")

GERP_bins_plot
# Now save the plot as a pdf of the desired dimensions for a publication quality figure





##### Linear models of VEP bins
Average_VEP_bins <- VEP_bins %>% 
  group_by(VEP_bin, Landscape, Location) %>%
  #filter(frequency !=0) %>%
  summarise(frequency = mean(frequency))

Average_VEP_bins <- Average_VEP_bins %>%
  mutate(Block = case_when(
    Landscape >= 21 & Landscape <= 40 ~ 2,
    Landscape >= 41 & Landscape <= 60 ~ 3,
    Landscape >= 61 & Landscape <= 80 ~ 4
  ))



VEP_bins_model<- lmer(frequency ~ VEP_bin*Location + (1 | Block), data = Average_VEP_bins)

summary(VEP_bins_model)

## use em to look at comparisions between core and edge for each bin
em_vep_average <- emmeans(VEP_bins_model, ~ VEP_bin *Location)
pairs(em_vep_average, by = "VEP_bin")

# Now use the same process to generate confidence intervals
newdat <- expand.grid(VEP_bin = c("Low", "Moderate", "High"), Location = c("Core", "Edge"))
bootpreds <- bootMer(VEP_bins_model, preds, nsim=50)
conf_ints <- vector(mode='list', length=6)  # 6 combinations of fixed effects
for(i in 1:6){
  conf_ints[[i]] <- boot.ci(bootpreds, type='perc', index=i)
}

# Store everything together for plotting
VEP_results <- NULL
for(i in 1:6){
  VEP_results <- rbind(VEP_results, c(conf_ints[[i]]$t0, conf_ints[[i]]$percent[1,4:5]))
}
colnames(VEP_results) <- c('mean', 'lower', 'upper')
VEP_results <- cbind(newdat, VEP_results)

####### plot using the confidence intervals from the bootstrap

VEP_bins_plot <- ggplot(VEP_results, aes(x=VEP_bin, y=mean, color = Location)) + 
  geom_point(shape = 19, size = 4, position = position_dodge(0.9)) + 
  geom_point(data = Average_VEP_bins, aes(y = frequency), shape = 21, size = 2, 
             position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.9),
             stroke = 0.5) +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=0.5, position=position_dodge(0.9)) +
  theme_classic() + My_Theme + 
  xlab("VEP Bin")+  scale_color_manual(values=c("orangered4", "#56B4E9"))+
  ylab("")

VEP_bins_plot

combined_bins_plot <- GERP_bins_plot +VEP_bins_plot + plot_annotation(tag_levels = 'A') & 
  theme(plot.tag.position  = c(.18, 1), plot.tag = element_text(size = 25))

combined_bins_plot

pdf(file = "combined_bins_plot .pdf", width = 12, height = 6, onefile = FALSE, paper = "special")
combined_bins_plot
dev.off()

