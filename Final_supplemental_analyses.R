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
HF_GERP <- lmer(Fitness ~ HF_GERP + (1 | Block)  , data = Scaled_summary)
HF_VEP <- lmer(Fitness ~ HF_VEP + (1 | Block)  , data = Scaled_summary)
Filtered_Load <-  lmer(Fitness ~ Load + (1 | Block)  , data = Scaled_summary)

Homo_Load <- lmer(Fitness ~ Load_homozygous + (1 | Block)  , data = Scaled_summary)
Homo_load_filtered <- lmer(Fitness ~ Load_homo_filtered + (1 | Block)  , data = Scaled_summary)
Only_High_homo <- lmer(Fitness ~ High_effect_homo + (1 | Block)  , data = Scaled_summary)
Grantham_homo <- lmer(Fitness ~ Grantham_homo + (1 | Block)  , data = Scaled_summary)
Load_high_homo <- lmer(Fitness ~ Load_homozygous +High_effect_homo + (1 | Block)  , data = Scaled_summary)
Load_high_homo_interaction <- lmer(Fitness ~ Load_homozygous*High_effect_homo + (1 | Block)  , data = Scaled_summary)
Grantham_Load_homo <- lmer(Fitness ~ Load_homozygous + Grantham_homo + (1 | Block)  , data = Scaled_summary)
Grantham_Load_homo_interaction <- lmer(Fitness ~ Load_homozygous*Grantham_homo + (1 | Block)  , data = Scaled_summary)
Grantham_high_homo <- lmer(Fitness ~ High_effect_homo+ Grantham_homo + (1 | Block)  , data = Scaled_summary)
Grantham_high_homo_interaction <- lmer(Fitness ~ Load_homozygous*Grantham_homo + (1 | Block)  , data = Scaled_summary)
Homo_all_loads <-  lmer(Fitness ~ Load_homozygous+Grantham_homo+High_effect_homo + (1 | Block)  , data = Scaled_summary)


# Look at AICc of models

models <- list(all_interactions, Load_Grantham_interactions, Grantham_High_interactions,
               Load_high_interactions, No_interactions, Only_Load_high_interactions,
               Only_Load_high, Only_Load_Grantham_interactions, Only_Load_Grantham,
               Only_High_Grantham_interactions, Only_High_Grantham, Only_Load, Only_Grantham,
               Only_High, Homo_Load, Only_High_homo, Load_high_homo, Grantham_homo, Grantham_Load_homo,
               Load_high_homo_interaction, Grantham_high_homo, Grantham_high_homo_interaction,
               Homo_all_loads, Grantham_Load_homo_interaction,Only_Diversity, HF_GERP, HF_VEP, Homo_load_filtered, Filtered_Load)


mod.names <- c(
  "all_interactions", "Load_Grantham_interactions","Grantham_High_interactions","Load_high_interactions",
  "No_interactions","Only_Load_high_interactions","Only_Load_high",
  "Only_Load_Grantham_interactions","Only_Load_Grantham","Only_High_Grantham_interactions",
  "Only_High_Grantham"," Only_Load","Only_Grantham","Only_High", "Homo_load", "Only_High_homo", 
  "Load_high_homo", "Grantham_homo", "Grantham_Load_homo", "Load_high_homo_interaction",
  "Grantham_high_homo", "Grantham_high_homo_interaction", "Homo_all_loads",
  "Grantham_Load_homo_interaction", "Only_Diversity", "HF_GERP", "HF_VEP", "Homo_load_filtered", "Filtered_Load")


AICctab(models, mnames = mod.names)





## make a summary DF that has the different in load for each landscape between core and edge 
##with a positive value being an increase at the edge
Pairwise_summary <- Scaled_summary %>%
  filter(Landscape != 75) %>%
  group_by(Landscape) %>%
  summarise(Grantham_homo = Grantham_homo[1]-Grantham_homo[2], 
            High_effect_homo = High_effect_homo[1]- High_effect_homo[2],
            Load_homo_filtered = Load_homo_filtered[1]- Load_homo_filtered[2],
            Load_homozygous = Load_homozygous[1]- Load_homozygous[2],
            HF_GERP = HF_GERP[1]- HF_GERP[2])




# Calculate P values to include in the figures
Grantham_homo_T <- t.test(Pairwise_summary$Grantham_homo, mu = 0, alternative = "greater")
Load_homo_T <- t.test(Pairwise_summary$Load_homozygous, mu = 0, alternative = "greater")
High_effect_homo_T <- t.test(Pairwise_summary$High_effect_homo, mu = 0, alternative = "greater")
HF_GERP_T <- t.test(Pairwise_summary$HF_GERP, mu = 0, alternative = "greater")
Load_homo_filtered_T <- t.test(Pairwise_summary$Load_homo_filtered, mu = 0, alternative = "greater")



#Pairwise summary of GERP load
Pairwise_GERP_load_expressed <-Pairwise_summary %>%
  ggplot(aes(x = Load_homozygous)) +
  geom_histogram(bins = 10, fill = "black")+
  theme_classic() +My_Theme +
  xlab("Expressed GERP load")+
  ylab("")+
  #coord_cartesian(xlim=c(-2,4), ylim = c(0,5))+
  geom_text(aes(x = Inf, y = Inf, label = paste("p =", signif(Load_homo_T$p.value, 3))),
            hjust = 1, vjust = 1, size = 5) +
  theme(plot.margin = margin(0,0,1.5,0, "cm"))



#Pairwise summary of GERP load
Filtered_GERP_load_expressed <-Pairwise_summary %>%
  ggplot(aes(x = Load_homo_filtered)) +
  geom_histogram(bins = 10, fill = "black")+
  theme_classic() +My_Theme +
  xlab("Expressed Filtered GERP load")+
  ylab("")+
  #coord_cartesian(xlim=c(-2,4), ylim = c(0,5))+
  geom_text(aes(x = Inf, y = Inf, label = paste("p =", signif(Load_homo_filtered_T$p.value, 3))),
            hjust = 1, vjust = 1, size = 5) +
  theme(plot.margin = margin(0,0,1.5,0, "cm"))


#Pairwise summary of high effect GERP
HF_GERP_load <-Pairwise_summary %>%
  ggplot(aes(x = HF_GERP)) +
  geom_histogram(bins = 10, fill = "black")+
  theme_classic() +My_Theme +
  xlab("High Frequency GERP sites")+
  ylab("")+
  #coord_cartesian(xlim=c(-2,4), ylim = c(0,5))+
  geom_text(aes(x = Inf, y = Inf, label = paste("p =", signif(HF_GERP_T$p.value, 3))),
            hjust = 1, vjust = 1, size = 5) +
  theme(plot.margin = margin(0,0,1.5,0, "cm"))



#Pairwise summary of High Effect load
Pairwise_Load_High_expressed <-Pairwise_summary %>%
  ggplot(aes(x = High_effect_homo)) +
  geom_histogram(bins = 10 , fill = "black")+
  theme_classic() +My_Theme +
  xlab(" Expressed High Effect Load")+
  ylab("")+
  #coord_cartesian(xlim=c(-2,4), ylim = c(0,5))+
  geom_text(aes(x = Inf, y = Inf, label = paste("p =", signif(High_effect_homo_T$p.value, 3))),
            hjust = 1, vjust = 1, size = 5)+
  theme(plot.margin = margin(0,0,1.5,0, "cm"))


#Pairwise summary of Grantham load
Pairwise_Grantham_Load_expressed <-Pairwise_summary %>%
  ggplot(aes(x = Grantham_homo)) +
  geom_histogram(bins = 15 , fill = "black")+
  theme_classic() +My_Theme +
  xlab(" Expressed Grantham Load")+
  ylab("")+
  #coord_cartesian(xlim=c(-2,4), ylim = c(0,5))+
  geom_text(aes(x = Inf, y = Inf, label = paste("p =", signif(Grantham_homo_T$p.value, 3))),
            hjust = 1, vjust = 1, size = 5) +
  theme(plot.margin = margin(0,0,1.5,0, "cm"))


## combine pairwise plots into a single figure
Pairwise_plots_combined_supp <-( Pairwise_GERP_load_expressed + Pairwise_Load_High_expressed + Pairwise_Grantham_Load_expressed)/
                                   (HF_GERP_load+Filtered_GERP_load_expressed) + plot_annotation(tag_levels = 'A') & 
  theme(plot.tag.position  = c(.16, .96), plot.tag = element_text(size = 25))

Pairwise_plots_combined_supp

pdf(file = "Pairwise_plots_supplemental_homo.pdf", width = 12, height = 6, onefile = FALSE, paper = "special")
Pairwise_plots_combined_supp
dev.off()




# calculate R squared to use add to figures
#de we want to do load filtered and homo load filtered in the supplement?
R2_Homo_Load <- r.squaredGLMM(Homo_Load)
R2_Load_filtered <- r.squaredGLMM(Filtered_Load)
R2_Only_High_homo <- r.squaredGLMM(Only_High_homo)
R2_Grantham_homo <- r.squaredGLMM(Grantham_homo)
R2_HF_GERP <-  r.squaredGLMM(HF_GERP)

## high frequency gerp load vs fitness
HF_GERP_linear_plot <- Scaled_summary %>%
  ggplot(aes(x= HF_GERP, y= Fitness)) +
  geom_point()+
  geom_smooth(method = "lm", fill = "springgreen3") +
  xlab("HF GERP sites") + ylim(3.5, 8.5)+
  theme_classic() +My_Theme +
  geom_text(aes(x = Inf, y = Inf, label = paste("R²m =", signif(digits = 4,R2_HF_GERP[1,1]))),
            hjust = 1, vjust = 1.5, size = 6) +
  geom_text(aes(x = Inf, y = Inf, label = paste("R²c =", signif(digits = 4, R2_HF_GERP[1,2]))),
            hjust = 1, vjust = 3, size = 6)


## gerp load with coding sites filtered out vs fitness
Filtered_GERP_linear_plot <- Scaled_summary %>%
  ggplot(aes(x= Load, y= Fitness)) +
  geom_point()+
  geom_smooth(method = "lm", fill = "springgreen3") +
  xlab("Filtered GERP sites") + ylim(3.5, 8.5)+
  theme_classic() +My_Theme +
  geom_text(aes(x = Inf, y = Inf, label = paste("R²m =", signif(digits = 4,R2_Load_filtered[1,1]))),
            hjust = 1, vjust = 1.5, size = 6) +
  geom_text(aes(x = Inf, y = Inf, label = paste("R²c =", signif(digits = 4, R2_Load_filtered[1,2]))),
            hjust = 1, vjust = 3, size = 6)

#GERP load vs fitness
GERP_expressed_linear_plot <- Scaled_summary %>%
  ggplot(aes(x= Load_homozygous, y= Fitness)) +
  geom_point()+
  geom_smooth(method = "lm") +
  xlab("Expressed GERP Load") + ylim(3.5, 8.5)+
  theme_classic() +My_Theme +
  geom_text(aes(x = Inf, y = Inf, label = paste("R²m =", signif(digits = 4,R2_Homo_Load[1,1]))),
            hjust = 1, vjust = 1.5, size = 6) +
  geom_text(aes(x = Inf, y = Inf, label = paste("R²c =", signif(digits = 4, R2_Homo_Load[1,2]))),
            hjust = 1, vjust = 3, size = 6)



#Grantham load vs fitness
Grantham_expressed_linear_plot <- Scaled_summary %>%
  ggplot(aes(x= Grantham_homo, y= Fitness)) +
  geom_point()+
  geom_smooth(method = "lm") +
  xlab("Expressed Grantham Load") + ylab("Fitness")+
  theme_classic() +My_Theme + ylim(3.5, 8.5)+
  geom_text(aes(x = Inf, y = Inf, label = paste("R²m =",  signif(digits = 4, R2_Grantham_homo[1,1]))),
            hjust = 1, vjust = 1.5, size = 6) +
  geom_text(aes(x = Inf, y = Inf, label = paste("R²c =", signif(digits = 4,R2_Grantham_homo[1,2]))),
            hjust = 1, vjust = 3, size = 6)


#High Effect load vs fitness
HighE_Expressed_linear_plot <- Scaled_summary %>%
  ggplot(aes(x= High_effect_homo, y= Fitness)) +
  geom_point()+
  geom_smooth(method = "lm") +
  xlab("Expressed High Effect Load") + ylab("Fitness")+
  theme_classic() +My_Theme + ylim(3.5, 8.5)+
  geom_text(aes(x = Inf, y = Inf, label = paste("R²m =", signif(digits = 4,R2_Only_High_homo[1,1]))),
            hjust = 1, vjust = 1.5, size = 6) +
  geom_text(aes(x = Inf, y = Inf, label = paste("R²c =", signif(digits = 4,R2_Only_High_homo[1,2]))),
            hjust = 1, vjust = 3, size = 6) 

## combine all plots together to make the final figure

Combined_linear_plot_supplemental <- (HighE_Expressed_linear_plot+Grantham_expressed_linear_plot + GERP_expressed_linear_plot)/
(Filtered_GERP_linear_plot+HF_GERP_linear_plot) &
  plot_annotation(tag_levels = 'A') & 
  theme(plot.tag.position  = c(.25, .96), plot.tag = element_text(size = 25))

Combined_linear_plot_supplemental


pdf(file = "Supplement_linear_plot.pdf", width = 12, height = 6, onefile = FALSE, paper = "special")
Combined_linear_plot_supplemental
dev.off()


####################### do again for homozygosity #############################################################
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



######Linear models on GERP bins ######################################

Average_GERP_bins <- GERP_bins %>% 
  group_by(GERP_bin, Landscape) %>%
  #filter(frequency !=0) %>%
  summarise(frequency = mean(frequency^2))

Average_GERP_bins <- Average_GERP_bins %>%
  mutate(Block = case_when(
    Landscape >= 21 & Landscape <= 40 ~ 2,
    Landscape >= 41 & Landscape <= 60 ~ 3,
    Landscape >= 61 & Landscape <= 80 ~ 4
  ))



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


# Now save the plot as a pdf of the desired dimensions for a publication quality figure





##### Linear models of VEP bins
Average_VEP_bins <- VEP_bins %>% 
  group_by(VEP_bin, Landscape, Location) %>%
  #filter(frequency !=0) %>%
  summarise(frequency = mean(frequency^2))

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



########## combine plots 

GERP_bins_plot_homo <- ggplot(GERP_results, aes(x=Bin, y=mean, color = Location)) + 
  geom_point(shape = 19, size = 4, position = position_dodge(0.9)) + 
  geom_point(data = Average_GERP_bins, aes(y = frequency), shape = 21, size = 2, 
             position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.9),
             stroke = 0.5) +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=0.5, position=position_dodge(0.9)) +
  theme_classic() + My_Theme + theme(legend.position = "none") +
  xlab("GERP Score")+  scale_color_manual(values=c("orangered4", "#56B4E9"))+
  ylab("Average Proportion of Homozygous Individuals")

VEP_bins_plot_homo <- ggplot(VEP_results, aes(x=VEP_bin, y=mean, color = Location)) + 
  geom_point(shape = 19, size = 4, position = position_dodge(0.9)) + 
  geom_point(data = Average_VEP_bins, aes(y = frequency), shape = 21, size = 2, 
             position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.9),
             stroke = 0.5) +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=0.5, position=position_dodge(0.9)) +
  theme_classic() + My_Theme +
  xlab("VEP Category")+   scale_color_manual(values=c("orangered4", "#56B4E9")) +
  ylab("")

Combined_averages_homo_plots <- (GERP_bins_plot_homo + VEP_bins_plot_homo) + plot_annotation(tag_levels = 'A') & 
  theme(plot.tag.position  = c(.3, .96), plot.tag = element_text(size = 25))
Combined_averages_homo_plots

pdf(file = "supp_Combined_averages_homo_plots.pdf", width = 12, height = 6, onefile = FALSE, paper = "special")
Combined_averages_homo_plots
dev.off()



##################################### fixed sites for benign alleles#######################


## get the total number of sites fixed for the derived allele in each landscape and bin

GERP_bins_fixed_benign <- GERP_bins %>%
  group_by(Landscape, GERP_bin) %>%
  dplyr::summarise(fixed = sum(frequency == 0)/n()) %>%
  mutate(Block = case_when(
    Landscape >= 21 & Landscape <= 40 ~ 2,
    Landscape >= 41 & Landscape <= 60 ~ 3,
    Landscape >= 61 & Landscape <= 80 ~ 4
  ))


## split the column so that the Bin identifier and the core vs edge identifier are in separate columns
GERP_bins_fixed_benign[c('Bin', 'Location')] <- str_split_fixed(GERP_bins_fixed_benign$GERP_bin, " ", n = 2)


## graph the data so that each point represents a landscape and bin 
##and use the mean and 95% confidence intervals from the model
GERP_fixed_plot_benign <- ggplot(GERP_bins_fixed_benign, aes(x=Bin, y=fixed, color = Location)) + 
  geom_boxplot()+
  theme_classic() + My_Theme +
  xlab("GERP Score")+  theme(legend.position = "none") +
  ylab("Proportion of Fixed sites") + scale_color_manual(values=c("orangered4", "#56B4E9"))+ylim(.1,.8)
GERP_fixed_plot_benign

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
VEP_bins_fixed_benign <- VEP_bins %>%
  group_by(Location, Landscape, VEP_bin) %>%
  dplyr::summarise(fixed = sum(frequency == 0)/n()) %>%
  mutate(Block = case_when(
    Landscape >= 21 & Landscape <= 40 ~ 2,
    Landscape >= 41 & Landscape <= 60 ~ 3,
    Landscape >= 61 & Landscape <= 80 ~ 4
  ))


VEP_bins_fixed_benign$VEP_bin <- factor(VEP_bins_fixed_benign$VEP_bin, levels = c("Low", "Moderate", "High"))

VEP_fixed_plot_benign <- ggplot(VEP_bins_fixed_benign, aes(x=VEP_bin, y=fixed, color = Location)) + 
  geom_boxplot()+
  theme_classic() + My_Theme + 
  xlab("VEP Category")+
  ylab("") + scale_color_manual(values=c("orangered4", "#56B4E9")) +ylim(.1,.8)

VEP_fixed_plot_benign
##### Combine fixed plots into 1


Combined_fixed_plots_benign <- GERP_fixed_plot_benign + VEP_fixed_plot_benign+ plot_annotation(tag_levels = 'A') & 
  theme(plot.tag.position  = c(.25, .96), plot.tag = element_text(size = 25))
Combined_fixed_plots_benign


pdf(file = "Supplemental_fixed_benign.pdf", width = 12, height = 6, onefile = FALSE, paper = "special")
Combined_fixed_plots_benign
dev.off()



############ correlation plots ####################


## correlation between Grantham Laod and High Effect Load
cor1 <-paste("p =", signif(cor(Final_summary$Grantham_load, Final_summary$High_effect_alleles), 3))

cor1_plot <- Final_summary %>%
  ggplot(aes(x= Grantham_load, y= High_effect_alleles))+
  geom_point()+My_Theme + theme_bw() +
  xlab("Grantham Load")+ ylab("High Effect Load") +
  geom_text(aes(x = Inf, y = Inf, label = cor1),
            hjust = 1.2, vjust = 1.2, size = 5)+ My_Theme

## correlation between Grantham Laod and GERP Load
cor2 <-paste("p =", signif(cor(Final_summary$Grantham_load, Final_summary$Load), 3))

cor2_plot <- Final_summary %>%
  ggplot(aes(x= Grantham_load, y= Load))+
  geom_point()+My_Theme + theme_bw() +
  xlab("Grantham Load")+ ylab("GERP Load") +
  geom_text(aes(x = Inf, y = Inf, label = cor2),
            hjust = 1.2, vjust = 1.2, size = 5) + My_Theme

## correlation between GERP Load and High Effect Load
cor3 <-paste("p =", signif(cor(Final_summary$High_effect_alleles, Final_summary$Load), 3))

cor3_plot <- Final_summary %>%
  ggplot(aes(x= High_effect_alleles, y= Load))+
  geom_point()+My_Theme + theme_bw() +
  xlab("High Effect Load")+ ylab("GERP Load") +
  geom_text(aes(x = Inf, y = Inf, label = cor3),
            hjust = 1.2, vjust = 1.2, size = 5) + My_Theme

## combine correlatuion plots into a single figure
cor_plots <-(cor1_plot + cor2_plot+ cor3_plot) + plot_annotation(tag_levels = 'A') & 
  theme(plot.tag.position  = c(.25, .96), plot.tag = element_text(size = 25))

cor_plots
pdf(file = "cor_plots.pdf", width = 12, height = 6, onefile = FALSE, paper = "special")
cor_plots
dev.off()


### diversity correlations

## correlation between GERP Load and High Effect Load
cor4 <-paste("p =", signif(cor(Final_summary$HF_GERP, Final_summary$Diversity), 3))

Final_summary %>%
  ggplot(aes(x= HF_GERP, y= Diversity))+
  geom_point()+My_Theme + theme_bw() +
  xlab("High Frequency GERP")+ ylab("Diversity") +
  geom_text(aes(x = Inf, y = Inf, label = cor4),
            hjust = 1, vjust = 1, size = 5)



########## diversity plot for core vs edge



#diversity vs fitness
Diversity_linear_plot_bytreatmnet <- Scaled_summary %>%
  ggplot(aes(x= Diversity, y= Fitness, color = treatment)) +
  geom_point()+
  geom_smooth(method = "lm", fill = "springgreen3") +
  xlab("Diversity") + ylab("")+
  theme_classic() +My_Theme + ylim(3.5, 8.5)


Only_Diversity_edge <- lmer(Fitness ~ Diversity + treatment+ (1 | Block)  , data = Scaled_summary)
summary(Only_Diversity_edge)
