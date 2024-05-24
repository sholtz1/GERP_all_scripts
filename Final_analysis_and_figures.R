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



### THeme for GG plots

My_Theme = theme(   axis.title.x = element_text(size = 16),   
                    axis.title.y = element_text(size = 16),
                    legend.text = element_text(size=17), 
                    axis.text= element_text(size=15),
                    legend.title=element_blank())



### Read in final data
Final_summary <- read.delim("Final_summary.delim", sep = "")

##Scale all data for comparison



Scaled_summary <- Final_summary %>%
  mutate(Load = c(scale(Load)), 
         Load_full = c(scale(Load_full)),
         Load_full_filtered = c(scale(Load_full_filtered)),
         High_effect_alleles = c(scale(High_effect_alleles)),
         change_squared = c(scale(change_squared)),
         Grantham_load = c(scale(Grantham_load)),
         Global_change = c(scale(change_global)),
         Grantham_global = c(scale(Grantham_global))) %>%
  select(treatment, Fitness, Load_full, Load_full_filtered, High_effect_alleles, change_squared, 
         Grantham_load, Load, Global_change, Grantham_global, Landscape)

Scaled_summary <- Scaled_summary %>%
  mutate(All_loads = High_effect_alleles+Load+Grantham_load) %>%
  mutate(GERP_high_load = High_effect_alleles+Load)
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
##High_effect_averages = VEP high effect alleles times the allele frequency averaged
## Change squared = all allele frequency changes across the genome.
## Grantham_load = grantham score X allele frequency summed across the population
## Grantham_average = grantham score X allele frequency  averaged across the population
## major_0fold_frq = average of the major allele for all 0 fold sites in a population


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
Only_Load_filtered <- lmer(Fitness ~ Load_full_filtered + (1 | Block)  , data = Scaled_summary)
Only_Grantham <- lmer(Fitness ~ Grantham_load + (1 | Block)  , data = Scaled_summary)
Only_High <- lmer(Fitness ~ High_effect_alleles + (1 | Block)  , data = Scaled_summary)
Combined_loads <- lmer(Fitness ~ All_loads + (1 | Block)  , data = Scaled_summary)
Combined_GERP_high <- lmer(Fitness ~ GERP_high_load + (1 | Block)  , data = Scaled_summary)
  
#Global_change <- lmer(Fitness ~ Global_change + (1 | Block)  , data = Scaled_summary)

summary(Only_High)


# Look at AIC of models

models <- list(all_interactions, Load_Grantham_interactions, Grantham_High_interactions,
               Load_high_interactions, No_interactions, Only_Load_high_interactions,
               Only_Load_high, Only_Load_Grantham_interactions, Only_Load_Grantham,
               Only_High_Grantham_interactions, Only_High_Grantham, Only_Load, Only_Grantham,
               Only_High)


mod.names <- c(
  "all_interactions", "Load_Grantham_interactions","Grantham_High_interactions","Load_high_interactions",
  "No_interactions","Only_Load_high_interactions","Only_Load_high",
  "Only_Load_Grantham_interactions","Only_Load_Grantham","Only_High_Grantham_interactions",
  "Only_High_Grantham"," Only_Load","Only_Grantham","Only_High"
)

AICctab(models, mnames = mod.names)

summary(Only_Load_Grantham_interactions)
## Get R squared values to add to plots 
R2_Only_Load <- r.squaredGLMM(Only_Load)
R2_Only_High <- r.squaredGLMM(Only_High)
R2_Only_Grantham <- r.squaredGLMM(Only_Grantham)
r.squaredGLMM(Only_Load_high)
r.squaredGLMM(No_interactions)
r.squaredGLMM(all_interactions)
r.squaredGLMM(Only_Load_Grantham)

### Read in full data
Full_filtered_data <- read.table("Grantham_GERP_Full.delim", header = TRUE)


##################################Pairwise comparisons ############################################
Pairwise_summary <- Scaled_summary %>%
  filter(Landscape != 75) %>%
  group_by(Landscape) %>%
  summarise(Grantham_load = Grantham_load[1]-Grantham_load[2], 
            Grantham_global = Grantham_global[1]-Grantham_global[2], 
            Load = Load[1]- Load[2], 
            Global_change = Global_change[1] - Global_change[2],
            High_effect_alleles = High_effect_alleles[1]- High_effect_alleles[2],
            Load_full = Load_full[1]- Load_full[2]) %>%
  mutate(All_loads = High_effect_alleles+Load+Grantham_load) %>%
  mutate(GERP_high_load = High_effect_alleles+Load) %>%
  mutate(Grantham_high_load = High_effect_alleles+Grantham_load) %>%
  mutate(Grantham_GERP_load = Load+Grantham_load)

# Calculate P vaues to include in the figures
Grantham_load_T <- t.test(Pairwise_summary$Grantham_load, mu = 0, alternative = "greater")
Load_full_T <- t.test(Pairwise_summary$Load_full, mu = 0, alternative = "greater")
Global_change_T <- t.test(Pairwise_summary$Global_change, mu = 0, alternative = "greater")
High_effect_alleles_T <- t.test(Pairwise_summary$High_effect_alleles, mu = 0, alternative = "greater")
All_loads_T <- t.test(Pairwise_summary$All_loads, mu = 0, alternative = "greater")
GERP_high_loads_T <- t.test(Pairwise_summary$GERP_high_load, mu = 0, alternative = "greater")
Grantham_GERP_load_T <- t.test(Pairwise_summary$Grantham_GERP_load, mu = 0, alternative = "greater")
Grantham_high_load_T <- t.test(Pairwise_summary$Grantham_high_load, mu = 0, alternative = "greater")


#Pairwise summary for Grantham_high_load
Pairwise_Grantham_high_load <- Pairwise_summary %>%
  ggplot(aes(x = Grantham_high_load)) +
  geom_histogram(bins= 20) +
  theme_classic() + My_Theme + 
  xlab("High Effect Load+Grantham load")+
  ylab("")+
  ylim(0,4) +xlim(-4,4)+
  geom_text(aes(x = Inf, y = Inf, label = paste("p =", signif(Grantham_high_load_T$p.value, 3))),
            hjust = 1, vjust = 1, size = 5) + theme(axis.title.x = element_text(size = 12))


#Pairwise summary for Grantham_GERP_load
Pairwise_Grantham_GERP_load <- Pairwise_summary %>%
  ggplot(aes(x = Grantham_GERP_load)) +
  geom_histogram(bins= 20) +
  theme_classic() + My_Theme + 
  xlab("GERP Load+Grantham load")+
  ylab("")+
  ylim(0,4) + xlim(-4,4)+
  geom_text(aes(x = Inf, y = Inf, label = paste("p =", signif(Grantham_GERP_load_T$p.value, 3))),
            hjust = 1, vjust = 1, size = 5)  + theme(axis.title.x = element_text(size = 12))

Pairwise_Grantham_GERP_load
#Pairwise summary for all loads
Pairwise_All_loads <- Pairwise_summary %>%
  ggplot(aes(x = All_loads)) +
  geom_histogram(bins= 20) +
  theme_classic() + My_Theme + 
  xlab("All Loads")+
  ylab("")+
  ylim(0,4) + xlim(-4,4)+
  geom_text(aes(x = Inf, y = Inf, label = paste("p =", signif(All_loads_T$p.value, 3))),
hjust = 1, vjust = 1, size = 5)  + theme(axis.title.x = element_text(size = 12))



#Pairwise summary for GERP and highE loads
Pairwise_GERP_High <- Pairwise_summary %>%
  ggplot(aes(x = GERP_high_load)) +
  geom_histogram(bins= 20) +
  theme_classic() + My_Theme + 
  xlab("GERP load+High Effect Load")+
  ylab("")+
  ylim(0,4) + xlim(-4,4)+
  geom_text(aes(x = Inf, y = Inf, label = paste("*p =", signif(GERP_high_loads_T$p.value, 3))),
            hjust = 1, vjust = 1, size = 5)  + theme(axis.title.x = element_text(size = 12))



#Pairwise summary of GERP load
Pairwise_GERP_load <-Pairwise_summary %>%
  ggplot(aes(x = Load_full)) +
  geom_histogram(bins = 20, fill = "blue3")+
  theme_classic() +My_Theme +
  xlab("GERP load")+
  ylab("")+
  ylim(0,4) + xlim(-4,4)+
  geom_text(aes(x = Inf, y = Inf, label = paste("p =", signif(Load_full_T$p.value, 3))),
            hjust = 1, vjust = 1, size = 5) +
  theme(plot.margin = margin(0,0,1.5,0, "cm"))


#Pairwise summary of HIgh Effect load
Pairwise_Load_High <-Pairwise_summary %>%
  ggplot(aes(x = High_effect_alleles)) +
  geom_histogram(bins = 20 , fill = "springgreen3")+
  theme_classic() +My_Theme +
  xlab("High Effect Load")+
  ylab("")+
  ylim(0,4) + 
  geom_text(aes(x = Inf, y = Inf, label = paste("p =", signif(High_effect_alleles_T$p.value, 3))),
            hjust = 1, vjust = 1, size = 5)+
  theme(plot.margin = margin(0,0,1.5,0, "cm"))


#Pairwise summary of Grantham load
Pairwise_Grantham_Load <-Pairwise_summary %>%
  ggplot(aes(x = Grantham_load)) +
  geom_histogram(bins = 20 , fill = "plum")+
  theme_classic() +My_Theme +
  xlab("Grantham Load")+
  ylab("")+
  ylim(0,4) + xlim(-4,4)+
  geom_text(aes(x = Inf, y = Inf, label = paste("p =", signif(Grantham_load_T$p.value, 3))),
            hjust = 1, vjust = 1, size = 5) +
  theme(plot.margin = margin(0,0,1.5,0, "cm"))


### combine plots and save as PDF

Pairwise_plots_combined <-( Pairwise_GERP_load + Pairwise_Load_High + Pairwise_Grantham_Load)/
 (Pairwise_GERP_High | Pairwise_All_loads | Pairwise_Grantham_GERP_load | Pairwise_Grantham_high_load
  ) + plot_annotation(tag_levels = 'A') & 
  theme(plot.tag.position  = c(.25, .96), plot.tag = element_text(size = 25))

Pairwise_plots_combined

pdf(file = "Pairwise_plots_combined.pdf", width = 12, height = 6, onefile = FALSE, paper = "special")
Pairwise_plots_combined
dev.off()


############################################# Linear comparisons##################### ##################

Scaled_summary <- Scaled_summary %>%
  mutate(All_loads = High_effect_alleles+Load+Grantham_load) %>%
  mutate(GERP_high_load = High_effect_alleles+Load)





###### Linear comparisons between fitness and load

#GERP
GERP_linear_plot <- Scaled_summary %>%
  ggplot(aes(x= Load_full, y= Fitness)) +
  geom_point()+
  geom_smooth(method = "lm", fill = "blue3") +
  xlab("GERP Load") + ylim(3.5, 8.5)+
  theme_classic() +My_Theme +
  geom_text(aes(x = Inf, y = Inf, label = paste("R²m =", signif(digits = 4,R2_Only_Load[1,1]))),
            hjust = 1, vjust = 1.5, size = 6) +
  geom_text(aes(x = Inf, y = Inf, label = paste("R²c =", signif(digits = 4, R2_Only_Load[1,2]))),
            hjust = 1, vjust = 3, size = 6)

GERP_linear_plot


#Grantham
Grantham_linear_plot <- Scaled_summary %>%
  ggplot(aes(x= Grantham_load, y= Fitness)) +
  geom_point()+
  geom_smooth(method = "lm",  fill = "plum") +
  xlab("Grantham Load") + ylab("")+
  theme_classic() +My_Theme + ylim(3.5, 8.5)+
  geom_text(aes(x = Inf, y = Inf, label = paste("R²m =",  signif(digits = 4, R2_Only_Grantham[1,1]))),
            hjust = 1, vjust = 1.5, size = 6) +
  geom_text(aes(x = Inf, y = Inf, label = paste("R²c =", signif(digits = 4,R2_Only_Grantham[1,2]))),
            hjust = 1, vjust = 3, size = 6)

Grantham_linear_plot

#High Effect
HighE_linear_plot <- Scaled_summary %>%
  ggplot(aes(x= High_effect_alleles, y= Fitness)) +
  geom_point()+
  geom_smooth(method = "lm", fill = "springgreen3") +
  xlab("High Effect Load") + ylab("")+
  theme_classic() +My_Theme + ylim(3.5, 8.5)+
  geom_text(aes(x = Inf, y = Inf, label = paste("R²m =", signif(digits = 4,R2_Only_High[1,1]))),
            hjust = 1, vjust = 1.5, size = 6) +
  geom_text(aes(x = Inf, y = Inf, label = paste("R²c =", signif(digits = 4,R2_Only_High[1,2]))),
            hjust = 1, vjust = 3, size = 6) 
  

HighE_linear_plot

Combined_linear_plot <-  GERP_linear_plot + HighE_linear_plot + Grantham_linear_plot &
  plot_annotation(tag_levels = 'A') & 
  theme(plot.tag.position  = c(.25, .96), plot.tag = element_text(size = 25))

Combined_linear_plot


pdf(file = "Combined_linear_plot.pdf", width = 12, height = 6, onefile = FALSE, paper = "special")
Combined_linear_plot
dev.off()

#Scaled_summary %>%
 # ggplot(aes(x= Global_change, y= Fitness)) +
  #geom_point()+
  #geom_smooth(method = "lm") +
  #xlab("All change assuming stabalizing selection") +
  #theme_classic() +My_Theme

#Scaled_summary %>%
 # mutate( All_loads = High_effect_alleles+Load+Grantham_load) %>%
  #ggplot(aes(x= All_loads, y= Fitness)) +
  #geom_point()+
  #geom_smooth(method = "lm") +
  #xlab("All Loads") +
  #theme_classic() +My_Theme










################# Split data into bins for analysis  ################################################

### Make GERP bins####

GERP_bins <- Full_filtered_data %>%
  filter(RS > 0) %>%
  filter(!is.na(X8_C)) %>%
  filter(!is.na(X8_E)) %>%
  select(Chromosome, Location, Landscape ,X8_C, X8_E, RS) %>%
  pivot_longer(cols = c(X8_C, X8_E), names_to = "Treatment", values_to = "frequency") %>%
  mutate(frequency = 1- frequency)

GERP_bins <- GERP_bins %>%
  mutate( GERP_bin = case_when(
    RS > 0 & RS <=2 & Treatment == "X8_E" ~ "0-2 Edge" ,
    RS > 0 & RS <=2 & Treatment == "X8_C" ~ "0-2 Core",
    RS > 2 & RS <=4 & Treatment == "X8_E" ~ "2-4 Edge",
    RS > 2 & RS <=4 & Treatment == "X8_C" ~ "2-4 Core",
    RS > 4 & Treatment == "X8_E" ~ "4+ Edge",
    RS > 4 & Treatment == "X8_C" ~ "4+ Core",
  ))


GERP_bins <- GERP_bins %>%
  unique()

GERP_bins$GERP_bin <- as.character(GERP_bins$GERP_bin)


######################### Split VEP results into bins

VEP_bins <- Full_filtered_data %>%
  filter(!is.na(IMPACT)) %>%
  filter(!is.na(X8_C)) %>%
  filter(!is.na(X8_E)) %>%
  select(Chromosome, Location, Landscape ,X8_C, X8_E, IMPACT) %>%
  pivot_longer(cols = c(X8_C, X8_E), names_to = "Treatment", values_to = "frequency") %>%
  mutate(frequency = 1- frequency) 

VEP_bins <- VEP_bins %>%
  mutate( VEP_bin = case_when(
    IMPACT == "MODIFIER" & Treatment == "X8_E" ~ "Modifier Edge",
    IMPACT == "MODIFIER" & Treatment == "X8_C" ~ "Modifier Core",
    IMPACT == "LOW" & Treatment == "X8_C" ~ "Low Core",
    IMPACT == "LOW" & Treatment == "X8_E" ~ "Low Edge",
    IMPACT == "MODERATE" & Treatment == "X8_C" ~ "Moderate Core",
    IMPACT == "MODERATE" & Treatment == "X8_E" ~ "Moderate Edge",
    IMPACT == "HIGH" & Treatment == "X8_C" ~ "High Core",
    IMPACT == "HIGH" & Treatment == "X8_E" ~ "High Edge",
  ))

VEP_bins$VEP_bin <- factor(VEP_bins$VEP_bin,levels = c("Modifier Core","Modifier Edge",
                                                       "Low Core", "Low Edge", "Moderate Core", 
                                                       "Moderate Edge", "High Core", "High Edge"))




#############################   Linear models of bins and edge effects ##################################

############################# Fixed sites GERP models

GERP_bins_fixed <- GERP_bins %>%
  group_by(Landscape, GERP_bin) %>%
  dplyr::summarise(fixed = sum(frequency == 1)) %>%
  mutate(Block = case_when(
    Landscape >= 21 & Landscape <= 40 ~ 2,
    Landscape >= 41 & Landscape <= 60 ~ 3,
    Landscape >= 61 & Landscape <= 80 ~ 4
  ))

GERP_bins_fixed[c('Bin', 'Location')] <- str_split_fixed(GERP_bins_fixed$GERP_bin, " ", n = 2)

# Fit a generalized linear model with a Poisson distribution to the data
GERP_GLMM <- glmer(fixed ~ Bin*Location + (1 | Block), data = GERP_bins_fixed,
                   family = poisson(link = "log"))
summary(GERP_GLMM)

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
bootpreds <- bootMer(GERP_GLMM, preds, nsim=5000)

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
VEP_bins <- Full_filtered_data %>%
  filter(!is.na(IMPACT)) %>%
  filter(!is.na(X8_C)) %>%
  filter(!is.na(X8_E)) %>%
  select(Chromosome, Location, Landscape ,X8_C, X8_E, IMPACT) %>%
  pivot_longer(cols = c(X8_C, X8_E), names_to = "Treatment", values_to = "frequency") %>%
  mutate(frequency = 1- frequency) 

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

VEP_bins <- VEP_bins%>%
  filter(VEP_bin != "Modifier")

# Updated code to capture the 0 data
VEP_bins_fixed <- VEP_bins %>%
  group_by(Location, Landscape, VEP_bin) %>%
  dplyr::summarise(fixed = sum(frequency == 1)) %>%
  mutate(Block = case_when(
    Landscape >= 21 & Landscape <= 40 ~ 2,
    Landscape >= 41 & Landscape <= 60 ~ 3,
    Landscape >= 61 & Landscape <= 80 ~ 4
  ))



VEP_GLMM <- glmer(fixed ~ VEP_bin*Location + (1 | Block), data = VEP_bins_fixed,
                  family = poisson(link = "log"))
summary(VEP_GLMM)

# Now use the same process to generate confidence intervals
newdat <- expand.grid(VEP_bin = c("Low", "Moderate", "High"), Location = c("Core", "Edge"))
bootpreds <- bootMer(VEP_GLMM, preds, nsim=5000)
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
  theme_classic() + My_Theme + ylim(0,38) +
  xlab("VEP Category")+
  ylab("") + scale_color_manual(values=c("orangered4", "#56B4E9"))

VEP_fixed_plot 
##### Combine fixed plots into 1

 


Combined_fixed_plots <- GERP_fixed_plot + VEP_fixed_plot + plot_annotation(tag_levels = 'A') & 
  theme(plot.tag.position  = c(.25, .96), plot.tag = element_text(size = 25))
Combined_fixed_plots 

pdf(file = "Combined_fixed_plots.pdf", width = 9, height = 6, onefile = FALSE, paper = "special")
Combined_fixed_plots
dev.off()



####################### Linear models on GERP bins ######################################

Average_GERP_bins <- GERP_bins %>% 
  group_by(GERP_bin, Landscape) %>%
  #filter(frequency !=0) %>%
  summarise(frequency = mean(frequency))

Average_GERP_bins <- Average_GERP_bins %>%
  mutate(Block = case_when(
    Landscape >= 21 & Landscape <= 40 ~ 2,
    Landscape >= 41 & Landscape <= 60 ~ 3,
    Landscape >= 61 & Landscape <= 80 ~ 4
  ))



Average_GERP_bins[c('Bin', 'Location')] <- str_split_fixed(Average_GERP_bins$GERP_bin, " ", n = 2)



GERP_bins_model<- lmer(frequency ~ Bin*Location + (1 | Block), data = Average_GERP_bins)

summary(GERP_bins_model)

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
bootpreds <- bootMer(GERP_bins_model, preds, nsim=5000)

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
  ylab("Average Minor Allele Frequency")


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

# Now use the same process to generate confidence intervals
newdat <- expand.grid(VEP_bin = c("Low", "Moderate", "High"), Location = c("Core", "Edge"))
bootpreds <- bootMer(VEP_bins_model, preds, nsim=5000)
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

GERP_bins_plot <- ggplot(GERP_results, aes(x=Bin, y=mean, color = Location)) + 
  geom_point(shape = 19, size = 4, position = position_dodge(0.9)) + 
  geom_point(data = Average_GERP_bins, aes(y = frequency), shape = 21, size = 2, 
             position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.9),
             stroke = 0.5) +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=0.5, position=position_dodge(0.9)) +
  theme_classic() + My_Theme + theme(legend.position = "none") +
  xlab("GERP Score")+  scale_color_manual(values=c("orangered4", "#56B4E9"))+
  ylab("Average Minor Allele Frequency") + ylim(0.07, 0.24)

VEP_bins_plot <- ggplot(VEP_results, aes(x=VEP_bin, y=mean, color = Location)) + 
  geom_point(shape = 19, size = 4, position = position_dodge(0.9)) + 
  geom_point(data = Average_VEP_bins, aes(y = frequency), shape = 21, size = 2, 
             position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.9),
             stroke = 0.5) +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=0.5, position=position_dodge(0.9)) +
  theme_classic() + My_Theme + ylim(0.07, 0.24)+
  xlab("VEP Category")+   scale_color_manual(values=c("orangered4", "#56B4E9")) +
  ylab("")


Combined_averages_plots <- GERP_bins_plot + VEP_bins_plot + plot_annotation(tag_levels = 'A') & 
  theme(plot.tag.position  = c(.3, .96), plot.tag = element_text(size = 25))
Combined_averages_plots

# Now save the plot as a pdf of the desired dimensions for a publication quality figure
pdf(file = "Combined_averages_plots.pdf", width = 10, height = 6, onefile = FALSE, paper = "special")
Combined_averages_plots
dev.off()



#############Fixed for major allele plots#################################################

GERP_bins <- Full_filtered_data %>%
  filter(RS > 0) %>%
  filter(!is.na(X8_C)) %>%
  filter(!is.na(X8_E)) %>%
  select(Chromosome, Location, Landscape ,X8_C, X8_E, RS) %>%
  pivot_longer(cols = c(X8_C, X8_E), names_to = "Treatment", values_to = "frequency") %>%
  mutate(frequency = 1- frequency)

GERP_bins <- GERP_bins %>%
  mutate( GERP_bin = case_when(
    RS > 0 & RS <=2 & Treatment == "X8_E" ~ "0-2 Edge" ,
    RS > 0 & RS <=2 & Treatment == "X8_C" ~ "0-2 Core",
    RS > 2 & RS <=4 & Treatment == "X8_E" ~ "2-4 Edge",
    RS > 2 & RS <=4 & Treatment == "X8_C" ~ "2-4 Core",
    RS > 4 & Treatment == "X8_E" ~ "4+ Edge",
    RS > 4 & Treatment == "X8_C" ~ "4+ Core",
  ))


GERP_bins <- GERP_bins %>%
  unique()

GERP_bins$GERP_bin <- as.character(GERP_bins$GERP_bin)


######################### Split VEP results into bins

VEP_bins <- Full_filtered_data %>%
  filter(!is.na(IMPACT)) %>%
  filter(!is.na(X8_C)) %>%
  filter(!is.na(X8_E)) %>%
  select(Chromosome, Location, Landscape ,X8_C, X8_E, IMPACT) %>%
  pivot_longer(cols = c(X8_C, X8_E), names_to = "Treatment", values_to = "frequency") %>%
  mutate(frequency = 1- frequency) 

VEP_bins <- VEP_bins %>%
  mutate( VEP_bin = case_when(
    IMPACT == "MODIFIER" & Treatment == "X8_E" ~ "Modifier Edge",
    IMPACT == "MODIFIER" & Treatment == "X8_C" ~ "Modifier Core",
    IMPACT == "LOW" & Treatment == "X8_C" ~ "Low Core",
    IMPACT == "LOW" & Treatment == "X8_E" ~ "Low Edge",
    IMPACT == "MODERATE" & Treatment == "X8_C" ~ "Moderate Core",
    IMPACT == "MODERATE" & Treatment == "X8_E" ~ "Moderate Edge",
    IMPACT == "HIGH" & Treatment == "X8_C" ~ "High Core",
    IMPACT == "HIGH" & Treatment == "X8_E" ~ "High Edge",
  ))

VEP_bins$VEP_bin <- factor(VEP_bins$VEP_bin,levels = c("Modifier Core","Modifier Edge",
                                                       "Low Core", "Low Edge", "Moderate Core", 
                                                       "Moderate Edge", "High Core", "High Edge"))




#############################   Linear models of bins and edge effects ##################################

############################# Fixed sites GERP models

GERP_bins_fixed <- GERP_bins %>%
  group_by(Landscape, GERP_bin) %>%
  dplyr::summarise(fixed = (sum(frequency == 0)/(sum(frequency > 0)+ sum(frequency == 0)))) %>%
  mutate(Block = case_when(
    Landscape >= 21 & Landscape <= 40 ~ 2,
    Landscape >= 41 & Landscape <= 60 ~ 3,
    Landscape >= 61 & Landscape <= 80 ~ 4
  ))

GERP_bins_fixed[c('Bin', 'Location')] <- str_split_fixed(GERP_bins_fixed$GERP_bin, " ", n = 2)



GERP_fixed_major_plot <- ggplot(GERP_bins_fixed, aes(x=Bin, y=fixed, color = Location)) + 
  geom_boxplot() +
  theme_classic() + My_Theme +
  xlab("GERP Score")+  theme(legend.position = "none") +
  ylab("Percent of sites fixed for benign allele") + scale_color_manual(values=c("orangered4", "#56B4E9"))

GERP_fixed_major_plot





###############vep fixed major frequency ###################################

Average_VEP_bins <- VEP_bins %>% 
  group_by(VEP_bin, Landscape) %>%
  #filter(frequency !=0) %>%
  summarise(frequency = (sum(frequency == 0)/(sum(frequency > 0)+ sum(frequency == 0))))

Average_VEP_bins <- Average_VEP_bins %>%
  mutate(Block = case_when(
    Landscape >= 21 & Landscape <= 40 ~ 2,
    Landscape >= 41 & Landscape <= 60 ~ 3,
    Landscape >= 61 & Landscape <= 80 ~ 4
  ))

Average_VEP_bins[c('Bin', 'Location')] <- str_split_fixed(Average_VEP_bins$VEP_bin, " ", n = 2)

VEP_Major_fixed_plot <- ggplot(Average_VEP_bins, aes(x=Bin, y=frequency, color = Location)) + 
  geom_boxplot() +
  theme_classic() + My_Theme + 
  scale_x_discrete(limits = c("Low", "Moderate", "High")) +
  xlab("VEP Category")+   scale_color_manual(values=c("orangered4", "#56B4E9")) +
  ylab("")


VEP_Major_fixed_plot 
########## combine plots 

Combined_major_fixed_plots <- GERP_fixed_major_plot + VEP_Major_fixed_plot + plot_annotation(tag_levels = 'A') & 
  theme(plot.tag.position  = c(.3, .96), plot.tag = element_text(size = 20))
Combined_major_fixed_plots

pdf(file = "Combined_major_fixed_plots.pdf", width = 10, height = 6, onefile = FALSE, paper = "special")
Combined_major_fixed_plots
dev.off()
