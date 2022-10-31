## Flnal population summary data with analyses

library(RColorBrewer)
library(tidyverse) 
library(AICcmodavg)
library(gridExtra)




### Read in final data
Final_summary <- read.delim("Final_summary.delim", sep = "")
##cahnge to minor 0 fold frq to represent a potential load


final_model_directional <- lm(Fitness ~ Load*directional_0fold +treatment, data = Final_summary)
final_model_squared <- lm(Fitness ~ change_0fold_squared*Load + treatment, data = Final_summary)
final_model <- lm(Fitness ~ major_0fold_frq*Load + treatment, data = Final_summary)
change_model <- lm(Fitness ~ change_squared + treatment, data = Final_summary)
final_model_grantham <- lm(Fitness ~ Grantham_load*Load + treatment, data = Final_summary)

summary(change_model)
summary(final_model)
summary(final_model_squared)
summary(final_model_directional)
summary(final_model_grantham)
## check residuals ec
res <- resid(final_model)
plot(fitted(final_model), res)


library(AICcmodavg)
models <- list(final_model, final_model_squared, final_model_directional, final_model_grantham)
mod.names <- c("no_square", "squared", "directional", "grantham")
aictab(cand.set = models, modnames = mod.names)

mean(Scaled_summary$major_0fold_frq)
?scale

# Will be worth making some plots with 95% intervals and non-linear predictions
### for non major frq
NewData <- subset(Final_summary, select = c("treatment", "major_0fold_frq", "Load"))
preds <- as.data.frame(predict(final_model, newdata = NewData, interval = "conf"))

plot(x = Final_summary$Fitness, y = preds$fit, pch = 16, col = as.factor(Final_summary$treatment), ylab = "Predicted values", xlab = "Observed values", main = "")
abline(a = 0, b = 1, lty = 2)
segments(x0 = Final_summary$Fitness, y0 = preds$lwr, y1 = preds$upr, col = as.factor(Final_summary$treatment))
legend("topleft", legend = c("Core", "Edge"), col = c(1,2), pch = 16, lty = 1)



## for directional data


NewData <- subset(Final_summary, select = c("treatment", "directional_0fold", "Load"))
preds <- as.data.frame(predict(final_model_directional, newdata = NewData, interval = "conf"))

plot(x = Final_summary$Fitness, y = preds$fit, pch = 16, col = as.factor(Final_summary$treatment), ylab = "Predicted values", xlab = "Observed values", main = "")
abline(a = 0, b = 1, lty = 2)
segments(x0 = Final_summary$Fitness, y0 = preds$lwr, y1 = preds$upr, col = as.factor(Final_summary$treatment))
legend("topleft", legend = c("Core", "Edge"), col = c(1,2), pch = 16, lty = 1)







### Visualize interactions
interaction_plot <- Final_summary %>%
  ggplot(aes(x= Load, y = major_0fold_frq, color = Fitness))+
  geom_point(size = 6)


interaction_plot + scale_color_gradient(low="red", high="blue")








#### Standardize data so effecct sizes can be dierctly compared


Scaled_summary <- Final_summary %>%
  mutate(Load = c(scale(Load)), 
         change_0fold_squared = c(scale(change_0fold_squared)), 
          major_0fold_frq = c(scale(major_0fold_frq)),
         directional_0fold = c(scale(directional_0fold)),
         change_squared = c(scale(change_squared)),
         Grantham_load = c(scale(Grantham_load)))


Scaled_final_model_squared <- lm(Fitness ~ Load*change_0fold_squared+ treatment, data = Scaled_summary)
Scaled_final_model <- lm(Fitness ~ Load*major_0fold_frq+ treatment, data = Scaled_summary)
Scaled_model_directional <- lm(Fitness ~ Load*directional_0fold +treatment, data = Scaled_summary)
Scaled_model_allchange <- lm(Fitness ~ change_squared, data = Scaled_summary)
Scaled_model_allchange_treatment <- lm(Fitness ~ change_squared + treatment, data = Scaled_summary)
Scaled_grantham_model <- lm(Fitness ~ Load*Grantham_load+ treatment, data = Scaled_summary)

# Look at AIC of models

models <- list(Scaled_final_model, Scaled_final_model_squared, Scaled_model_directional, Scaled_model_allchange, Scaled_model_allchange_treatment, Scaled_grantham_model)
mod.names <- c( "no_square", "squared", "directional", "allchange", "allchange_treatment", "Grantham")
aictab(cand.set = models, modnames = mod.names)


### look at fina model
summary(Scaled_final_model)
summary(Scaled_final_model_squared)
summary(Scaled_model_directional)
summary(Scaled_model_allchange)
summary(Scaled_model_allchange_treatment)
summary(Scaled_grantham_model)


res <- resid(Scaled_final_model)
plot(fitted(Scaled_final_model), res)




### use scaled range to create the dataframe to generate a heatmap


load_range <- seq(from = -2.2, to = 2.2, by = .1)
fold_change_range <- seq(from = -2.2, to = 2.2, by = .1)
major_frq_range <- (seq(from = -2.2, to = 2.2, by = .1))


## make heatmap of 0 fold change squared data vs fitness
data_test <- expand.grid(Load=load_range, major_0fold_frq=fold_change_range)
data_test$treatment <- rep("Change_edge", nrow(data_test))


data_test <- data_test %>%
  mutate(fitness = predict(Scaled_final_model, newdata = data_test))


Scaled_summary_edge <- Scaled_summary %>%
  mutate(scaled_fitness = case_when( treatment =="Change_core" ~ Fitness + .37,
                                     treatment == "Change_edge" ~ Fitness))

grid_point_plot <- ggplot(data_test, aes(Load, major_0fold_frq, fill= fitness)) + 
  geom_tile()+
  geom_point(data = Scaled_summary_edge, aes(x= Load, y = major_0fold_frq, fill = Fitness), shape = 22, size = 5)+
  xlab("GERP Load")+
  ylab("Major allele frequency at 0-fold loci")+
  labs(fill='Fitness')+
  My_Theme


grid_point_plot
### make heatmap of 0 fold major allele frq vs fitness


data_test_2 <- expand.grid(Load=load_range, directional_0fold=major_frq_range)
data_test_2$treatment <- rep("Change_edge", nrow(data_test_2))


data_test_2 <- data_test_2 %>%
  mutate(fitness = predict(Scaled_final_model, newdata = data_test_2))




ggplot(data_test_2, aes(Load,directional_0fold , fill= fitness)) + 
  geom_tile()

### point plot
interaction_plot <- Scaled_summary %>%
  ggplot(aes(x= Load, y = directional_0fold, color = Fitness))+
  geom_point(size = 6)


interaction_plot + scale_color_gradient(low="red", high="blue")


# Will be worth making some plots with 95% intervals and non-linear predictions
### for non squared data
NewData <- subset(Scaled_summary, select = c("treatment", "major_0fold_frq", "Load"))
preds <- as.data.frame(predict(Scaled_final_model, newdata = NewData, interval = "conf"))

plot(x = Scaled_summary$Fitness, y = preds$fit, pch = 16, col = as.factor(Scaled_summary$treatment), ylab = "Predicted values", xlab = "Observed values", main = "")
abline(a = 0, b = 1, lty = 2)
segments(x0 = Scaled_summary$Fitness, y0 = preds$lwr, y1 = preds$upr, col = as.factor(Scaled_summary$treatment))
legend("topleft", legend = c("Core", "Edge"), col = c(1,2), pch = 16, lty = 1)



## for squared data


NewData <- subset(Scaled_summary, select = c("treatment", "change_0fold_squared", "Load"))
preds <- as.data.frame(predict(Scaled_final_model_squared, newdata = NewData, interval = "conf"))

plot(x = Scaled_summary$Fitness, y = preds$fit, pch = 16, col = as.factor(Scaled_summary$treatment), ylab = "Predicted values", xlab = "Observed values", main = "")
abline(a = 0, b = 1, lty = 2)
segments(x0 = Scaled_summary$Fitness, y0 = preds$lwr, y1 = preds$upr, col = as.factor(Scaled_summary$treatment))
legend("topleft", legend = c("Core", "Edge"), col = c(1,2), pch = 16, lty = 1)



## pairwise core vs edge differences in metrics

## edge vs core within a replicate genetic loads


pop_genetic_loads_diff <- pop_genetic_loads %>%
  filter(edge != 0) %>%
  mutate(edgeVcore = edge - core ) %>%
  pivot_longer(cols = c("edge", "core"), names_to = "Location" , values_to = "Load_change") %>%
  select(Landscape, Location, edgeVcore) %>%
  filter(Location == "core")

genetic_loads_diff_plot <- pop_genetic_loads_diff %>%
  ggplot(aes(x= edgeVcore))+
  geom_histogram(bins = 15) +
  xlab("GERP Load Difference")+
  ylab("")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        plot.title = element_text(size=16,face="bold"))


t.test(pop_genetic_loads_diff$edgeVcore, mu = 0, alternative = "greater")


### 0-fold allele frequency pairwise diffs

genetic_loads_diff_0fold <- allele_fold_change_squared %>%
  select( -Fitness)%>%
  pivot_wider(values_from = change_0fold_squared, names_from = treatment)%>%
  mutate(edgeVcore = Change_edge - Change_core )

diff_plot_0fold <- genetic_loads_diff_0fold %>%
  ggplot(aes(x= edgeVcore))+
  geom_histogram(bins = 15) +
  xlab("Squred change for 0-fold loci")+
  ylab("")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        plot.title = element_text(size=16,face="bold"))

diff_plot_0fold

t.test(genetic_loads_diff_0fold$edgeVcore, mu = 0, alternative = "greater")

### all change pairwise diffs

all_change_pairwise_diff <- Final_summary %>%
  select("treatment", "Landscape", "change_squared")%>%
  pivot_wider(values_from = change_squared, names_from = treatment)%>%
  mutate(edgeVcore = Change_edge - Change_core )

allchange_diff_plot <- all_change_pairwise_diff %>%
  ggplot(aes(x= edgeVcore))+
  geom_histogram(bins = 15) +
  ylab("Number of Landscapes") +
  xlab("Squred change for all loci")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        plot.title = element_text(size=16,face="bold"))

allchange_diff_plot

t.test(all_change_pairwise_diff$edgeVcore, mu = 0, alternative = "greater")


### combine all plots 
ggarange
ggarrange(allchange_diff_plot, diff_plot_0fold, genetic_loads_diff_plot, ncol = 3, 
          labels = c("A", "B", "C"), heights = c(.5,.5,.5))


## theme for graphs
My_Theme = theme(   axis.title.x = element_text(size = 16),   
                    axis.title.y = element_text(size = 16),
                    legend.text = element_text(size=15), 
                    axis.text= element_text(size=12))
## 
Final_summary %>%
ggplot(aes(x = major_0fold_frq, y= Fitness, color = treatment))+
  geom_point()+
  geom_smooth(method = "lm") +
  xlab("Major allele frequency at 0-fold loci")+
  theme (legend.position="none")+
   My_Theme
  



### Effect sizes of coefficients with confidence intervals

conf_ints <- confint(Scaled_final_model, c('Load','major_0fold_frq', 'treatmentChange_edge', 'Load:major_0fold_frq'), level=0.95)
model_coefficents <- Scaled_final_model$coefficients[2:5]

effect_sizes <- data.frame(lower = conf_ints[,1],
           upper = conf_ints[,2],
           Slope = model_coefficents,
           Variable = c('Gerp Load','Major 0-fold allele frequency', '0-fold/GERP Interaction', 'Edge Effect'))


ggplot(effect_sizes ,aes(x = Variable,  y=Slope,ymin=lower,ymax=upper))+
  geom_pointrange()  +
  xlab("")+
  geom_hline(yintercept = 0, linetype = "dashed")+
  My_Theme

summary(Scaled_final_model)

model_coefficents


### individual models and coefficient confidence intervals 

GERP_model <- lm(Fitness ~ Load +treatment, data = Final_summary)

summary(GERP_model)

confint(GERP_model, "treatmentChange_edge", level=0.95)

##

Zerof_model <- lm(Fitness ~ major_0fold_frq +treatment, data = Final_summary)

summary(Zerof_model)

confint(Zerof_model, "treatmentChange_edge", level=0.95)

summary(Scaled_final_model)

