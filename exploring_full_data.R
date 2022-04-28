#### Lets look and see if RS score is a ssociated with less allele frequency chagne.
GERPvsallele <-read.table("Full_filtered_data.delim", header = TRUE, sep = " ")
GERPvsallele <- GERPvsallele %>%
  select(-Ref)


useful_sites <- GERPvsallele %>%
  filter(RS > 2) 

ggplot(GERPvsallele, aes(x= RS, y= Frequency_change)) +
  geom_point()+
  geom_smooth(method = lm) +
  ylab("Allele Frequency Change") +
  xlab("RS")



conserved_sites <- GERPvsallele %>%
  filter(RS > 2) %>%
  group_by(Chromosome , Location, Treatment) %>%
  summarise(Loc_FRQ_change = mean(Frequncy_change, na.rm = TRUE), RS = mean(RS))
conserved_sites$CONS <- rep(0 , nrow(conserved_sites))

unconserved_sites <- GERPvsallele %>%
  filter(RS < .01) %>%
  filter(Max_RS > 0) %>%
  group_by(Chromosome , Location, Treatment) %>%
  summarise(Loc_FRQ_change = mean(Frequncy_change, na.rm = TRUE), RS = mean(RS))
unconserved_sites$CONS <- rep(1 , nrow(unconserved_sites))

both_sites <- rbind(conserved_sites, unconserved_sites)
both_sites$CONS <- as.character(both_sites$CONS)

ggplot(both_sites, aes(y = Loc_FRQ_change, x = CONS)) +
  geom_boxplot() +
theme_classic() +
  stat_compare_means(label.x = 1.5, label.y = 1)+
  stat_compare_means(method = "t.test",label.x = 2, label.y = 1 )



###Lets look at allele frequency changes for high GERP score sites and see if
##they are associated with a reduction in fitness. 

test_fitness <- GERPvsallele %>%
  filter(Element_score > 4) %>%
  group_by(Landscape, Treatment) %>%
  summarise(mean_change = mean(Frequency_change), Fitness = mean(Fitness))

test_fitness_unconserved <- GERPvsallele %>%
  filter(is.na(Element_score)) %>%
  group_by(Landscape, Treatment) %>%
  summarise(mean_change = mean(Frequency_change), Fitness = mean(Fitness))



ggplot(test_fitness, aes(y = Fitness, x = mean_change, color = Treatment)) +
  geom_point()+
  geom_smooth(method = lm) +
  ylab("Fitness") +
  xlab("Frq Change")




### lets look at starting allele frequency and see how that effects change
## We can also see how GERP scores affect starting frequency.
## Purifying selection should reduce starting allele frequencys.

#startin_frqs <- read.table("GERP_Allele_combined.delim", header = TRUE, sep = " ")

startin_frqs <- Allele_frequencies

startin_frqs <- startin_frqs %>%
  filter(X0_NA != 1)

startin_frqs %>%
  filter(Max_RS > 0)%>%
  filter(RS > 2) %>%
  ggplot(aes(y = X0_NA, x= RS))+
  geom_point() +
  geom_smooth(method = lm)

## compare staring FRQs of RS sites >2 to sites around 0

## how does starting FRQ affect change

startin_frqs %>%
  filter(`0_NA` != 1) %>%
  ggplot(aes(x = `0_NA`, y= Change_stationary))+
  geom_point()+ geom_smooth()



## Join with fitness
fit_allele_test <- right_join(beetle_fitness , Allele_frequencies_test)

fit_allele_test <- fit_allele_test %>%
  group_by(Landscape, Treatment) %>%
  summarise(mean_change = mean(Frequency_change), Fitness = mean(Fitness))

fit_allele_test_long <- fit_allele_test %>%
  pivot_longer(cols = c(Change_core, Change_stationary, Change_edge),
  names_to = "type") %>%
  filter(!is.na(value))

fit_allele_test_sum <- fit_allele_test_long %>%
  group_by(Landscape, type) %>%
  summarise(fit_mean = mean(Fitness), change_mean = mean(value))

ggplot(fit_allele_test_sum, aes(x= change_mean, y= fit_mean, color = type)) +
  geom_point()+
  geom_smooth(method = lm)




## Lets look at the starting allele frequencies by type this should be the same
fit_allele_test_long %>%
  filter(`0_NA` != 1) %>%
  group_by(Landscape, type) %>%
  summarise(start_mean = mean(`0_NA`), change_mean = mean(value)) %>%
ggplot(aes(y = change_mean, x = start_mean, color = type))+ 
  geom_point()+
  geom_smooth(method = lm)
