library(ape)
library(tidyverse)
## Read in annotation data for tribolium
Trib_annotation <- read.gff("Tribolium_annotation.gff")
#Trib_annotation <- head(Trib_annotation, 10000)
## Read in final filtered data to compare 
Full_filtered_data <- read.table("Full_filtered_data.delim", header = T)
## Remove refseq rows since they contain the whole ganome
Trib_annotation <- Trib_annotation %>%
  filter(type != "region")


## We need to find out for each of our GERP sites whether they fall inside an exon or not

## Remove duplicates from the annotation
Trib_annotation_unique <- Trib_annotation[!duplicated(Trib_annotation[c('seqid', 'start', 'end')]),]

## Change to chacacter for comparison
Trib_annotation_unique$seqid <- as.character(Trib_annotation_unique$seqid)

## Take out all locations in the annotation that arent exons
Trib_annotation_unique <- Trib_annotation_unique %>%
  filter(type == "exon")
## make a new row to put results into
Full_filtered_data$Exon <- rep(NA, nrow(Full_filtered_data))

for(i in 1:nrow(Full_filtered_data)){
  # Find if a row in our data matches with an exon from the refrence.
  GERP_row <- which(Trib_annotation_unique$seqid == Full_filtered_data$Chromosome[i] &
                      Trib_annotation_unique$start <= Full_filtered_data$Location[i] &
                      Trib_annotation_unique$end >= Full_filtered_data$Location[i])
# If its matches an exon put yes if it doesnt put no into the new column
  if(is_empty(GERP_row)){Yes_No <- "No"} else{Yes_No <- "Yes"}

  Full_filtered_data$Exon[i] <- Yes_No

}




## Lets visualize

## remove extra replicates since we are just looking at the RS score for specific sites
Full_filtered_data_unique <- Full_filtered_data[!duplicated(Full_filtered_data[c('Chromosome', 'Location')]),]

## DO we see diffierent ratios of sites appearing in exons when we filter for higher RS scores?
table(Full_filtered_data_unique$Exon)

Full_filtered_data_unique2 <- Full_filtered_data_unique %>%
  filter(RS >= 4)
table(Full_filtered_data_unique2$Exon)

## Do sites that appear in exons have higher RS on average like we would expect?
library(ggpubr)
Full_filtered_data_unique %>%
  filter(Element_score >= 0) %>%
ggplot(aes(x=Exon, y= RS )) +
  geom_boxplot() +
  theme_classic() +
  stat_compare_means(label.x = 1.5, label.y = 8)+
  stat_compare_means(method = "t.test",label.x = 2, label.y = 8 )

detach(package:ggpubr,unload=TRUE)
