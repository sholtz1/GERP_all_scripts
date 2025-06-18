library(caroline)
library(tidyverse)
library(ape)

##Read in formatted allele frequencies
allele <- read.table("Derived_FRQs_test.delim", header = TRUE)


##Make names match
Full_filtered_data  <- allele

## read in all gerp information

NC_007416.3_rm.rates.bed <-read.delim("/project/beetlegenes/sholtz1/GERP/analyses/gerp/NC_007416.3_rm_gapped.rates", header = FALSE)
NC_007417.3_rm.rates.bed <-read.delim("/project/beetlegenes/sholtz1/GERP/analyses/gerp/NC_007417.3_rm_gapped.rates", header = FALSE)
NC_007418.3_rm.rates.bed <-read.delim("/project/beetlegenes/sholtz1/GERP/analyses/gerp/NC_007418.3_rm_gapped.rates", header = FALSE)
NC_007419.2_rm.rates.bed <-read.delim("/project/beetlegenes/sholtz1/GERP/analyses/gerp/NC_007419.2_rm_gapped.rates", header = FALSE)
NC_007420.3_rm.rates.bed <-read.delim("/project/beetlegenes/sholtz1/GERP/analyses/gerp/NC_007420.3_rm_gapped.rates", header = FALSE)
NC_007421.3_rm.rates.bed <-read.delim("/project/beetlegenes/sholtz1/GERP/analyses/gerp/NC_007421.3_rm_gapped.rates", header = FALSE)
NC_007422.5_rm.rates.bed <-read.delim("/project/beetlegenes/sholtz1/GERP/analyses/gerp/NC_007422.5_rm_gapped.rates", header = FALSE)
NC_007423.3_rm.rates.bed <-read.delim("/project/beetlegenes/sholtz1/GERP/analyses/gerp/NC_007423.3_rm_gapped.rates", header = FALSE)
NC_007424.3_rm.rates.bed <-read.delim("/project/beetlegenes/sholtz1/GERP/analyses/gerp/NC_007424.3_rm_gapped.rates", header = FALSE)
NC_007425.3_rm.rates.bed <-read.delim("/project/beetlegenes/sholtz1/GERP/analyses/gerp/NC_007425.3_rm_gapped.rates", header = FALSE)




#Join each together to prevent duplicating columns
rates_bed <- rbind(NC_007416.3_rm.rates.bed, NC_007417.3_rm.rates.bed, NC_007418.3_rm.rates.bed,
                   NC_007419.2_rm.rates.bed, NC_007420.3_rm.rates.bed, NC_007421.3_rm.rates.bed,
                   NC_007422.5_rm.rates.bed, NC_007423.3_rm.rates.bed, NC_007424.3_rm.rates.bed,
                   NC_007425.3_rm.rates.bed)



##Change column names for joining
colnames(rates_bed) <- c("Chromosome", "Location", "RS_MAX_RS")

## seperate RS and max RS into two columns
rates_bed <- rates_bed %>%
  separate(RS_MAX_RS, into = c("RS", "MAX_RS"), sep = " ")

## change chromosome names for joining
rates_bed <- rates_bed %>%
  mutate(Chromosome = case_when(
    Chromosome == 16 ~ "NC_007416.3",
    Chromosome == 17 ~ "NC_007417.3",
    Chromosome == 18 ~ "NC_007418.3",
    Chromosome == 19 ~ "NC_007419.2",
    Chromosome == 20 ~ "NC_007420.3",
    Chromosome == 21 ~ "NC_007421.3",
    Chromosome == 22 ~ "NC_007422.5",
    Chromosome == 23 ~ "NC_007423.3",
    Chromosome == 24 ~ "NC_007424.3",
    Chromosome == 25 ~ "NC_007425.3"
  ))

## add 1 to location to start at base 1 currently it is coded to start at base 0


rates_bed$Location <- rates_bed$Location+1

## Join to add GERP information for each locus we have information for
AlleleFreqs <- left_join(Full_filtered_data, rates_bed, by = c("Chromosome", "Location"))


write_delim(AlleleFreqs, "/project/beetlegenes/sholtz1/GERP/Full_data_newgerp.delim")
