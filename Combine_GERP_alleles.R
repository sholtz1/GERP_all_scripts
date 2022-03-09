library(caroline)
library(tidyverse)

# This script will read in both tables and fill in RS_score and P_values from the
#   GERP value to the corresponding row of the allele frequency data.

GERP_7416 <- read.tab("/project/beetlegenes/sholtz1/GERP/analyses/gerp/NC_007416.3_rm.rates.elems.bed", head = FALSE)
GERP_7417 <- read.tab("/project/beetlegenes/sholtz1/GERP/analyses/gerp/NC_007417.3_rm.rates.elems.bed", head = FALSE)
GERP_7418 <- read.tab("/project/beetlegenes/sholtz1/GERP/analyses/gerp/NC_007418.3_rm.rates.elems.bed", head = FALSE)
GERP_7419 <- read.tab("/project/beetlegenes/sholtz1/GERP/analyses/gerp/NC_007419.2_rm.rates.elems.bed", head = FALSE)
GERP_7420 <- read.tab("/project/beetlegenes/sholtz1/GERP/analyses/gerp/NC_007420.3_rm.rates.elems.bed", head = FALSE)
GERP_7421 <- read.tab("/project/beetlegenes/sholtz1/GERP/analyses/gerp/NC_007421.3_rm.rates.elems.bed", head = FALSE)
GERP_7422 <- read.tab("/project/beetlegenes/sholtz1/GERP/analyses/gerp/NC_007422.5_rm.rates.elems.bed", head = FALSE)
GERP_7423 <- read.tab("/project/beetlegenes/sholtz1/GERP/analyses/gerp/NC_007423.3_rm.rates.elems.bed", head = FALSE)
GERP_7424 <- read.tab("/project/beetlegenes/sholtz1/GERP/analyses/gerp/NC_007424.3_rm.rates.elems.bed", head = FALSE)
GERP_7425 <- read.tab("/project/beetlegenes/sholtz1/GERP/analyses/gerp/NC_007425.3_rm.rates.elems.bed", head = FALSE)


colnames(GERP_7416) <- c("Chromosome", "start", "end", "Element_score", "P_value", "","","")
colnames(GERP_7417) <- c("Chromosome", "start", "end", "Element_score", "P_value", "","","")
colnames(GERP_7418) <- c("Chromosome", "start", "end", "Element_score", "P_value", "","","")
colnames(GERP_7419) <- c("Chromosome", "start", "end", "Element_score", "P_value", "","","")
colnames(GERP_7420) <- c("Chromosome", "start", "end", "Element_score", "P_value", "","","")
colnames(GERP_7421) <- c("Chromosome", "start", "end", "Element_score", "P_value", "","","")
colnames(GERP_7422) <- c("Chromosome", "start", "end", "Element_score", "P_value", "","","")
colnames(GERP_7423) <- c("Chromosome", "start", "end", "Element_score", "P_value", "","","")
colnames(GERP_7424) <- c("Chromosome", "start", "end", "Element_score", "P_value", "","","")
colnames(GERP_7425) <- c("Chromosome", "start", "end", "Element_score", "P_value", "","","")

Gerp <- rbind(GERP_7416, GERP_7417, GERP_7418, GERP_7419, GERP_7420, GERP_7421,
              GERP_7422, GERP_7423, GERP_7424, GERP_7425)

# load in the data
AlleleFreqs <- read.tab("/project/beetlegenes/sholtz1/beetle_test/Filtered_FRQs.delim", header = TRUE)

Chromosomes <- c("NC_007416.3",    "NC_007417.3" ,   "NC_007418.3",    "NC_007419.2",
                 "NC_007420.3",   "NC_007421.3",    "NC_007422.5",    "NC_007423.3",
                 "NC_007424.3",    "NC_007425.3")

AlleleFreqs <- AlleleFreqs %>%
  filter(Chromosome %in% Chromosomes)


# start a timer for reference
start_time <- proc.time()

# First create a "small" data frame that has only one entry per chromosome by
#   Location combination
UniqueLocations <- unique(AlleleFreqs[,c("Chromosome", "Location")])

# Add columns to the allele frequency data to hold the RS scores and p values
AlleleFreqs$RS_score <- rep(NA, nrow(AlleleFreqs))
AlleleFreqs$P_value <- rep(NA, nrow(AlleleFreqs))

# Now loop through the unique combinations, find the GERP data, and add it to all
# corresponding entries for the allele frequency data
for(i in 1:nrow(UniqueLocations)){
  # Find the row in the GERP data that corresponds to this chromosome and location combination
  GERP_row <- which(Gerp$Chromosome == UniqueLocations$Chromosome[i] &
                      Gerp$start <= UniqueLocations$Location[i] &
                      Gerp$end _007416.3_rm.rates.bed
                    # It seems like many of the allele locations don't have a corresponding GERP score,
                    #   so double check that each one has an entry in the GERP data
                    if(length(GERP_row) == 1){
                      # Find all the allele frequency rows that correspond to this unique entry
                      AlleleRows <- which(AlleleFreqs$Chromosome == UniqueLocations$Chromosome[i] &
                                            AlleleFreqs$Location == UniqueLocations$Location[i])
                      # Add the RS_score and P_value values for this unique combination to each of
                      # those entries
                      AlleleFreqs$RS_score[AlleleRows] <- Gerp$RS_score[GERP_row]
                      AlleleFreqs$P_value[AlleleRows] <- Gerp$P_value[GERP_row]
                    }
}


##Code for combining the RS scores to each location

##read in Files bed files containing max RS and RS for each location

NC_007416.3_rm.rates.bed <-read.delim("/project/beetlegenes/sholtz1/GERP/analyses/gerp/NC_007416.3_rm.rates.bed", header = FALSE)
NC_007417.3_rm.rates.bed <-read.delim("/project/beetlegenes/sholtz1/GERP/analyses/gerp/NC_007417.3_rm.rates.bed", header = FALSE)
NC_007418.3_rm.rates.bed <-read.delim("/project/beetlegenes/sholtz1/GERP/analyses/gerp/NC_007418.3_rm.rates.bed", header = FALSE)
NC_007419.2_rm.rates.bed <-read.delim("/project/beetlegenes/sholtz1/GERP/analyses/gerp/NC_007419.2_rm.rates.bed", header = FALSE)
NC_007420.3_rm.rates.bed <-read.delim("/project/beetlegenes/sholtz1/GERP/analyses/gerp/NC_007420.3_rm.rates.bed", header = FALSE)
NC_007421.3_rm.rates.bed <-read.delim("/project/beetlegenes/sholtz1/GERP/analyses/gerp/NC_007421.3_rm.rates.bed", header = FALSE)
NC_007422.5_rm.rates.bed <-read.delim("/project/beetlegenes/sholtz1/GERP/analyses/gerp/NC_007422.5_rm.rates.bed", header = FALSE)
NC_007423.3_rm.rates.bed <-read.delim("/project/beetlegenes/sholtz1/GERP/analyses/gerp/NC_007423.3_rm.rates.bed", header = FALSE)
NC_007424.3_rm.rates.bed <-read.delim("/project/beetlegenes/sholtz1/GERP/analyses/gerp/NC_007424.3_rm.rates.bed", header = FALSE)
NC_007425.3_rm.rates.bed <-read.delim("/project/beetlegenes/sholtz1/GERP/analyses/gerp/NC_007425.3_rm.rates.bed", header = FALSE)

#Change column names so they match for joining

colnames(NC_007416.3_rm.rates.bed) <- c("Chromosome","not_important", "Location", "Max_RS", "RS")
NC_007417.3_rm.rates.bed <-read.delim("/project/beetlegenes/sholtz1/GERP/analyses/gerp/NC_007417.3_rm.rates.bed", header = FALSE)
NC_007418.3_rm.rates.bed <-read.delim("/project/beetlegenes/sholtz1/GERP/analyses/gerp/NC_007418.3_rm.rates.bed", header = FALSE)
NC_007419.2_rm.rates.bed <-read.delim("/project/beetlegenes/sholtz1/GERP/analyses/gerp/NC_007419.2_rm.rates.bed", header = FALSE)
NC_007420.3_rm.rates.bed <-read.delim("/project/beetlegenes/sholtz1/GERP/analyses/gerp/NC_007420.3_rm.rates.bed", header = FALSE)
NC_007421.3_rm.rates.bed <-read.delim("/project/beetlegenes/sholtz1/GERP/analyses/gerp/NC_007421.3_rm.rates.bed", header = FALSE)
NC_007422.5_rm.rates.bed <-read.delim("/project/beetlegenes/sholtz1/GERP/analyses/gerp/NC_007422.5_rm.rates.bed", header = FALSE)
NC_007423.3_rm.rates.bed <-read.delim("/project/beetlegenes/sholtz1/GERP/analyses/gerp/NC_007423.3_rm.rates.bed", header = FALSE)
NC_007424.3_rm.rates.bed <-read.delim("/project/beetlegenes/sholtz1/GERP/analyses/gerp/NC_007424.3_rm.rates.bed", header = FALSE)
NC_007425.3_rm.rates.bed <-read.delim("/project/beetlegenes/sholtz1/GERP/analyses/gerp/NC_007425.3_rm.rates.bed", header = FALSE)

#Change column names so they match for joining

colnames(NC_007416.3_rm.rates.bed) <- c("Chromosome","not_important", "Location", "Max_RS", "RS")
colnames(NC_007417.3_rm.rates.bed) <- c("Chromosome","not_important", "Location", "Max_RS", "RS")
colnames(NC_007418.3_rm.rates.bed) <- c("Chromosome","not_important", "Location", "Max_RS", "RS")
colnames(NC_007419.2_rm.rates.bed) <- c("Chromosome","not_important", "Location", "Max_RS", "RS")
colnames(NC_007420.3_rm.rates.bed) <- c("Chromosome","not_important", "Location", "Max_RS", "RS")
colnames(NC_007421.3_rm.rates.bed) <- c("Chromosome","not_important", "Location", "Max_RS", "RS")
colnames(NC_007422.5_rm.rates.bed) <- c("Chromosome","not_important", "Location", "Max_RS", "RS")
colnames(NC_007423.3_rm.rates.bed) <- c("Chromosome","not_important", "Location", "Max_RS", "RS")
colnames(NC_007424.3_rm.rates.bed) <- c("Chromosome","not_important", "Location", "Max_RS", "RS")
colnames(NC_007425.3_rm.rates.bed) <- c("Chromosome","not_important", "Location", "Max_RS", "RS")


#Join each object 1 by 1, this is done to help with memory issues

AlleleFreqs <- left_join(AlleleFreqs, NC_007416.3_rm.rates.bed, by = c("Chromosome", "Location"))
AlleleFreqs <- left_join(AlleleFreqs, NC_007417.3_rm.rates.bed, by = c("Chromosome", "Location"))
AlleleFreqs <- left_join(AlleleFreqs, NC_007418.3_rm.rates.bed, by = c("Chromosome", "Location"))
AlleleFreqs <- left_join(AlleleFreqs, NC_007419.2_rm.rates.bed, by = c("Chromosome", "Location"))
AlleleFreqs <- left_join(AlleleFreqs, NC_007420.3_rm.rates.bed, by = c("Chromosome", "Location"))
AlleleFreqs <- left_join(AlleleFreqs, NC_007421.3_rm.rates.bed, by = c("Chromosome", "Location"))
AlleleFreqs <- left_join(AlleleFreqs, NC_007422.5_rm.rates.bed, by = c("Chromosome", "Location"))
AlleleFreqs <- left_join(AlleleFreqs, NC_007423.3_rm.rates.bed, by = c("Chromosome", "Location"))
AlleleFreqs <- left_join(AlleleFreqs, NC_007424.3_rm.rates.bed, by = c("Chromosome", "Location"))
AlleleFreqs <- left_join(AlleleFreqs, NC_007425.3_rm.rates.bed, by = c("Chromosome", "Location"))

#Remove unnecessary column

AlleleFreqs <- AlleleFreqs %>%
  select(-not_important)


write_delim(AlleleFreqs, "/project/beetlegenes/sholtz1/GERP/GERP_Allele_combined.delim")
# end the timer. On my system, this took 0.3 seconds with the abbridged data you sent
proc.time() - start_time


