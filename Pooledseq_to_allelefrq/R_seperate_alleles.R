library("tidyverse")


pairwise_frq <- read.delim( "AllPops_allelefrq_pwc" ,header = TRUE)

test_frqs<-read.delim("AllPops.sync"  ,header = TRUE)

pool_names <- c("Chromosome","Location","Ref","22_0_NA_S31_R1_001.fastq.gz","22_8_NA_S38_R1_001.fastq.gz","26_0_NA_S20_R1_001.fastq.gz","26_8_NA_S30_R1_001.fastq.gz","30_0_NA_S41_R1_001.fastq.gz","30_8_NA_S36_R1_001.fastq.gz","31_0_NA_S34_R1_001.fastq.gz","31_8_C_S21_R1_001.fastq.gz","31_8_E_S25_R1_001.fastq.gz","32_0_NA_S35_R1_001.fastq.gz","32_8_C_S2_R1_001.fastq.gz","32_8_E_S17_R1_001.fastq.gz","33_0_NA_S16_R1_001.fastq.gz","33_8_C_S29_R1_001.fastq.gz","33_8_E_S16_R1_001.fastq.gz","35_0_NA_S39_R1_001.fastq.gz","35_8_C_S3_R1_001.fastq.gz","35_8_E_S1_R1_001.fastq.gz","36_0_NA_S23_R1_001.fastq.gz","36_8_C_S28_R1_001.fastq.gz","36_8_E_S26_R1_001.fastq.gz","40_0_NA_S6_R1_001.fastq.gz","40_8_C_S13_R1_001.fastq.gz","40_8_E_S24_R1_001.fastq.gz","42_0_NA_S26_R1_001.fastq.gz","42_8_NA_S48_R1_001.fastq.gz","44_0_NA_S4_R1_001.fastq.gz","44_8_NA_S25_R1_001.fastq.gz","45_0_NA_S28_R1_001.fastq.gz","45_8_NA_S33_R1_001.fastq.gz","47_0_NA_S22_R1_001.fastq.gz","47_8_NA_S22_R1_001.fastq.gz","48_0_NA_S2_R1_001.fastq.gz","48_8_NA_S18_R1_001.fastq.gz","49_0_NA_S44_R1_001.fastq.gz","49_8_NA_S21_R1_001.fastq.gz","51_0_NA_S15_R1_001.fastq.gz","51_8_C_S12_R1_001.fastq.gz","51_8_E_S15_R1_001.fastq.gz","52_0_NA_S14_R1_001.fastq.gz","52_8_C_S36_R1_001.fastq.gz","52_8_E_S10_R1_001.fastq.gz","53_0_NA_S1_R1_001.fastq.gz","53_8_C_S47_R1_001.fastq.gz","53_8_E_S40_R1_001.fastq.gz","54_0_NA_S29_R1_001.fastq.gz","54_8_C_S20_R1_001.fastq.gz","54_8_E_S30_R1_001.fastq.gz","55_0_NA_S45_R1_001.fastq.gz","55_8_C_S34_R1_001.fastq.gz","55_8_E_S14_R1_001.fastq.gz","56_0_NA_S45_R1_001.fastq.gz","56_8_C_S18_R1_001.fastq.gz","56_8_E_S38_R1_001.fastq.gz","57_0_NA_S27_R1_001.fastq.gz","57_8_C_S7_R1_001.fastq.gz","57_8_E_S3_R1_001.fastq.gz","58_0_NA_S8_R1_001.fastq.gz","58_8_C_S32_R1_001.fastq.gz","58_8_E_S8_R1_001.fastq.gz","59_0_NA_S47_R1_001.fastq.gz","59_8_C_S4_R1_001.fastq.gz","59_8_E_S32_R1_001.fastq.gz","61_0_NA_S27_R1_001.fastq.gz","61_8_NA_S5_R1_001.fastq.gz","63_0_NA_S37_R1_001.fastq.gz","63_8_NA_S40_R1_001.fastq.gz","64_0_NA_S19_R1_001.fastq.gz","64_8_NA_S9_R1_001.fastq.gz","65_0_NA_S46_R1_001.fastq.gz","65_8_NA_S37_R1_001.fastq.gz","66_0_NA_S43_R1_001.fastq.gz","66_8_NA_S7_R1_001.fastq.gz","67_0_NA_S9_R1_001.fastq.gz","67_8_NA_S48_R1_001.fastq.gz","73_0_NA_S31_R1_001.fastq.gz","73_8_C_S46_R1_001.fastq.gz","73_8_E_S19_R1_001.fastq.gz","74_0_NA_S39_R1_001.fastq.gz","74_8_C_S41_R1_001.fastq.gz","74_8_E_S11_R1_001.fastq.gz","75_0_NA_S6_R1_001.fastq.gz","75_8_C_S17_R1_001.fastq.gz","75_8_E_S35_R1_001.fastq.gz","77_0_NA_S44_R1_001.fastq.gz","77_8_C_S24_R1_001.fastq.gz","77_8_E_S43_R1_001.fastq.gz","78_0_NA_S10_R1_001.fastq.gz","78_8_C_S11_R1_001.fastq.gz","78_8_E_S42_R1_001.fastq.gz","79_0_NA_S23_R1_001.fastq.gz","79_8_C_S42_R1_001.fastq.gz","79_8_E_S12_R1_001.fastq.gz","80_0_NA_S5_R1_001.fastq.gz","80_8_C_S13_R1_001.fastq.gz","80_8_E_S33_R1_001.fastq.gz")



colnames(test_frqs) <-c(pool_names)

test_frqs <- test_frqs %>%
  select(1:99)

for (i in 4:length(pool_names)) { 
  
  good_names_A <-paste(pool_names[i], "A", sep = "_")
  good_names_T <- paste(pool_names[i], "T", sep = "_")
  good_names_C <- paste(pool_names[i], "C", sep = "_")
  good_names_G <- paste(pool_names[i], "G", sep = "_")
  good_names_N <- paste(pool_names[i], "N", sep = "_")
  good_names_del <- paste(pool_names[i], "del", sep = "_")
  test_frqs <- separate(test_frqs, pool_names[i], 
  paste(c(good_names_A,good_names_T,good_names_C, good_names_G, good_names_N, good_names_del  ) , sep = ":"))
  


}


pairwise_frq <- pairwise_frq %>%
  mutate(SNP_LOC = paste(X..chr, pos, sep = "_"))


SNP_locations <- c(pairwise_frq$SNP_LOC)

test_frqs <- test_frqs %>%
  mutate(LOC_SNP = paste(Chromosome, Location, sep = "_")) %>%
filter(LOC_SNP %in% SNP_locations)

#Get the data organized by landscape

test_frqs <- test_frqs %>%
  pivot_longer(cols = 4:ncol(test_frqs), names_to = "Identity", values_to = "Base_frq") %>%
  separate(col = Identity, into = c("Landscape_Loc", "Base" ), sep = "_R1_001.fastq.gz_", ) %>%
  filter(!is.na(Base))%>%
  pivot_wider(names_from = Base, values_from = Base_frq)

test_frqs <- test_frqs %>%
  separate(col = Landscape_Loc, into = c("Landscape", "Generation", "Core_Edge", "S" ), sep = "_", )

write_delim(test_frqs, "/project/beetlegenes/sholtz1/beetle_test/Allele_frequencies.delim")
