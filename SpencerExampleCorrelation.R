library(nlme)

SelData <- read.csv("SelData.csv")
# Set Chromosome to a factor and then make a useful interaction column
SelData$Chrom <- as.factor(SelData$Chrom)
SelData$ChromLoc <- interaction(SelData$Chrom, SelData$location)


SigEstModel <- lme(fixed = sigma ~  -1 + ChromLoc, data = SelData, random = ~ 1|window,
                   correlation = corExp(), method = "ML")
# nested random effects: chromosome/location

SigEstModel
# Generate approximate 95% confidence intervals
SigInts <- intervals(SigEstModel, which = "fixed")
SigInts

Test_LME_data <- read.table("Full_filtered_data.delim", header = TRUE)
?lme
EstModel <- lme(fixed = Frequency_change ~ RS, 
                data = Test_LME_data,
                random = ~1|Chromosome/Location,  
                method = "ML",
                correlation = corExp())

EstModel

mod.corExp <- update(EstModel, correlation = nlme::corExp(form = ~ RS + Location, nugget=T))
