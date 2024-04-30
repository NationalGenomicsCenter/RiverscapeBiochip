### Amplification versus eDNAssay
### Test of how specificity tests relate to eDNAssay assignment probabilities

#data file has non-targets (species column), maximum assignment probablity (Max AP column) for respective assay (assay column), and number of wells that amplified on HT-qPCR (columns: chip1-chip6)

#Columns:
#- assay
#- species
#- chip1 ... chip6 are number of wells with amplification with 1 - 3 individuals each run in duplicate on the sample plate (thus 2 - 6 entries per taxon)
#- chip_binary = was there apparent cross-amplification? requires that there was at least 2/4 replicates for a sample amplifying
#- chip_binary_2 = was there apparent cross-amplification? any >1/4 replicates amplifying


data <- read.csv(file.choose())
data$chip_binary # stringent filtering - at least 2/4 replicates amp to call positive
data$chip_binary_2 # lax filtering - any amplification (including 1/4 for a single tissue) is positive

par(mfrow = c(1,2))
plot(data$chip_binary ~ data$Max.AP,
     xlab = "eDNAssay AP", ylab = "Amplification?", main = "Stringent filtering")
summary(glm(data$chip_binary ~ data$Max.AP, data = data, family = "binomial"))

range(subset(data, data$chip_binary == 1)$Max.AP)
range(subset(data, data$chip_binary == 0)$Max.AP, na.rm = T)

plot(data$chip_binary_2 ~ data$Max.AP,
     xlab = "eDNAssay AP", ylab = "Amplification?", main = "Non-stringent filtering")
summary(glm(data$chip_binary_2 ~ data$Max.AP, data = data, family = "binomial"))


range(subset(data, data$chip_binary_2 == 1)$Max.AP, na.rm = T)
range(subset(data, data$chip_binary_2 == 0)$Max.AP, na.rm = T)
