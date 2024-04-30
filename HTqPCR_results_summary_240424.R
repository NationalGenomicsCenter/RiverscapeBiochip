##################################################################################################
##################################################################################################
### Script parses HTqPCR data exported by WaferGen software; updated 10/12/23, JAK
##################################################################################################
##################################################################################################
library(plyr)
library(tidyverse)

### Load results files
welldat <-
  read.delim(file.choose()) # TXT file with well data
curvedat <-
  read.delim(file.choose(), header = FALSE) # TXT file with curve data

### Update the following fields <-----------------!!!
negative_names <-
  "NEG1|NEG2|NEG3|NEG4|NEG5|NEG6" # Negative control names (delimited with "|")
positive_names <-
  "POS1|POS2|POS3|POS4|POS5|POS6" # Positive control names (delimited with "|")
results_file <-
  "HTqPCR_results_summary" # Desired name for the results file

##################################################################################################
### Parse amplification count data from well data output file
##################################################################################################
### Count amplifications
names(welldat) <-
  tolower(names(welldat)) # Make column names lowercase for consistency
welldat <- welldat[order(welldat[, 1], welldat[, 2]), ]
welldat <- welldat[, c(1:4, 6, 17)]
welldat$amp <- welldat$flags
welldat$amp <-
  ifelse(welldat$amp == "CurveFitFailed", FALSE, TRUE)

welldat$sample <- gsub(positive_names, "POS", welldat$sample)
welldat$sample <- gsub(negative_names, "NEG", welldat$sample)

amp_counts <- welldat[(welldat$amp == TRUE),]
names(amp_counts)[7] <- "amps"

assays_total <- length(unique(welldat$assay))
assays_that_amped <- length(unique(amp_counts$assay))
samples_total <- length(unique(welldat$sample))
samples_that_amped <- length(unique(amp_counts$sample))

amp_counts <-
  ddply(amp_counts, .(amp_counts$sample, amp_counts$assay), nrow)
amp_counts <- amp_counts[, c(2, 1, 3)]
names(amp_counts) <- c("assay", "sample", "amps")

### Account for any assays or samples without amplifications
if (assays_that_amped != assays_total | samples_that_amped != samples_total) {
  amp_summary_wo_missing <- amp_counts %>% spread(assay, amps)
  amp_summary_wo_missing[is.na(amp_summary_wo_missing)] <- "0"
  
  amp_counts_missing_assays <-
    subset(welldat,!(assay %in% amp_counts$assay))
  amp_counts_missing_assays <- unique(amp_counts_missing_assays$assay)
  amp_summary_wo_missing[, amp_counts_missing_assays] <- 0
  
  amp_counts_missing_samples <-
    subset(welldat,!(sample %in% amp_counts$sample))
  amp_counts_missing_samples <-
    unique(amp_counts_missing_samples$sample)
  amp_counts_missing_samples <-
    data.frame(samples = amp_counts_missing_samples, row.names = NULL)
  amp_counts_missing_samples[, c(2:ncol(amp_summary_wo_missing))] <- 0
  
  names(amp_counts_missing_samples) = names(amp_summary_wo_missing)
  
  amp_summary <-
    rbind(amp_summary_wo_missing, amp_counts_missing_samples)
} else{
  amp_summary <- amp_counts %>% spread(assay, amps)
  amp_summary[is.na(amp_summary)] <- "0"
}

### Prepare amplification summary dataframe
amp_summary <- amp_summary[order(amp_summary$sample),]
amp_summary <- select(amp_summary, order(colnames(amp_summary)))
amp_summary <- select(amp_summary, sample, everything())

amp_summary <-
  amp_summary[which(amp_summary$sample != "NEG" &
                      amp_summary$sample != "POS"), ]

names(amp_summary)[-1] <-
  paste(colnames(amp_summary)[-1], "amps", sep = "_")

##################################################################################################
### Parse Ct data from well data output file
##################################################################################################
ctdat <- welldat
ctdat[which(ctdat[,6]=="CtIsOutlier"),5] <- "Outlier"
ctdat[which(ctdat[,6]=="CtIsOutlier, CtIsLarge"),5] <- "Outlier"
ctdat <- ctdat[, 3:5]
ctdat <- ctdat[order(ctdat$assay, ctdat$sample),]

ctdat$replicate <- rep(seq(1:4), length(unique(ctdat$assay)))
ctdat$assayreplicate <- paste(ctdat$replicate, ctdat$assay)
ctdat <- ctdat[,-c(1, 4)]
ctdat <-
  ctdat[which(ctdat$sample != "NEG" & ctdat$sample != "POS"), ]
                          
ctdat_summary <- reshape(
  data = ctdat,
  idvar = "sample",
  v.names = "ct",
  timevar = c("assayreplicate"),
  direction = "wide"
)

names(ctdat_summary) <-
  gsub(
    x = names(ctdat_summary),
    pattern = "ct.1 ",
    replacement = ""
  )
names(ctdat_summary) <-
  gsub(
    x = names(ctdat_summary),
    pattern = "ct.2 ",
    replacement = ""
  )
names(ctdat_summary) <-
  gsub(
    x = names(ctdat_summary),
    pattern = "ct.3 ",
    replacement = ""
  )
names(ctdat_summary) <-
  gsub(
    x = names(ctdat_summary),
    pattern = "ct.4 ",
    replacement = ""
  )

names(ctdat_summary)[seq(2, ncol(ctdat_summary), 4)] <-
  paste(names(ctdat_summary)[seq(2, ncol(ctdat_summary), 4)], "ct1", sep = "_")
names(ctdat_summary)[seq(3, ncol(ctdat_summary), 4)] <-
  paste(names(ctdat_summary)[seq(3, ncol(ctdat_summary), 4)], "ct2", sep = "_")
names(ctdat_summary)[seq(4, ncol(ctdat_summary), 4)] <-
  paste(names(ctdat_summary)[seq(4, ncol(ctdat_summary), 4)], "ct3", sep = "_")
names(ctdat_summary)[seq(5, ncol(ctdat_summary), 4)] <-
  paste(names(ctdat_summary)[seq(5, ncol(ctdat_summary), 4)], "ct4", sep = "_")

### Calculate mean Ct values
ctdat_mean <- welldat
ctdat_mean <- ctdat_mean[, 3:5]
ctdat_mean <- ctdat_mean[order(ctdat_mean$assay, ctdat_mean$sample),]

ctdat_mean <-
  aggregate(
    as.numeric(ctdat_mean$ct),
    by = list(ctdat_mean$assay, ctdat_mean$sample),
    FUN = mean,
    na.rm = TRUE
  )

ctdat_mean <- reshape(
  data = ctdat_mean,
  idvar = "Group.2",
  v.names = "x",
  timevar = c("Group.1"),
  direction = "wide"
)

ctdat_mean <-
  ctdat_mean[which(ctdat_mean$Group.2 != "NEG" &
                     ctdat_mean$Group.2 != "POS"), ]
names(ctdat_mean) <-
  gsub(x = names(ctdat_mean),
       pattern = "x.",
       replacement = "")
names(ctdat_mean)[2:ncol(ctdat_mean)] <-
  paste(names(ctdat_mean)[2:ncol(ctdat_mean)], "ctmean", sep = "_")

### Combine Ct and mean Ct dataframes
ctdat_summary <- cbind(ctdat_summary, ctdat_mean[-1])
ctdat_summary[is.na(ctdat_summary)] <- "Undetermined"

##################################################################################################
### Parse fluorescence data from curve data output file
##################################################################################################
### Extract final (i.e., maximum) fluorescence value from each curve
names(curvedat) <-
  tolower(names(curvedat)) # Make column names lowercase for consistency
fluorodat <- curvedat[curvedat$v3 == "FAM",]
fluorodat <- fluorodat[, c(1, 2, (ncol(fluorodat) - 1))]
fluorodat <- cbind(welldat[, c(3, 4, 6)], fluorodat[, 3])
names(fluorodat) <- c("assay", "sample", "flags", "fluorescence")

### Determine mean positive control fluorescence among replicates for each assay, first removing
### outlier fluorescence values and positive controls that did not amplify
fluorodat_wo_out <-
  fluorodat[!(fluorodat$fluorescence > 20000), ] # Change value deemed an outlier as appropriate
fluorodat_wo_out_amp <-
  fluorodat_wo_out[which(fluorodat_wo_out$flags != "CurveFitFailed"),]

mean_fluoro <-
  aggregate(
    fluorodat_wo_out_amp$fluorescence,
    by = list(fluorodat_wo_out_amp$assay, fluorodat_wo_out_amp$sample),
    FUN = mean
  )
names(mean_fluoro) <- c("assay", "sample", "sample_fluoro")

pos_fluoro <- mean_fluoro[(mean_fluoro$sample == "POS"),]

### Add in assays that were dropped because their positive controls did not amplify
pos_fluoro_missing_assays <-
  subset(welldat,!(assay %in% pos_fluoro$assay))
pos_fluoro_missing_assays <-
  unique(pos_fluoro_missing_assays$assay)
pos_fluoro_missing_assays <-
  data.frame(
    assay = pos_fluoro_missing_assays,
    sample = rep("POS", length(pos_fluoro_missing_assays)),
    sample_fluoro = rep(0, length(pos_fluoro_missing_assays))
  )

pos_fluoro <-
  rbind(pos_fluoro, pos_fluoro_missing_assays)
pos_fluoro <-
  pos_fluoro[order(pos_fluoro$assay),]

pos_fluoro <-
  pos_fluoro %>% slice(rep(1:n(), each = nrow(ctdat_summary) * 4))

### Determine the fluorescence of each amplification curve
fluorodat <- fluorodat[-3]
fluorodat <- fluorodat[order(fluorodat$assay, fluorodat$sample),]
fluorodat <-
  fluorodat[which(fluorodat$sample != "NEG" &
                    fluorodat$sample != "POS"), ]
fluorodat$replicate <-
  rep(seq(1:4), length(unique(fluorodat$assay)))
fluorodat <- cbind(fluorodat, pos_fluoro[3])
fluorodat$fluorodiff <-
  fluorodat$fluorescence - fluorodat$sample_fluoro

fluorodat$assayreplicate <-
  paste(fluorodat$replicate, fluorodat$assay)
fluorodat <- fluorodat[,-c(1, 3:5)]

fluorodat_summary <- reshape(
  data = fluorodat,
  idvar = "sample",
  v.names = "fluorodiff",
  timevar = c("assayreplicate"),
  direction = "wide"
)

names(fluorodat_summary) <-
  gsub(
    x = names(fluorodat_summary),
    pattern = "fluorodiff.1 ",
    replacement = ""
  )
names(fluorodat_summary) <-
  gsub(
    x = names(fluorodat_summary),
    pattern = "fluorodiff.2 ",
    replacement = ""
  )
names(fluorodat_summary) <-
  gsub(
    x = names(fluorodat_summary),
    pattern = "fluorodiff.3 ",
    replacement = ""
  )
names(fluorodat_summary) <-
  gsub(
    x = names(fluorodat_summary),
    pattern = "fluorodiff.4 ",
    replacement = ""
  )

names(fluorodat_summary)[seq(2, ncol(fluorodat_summary), 4)] <-
  paste(names(fluorodat_summary)[seq(2, ncol(fluorodat_summary), 4)], "fluorodiff1", sep = "_")
names(fluorodat_summary)[seq(3, ncol(fluorodat_summary), 4)] <-
  paste(names(fluorodat_summary)[seq(3, ncol(fluorodat_summary), 4)], "fluorodiff2", sep = "_")
names(fluorodat_summary)[seq(4, ncol(fluorodat_summary), 4)] <-
  paste(names(fluorodat_summary)[seq(4, ncol(fluorodat_summary), 4)], "fluorodiff3", sep = "_")
names(fluorodat_summary)[seq(5, ncol(fluorodat_summary), 4)] <-
  paste(names(fluorodat_summary)[seq(5, ncol(fluorodat_summary), 4)], "fluorodiff4", sep = "_")

##################################################################################################
### Combine amplification, Ct, and fluorescence dataframes and save as CSV files
##################################################################################################
results_summary <-
  cbind(amp_summary, ctdat_summary, fluorodat_summary)
results_summary <-
  results_summary[, !(names(results_summary) %in% "sample")]
results_summary <-
  results_summary[, order(colnames(results_summary))]
results_summary <- cbind(amp_summary[1], results_summary)
results_summary[] <- lapply(results_summary, as.character)

### Prep dataframe for long format
welldat <- welldat[order(welldat$assay, welldat$sample),]
assays <- unique(welldat$assay)

amps <-
  results_summary[, c("sample", paste(assays, "amps", sep = "_"))]
names(amps)[-1] <- assays
amps <- amps %>% pivot_longer(cols = names(amps[-1]),
                              names_to = "assay",
                              values_to = "amps")

ct1 <-
  results_summary[, c("sample", paste(assays, "ct1", sep = "_"))]
names(ct1)[-1] <- assays
ct1 <- ct1 %>% pivot_longer(cols = names(ct1[-1]),
                            names_to = "assay",
                            values_to = "ct1")

ct2 <-
  results_summary[, c("sample", paste(assays, "ct2", sep = "_"))]
names(ct2)[-1] <- assays
ct2 <- ct2 %>% pivot_longer(cols = names(ct2[-1]),
                            names_to = "assay",
                            values_to = "ct2")

ct3 <-
  results_summary[, c("sample", paste(assays, "ct3", sep = "_"))]
names(ct3)[-1] <- assays
ct3 <- ct3 %>% pivot_longer(cols = names(ct3[-1]),
                            names_to = "assay",
                            values_to = "ct3")

ct4 <-
  results_summary[, c("sample", paste(assays, "ct4", sep = "_"))]
names(ct4)[-1] <- assays
ct4 <- ct4 %>% pivot_longer(cols = names(ct4[-1]),
                            names_to = "assay",
                            values_to = "ct4")

ctmean <-
  results_summary[, c("sample", paste(assays, "ctmean", sep = "_"))]
names(ctmean)[-1] <- assays
ctmean <- ctmean %>% pivot_longer(cols = names(ctmean[-1]),
                                  names_to = "assay",
                                  values_to = "ctmean")

fluorodiff1 <-
  results_summary[, c("sample", paste(assays, "fluorodiff1", sep = "_"))]
names(fluorodiff1)[-1] <- assays
fluorodiff1 <-
  fluorodiff1 %>% pivot_longer(
    cols = names(fluorodiff1[-1]),
    names_to = "assay",
    values_to = "fluorodiff1"
  )

fluorodiff2 <-
  results_summary[, c("sample", paste(assays, "fluorodiff2", sep = "_"))]
names(fluorodiff2)[-1] <- assays
fluorodiff2 <-
  fluorodiff2 %>% pivot_longer(
    cols = names(fluorodiff2[-1]),
    names_to = "assay",
    values_to = "fluorodiff2"
  )

fluorodiff3 <-
  results_summary[, c("sample", paste(assays, "fluorodiff3", sep = "_"))]
names(fluorodiff3)[-1] <- assays
fluorodiff3 <-
  fluorodiff3 %>% pivot_longer(
    cols = names(fluorodiff3[-1]),
    names_to = "assay",
    values_to = "fluorodiff3"
  )

fluorodiff4 <-
  results_summary[, c("sample", paste(assays, "fluorodiff4", sep = "_"))]
names(fluorodiff4)[-1] <- assays
fluorodiff4 <-
  fluorodiff4 %>% pivot_longer(
    cols = names(fluorodiff4[-1]),
    names_to = "assay",
    values_to = "fluorodiff4"
  )

results_summary_long <-
  cbind(
    amps,
    ct1[3],
    ct2[3],
    ct3[3],
    ct4[3],
    ctmean[3],
    fluorodiff1[3],
    fluorodiff2[3],
    fluorodiff3[3],
    fluorodiff4[3]
  )

### Save wide and long dataframes as CSV files
write.csv(
  results_summary,
  file = paste(results_file, "wide.csv", sep = "_"),
  row.names = FALSE
)
write.csv(
  results_summary_long,
  file = paste(results_file, "long.csv", sep = "_"),
  row.names = FALSE
)
print("Finished!")
