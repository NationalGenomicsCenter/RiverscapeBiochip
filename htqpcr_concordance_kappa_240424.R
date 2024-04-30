# JE: need long data frame with 6 columns:
#one column with sample ID
#one column for "Assay"
#two columns for HT results (# wells amped named "ht" and binary values named "ht_bin"),
#two columns for qPCR results (# wells amped named "qPCR" and binary values named "qPCR_bin")

data<- read.csv(file.choose())
head(data)

sum(is.na(data$ht_bin))
sum(is.na(data$qPCR_bin)) # rows to remove from consideration here

data_concord <- subset(data, is.na(data$qPCR_bin) == F)
nrow(data_concord) # 696


# JE: fixed error in code for function
require("irr")
### Build confusion matrices 

confusion_matrix <- function(HTqpcr, qPCR){
  ht_binary <- ifelse(HTqpcr > 1, 1, 0)
  qPCR_binary <- ifelse(as.numeric(qPCR) > 0, 1, 0)
  
  PP <- sum(ht_binary == 1 & qPCR_binary == 1, na.rm = T) #++
  NN <- sum(ht_binary == 0 & qPCR_binary == 0, na.rm = T) #--
  PN <- sum(ht_binary == 1 & qPCR_binary == 0, na.rm = T) #+-
  NP <- sum(ht_binary == 0 & qPCR_binary == 1, na.rm = T) #-+
  
  print(matrix(c(PN,NN,PP,NP), nrow = 2, ncol = 2))
  
}


assays <- levels(as.factor(data_concord$Assay))
confusion_matrix_list <- list()
concordance <- c()
concordance_interval_lwr <- c()
concordance_interval_upr <- c()
kappa <- c()
for(i in 1:length(assays)){
  
  ### Confusion matrices
  data_assay <- subset(data_concord, data_concord$Assay == assays[i])
  confusion_matrix_list[[i]] <- confusion_matrix(HTqpcr = data_assay$ht, qPCR = data_assay$qPCR)

  ### Raw concordance
  concordance[i] <- (confusion_matrix_list[[i]][2,1] + confusion_matrix_list[[i]][1,2])/sum(confusion_matrix_list[[i]])
  
  ### Bootstrap
  boot <- c()
  for(j in 1:500){
    data_assay_boot <- data_assay[sample(1:nrow(data_assay), nrow(data_assay), replace = T),]
    data_assay_confusion_matrix <- confusion_matrix(HTqpcr = data_assay_boot $ht, qPCR = data_assay_boot$qPCR)
    boot[j] <- (data_assay_confusion_matrix[2,1] + data_assay_confusion_matrix[1,2])/sum(data_assay_confusion_matrix)
  }
  
  concordance_interval_lwr[[i]] <- quantile(boot, 0.025)
  concordance_interval_upr[[i]] <- quantile(boot, 0.975)
  
  ### Kappa
  kappa2
  kappa[i] <- kappa2(data_assay[, colnames(data_assay) == "ht_bin" | colnames(data_assay) == "qPCR_bin"])$value
}

out_table <- cbind(assays,
                   concordance,
                   concordance_interval_lwr,
                   concordance_interval_upr,
                   kappa)

out_table

### Figure
concord_plot <- barplot(concordance, ylim = c(0,1))
concord_plot[,1]
segments(concord_plot[,1],
         sapply(concordance_interval_lwr, "[[", 1),
         concord_plot[,1],
         sapply(concordance_interval_upr, "[[", 1),
         lwd = 3)

