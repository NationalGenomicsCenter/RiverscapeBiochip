# RiverscapeBiochip

Datasets and scripts used for the development and validation of a HT-qPCR biochip (i.e., biochip) protocol for simultaneous detection of 10 freshwater fishes.

## File guide

•	Master_Long_PNW10242023.csv – Dataset containing paired HT-qPCR and ST-qPCR results for 367 environmental samples, required for htqpcr_concordance_kappa_240424.R script.

•	eDNAssay_cross_amp_HTqPCR.csv - Dataset containing information on putative cross-amplifications of taxon-specific assays with tissue samples on the HT-qPCR platform described in Appendix 2, required for eDNAssay_cross_amp_HTqPCR_240424.R script. 

•	htqpcr_concordance_kappa_240424.R - Script used to generate raw concordance (the proportion of results in agreement between ST-qPCR and HT-qPCR), Cohen’s kappa statistic, bootstrapped 95% confidence intervals around raw concordance point estimates for each assay.

•	eDNAssay_cross_amp_HTqPCR_240424.R - Script used to run logistic regression models (eDNAssay assignment probabilities for each assay/non-target combination against non-target tissue amplification results) and generate figures in Appendix 2. 

•	HTqPCR_results_summary_240424.R - Script used to summarize HT-qPCR data from RAW files generated from each HT-qPCR run using the WaferGen qPCR Software (version 2.8.33.0, Wafergen Biosystems, Inc.). 

## Contact information
Please reach out to us at the National Genomics Center for Wildlife and Fish Conservation with any questions or comments. Scripts and data summaries were created by Joanna Elmore at joanna.elmore@usda.gov, John Kronenberger at john.kronenberger@usda.gov, and Taylor Wilcox at taylor.wilcox@usda.gov.
