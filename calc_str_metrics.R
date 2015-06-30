

## Description: Creates new data frame with only majority peaks  
## Input
## Output
cov_maj_peaks <- function (allele_counts_df) {
    allele_counts_df[is.na(allele_counts_df)] <- 0
    allele_counts_df$Coverage_of_Majority_Peaks <- (allele_counts_df$D2_R1 + allele_counts_df$D2_R2 + allele_counts_df$D1_R1 + allele_counts_df$D1_R2)
    allele_counts_df  %>% 
                    group_by(Locus)  %>%  
                    top_n(2)
}

#preliminary code for genotype calls (.2 can be changed to any other ratio if needed)
genotype_call <- function (maj_peak_df) {
    maj_peak_df %>% 
      group_by(Locus) %>% 
      mutate(Genotype = ifelse(min(Coverage_of_Majority_Peaks)/max(Coverage_of_Majority_Peaks) > 0.2, "Heterogyzous", "Homozygous"))
}

#preliminary code for PHR

PHR <- function (genotype_call_df) {
      genotype_call_df %>% 
      group_by(Locus) %>% 
      mutate(PHR =  ifelse((min(Coverage_of_Majority_Peaks)/max(Coverage_of_Majority_Peaks))> .2,
                           round((min(Coverage_of_Majority_Peaks)/max(Coverage_of_Majority_Peaks)),4), NA_real_))
}


#preliminary code for Read Bias
Read_bias <- function(PHR_df) {
    ungroup(PHR_df) %>% 
      mutate(Sum_R1 = (D2_R1 + D1_R1)) %>% 
      mutate(Sum_R2 = (D2_R2 + D1_R2)) %>% 
      rowwise() %>%
      mutate(Read_Bias = min(c(Sum_R1, Sum_R2))/Coverage_of_Majority_Peaks)
}


#preliminary code for Strand Bias
Strand_bias <- function (Read_bias_df) {
  Read_bias_df %>% 
      mutate(Sum_D1 = (D1_R1 + D1_R2),Sum_D2 = (D2_R1 + D2_R2)) %>% 
      rowwise() %>%
      mutate(Strand_Bias = min(c(Sum_D1, Sum_D2))/Coverage_of_Majority_Peaks)
}

#Prelimiary coverage of non majority peaks (need to figure out how to do this all in one data frame, and not include the stutter alleles)
# new_seq3 <- new_seq  %>% 
# group_by(Locus, Allele)  %>%  
# mutate(Percentage_of_non_Majority_peaks = (sum(Coverage_of_Majority_Peaks)-max(Coverage_of_Majority_Peaks))/(sum(Coverage_of_Majority_Peaks)))
#                

calc_allele_metrics <- function(allele_counts_df){
    cov_maj_peaks(allele_counts_df) %>% 
        genotype_call() %>% 
        PHR() %>% 
        Read_bias() %>% 
        Strand_bias()
}