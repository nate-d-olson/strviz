

## Description: Creates new data frame with only majority peaks  
## Input: A data frame of the counts split up by read and direction.
## Output: Only the top two allele counts for each Locus, disregards stutter for heterozygous genotypes.

cov_maj_peaks <- function (allele_counts_df) {
    allele_counts_df[is.na(allele_counts_df)] <- 0
    allele_counts_df$Coverage_of_Majority_Peaks <- (allele_counts_df$D2_R1 + allele_counts_df$D2_R2 + allele_counts_df$D1_R1 + allele_counts_df$D1_R2)
    allele_counts_df  %>% 
                    group_by(Locus)  %>%  
                    top_n(2)
}
#Want to NA out all homozygous stutter counts



## Description: Assigns a genotype call to each locus based on the peak height ratio (20%).
## Input: Data frame with only top two allele counts per locus.
## Output: Data frame with extra column containing genotype call per locus.

genotype_call <- function (maj_peak_df) {
    maj_peak_df %>% 
      group_by(Locus) %>% 
      mutate(Genotype = ifelse(min(Coverage_of_Majority_Peaks)/max(Coverage_of_Majority_Peaks) > 0.2, "Heterogyzous", "Homozygous"))
}

## Description: Assigns a peak height ratio to each heterozygous locus based on the top two peaks.
## Input: Data frame with only top two allele counts per locus, and genotype call.
## Output: Data frame with extra column containing the peak heigh ratio for heterzygous loci, and NA for homozygous loci. One PHR per locus.

PHR <- function (genotype_call_df) {
      genotype_call_df %>% 
      group_by(Locus) %>% 
      mutate(PHR =  ifelse((min(Coverage_of_Majority_Peaks)/max(Coverage_of_Majority_Peaks))> .2,
                           round((min(Coverage_of_Majority_Peaks)/max(Coverage_of_Majority_Peaks)),4), NA_real_))
}

## Description: Assigns a read bias to each allele (two for heterzygotes and only the majority peak for homozygotes).
## Input: Data frame with top two allele counts per locus, genotype call, and PHR per loci.
## Output: Data frame with extra column containing the read bias for each allele. 

Read_bias <- function(PHR_df) {
    ungroup(PHR_df) %>% 
      mutate(Sum_R1 = (D2_R1 + D1_R1)) %>% 
      mutate(Sum_R2 = (D2_R2 + D1_R2)) %>% 
      rowwise() %>%
      mutate(Read_Bias = min(c(Sum_R1, Sum_R2))/Coverage_of_Majority_Peaks)
}
#Want to NA out all homozygous stutter counts


## Description: Assigns a strand bias to each allele (two for heterzygotes and only the majority peak for homozygotes).
## Input: Data frame with top two allele counts per locus, genotype call, PHR per loci, and read bias per allele.
## Output: Data frame with extra column containing the strand bias for each allele. 

Strand_bias <- function (Read_bias_df) {
  Read_bias_df %>% 
      mutate(Sum_D1 = (D1_R1 + D1_R2),Sum_D2 = (D2_R1 + D2_R2)) %>% 
      rowwise() %>%
      mutate(Strand_Bias = min(c(Sum_D1, Sum_D2))/Coverage_of_Majority_Peaks)
}
#Want to NA out all homozygous stutter counts


#Prelimiary coverage of non majority peaks (need to figure out how to do this all in one data frame, and not include the stutter alleles)
non_maj_peaks <- function (allele_counts_df) {
    group_by(Locus, Allele) %>% 
    mutate(Percentage_of_non_Majority_peaks = (sum(Coverage_of_Majority_Peaks)-max(Coverage_of_Majority_Peaks))/(sum(Coverage_of_Majority_Peaks)))
   #I want to make a df with all the counts only for the top 2 alleles, after I have that the mutate function properly calcuates the non maj peaks ratio  
}
             

calc_allele_metrics <- function(allele_counts_df){
    cov_maj_peaks(allele_counts_df) %>% 
        genotype_call() %>% 
        PHR() %>% 
        Read_bias() %>% 
        Strand_bias() %>% 
        non_maj_peaks ()
}