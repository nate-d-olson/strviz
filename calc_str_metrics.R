

## Description: Creates new data frame with only majority peaks  
## Input: A data frame of the counts split up by read and direction.
## Output: Only the top two allele counts for each Locus, disregards stutter for heterozygous genotypes.

peak_coverage <- function (allele_counts_df) {
    allele_counts_df[is.na(allele_counts_df)] <- 0
    allele_counts_df %>% 
        rowwise() %>% 
        mutate(Sum_R1 = (D2_R1 + D1_R1),
               Sum_R2 = (D2_R2 + D1_R2),
               Sum_D1 = (D1_R1 + D1_R2),
               Sum_D2 = (D2_R1 + D2_R2),
               Seq_Coverage = sum(D2_R1, D2_R2, D1_R1, D1_R2)) %>% 
        group_by(Allele,Locus) %>% 
        mutate(Allele_Coverage = sum(Seq_Coverage)) %>% 
        # changed to peak coverage as it represents the coverage for the 
        # peak and not the coverage of the majority peaks coverage
        group_by(Locus) %>% 
        mutate(Locus_Coverage = sum(Seq_Coverage)) %>% 
        top_n(n = 2, wt = Seq_Coverage)
}
#Want to NA out all homozygous stutter counts



## Description: Assigns a genotype call to each locus based on the peak height 
##      ratio threshold and uses call to calculate the majority peaks coverage.
## Input: Data frame with only top two allele counts per locus. 
##          gt_threshold - cutoff for the ratio between top two allele calls,
##          default value 0.2
## Output: Data frame with extra column containing genotype call per locus.

genotype_call <- function (peak_cov_df, gt_threshold = 0.2) {
    peak_cov_df %>% 
      group_by(Locus) %>% 
      mutate(Genotype = ifelse(min(Seq_Coverage)/max(Seq_Coverage) > gt_threshold, 
                               "Heterogyzous", "Homozygous"))
#              ,Coverage_of_Majority_Peaks = ifelse(Genotype == "Heterozygous", 
#                                             sum(Seq_Coverage), 
#                                             max(Seq_Coverage)))
}



## Description: Assigns a peak height ratio to each heterozygous locus based on the top two peaks.
## Input: Data frame with only top two allele counts per locus, and genotype call.
## Output: Data frame with extra column containing the peak heigh ratio for heterzygous loci, and NA for homozygous loci. One PHR per locus.
peak_height_ratio <- function (genotype_call_df) {
      genotype_call_df %>% 
        group_by(Locus) %>% 
        mutate(peak_height_ratio =  ifelse((min(Seq_Coverage)/max(Seq_Coverage))> .2,
                           round((min(Seq_Coverage)/max(Seq_Coverage)),4), NA_real_))
}

## Description: Assigns a read bias to each allele (two for heterzygotes and
## only the majority peak for homozygotes). 
## Input: Data frame with top two allele counts per locus, genotype call, 
##          and PHR per loci. 
## Output: Data frame with extra column containing the read bias for each allele.

read_bias <- function(peak_height_ratio_df) {
    ungroup(peak_height_ratio_df) %>% 
      rowwise() %>%
      mutate(Read_Bias = min(c(Sum_R1, Sum_R2))/Seq_Coverage)
}
#Want to NA out all homozygous stutter counts


## Description: Assigns a strand bias to each allele (two for heterzygotes and only the majority peak for homozygotes).
## Input: Data frame with top two allele counts per locus, genotype call, PHR per loci, and read bias per allele.
## Output: Data frame with extra column containing the strand bias for each allele. 

strand_bias <- function (read_bias_df) {
  read_bias_df %>% 
      rowwise() %>%
      mutate(Strand_Bias = min(c(Sum_D1, Sum_D2))/Seq_Coverage)
}
#Want to NA out all homozygous stutter counts


#Prelimiary coverage of non majority peaks (need to figure out how to do this
#all in one data frame, and not include the stutter alleles)
non_maj_peaks <- function (allele_counts_df) {
    allele_counts_df %>% 
        rowwise() %>% 
        mutate(Percentage_of_non_Majority_peaks = 
                   (Allele_Coverage - Seq_Coverage)/Allele_Coverage)

}
             

calc_allele_metrics <- function(allele_counts_df){
    peak_coverage(allele_counts_df) %>% 
        genotype_call() %>% 
        peak_height_ratio() %>% 
        read_bias() %>% 
        strand_bias() %>% 
        non_maj_peaks ()
}