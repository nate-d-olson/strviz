

## Description: Creates new data frame with only majority peaks  
## Input: A data frame of the counts split up by read and direction.
## Output: Only the top two allele counts for each Locus, disregards stutter for heterozygous genotypes.

#Have function that has all the counting information in a df
#Create new df from count df with genotype (each locus and genotype) 
#Create a homozygous metrics df and heterzygous metrics df (filter by genotype(join-> filter))
#Then calculate rest of metrics

coverage_calc <- function (allele_counts_df) {
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
        mutate(Locus_Coverage = sum(Seq_Coverage)) 
    
}




## Description: Assigns a genotype call to each locus based on the peak height 
##      ratio threshold and uses call to calculate the majority peaks coverage.
## Input: Data frame with only top two allele counts per locus. 
##          gt_threshold - cutoff for the ratio between top two allele calls,
##          default value 0.2
## Output: Data frame with extra column containing genotype call per locus.

genotype_call <- function (peak_cov_df, gt_threshold = 0.2) {
    peak_cov_df %>% 
      top_n(n = 2, wt = Seq_Coverage) %>% 
      group_by(Locus) %>% 
      mutate(Genotype = ifelse(min(Seq_Coverage)/max(Seq_Coverage) > gt_threshold, 
                               "Heterozygous", "Homozygous")) %>% 
        select(Locus, Genotype) %>% 
        unique ()
}


## Description: 
## Input: 
## Output: 
calc_het_metrics <- function(geno_df, cov_df){
     het_cov_df <- filter(geno_df, Genotype == "Heterozygous") %>% 
     left_join(cov_df)  
     het_cov_df$Allele <- het_cov_df$Allele %>% str_replace("X", "100") %>% 
     str_replace("Y", "1010") %>%  
     as.numeric() 
     het_cov_df <- het_cov_df %>%      
     group_by(Locus) %>% 
     stutter_homo() %>% 
     as.character() %>% 
     str_replace("100", "X") %>% 
     str_replace("1010", "Y") %>% 
     top_n(2, wt= Seq_Coverage)
     peak_height_ratio() %>% 
     read_bias() %>% 
     strand_bias() %>% 
     non_maj_peaks ()
 }

## Description: 
## Input: 
## Output:
calc_homo_metrics <- function(geno_df, cov_df){
    homo_cov_df <- filter(geno_df, Genotype == "Homozygous") %>% 
        left_join(cov_df) %>% 
        mutate(Allele = as.numeric(Allele)) %>% 
        stutter_homo() %>%  
        top_n(1, wt = Seq_Coverage)  %>% 
        stutter_rat() %>% 
        peak_height_ratio() %>% 
        read_bias() %>% 
        strand_bias() %>% 
        non_maj_peaks ()
     
    homo_cov_df[homo_cov_df == 1] <- NA
    ## see note about function returning values in calc_allele_metric
    homo_cov_df     
}

## Description: 
## Input: 
## Output:
stutter_homo <- function(cov_df) {
    cov_df %>% 
    group_by(Locus) %>%    
        mutate(Stutter_1 = lag(Allele_Coverage)) 
}
 
## Description: 
## Input: 
## Output:
stutter_rat <- function(stutter_df) {
    stutter_df %>% 
        mutate(Stutter_Ratio = Stutter_1/Allele_Coverage ) 
        
}
               
         
         
     
## Not sure if you still need this code but it was presenting the script from being sourced    
#      peak_height_ratio(homo_cov_df) %>% 
#          read_bias() %>% 
#          strand_bias() %>% 
#          non_maj_peaks ()
#  }
#  


## Description: Assigns a peak height ratio to each heterozygous locus based on the top two peaks.
## Input: Data frame with only top two allele counts per locus, and genotype call.
## Output: Data frame with extra column containing the peak heigh ratio for heterzygous loci, and NA for homozygous loci. One PHR per locus.
peak_height_ratio <- function (genotype_call_df, gt_threshold = 0.2) {
      genotype_call_df %>% 
        group_by(Locus) %>% 
        mutate(Peak_Height_Ratio =  ifelse((min(Seq_Coverage)/max(Seq_Coverage))> gt_threshold,
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
 
            

# calc_allele_metrics <- function(allele_counts_df){
#     peak_coverage(allele_counts_df) %>% 
#         genotype_call() %>% 
#         peak_height_ratio() %>% 
#         read_bias() %>% 
#         strand_bias() %>% 
#         non_maj_peaks ()
# }


## Description: 
## Input: 
## Output:
calc_geno_metrics <- function(geno_df, cov_df){
    ## function names did not match to function calls, changed function names in
    ## definitions above
    ## Their are two issues with the code 
    ## 1. homo_df and het_df have different numbers of rows, for bind_rows to 
    ##    work the data frames need to have the same number of rows and same 
    ##    column names
    ## 2. het_df have 1592 row,  I think you skipped a step ;)
    homo_df <- calc_homo_metrics(geno_df, cov_df) 
    het_df <- calc_het_metrics(geno_df, cov_df)
    bind_rows(homo_df, het_df)
}

## Description: 
## Input: 
## Output:
calc_allele_metrics <- function(allele_counts_df){
    cov_df <- coverage_calc(allele_counts_df)  
    geno_df<- genotype_call(cov_df)
    ## for a funtion to return a value the last command executed needs to return a value
    ## when a command is assigned to a value the function returns NULL 
    calc_geno_metrics(geno_df, cov_df)
    ## can also use return()
    ## metric_df <- calc_geno_metrics(geno_df, cov_df)
    ## return(metric_df) 
    ## in this case return() is not necessary
    ## but return() can be useful in other situtations
    ## for example if you want return different values when using a if statement
}