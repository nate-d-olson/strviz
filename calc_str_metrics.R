library(stringr)


## Description: Creates new data frame with all the coverage values, including the sums of reads 
##              and direction as well as sequence, locus, and allele sums. 
## Input: A data frame of the counts split up by read and direction.
## Output: One dataframe will all various coverage values.

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
        group_by(Locus) %>% 
        mutate(Locus_Coverage = sum(Seq_Coverage)) 
}

## Description: Assigns a genotype call to each locus based on the peak height 
##              ratio threshold and uses call to calculate the majority peaks coverage.
## Input: Data frame with various coverage values - 
##        cutoff for the ratio between top two allele calls,
##        default value 0.2
## Output: Data frame with only the loci and their respective genotypes.

genotype_call <- function (cov_df, gt_threshold = 0.2) {
    cov_df %>% 
      top_n(n = 2, wt = Seq_Coverage) %>% 
      group_by(Locus)
    cov_df <- cov_df[!(cov_df$Locus_Coverage==0),]
    cov_df <- mutate(Genotype = ifelse(min(Seq_Coverage)/max(Seq_Coverage) > gt_threshold, 
                               "Heterozygous", "Homozygous")) %>% 
    ##Need if statement if seq_coverage for both are equal   
        select(Locus, Genotype) %>% 
        unique ()
}

## Description: Calculates stutter
## Input: Data frame with coverage values
## Output: Data frame with additional column containin the stutter,
##         the total allele coverage for the previous allele in the locus 

stutter <- function(cov_df) {
    cov_df %>% 
    group_by(Locus) %>%    
        mutate(Stutter_1 = lag(Allele_Coverage)) 
}

## Description: Adjusts stutter for heterzygous loci
## Input: Data frame with stutter calculated for hererozygous loci. 
## Output: Makes the higher of two consecutive alleles N/A, and makes all stutter for AMEL N/A
stutter_het_adj <- function(stutter_df) {
    stutter_df %>% 
        group_by(Locus) %>% 
           mutate(Stutter_1 = ifelse(max(Allele)-min(Allele)==1 & 
              Allele==max(Allele), NA_integer_, Stutter_1)) %>% 
           mutate(Stutter_1 = ifelse(Locus=="AMEL", NA_real_, Stutter_1) )
}
 
## Description: Cacluates sutter ratio
## Input: Data frame that includes the stutter as well as coverage metrics  
## Output: Data frame with added column for stutter ratio, 
##         the stutter divided be the allele coverage.

stutter_rat <- function(stutter_df) {
    stutter_df %>% 
        mutate(Stutter_Ratio = Stutter_1/Allele_Coverage ) 
}

## Description: Assigns a peak height ratio to each heterozygous locus based on the top two peaks.
## Input: Data frame with only top two allele counts per locus, and genotype call.
## Output: Data frame with extra column containing the peak heigh ratio for heterzygous loci, 
##         and NA for homozygous loci. One PHR per locus.

peak_height_ratio <- function (genotype_call_df, gt_threshold = 0.2) {
      genotype_call_df %>% 
        group_by(Locus) %>% 
        mutate(Peak_Height_Ratio =  ifelse((min(Seq_Coverage)/max(Seq_Coverage))> gt_threshold,
                           round((min(Seq_Coverage)/max(Seq_Coverage)),4), NA_real_))
}

## Description: Assigns a read bias to each allele (two for heterzygotes and
##              only the majority peak for homozygotes). 
## Input: Data frame with top two allele counts per locus, genotype call, 
##        and PHR per loci. 
## Output: Data frame with extra column containing the read bias for each allele.

read_bias <- function(peak_height_ratio_df) {
    ungroup(peak_height_ratio_df) %>% 
      rowwise() %>%
      mutate(Read_Bias = min(c(Sum_R1, Sum_R2))/Seq_Coverage)
}

## Description: Assigns a strand bias to each allele (two for heterzygotes and only the majority peak for homozygotes).
## Input: Data frame with top two allele counts per locus, genotype call, PHR per loci, and read bias per allele.
## Output: Data frame with extra column containing the strand bias for each allele. 

strand_bias <- function (read_bias_df) {
  read_bias_df %>% 
      rowwise() %>%
      mutate(Strand_Bias = min(c(Sum_D1, Sum_D2))/Seq_Coverage)
}

## Description: Calculates the ration of non majority peaks over majority peaks for largest sequence per allele
## Input: Data frame with coverage values
## Output: Data frame with additional column containin the stutter,
##         the total allele coverage for the previous allele in the locus.

non_maj_peaks <- function (cov_df) {
    cov_df %>% 
        rowwise() %>% 
        mutate(Percentage_of_non_Majority_peaks = 
                   (Allele_Coverage - Seq_Coverage)/Allele_Coverage)
}

## Description: This function is used to calculate all the various metris for all the heterozygous loci.
## Input: A data frame with the list of genotypes and a data frame with all the coverage information.
## Output: A data frame with all the metrics for all heterozygous loci.

calc_het_metrics <- function(geno_df, cov_df){
    het_cov_df <- filter(geno_df, Genotype == "Heterozygous") %>% 
        left_join(cov_df)  
    het_cov_df$Allele <- het_cov_df$Allele %>% 
        str_replace("X", "100") %>% 
        str_replace("Y", "1010") %>%  
        as.numeric() 
    het_cov_df <- het_cov_df %>% 
        stutter() %>% 
        top_n(2, wt= Seq_Coverage)  %>%
         stutter_het_adj() %>%         
         stutter_rat()  
    het_cov_df$Allele <- het_cov_df$Allele %>% 
        as.character() %>% 
        str_replace("100", "X") %>% 
        str_replace("1010", "Y") 
    het_cov_df <- het_cov_df %>% 
        peak_height_ratio() %>% 
        read_bias() %>% 
        strand_bias() %>% 
        non_maj_peaks () 
}

## Description: This function is used to calculate all the various metris for all the homozygous loci.
## Input: A data frame with the list of genotypes and a data frame with all the coverage information.
## Output: A data frame with all the metrics for all homozygous loci.

calc_homo_metrics <- function(geno_df, cov_df){
    homo_cov_df <- filter(geno_df, Genotype == "Homozygous") %>% 
        left_join(cov_df) 
    homo_cov_df$Allele <- homo_cov_df$Allele %>% 
        str_replace("X", "100") %>% 
        as.numeric()
     homo_cov_df <- homo_cov_df %>% 
        stutter() %>%  
        top_n(1, wt = Seq_Coverage)  %>% 
        stutter_rat()  
     homo_cov_df$Allele <- homo_cov_df$Allele %>% 
         as.character() %>% 
         str_replace("100", "X")
     homo_cov_df <- homo_cov_df %>%   
        peak_height_ratio() %>% 
        read_bias() %>% 
        strand_bias() %>% 
        non_maj_peaks () %>% 
        mutate(Allele = as.character(Allele))  
     homo_cov_df <- homo_cov_df %>% 
        bind_rows(homo_cov_df) %>% 
        arrange(Locus) %>% 
        mutate_each(funs(./2), D1_R1:Locus_Coverage)
    homo_cov_df$Peak_Height_Ratio[homo_cov_df$Peak_Height_Ratio == 1] <- NA
    homo_cov_df     
}

## Description: Calcualtes metrics for both heterozygous and homozygous loci.
## Input: Two data frames, one with all coverage values, and one with all loci and their respective genotypes.
## Output: Data frame with all coverage counts and summary metrics for all loci. 

calc_geno_metrics <- function(geno_df, cov_df){
    homo_df <- calc_homo_metrics(geno_df, cov_df) 
    het_df <- calc_het_metrics(geno_df, cov_df)
    bind_rows(homo_df, het_df)
}

## Description: Created two data frames, one with all coverage values, and one with all loci and their respective genotypes, 
##              and then calcuates all metrics using calc_geno_metrics.
## Input:  Data frame with all the allele counts.
## Output: Data frame with all coverage counts and summary metrics for all loci.

calc_allele_metrics <- function(allele_counts_df){
    cov_df <- coverage_calc(allele_counts_df)  
    geno_df<- genotype_call(cov_df)
    calc_geno_metrics(geno_df, cov_df)
}