## Workflow for calculating sequece metrics
source("process_straitrazor_output.R")
source("get_counts.R")

seq_dir <- "H12-H-ZT80925_S96_L001_R1_001.fastq.STRaitRazor/PE_MiSeq_ALL/"

seq_df <- process_sequences_directory(seq_dir)
allele_counts_df <- str_allele_counts(seq_df)

str_allele_counts_test(seq_df, allele_counts_df)

#Creates new DF with only majority peaks 
new_seq <- seq_df
new_seq[is.na(new_seq)] <- 0
new_seq$Coverage_of_Majority_Peaks <- (new_seq$D2_R1 + new_seq$D2_R2 + new_seq$D1_R1 + new_seq$D1_R2)
new_seq2 <- new_seq  %>% 
    group_by(Locus)  %>%  
        top_n(2)

#preliminary code for genotype calls (.2 can be changed to any other ratio if needed)
new_seq2 <- new_seq2 %>% 
    group_by(Locus) %>% 
    mutate(Genotype = ifelse(min(Coverage_of_Majority_Peaks)/max(Coverage_of_Majority_Peaks) > 0.2, "Heterogyzous", "Homozygous"))

#preliminary code for PHR

new_seq2 <- new_seq2 %>% 
    group_by(Locus) %>% 
    mutate(PHR =  ifelse((min(Coverage_of_Majority_Peaks)/max(Coverage_of_Majority_Peaks))> .2,
                         round((min(Coverage_of_Majority_Peaks)/max(Coverage_of_Majority_Peaks)),4), NA_real_))

#preliminary code for Read Bias
ungroup(new_seq2)
new_seq2 <- new_seq2  %>%  
    mutate(Sum_R1 = (D2_R1 + D1_R1))
new_seq2 <- new_seq2  %>%  
    mutate(Sum_R2 = (D2_R2 + D1_R2))
new_seq2 <- new_seq2 %>%
    rowwise() %>%
    mutate(Read_Bias = min(c(Sum_R1, Sum_R2))/Coverage_of_Majority_Peaks)


#preliminary code for Strand Bias
ungroup(new_seq2)
new_seq2 <- new_seq2  %>%
    mutate(Sum_D1 = (D1_R1 + D1_R2))
new_seq2 <- new_seq2  %>%  
    mutate(Sum_D2 = (D2_R1 + D2_R2))
new_seq2 <- new_seq2 %>%
    rowwise() %>%
    mutate(Strand_Bias = min(c(Sum_D1, Sum_D2))/Coverage_of_Majority_Peaks)

#Prelimiary coverage of non majority peaks (need to figure out how to do this all in one data frame, and not include the stutter alleles)
new_seq3 <- new_seq  %>% 
group_by(Locus, Allele)  %>%  
mutate(Percentage_of_non_Majority_peaks = (sum(Coverage_of_Majority_Peaks)-max(Coverage_of_Majority_Peaks))/(sum(Coverage_of_Majority_Peaks)))
               