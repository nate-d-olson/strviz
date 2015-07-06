## Workflow for calculating sequece metrics
source("process_straitrazor_output.R")
source("get_counts.R")
source("calc_str_metrics.R")

seq_dir <- "H12-H-ZT80925_S96_L001_R1_001.fastq.STRaitRazor/PE_MiSeq_ALL/"

seq_df <- process_sequences_directory(seq_dir)
allele_counts_df <- str_allele_counts(seq_df)


str_allele_counts_test(seq_df, allele_counts_df)

allele_metrics_df <- calc_allele_metrics(allele_counts_df)

allele_metrics <- allele_metrics_df[ -c(6:13) ]


allele_metrics <- (allele_metrics %>%
    group_by(Locus) %>% 
    mutate(make_na = ifelse(Genotype == "Homozygous" & Allele_Coverage == min(Allele_Coverage), 1, 0),
           Read_Bias = ifelse(make_na == 1, NA, Read_Bias),
           Strand_Bias = ifelse(make_na == 1, NA, Strand_Bias),
           Percentage_of_non_Majority_peaks = ifelse(make_na == 1, NA, Percentage_of_non_Majority_peaks)) %>% 
select(-make_na))

## next step evaluation allele_metrics_df

