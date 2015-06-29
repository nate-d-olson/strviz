## Workflow for calculating sequece metrics
source("process_straitrazor_output.R")
source("get_counts.R")
source("calc_str_metrics.R")

seq_dir <- "H12-H-ZT80925_S96_L001_R1_001.fastq.STRaitRazor/PE_MiSeq_ALL/"

seq_df <- process_sequences_directory(seq_dir)
allele_counts_df <- str_allele_counts(seq_df)

str_allele_counts_test(seq_df, allele_counts_df)

allele_metrics_df <- calc_allele_metrics(allele_counts_df)

## next step evaluation allele_metrics_df