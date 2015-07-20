## Workflow for calculating sequence metrics
source("process_straitrazor_output.R")
source("get_counts.R")
source("calc_str_metrics.R")


# seq_dir <- "H12-H-ZT80925_S96_L001_R1_001.fastq.STRaitRazor/PE_MiSeq_ALL/"

process_sample_dir <- function (seq_dir) {
  seq_df <- process_sequences_directory(seq_dir)
  allele_counts_df <- str_allele_counts(seq_df)
  calc_allele_metrics(allele_counts_df) %>% mutate(sample = seq_dir)
}
