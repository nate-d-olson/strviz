## Testing out batch processing 
source("process_straitrazor_output.R")

sample_dirs <- list.dirs("STRaitRazor2.0output", recursive = F) %>% 
                    paste0("/PE_MiSeq_ALL")

test_batch <- function (sample_dirs) {
    if (file.exists("summary.csv")){
        unlink("summary.csv")
    }
    for (dirs in sample_dirs) {    
        summary_met <- process_sample(dirs)
        write.table(summary_met, "summary.csv", append = TRUE, sep = ",", row.names = FALSE, col.names=!file.exists("summary.csv")) 
    }
}


test_batch(sample_dirs)

## Possible way to output specific summaries 
summary_read_bias <- function (sample_dirs) {
    for (dirs in sample_dirs) {    
        summary_met <- process_sample(dirs) 
        summary_met <- select(summary_met, Sample, Locus, Read_Bias)
        write.table(summary_met, "summary_read_bias.csv", append = TRUE, sep = ",", col.names=!file.exists("summary_read_bias.csv"), row.names = FALSE) 
    }
}

