## Testing out batch processing 
source("process_straitrazor_output.R")

sample_dirs <- list.dirs("STRaitRazor2.0output", recursive = F) %>% 
                    paste0("/PE_MiSeq_ALL")

test_batch <- function (sample_dirs) {
    for (dirs in sample_dirs) {    
        summary_met <- process_sample(dirs)
        ## this writes to the same file overwriting existing data 
        ## use append = TRUE to add to the file - need to be careful with this
        ## as it will add to the file even if it was from a previous run
        ## think about how to avoid this issue
        write.table(summary_met, "summary.csv", append = TRUE, sep = ",", row.names = FALSE, col.names=!file.exists("summary.csv")) 
    }
}


## batch process fails
test_batch(sample_dirs)


summary_read_bias <- function (sample_dirs) {
    for (dirs in sample_dirs) {    
        summary_met <- process_sample(dirs) 
        summary_met <- select(summary_met, Sample, Locus, Read_Bias)
        ## this writes to the same file overwriting existing data 
        ## use append = TRUE to add to the file - need to be careful with this
        ## as it will add to the file even if it was from a previous run
        ## think about how to avoid this issue
        write.table(summary_met, "summary_read_bias.csv", append = TRUE, sep = ",", col.names=!file.exists("summary_read_bias.csv"), row.names = FALSE) 
    }
}

