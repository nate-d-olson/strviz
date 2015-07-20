## Testing out batch processing 
source("process_straitrazor_output.R")

sample_dirs <- list.dirs("STRaitRazor2.0output", recursive = F) %>% 
                    paste0("/PE_MiSeq_ALL")

## batch process fails
batch_process_samples(sample_dirs)
