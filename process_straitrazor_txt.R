library(dplyr)
library(tidyr)
library(magrittr)

process_allelecalls_file <- function(txt_file, read){
    # no need to name differently from process_seqeunces as they are in 
    # different environments, they can't see eachother
    df <- try(read.table(txt_file,sep =",", stringsAsFactors=F))
    if(class(df)=='try-error') {
        warning(paste0(txt_file, "is empty"))
        return()
    }
    read.table(txt_file, sep=",", stringsAsFactors=F) %>% 
        head(-1) %>% 
        separate(V1, into = c("Locus", "Read", "Locus", "Allele", "Allele_Calls"), sep = ":")  
         
}

## changed name to be more descriptive
process_allelecalls_directory <- function(root_dir, paired = TRUE){
    df <- data_frame()
    if(paired){
        for(read in c("R1", "R2")){
            directory <- paste0(root_dir,"/", read,"/")
            for(txt_file in list.files(directory, "txt",full.names = TRUE)){
                df %<>% bind_rows(process_allelecalls_file(txt_file,read))
            }
            
        }
    }else{
        #%%TODO%%
        warning("No code for unpaired read data")
    }
    df
}

process_allelecalls_directory("H12-H-ZT80925_S96_L001_R1_001.fastq.STRaitRazor/PE_MiSeq_ALL/")


#mutate a column names sample
