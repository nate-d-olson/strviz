
process_txtfile <- function(txt_file, read){
    df2 <- try(read.table(txt_file,sep =",", stringsAsFactors=F))
    if(class(df2)=='try-error') {
        warning(paste0(txt_file, "is empty"))
        return()
    }
    read.table(txt_file, sep=",", stringsAsFactors=F) %>% 
        head(-1) %>% 
        separate(V1, into = c("Locus", "Read", "Locus", "Allele", "Allele_Calls"), sep = ":")  
         
}


process_directory2 <- function(root_dir, paired = TRUE){
    df2 <- data_frame()
    if(paired){
        for(read in c("R1", "R2")){
            directory <- paste0(root_dir,"/", read,"/")
            for(txt_file in list.files(directory, "txt",full.names = TRUE)){
         df2 %<>% bind_rows(process_txtfile(txt_file,read))
            }
            
        }
    }else{
        #%%TODO%%
        warning("No code for unpaired read data")
    }
    df
}





