## geting seq counts for R1, R2 and directions 1 and 2

## Find the reverse complement of a DNA sequence
revcom <- function(sequence) {
    reverse <- sapply(lapply(strsplit(sequence, NULL), rev), paste, collapse="") 
    complement <- chartr("ATGC", "TACG", reverse)
    return(complement)
    
}

## Function to remove duplicate allele entries from count data frame
remove_duplicates <- function (str_count_df) {
    seq_count <- data_frame()
    for(i in unique(str_count_df$Locus)){
        locus_D1D2 <- str_count_df %>% filter(Locus == i)
        locus_count <- locus_D1D2[1,]
        for(j in 2:nrow(locus_D1D2)){
            if(!(locus_D1D2[j,]$Sequence %in% locus_count$rev_seq)){
                locus_count <- bind_rows(locus_count, locus_D1D2[j,])
            }
        }
        seq_count <- bind_rows(seq_count, locus_count)
    }
    return(seq_count)
}

## Calculate Allele counts for R1 and R2 forward and reverse sequences, input from 
str_allele_counts <- function (seq_df) {
  seq_df <- seq_df %>%  
              mutate(rev_seq = revcom(Sequence)) %>% 
              spread(Read, Allele_Calls)
  
  seq_D2 <- seq_df %>% 
                  rename(rev_seq = Sequence, 
                         Sequence = rev_seq, 
                         D2_R1 = R1, D2_R2 = R2)
  
  seq_D1 <- seq_df %>% rename(D1_R1 = R1, D1_R2 = R2)
  
  left_join(seq_D1, seq_D2) %>% 
      remove_duplicates() %>% 
      return()
}


### Testing counts
str_allele_counts_test <- function(seq_df, allele_counts_df){
    seq_sum <- seq_df  %>%  
                group_by(Locus, Allele)  %>%  
                summarise(count_seq = sum(Allele_Calls))
    allele_count_sum <- allele_counts_df %>% 
                            gather(RD, Allele_Calls, 6:9) %>% 
                            group_by(Locus, Allele)  %>%  
                            summarise(count_allele = sum(Allele_Calls, na.rm = TRUE))
    test_df <- full_join(seq_sum, allele_count_sum) %>% 
        mutate(test = count_allele != count_seq)
    if(sum(test_df$test) > 0){
        return(test)
    }else{
        print("Counts Match")
    }
}





