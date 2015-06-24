## geting seq counts for R1, R2 and directions 1 and 2
source("process_straitrazor_output.R")
source("revcom.R")
seq_dir <- "H12-H-ZT80925_S96_L001_R1_001.fastq.STRaitRazor/PE_MiSeq_ALL/"
seq_df <- process_sequences_directory(seq_dir) %>%  
            mutate(rev_seq = revcom(Sequence)) %>% 
            spread(Read, Allele_Calls)

seq_D2 <- seq_df %>% 
                rename(rev_seq = Sequence, 
                       Sequence = rev_seq, 
                       D2_R1 = R1, D2_R2 = R2)

seq_D1<- seq_df %>% rename(D1_R1 = R1, D1_R2 = R2)

seq_D1D2 <- left_join(seq_D1, seq_D2)


## Function to remove duplicate allele entries from count data frame
remove_duplicates <- function (str_count_df) {
  seq_count <- data_frame()
  for(i in unique(str_count_df$Locus)){
      locus_D1D2 <- str_count_df %>% filter(Locus == i)
      locus_count <- locus_D1D2[1,]
      for(j in 2:nrow(locus_count)){
          if(!(locus_D1D2[j,]$rev_seq %in% locus_count$Sequence)){
              locus_count <- bind_rows(locus_count, locus_D1D2[j,])
          }
      }
      seq_count <- bind_rows(seq_count, locus_count)
  }
}


#This is just how I made sure that the allele calls were all equal...I will eventually turn this into a function
seq_df[is.na(seq_df)] <- 0
seq_df <- mutate(seq_df, sum = R1+R2)
test1 <- seq_df  %>%  group_by(Locus, Allele)  %>%  summarise(sum(sum))

new_seq[is.na(new_seq)] <- 0
new_seq$sum <- (new_seq$D2_R1 + new_seq$D2_R2 + new_seq$D1_R1 + new_seq$D1_R2)
test2 <- (new_seq  %>%  group_by(Locus, Allele)  %>%  summarise(sum(sum)))

test3 <- bind_cols(test1, test2)




