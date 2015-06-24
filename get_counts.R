## geting seq counts for R1, R2 and directions 1 and 2
source("process_straitrazor_output.R")

seq_df <- process_sequences_directory("H12-H-ZT80925_S96_L001_R1_001.fastq.STRaitRazor/PE_MiSeq_ALL/") 
seq_df <- seq_df %>%  mutate(rev_seq = revcom(seq_df$Sequence)) %>% 


seq_comp2 <- seq_df %>% rename(rev_seq = Sequence, Sequence = rev_seq, D2_R1 = R1, D2_R2 = R2)
seq_comp<- seq_df %>% rename(D1_R1 = R1, D1_R2 = R2)

seq_count <- left_join(seq_comp2, seq_comp)



new_seq <- data_frame()
for(i in unique(seq_count$Locus)){
    locus_count <- seq_count %>% filter(Locus == i)
    locus_seq <- locus_count[1,]
    for(j in 2:nrow(locus_count)){
        if(!(locus_count[j,]$rev_seq %in% locus_seq$Sequence)){
            locus_seq <- bind_rows(locus_seq, locus_count[j,])
        }
    }
    new_seq <- bind_rows(new_seq, locus_seq)
}


#This is just how I made sure that the allele calls were all equal... will eventually turn this into a function

seq_df[is.na(seq_df)] <- 0
seq_df <- mutate(seq_df, sum = R1+R2)
test1 <- seq_df  %>%  group_by(Locus, Allele)  %>%  summarise(sum(sum))

new_seq[is.na(new_seq)] <- 0
new_seq$sum <- (new_seq$D2_R1 + new_seq$D2_R2 + new_seq$D1_R1 + new_seq$D1_R2)
test2 <- (new_seq  %>%  group_by(Locus, Allele)  %>%  summarise(sum(sum)))

test3 <- bind_cols(test1, test2)




