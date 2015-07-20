
peak_coverage2 <- function (allele_counts_df) {
    allele_counts_df[is.na(allele_counts_df)] <- 0
    allele_counts_df %>% 
        rowwise() %>% 
        mutate(Sum_R1 = (D2_R1 + D1_R1),
               Sum_R2 = (D2_R2 + D1_R2),
               Sum_D1 = (D1_R1 + D1_R2),
               Sum_D2 = (D2_R1 + D2_R2),
               Seq_Coverage = sum(D2_R1, D2_R2, D1_R1, D1_R2)) %>% 
        group_by(Allele,Locus) %>% 
        mutate(Allele_Coverage = sum(Seq_Coverage)) %>% 
        # changed to peak coverage as it represents the coverage for the 
        # peak and not the coverage of the majority peaks coverage
        group_by(Locus) %>% 
        mutate(Locus_Coverage = sum(Seq_Coverage)) %>% 
        top_n(n = 4, wt = Seq_Coverage) 
}

        
        # changed to peak coverage as it represents the coverage for the 
        # peak and not the coverage of the majority peaks coverage



library(data.table)
test <- group_by(test, Locus)
tes34 <- setDT(test)[order(-Seq_Coverage), list(Seq_Coverage=Seq_Coverage[1:2], 
                                        Stutter_Coverage=Seq_Coverage[3:4],
                                        Allele=Allele[1:2], 
                                        Sequence=Sequence[1:2],
                                        rev_seq=rev_seq[1:2],
                                        D1_R1 = D1_R1[1:2],
                                        D1_R2 = D1_R2[1:2],
                                        D2_R2 = D2_R2[1:2],
                                        D2_R1 = D2_R1[1:2],
                                        Sum_D1 = Sum_D1[1:2],
                                        Sum_D2 = Sum_D2[1:2],
                                        Sum_R1 = Sum_R1[1:2],
                                        Sum_R2 = Sum_R2[1:2],
                                        Size = Size[1:2],
                                        Allele_Coverage=Allele_Coverage[1:2],
                                        Locus_Coverage=Locus_Coverage[1:2]),
                      keyby = Locus]


