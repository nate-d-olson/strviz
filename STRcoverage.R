setwd("/Users/sew/Downloads/H12-H-ZT80925_S96_L001_R1_001.fastq.STRaitRazor copy/PE_MiSeq_ALL/R1")v
loc_dat <- read.table(file="Amel.sequences", header=FALSE, sep="\t")
library(dplyr)
library(tidyr)

loc_dat2 <- loc_dat %>% separate(V1, into = c("Locus", "Allele"), sep = ":")
loc_dat2 <- rename(loc_dat2, "Allele Calls" = V4)
loc_dat2 <- loc_dat2[,-c(3,4)]

loc2_dat <- read.table(file="Amel.Allelecalls.txt", header=FALSE, sep="\t")
loc2_dat2 <- head(loc2_dat, -1) 
loc2_dat2 <- loc2_dat2 %>% separate(V1, into = c("Locus", "Read", "Locus", "Allele", "Allele Calls"), sep = ":")
loc2_dat2 <- loc2_dat2[,-c(2:3)]
loc2_dat2[loc2_dat2==""]  <- 0 
loc2_dat2$"Allele Calls"<-as.integer(loc2_dat2$"Allele Calls")
loc2_dat2 <- loc2_dat2[loc2_dat2$"Allele Calls" == 0,]

loc_dat3 <- bind_rows(loc_dat2, loc2_dat2)


