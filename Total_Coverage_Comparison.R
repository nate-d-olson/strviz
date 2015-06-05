str_dat <- read.csv(file="BracketNotation_A4-C-WA29594_S25_L001_R1_001.fastq.STRaitRazor_output.csv", header=TRUE, sep=",")

library(dplyr)
library(ggplot2)

tbl_df(str_dat)

str_dat2 <- slice(str_dat, 1:44)

str_dat2$Locus <- as.character(str_dat2$Locus)

str_dat2$Coverage.of.majority.peaks <- as.numeric(as.character(str_dat2$Coverage.of.majority.peaks))

str_dat2$Coverage.of.majority.peaks[is.na(str_dat2$Coverage.of.majority.peaks)] <- 0

str_dat2 <- str_dat2  %>%  group_by(Locus)  %>%  mutate(sum=sum(Coverage.of.majority.peaks))

ggplot(data=str_dat2, aes(x=Locus, y=sum, colour=Locus)) + 
  geom_point(position = "identity", size=3) + 
  guides(colour=FALSE) +
  ylab("Total Coverage") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  xlab("Locus") +
  ggtitle("Total Coverage Comparison") +
  theme(aspect.ratio = .125) +
  geom_hline(aes(yintercept=median(str_dat2$sum)))
