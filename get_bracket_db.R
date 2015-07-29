library(dplyr)
setwd("/Users/sew/Desktop/strviz/Bracket_Notation")
file_names <- dir()
Bracket_Database <- do.call(rbind,lapply(file_names,read.csv))
Bracket_Database <- select(Bracket_Database, 3:6)
Bracket_Database["Forward"] <- "Forward"
Bracket_Database["Reverse"] <- "Reverse"
Bracket_Database[Bracket_Database==""]  <- NA 
Bracket_Database <- na.omit(Bracket_Database)
rownames(Bracket_Database) <- NULL
Bracket_Database <- Bracket_Database[,c(3,1,5,4,2,6)]
df <- Bracket_Database[,c(1:3)]
df2 <- Bracket_Database[,c(4:6)]
names(df)[1] <- "Sequence"
names(df2)[1] <- "Sequence"
names(df)[2] <- "Bracket Notation"
names(df2)[2] <- "Bracket Notation"
names(df)[3] <- "Direction"
names(df2)[3] <- "Direction"
Bracket_Database <- bind_rows(df, df2)
Bracket_Database <- unique(Bracket_Database)