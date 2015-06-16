revcom <- function(sequence) {
    sapply(lapply(strsplit(sequence, NULL), rev), paste, collapse="")
    chartr("ATGC", "TACG", sequence)
    
}
