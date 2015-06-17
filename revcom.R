revcom <- function(sequence) {
    reverse <- sapply(lapply(strsplit(sequence, NULL), rev), paste, collapse="") 
    complement <- chartr("ATGC", "TACG", reverse)
    return(complement)
    
}
