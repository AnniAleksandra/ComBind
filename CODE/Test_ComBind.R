##################################################### Required packages
suppressMessages(library(ShortRead))
suppressMessages(library(randomForest))
suppressMessages(library(reshape2))
suppressMessages(library(pROC))

##################################################### Choose datasets here
args <- commandArgs(TRUE)
dataset_name <- args[1]
RF_name <- args[2]
print(paste0("Data set name: ", dataset_name))
print(paste0("RF names: ", RF_name))

##################################################### Place the data set here
sequence <- sread(readFasta(paste0('./DNA/SCORE-DNA/', dataset_name, '.fasta.gz')))
reverse <- reverseComplement(sequence)

if(width(sequence)[1] < 25){
  stop("Sequence length should be >25")
} else if(length(unique(width(sequence))) != 1){
  stop("Sequences must have the same length")
} else if(sum(rowSums(letterFrequency(sequence, c("A", "T", "C", "G"))) != width(sequence)) > 0) {
  stop("Ambigious DNA is not allowed -> Make sure that sequences entail only A, T, C and G")
} else {
  print("Sequences are >25 nucleic acids long, contain no ambigious letters and have the same length")
  print("Continue to sequence scoring")
}

########################################################################################################### Functions
test_all_sites <- function(sequence, reverse_sequence, rf, l){
  div_pos1 <- colsplit(string=as.character(sequence), pattern="", names = c(1:(width(sequence)[1])))
  div_pos2 <- colsplit(string=as.character(reverse_sequence), pattern="", names = c(1:(width(reverse_sequence)[1])))
  
  div_pos1[div_pos1 == TRUE] <- "T"
  div_pos2[div_pos2 == TRUE] <- "T"
  
  maximum_index <- width(sequence)[1] - l + 1
  probabilities <- matrix(0, length(sequence), maximum_index)
  
  for(i in 1:maximum_index){
    s1 <- div_pos1[,i:(i+l-1)]
    s2 <- div_pos2[,i:(i+l-1)]
    colnames(s1) <- c(1:l)
    colnames(s2) <- c(1:l)
    
    s1[] <- lapply(s1, factor, levels=c("A", "T", "G", "C"))
    s2[] <- lapply(s2, factor, levels=c("A", "T", "G", "C"))
    
    pp1 <- predict(rf, newdata = s1, type = 'prob')[, 2]
    pp2 <- predict(rf, newdata = s2, type = 'prob')[, 2]
    
    ind <- which(pp1 < pp2)
    pp <- pp1
    pp[ind] <- pp2[ind]
    probabilities[, i] <- pp
  }
  
  prob <- apply(probabilities, 1, max)
  return(prob)
}

########################################################################################################### Test with ComBind
rf_1 <- readRDS(paste0("./RF/", RF_name, "_1.rds.gz"))
rf_2 <- readRDS(paste0("./RF/", RF_name, "_2.rds.gz"))

rf.pred1 <- test_all_sites(sequence, reverse, rf_1, 25)
rf.pred2 <- test_all_sites(sequence, reverse, rf_2, 25)

score <- rowMeans(cbind(rf.pred1, rf.pred2))
predictions <- cbind(as.character(sequence), score)

#################################################################################### Save predictions
write.table(predictions, paste0("./SCORE/", dataset_name, ".txt"), sep = "\t", row.names = F, col.names = F, quote = F)
print("Sequences scored and downloaded to 'SCORE' folder")



