##################################################### Required packages
suppressMessages(library(ShortRead))
suppressMessages(library(randomForest))
suppressMessages(library(reshape2))
suppressMessages(library(pROC))

##################################################### Choose datasets here
args <- commandArgs(TRUE)
dataset_name <- args[1]
pwm_1_name <- args[2]
pwm_2_name <- args[3]

print(paste0("Data set name: ", dataset_name))
print(paste0("TF1: ", pwm_1_name))
print(paste0("TF2: ", pwm_2_name))

##################################################### Place the data set here
positive1 <- clean(sread(readFasta(paste0('./DNA/CAP-SELEX-Positive-Training/', dataset_name, '.fasta.gz'))))
positive2 <- reverseComplement(positive1)
negative1 <- clean(sread(readFasta(paste0('./DNA/CAP-SELEX-Negative-Training/', dataset_name, '.fasta.gz'))))
negative2 <- reverseComplement(negative1)

if(width(positive1)[1] < 25){
  stop("Sequence length should be >25")
} else if(length(unique(c(width(positive1), width(negative1)))) != 1){
  stop("Sequences must have the same length")
} else {
  print("Sequences are >25 nucleic acids long and have the same length")
  print("Continue to training sequence pre-search")
}

##################################################### Select TF1 here
pwm1 <- read.table(paste0('./PWM_HT-SELEX/', pwm_1_name, '.txt'), sep = "", row.names = 1)
colnames(pwm1) <- c(1:ncol(pwm1))

##################################################### Select TF2 here
pwm2 <- read.table(paste0('./PWM_HT-SELEX/', pwm_2_name, '.txt'), sep = "", row.names = 1)
colnames(pwm2) <- c(1:ncol(pwm2))

########################################################################################################### Functions
find_short_seq_pwm12 <- function(sequence, reverse_sequence, pwm1, pwm2){
  short_seqs <- character(0)
  
  for(i in 1:length(sequence)){ 
    seq1 <- strsplit(as.character(sequence[i]), "")[[1]]
    seq2 <- strsplit(as.character(reverse_sequence[i]), "")[[1]]
    
    match1 <- matrix(rownames(pwm1), 4, ncol(pwm1))
    match2 <- matrix(rownames(pwm2), 4, ncol(pwm2))
    
    values1 <- matrix(0, 1, length(seq1)-ncol(pwm1)+1)
    for(n in 1:(length(seq1)-ncol(pwm1)+1)){
      seq <- seq1[n:(n+ncol(pwm1)-1)]
      values1[n] <- prod(pwm1[rbind(seq, seq, seq, seq) == match1])
    }
    
    values2 <- matrix(0, 1, length(seq1)-ncol(pwm2)+1)
    for(n in 1:(length(seq1)-ncol(pwm2)+1)){
      seq <- seq1[n:(n+ncol(pwm2)-1)]
      values2[n] <- prod(pwm2[rbind(seq, seq, seq, seq) == match2])
    }
    
    values3 <- matrix(0, 1, length(seq2)-ncol(pwm1)+1)
    for(n in 1:(length(seq2)-ncol(pwm1)+1)){
      seq <- seq2[n:(n+ncol(pwm1)-1)]
      values3[n] <- prod(pwm1[rbind(seq, seq, seq, seq) == match1])
    }
    
    values4 <- matrix(0, 1, length(seq2)-ncol(pwm2)+1)
    for(n in 1:(length(seq2)-ncol(pwm2)+1)){
      seq <- seq2[n:(n+ncol(pwm2)-1)]
      values4[n] <- prod(pwm2[rbind(seq, seq, seq, seq) == match2])
    }
    
    ## Sequence score
    value <- which.max(c(max(values1), max(values2), max(values3), max(values4)))
    if(value == 1){
      n <- which.max(values1)
      chosen_sequence <- seq1
      chosen_pwm <- pwm1
    }
    if(value == 2){
      n <- which.max(values2)
      chosen_sequence <- seq1
      chosen_pwm <- pwm2
    }
    if(value == 3){
      n <- which.max(values3)
      chosen_sequence <- seq2
      chosen_pwm <- pwm1
    }
    if(value == 4){
      n <- which.max(values4)
      chosen_sequence <- seq2
      chosen_pwm <- pwm2
    }
    
    short_seqs[i] <- paste(chosen_sequence[n:(n+ncol(chosen_pwm)-1)], collapse = "")
  }
  
  return(short_seqs)
}

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
    
    s1[] <- lapply(s1, factor)
    s2[] <- lapply(s2, factor)
    
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

########################################################################################################### Create training and testing sets 
labels <- as.factor(c(rep(1, length(positive1)), rep(0, length(negative1))))

selex <- c(positive1, negative1)
reverse <- c(positive2, negative2)

########################################################################################################### Train ComBind
##################################################### Extend PWMs to 25
n1 <- 25-ncol(pwm1)
n2 <- 25-ncol(pwm2)

mat <- matrix(data=c(0.25,0.25,0.25,0.25), 4, 1)
ext1 <- t(apply(mat, MARGIN=1, function(x) rep(x, n1)))
ext2 <- t(apply(mat, MARGIN=1, function(x) rep(x, n2)))

pwm1_x <- cbind(pwm1, ext1)
pwmx_1 <- cbind(ext1, pwm1)

pwm2_x <- cbind(pwm2, ext2)
pwmx_2 <- cbind(ext2, pwm2)

colnames(pwm1_x) <- 1:25
colnames(pwmx_1) <- 1:25
colnames(pwm2_x) <- 1:25
colnames(pwmx_2) <- 1:25 

pwm1 <- NULL
pwm2 <- NULL

###################################################### Find sequences and represent as a data frame
position1 <- find_short_seq_pwm12(selex, reverse, pwm2_x, pwmx_1) ## Position PWM2 + PWM1
print("Position 1 sequences pre-searched")
position2 <- find_short_seq_pwm12(selex, reverse, pwm1_x, pwmx_2) ## Position PWM1 + PWM2
print("Position 2 sequences pre-searched")

train1 <- substr(position1, 1, 1)
for(i in 2:25){
  train1 <- cbind(train1, substr(position1, i, i))
}
colnames(train1) <- c(1:25)
train1 <- as.data.frame(train1)
train1[] <-  lapply(train1, factor, levels=c("A", "T", "G", "C"))

train2 <- substr(position2, 1, 1)
for(i in 2:25){
  train2 <- cbind(train2, substr(position2, i, i))
}
colnames(train2) <- c(1:25)
train2 <- as.data.frame(train2) 
train2[] <-  lapply(train2, factor, levels=c("A", "T", "G", "C"))

#################################################### Choose parameters with validation set using balanced classes
print("Select random forest parameters")
set.seed(1001)
tr.idx = c(sample(sum(labels==1), size=round(sum(labels==1) * 0.75)), sample((sum(labels==1) + 1:sum(labels==0)), size=round(sum(labels==0) * 0.75)))
te.idx = setdiff(1:nrow(train1), tr.idx)

train_position1_validation <- train1[tr.idx, ]
train_position2_validation <- train2[tr.idx, ]

labels_train_validation <- labels[tr.idx]
labels_test_validation <- labels[te.idx]

terminal_nodesize <- c(1, 5, 10, 15)
variables_at_split <- round(c(0.1, 0.2, 0.3, 0.4, 0.5)*ncol(pwm1_x))
grid <- matrix(0, nrow = 5, ncol = 4)

for(i in 1:5){
  for(n in 1:4){
    rf_1 = randomForest(train_position1_validation, labels_train_validation, ntree=200, mtry = variables_at_split[i], nodesize = terminal_nodesize[n])
    rf.pred1 = test_all_sites(selex[te.idx], reverse[te.idx], rf_1, (ncol(pwm1_x)))
    
    rf_2 = randomForest(train_position2_validation, labels_train_validation, ntree=200, mtry = variables_at_split[i], nodesize = terminal_nodesize[n])
    rf.pred2 = test_all_sites(selex[te.idx], reverse[te.idx], rf_2, (ncol(pwmx_1)))
    
    roc_obj <- roc(labels_test_validation, rowMeans(cbind(rf.pred1, rf.pred2)), quiet = T)
    grid[i, n] <- auc(roc_obj)
    
    rm(rf_1)
    rm(rf_2)
  }
}

parameters <- which(grid == max(grid), arr.ind = TRUE)
var <- variables_at_split[parameters[1]]
node <- terminal_nodesize[parameters[2]]

var_ComBind <- c(0.1, 0.2, 0.3, 0.4, 0.5)[parameters[1]]
node_ComBind <- node

## Remove validation data sets
rm(train_position1_validation)
rm(train_position2_validation)
print("Random forest parameters selected")
print(paste0("Varianbles at split: ", var_ComBind*10, "% of sequence length"))
print(paste0("Terminal node size: ", node_ComBind))

####################################################### Final model
rf_1 <- randomForest(train1, labels, ntree=200, mtry = var, nodesize = node)
rf_2 <- randomForest(train2, labels, ntree=200, mtry = var, nodesize = node)

#### Save random forests
saveRDS(rf_1, paste0("./RF/", dataset_name, "_1.rds.gz"), compress = "gzip")
saveRDS(rf_2, paste0("./RF/", dataset_name, "_2.rds.gz"), compress = "gzip")

print("Random forests trained and downloaded")