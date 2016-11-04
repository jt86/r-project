library(e1071)
set.seed(0)

# args = commandArgs(trailingOnly=TRUE)
# i <- commandArgs()[1]
# i <- 1


# number_selected <- commandArgs()[2]
# datasetnum <- commandArgs()[3]
# fold <- commandArgs()[4]
# seed <- commandArgs()[5]
# number_selected

number_selected <- 500
datasetnum <- 39
fold <- 3
seed <- 4
number_selected

tuneParametersVal <- function(X, Y, cost  =  exp(seq(-3, 7, length = 5)),  
                              gamma = exp(seq(-3, 3, length = 5)), percentage_of_data_for_validation = 0.1) {
  
  mat <- as.matrix(dist(X))^2 #
  med <- median(mat[ upper.tri(mat) ])
  
  gamma <- gamma / med
  
  i_pos <- which(Y == 1)
  i_neg <- which(Y == -1)
  half <- round(percentage_of_data_for_validation * length(Y))
  
  
  errors <- matrix(0, length(gamma), length(cost))
  
  for (fold in 1 : 5) {
    sample_pos = sample(i_pos, half)
    sample_neg = sample(i_neg, half)
    s = sample(c(sample_pos, sample_neg))
    Xval <- X[ s, ]
    Yval <- Y[ s ]
    
    s_remain <- (1 : nrow(X))[ -c(s) ]
    Xtrain <- X[ s_remain, ]
    Ytrain <- Y[ s_remain ]
    
    for (i in 1 : length(gamma))	
      for (j in 1 : length(cost)) {
        model <- svm(Xtrain, as.factor(Ytrain), kernel = "radial", gamma = gamma[ i ], cost = cost[ j ],
                     scale = FALSE)
        errorTest <- mean(predict(model, Xval) != Yval)
        
        errors[ i, j] <- errors[ i, j] + errorTest
        cat(".")
      }
    cat("\n")
  }
  
  errors <- errors/5.0
  bestCost <- NULL
  bestGamma <- NULL
  bestError <- +Inf
  
  
  for (i in 1 : length(gamma))	
    for (j in 1 : length(cost)) {
      
      if (bestError >= errors[ i, j ]) {
        bestError <- errors[ i, j]
        bestCost <- cost[ j ]
        bestGamma <- gamma[ i ]
      }
      
    }
  #retrain the model
  model <- svm(X, as.factor(Y), kernel = "radial", gamma = bestGamma, cost = bestCost, scale = FALSE)	
  list(model = model, bestError = bestError, bestGamma = bestGamma, bestCost = bestCost)
}

#number_selected <- 500
# datasetnum <- 39
# fold <- 3
# seed <- 4


#Loading the data
dataset <- (read.table(paste("/Volumes/LocalDataHD/j/jt/jt306/Documents/CVPR2016_Rcode/saved_datasets/tech",datasetnum,"data",sep = "")))
dim(dataset)
train_indices <- as.integer(read.table(paste('/home/j/jt/jt306/Documents/CVPR2016_Rcode/saved_indices/top',number_selected,'tech',datasetnum,fold,seed,'train_instances_indices', sep='-'))[[1]])+1
test_indices <- as.integer(read.table(paste('/home/j/jt/jt306/Documents/CVPR2016_Rcode/saved_indices/top',number_selected,'tech',datasetnum,fold,seed,'test_instances_indices', sep='-'))[[1]])+1
selected_indices <- as.integer(read.table(paste('/home/j/jt/jt306/Documents/CVPR2016_Rcode/saved_indices/top',number_selected,'tech',datasetnum,fold,seed,'selected_feat_indices', sep='-'))[[1]])+1
priv_indices <- as.integer(read.table(paste('/home/j/jt/jt306/Documents/CVPR2016_Rcode/saved_indices/top',number_selected,'tech',datasetnum,fold,seed,'unselected_feat_indices', sep='-'))[[1]])+1




Xtrain <- dataset[train_indices,selected_indices]
Xtest <- dataset[test_indices,selected_indices]

print(typeof(dataset))
print(typeof(Xtest))
print(typeof(Xtrain))


train_priv <- dataset[train_indices,priv_indices]
dim(train_priv)
test_priv <- dataset[test_indices,selected_indices]
dim(test_priv)

labels <-(as.double(unlist(read.table("/Volumes/LocalDataHD/j/jt/jt306/Documents/CVPR2016_Rcode/saved_datasets/tech39labels"))))
Ytrain <- labels[train_indices]
Ytest <- labels[test_indices]
dim(Ytrain)
length(Ytrain)
length(Ytest)





#zero mean unit variance normalization
meanTrain <- apply(Xtrain, 2, mean)
sdTrain <- apply(Xtrain, 2, sd)
sdTrain[ sdTrain == 0 ] <- 1
Xtrain <- (Xtrain - matrix(meanTrain, nrow(Xtrain), ncol(Xtrain), byrow = TRUE)) / 
  matrix(sdTrain, nrow(Xtrain), ncol(Xtrain), byrow = TRUE)
Xtest <- (Xtest - matrix(meanTrain, nrow(Xtest), ncol(Xtest), byrow = TRUE)) / 
  matrix(sdTrain, nrow(Xtest), ncol(Xtest), byrow = TRUE)

#SVM baseline
time <- system.time( ret <- tuneParametersVal(Xtrain, Ytrain))
print(ret)
resultRet <- ret
ret <- ret$model

errorTest <- mean(predict(ret, Xtest) != Ytest)

# write.table(errorTest, file = paste("./results/SVM/",i,"_errorTest_X_val.txt",sep = ""), row.names = F, col.names = F, append = FALSE)
# write.table(t(time), file = paste("./results/SVM/",i,"_time_X_val.txt",sep = ""), row.names = F, col.names = F, append = FALSE)

write.table(errorTest, file = paste("/home/j/jt/jt306/Documents/CVPR2016_Rcode/joe_results/SVM/error",number_selected,'tech',datasetnum,fold,seed,".txt",sep = "-"), row.names = F, col.names = F, append = FALSE)
write.table(t(time), file = paste("/home/j/jt/jt306/Documents/CVPR2016_Rcode/joe_results/SVM/time",number_selected,'tech',datasetnum,fold,seed,".txt",sep = "-"), row.names = F, col.names = F, append = FALSE)

