source("./Rcode/gFITC.R")
source("./Rcode/epGPC.R")
set.seed(0)

i <- commandArgs()[1]

number_selected <- commandArgs()[6]
datasetnum <- commandArgs()[7]
seed <- commandArgs()[8]
fold <- commandArgs()[9]
print(number_selected)
print(datasetnum)
print(fold)
print(seed)


#number of pseudo-inputs is set to 100
for (npi in c(100)) {
	
		print(i)

		#Loading the data
		# load(paste("./data/",i,"data.dat",sep = ""))
		# Xtrain <- data$x[ data$itrain, ]
		# Ytrain <- as.vector(c(data$y_train))
		# Xtest <- data$x[ data$itest, ]
		# Ytest <- as.vector(c(data$y_test))
    
    #Loading the data
    dataset <-(read.table(paste("../saved_datasets/tech",datasetnum,"data",sep = "")))
    dims=dim(dataset)
    dataset <-as.numeric(unlist(dataset))
    dataset<-matrix(dataset,dims[1],dims[2])
    train_indices <- as.integer(read.table(paste('../saved_indices/top',number_selected,'-tech',datasetnum,'-',fold,'-',seed,'-train_instances_indices', sep=''))[[1]])+1
    test_indices <- as.integer(read.table(paste('../saved_indices/top',number_selected,'-tech',datasetnum,'-',fold,'-',seed,'-test_instances_indices', sep=''))[[1]])+1
    selected_indices <- as.integer(read.table(paste('../saved_indices/top',number_selected,'-tech',datasetnum,'-',fold,'-',seed,'-selected_feat_indices', sep=''))[[1]])+1
    priv_indices <- as.integer(read.table(paste('../saved_indices/top',number_selected,'-tech',datasetnum,'-',fold,'-',seed,'-unselected_feat_indices', sep=''))[[1]])+1
    Xtrain <- dataset[train_indices,selected_indices]
    Xtest <- dataset[test_indices,selected_indices]
  
    labels <-(as.double(unlist(read.table(paste('../saved_datasets/tech',datasetnum,'labels', sep='')))))
    Ytrain <- labels[train_indices]
    Ytest <- labels[test_indices]

	
		#zero mean unit variance normalization
		meanTrain <- apply(Xtrain, 2, mean)
		sdTrain <- apply(Xtrain, 2, sd)
		sdTrain[ sdTrain == 0 ] <- 1
		Xtrain <- (Xtrain - matrix(meanTrain, nrow(Xtrain), ncol(Xtrain), byrow = TRUE)) / 
			matrix(sdTrain, nrow(Xtrain), ncol(Xtrain), byrow = TRUE)
		Xtest <- (Xtest - matrix(meanTrain, nrow(Xtest), ncol(Xtest), byrow = TRUE)) / 
			matrix(sdTrain, nrow(Xtest), ncol(Xtest), byrow = TRUE)

		#GPC baseline
		time <- system.time(
		ret <- epGPCExternal(Xtrain, Ytrain, npi, l = 5e-2, sigma0 = 1, optimize_flags = c(TRUE, TRUE, TRUE, TRUE))
		)
		
		
	
		errorTest <- mean(sign(predictGPC(ret, Xtest) - 0.5) != Ytest)
		
		print(errorTest)
		write.table(errorTest, file = paste("../results/GPC_selected/error" ,npi,number_selected,'tech',datasetnum,fold,seed,".txt",sep = "-"), row.names = F, col.names = F, append = FALSE)
		
		# write.table(errorTest, file = paste("./results/GPC/",i,"_errorTest_X_", npi, ".txt", sep = ""),
		# 	row.names = F, col.names = F, append = FALSE)
		# write.table(t(time), file = paste("./results/GPC/",i,"_time_X_", npi, ".txt", sep = ""),
		# 	row.names = F, col.names = F, append = FALSE)

}

