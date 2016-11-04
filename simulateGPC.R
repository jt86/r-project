source("./Rcode/gFITC.R")
source("./Rcode/epGPC.R")
set.seed(0)

i <- commandArgs()[1]

#number of pseudo-inputs is set to 100
for (npi in c(100)) {
	
		print(i)

		#Loading the data
		load(paste("./data/",i,"data.dat",sep = ""))
		Xtrain <- data$x[ data$itrain, ]
		Ytrain <- as.vector(c(data$y_train))
		Xtest <- data$x[ data$itest, ]
		Ytest <- as.vector(c(data$y_test))
	
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

		write.table(errorTest, file = paste("./results/GPC/",i,"_errorTest_X_", npi, ".txt", sep = ""),
			row.names = F, col.names = F, append = FALSE)
		write.table(t(time), file = paste("./results/GPC/",i,"_time_X_", npi, ".txt", sep = ""),
			row.names = F, col.names = F, append = FALSE)

}

