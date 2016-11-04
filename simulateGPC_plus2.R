source("./Rcode/gFITC.R")
source("./Rcode/epAGPC_het.R")
set.seed(0)

# i <- commandArgs()[1]
i <- 1

#number of pseudo-inputs is set to 100
for (npi in c(100)) {

		print(i)
		#Loading the data
		load(paste("./data/",i,"data.dat",sep = ""))
		Xtrain <- data$x[ data$itrain, ]
		Ytrain <- as.vector(c(data$y_train))
		Xtest <- data$x[ data$itest, ]
		Ytest <- as.vector(c(data$y_test))
		Xstar_train <- data$x_star[ data$itrain, ]
		Xstar_train <- matrix(Xstar_train, length(c(Xstar_train)), 1)

		#zero mean unit variance normalization
		meanTrain <- apply(Xtrain, 2, mean)
		sdTrain <- apply(Xtrain, 2, sd)
		sdTrain[ sdTrain == 0 ] <- 1
		Xtrain <- (Xtrain - matrix(meanTrain, nrow(Xtrain), ncol(Xtrain), byrow = TRUE)) / 
			matrix(sdTrain, nrow(Xtrain), ncol(Xtrain), byrow = TRUE)
		Xtest <- (Xtest - matrix(meanTrain, nrow(Xtest), ncol(Xtest), byrow = TRUE)) / 
			matrix(sdTrain, nrow(Xtest), ncol(Xtest), byrow = TRUE)

		#GPC+ baseline
		time <- system.time(
		ret <- epGPCExternal(Xtrain, Xstar_train, Ytrain, npi, sigmaF = 1, sigma0F = 1, 
			lF = 5e-2, sigmaG = 1, sigma0G = 1, lG = 1, optimize_flags_F = c(TRUE, TRUE, TRUE, TRUE),
			optimize_flags_G = c(TRUE, TRUE, TRUE, TRUE, TRUE))
		)

		errorTest <- mean(sign(predictGPC(ret, Xtest) - 0.5) != Ytest)	
		write.table(errorTest, file = paste("./results/GPC_plus/",i,"_errorTest_X_", npi, "NEW.txt", sep = ""),
			row.names = F, col.names = F, append = FALSE)
		write.table(t(time), file = paste("./results/GPC_plus/",i,"_time_X_", npi, "NEW.txt", sep = ""), 
			row.names = F, col.names = F, append = FALSE)
}
