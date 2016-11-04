library(e1071)
source("Rcode/SVMplus.R")
pathToExec <- "./Rcode"
set.seed(0)

i <- commandArgs()[1]

tune_svm_plus_validation_QP <- function(X, Xstar, Y, 
	C =  exp(seq(-3, 7, length = 5)), Cstar = exp(seq(-3, 7, length = 5)), 
	gamma = exp(seq(-3, 3, length = 5)), gammaStar = exp(seq(-3, 3, length = 5)), percentage_of_data_for_validation = 0.1) {

	mat <- as.matrix(dist(X))^2
	med <- median(mat[ upper.tri(mat) ])
	gamma <- gamma / med

        i_pos <- which(Y == 1)
        i_neg <- which(Y == -1)
        half <- round(percentage_of_data_for_validation * length(Y))

	errors <- array(0, dim = c(length(gamma), length(C), length(Cstar), length(gammaStar)))
	for (fold in 1 : 5) {
        	sample_pos = sample(i_pos, half)
        	sample_neg = sample(i_neg, half)
        	s = sample(c(sample_pos, sample_neg))
        	Xval <- X[ s, ]
        	Yval <- Y[ s ]

        	s_remain <- (1 : nrow(X))[ -c(s) ]
        	Xtrain <- X[ s_remain, ]
        	Ytrain <- Y[ s_remain ]
        	Xstar_train <- matrix(Xstar[s_remain,])

        	errors <- errors + svm_plus_performance_QP(Xtrain, Xstar_train, Ytrain, C, Cstar, gamma, gammaStar, Xval, Yval)
	}

	errors <- errors/5.0
	errors

        bestC <- NULL
        bestCstar <- NULL
        bestgammaStar <- NULL
        bestValue <- Inf

        for (l in 1 : length(gamma))
                for (k in 1 : length(C))
                        for (i in 1 : length(Cstar))
                                for (j in 1 : length(gammaStar))
                                        if (errors[ l, k, i, j ] < bestValue) {
                                                bestCstar <- Cstar[ i ]
                                                bestC <- C[ k ]
                                                bestgamma <- gamma[ l ]
                                                bestgammaStar <- gammaStar[ j ]
                                                bestValue <- errors[ l, k, i, j ]
                                        }

        model <- solve_svm_plus_QP(X, Xstar, Y, bestC, bestCstar, bestgamma, bestgammaStar)

        list(best.model = model, error = bestValue, bestgammaStar = bestgammaStar, bestCstar = bestCstar, bestC = bestC, bestgamma = bestgamma)
}


	print(i)

	#Loading the data
	load(paste("./data/",i,"data.dat",sep = ""))
	Xstar_train <- data$x_star[ data$itrain, ]
	Xstar_train <- matrix(Xstar_train, length(c(Xstar_train)), 1)
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

	#SVM+ baseline
	time <- system.time(result <- tune_svm_plus_validation_QP(Xtrain, Xstar_train, Ytrain))
	print(result$error)
	model <- result$best.model

	errorTest <- mean(predict_svm_plus_QP(model, Xtest)$classes != Ytest)
	write.table(errorTest, file = paste("./results/SVM_plus/",i,"_errorTest_X_QP_val.txt",sep = ""), row.names = F, col.names = F, append = FALSE)
	write.table(t(time), file = paste("./results/SVM_plus/",i,"_time_X_QP_val.txt",sep = ""), row.names = F, col.names = F, append = FALSE)
  
	errorTest
	

	