library(e1071)
set.seed(0)

# i <- commandArgs()[1]
i <- 1

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

	#SVM baseline
	time <- system.time( ret <- tuneParametersVal(Xtrain, Ytrain))
	print(ret)
	resultRet <- ret
	ret <- ret$model

	errorTest <- mean(predict(ret, Xtest) != Ytest)

	write.table(errorTest, file = paste("./results/SVM/",i,"_errorTest_X_val.txt",sep = ""), row.names = F, col.names = F, append = FALSE)
	write.table(t(time), file = paste("./results/SVM/",i,"_time_X_val.txt",sep = ""), row.names = F, col.names = F, append = FALSE)



