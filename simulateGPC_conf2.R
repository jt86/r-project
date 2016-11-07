source("./Rcode/gFITC.R")
source("./Rcode/epAGPC.R")
set.seed(0)

i <- 1
# i <- commandArgs()[1]


number_selected <- commandArgs()[7]
datasetnum <- commandArgs()[8]
seed <- commandArgs()[9]
fold <- commandArgs()[10]
print(number_selected)
print(datasetnum)
print(fold)
print(seed)

# We plot the values of the G function (Only works for uni-dimensional privileged data)
plot_function_g <- function(ret, Xstar_train) {
	grid_size <- 100

	par(oma = c(0,0,0,3))
	grid  <- matrix(seq(min(Xstar_train), max(Xstar_train), length = grid_size), grid_size, 1) 
	result <- predict_g_function(ret, grid)
	plot(grid, result$mean, type = "l", col = "blue", lwd = 5, 
		ylim = c(min(result$mean - sqrt(result$var)), max(result$mean + sqrt(result$var))), 
		xlab = "Confidence introduced in Xconf", ylab = "Value of the g(·) function")
	for (k in (1 : (length(grid) - 1))) {
		x <- c(grid[ k ], grid[ k ], grid[ k + 1 ], grid[ k + 1 ])
		y <- c(result$mean[ k ], result$mean[ k ] + sqrt(result$var[ k ]), result$mean[ k + 1 ] + 
			sqrt(result$var[ k + 1 ]), result$mean[ k + 1 ])
		polygon(x, y, col = rgb(0.8, 0.8, 1), border = rgb(0.8, 0.8, 1))
		x <- c(grid[ k ], grid[ k ], grid[ k + 1 ], grid[ k + 1 ])
		y <- c(result$mean[ k ], result$mean[ k ] - sqrt(result$var[ k ]), result$mean[ k + 1 ] - 
			sqrt(result$var[ k + 1 ]), result$mean[ k + 1 ])
		polygon(x, y, col = rgb(0.8, 0.8, 1), border = rgb(0.8, 0.8, 1))
	}

	lines(grid, result$mean + sqrt(result$var), col = "blue", lwd = 1, lty = 1)
	lines(grid, result$mean - sqrt(result$var), col = "blue", lwd = 1, lty = 1)
	prob_positive <- pnorm(result$mean / sqrt(result$var))
	trans_probs_for_plotting <- (prob_positive) * (max(result$mean + sqrt(result$var)) - 
		min(result$mean - sqrt(result$var))) + min(result$mean - sqrt(result$var))
	lines(grid, trans_probs_for_plotting, col = "red", lwd = 5)
	axis(4, seq(0, 1, length = 11) * (max(result$mean + sqrt(result$var)) - 
		min(result$mean - sqrt(result$var))) + min(result$mean - sqrt(result$var)), labels = seq(0, 1, length = 11))
	legend("topright", c("Posterior mean of g(·)", "Probability of g(·)>0"), col = c("blue", "red"), lwd = 5)
	lines(grid, result$mean, col = "blue", lwd = 5)
	mtext("Probability", side = 4, line = 3)

}

plot_training_confidence_vs_privileged_values <- function(conf_train, Xstar_train) {

	conf_train <- pmax(conf_train, 1 - conf_train)
	plot(Xstar_train, conf_train, col = "blue", pch = "x", lwd = 3, ylab = "Confidence in the Label by the Model",
		xlab = "Confidence Introduced to the Model in Xconf")

}

#number of pseudo-inputs is set to 100
for (npi in c(100)) {

		print(i)
		#Loading the data
		load(paste("./data/",i,"data.dat",sep = ""))
		# Xtrain <- data$x[ data$itrain, ]
		# Ytrain <- as.vector(c(data$y_train))
		# Xtest <- data$x[ data$itest, ]
		# Ytest <- as.vector(c(data$y_test))
		# Xstar_train <- data$x_star[ data$itrain, ]
		# Xstar_train <- matrix(Xstar_train, length(c(Xstar_train)), 1)
    
    #Loading the data
    dataset <-(read.table(paste("../saved_datasets/tech",datasetnum,"data",sep = "")))
    dims=dim(dataset)
    dataset <-as.numeric(unlist(dataset))
    print('dim of dataset')
    print (length(dataset))
    typeof(dataset)
    dataset<-matrix(dataset,dims[1],dims[2])
    typeof(dataset)
    dim(dataset)
    # /home/j/jt/jt306/Documents/CVPR2016_Rcode/
    train_indices <- as.integer(read.table(paste('../saved_indices/top',number_selected,'-tech',datasetnum,'-',fold,'-',seed,'-train_instances_indices', sep=''))[[1]])+1
    test_indices <- as.integer(read.table(paste('../saved_indices/top',number_selected,'-tech',datasetnum,'-',fold,'-',seed,'-test_instances_indices', sep=''))[[1]])+1
    selected_indices <- as.integer(read.table(paste('../saved_indices/top',number_selected,'-tech',datasetnum,'-',fold,'-',seed,'-selected_feat_indices', sep=''))[[1]])+1
    priv_indices <- as.integer(read.table(paste('../saved_indices/top',number_selected,'-tech',datasetnum,'-',fold,'-',seed,'-unselected_feat_indices', sep=''))[[1]])+1
    
    Xtrain <- dataset[train_indices,selected_indices]
    Xtest <- dataset[test_indices,selected_indices]
    
    Xstar_train <- dataset[train_indices,priv_indices]
    print('dim of train priv')
    print(dim(Xstar_train))

    labels <-(as.double(unlist(read.table(paste('../saved_datasets/tech',datasetnum,'labels', sep='')))))
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

		#GPC_conf baseline
		time <- system.time(
		ret <- epGPCExternal(Xtrain, Xstar_train, Ytrain, npi, sigmaF = 1, sigma0F = 1, 
			lF = 5e-2, sigmaG = 1, sigma0G = 1, lG = 1, optimize_flags_F = c(TRUE, TRUE, TRUE, TRUE),
			optimize_flags_G = c(TRUE, TRUE, TRUE, TRUE, TRUE))
		)

		errorTest <- mean(sign(predictGPC(ret, Xtest) - 0.5) != Ytest)		

		# write.table(errorTest, file = paste("./results/GPC_conf/",i,"_errorTest_X_", npi, ".txt", sep = ""),
		# 	row.names = F, col.names = F, append = FALSE)
		# 
		# write.table(t(time), file = paste("./results/GPC_conf/",i,"_time_X_", npi, ".txt", sep = ""), 
		# 	row.names = F, col.names = F, append = FALSE)


		write.table(errorTest, file = paste("../results/GPC_conf/error" ,npi,number_selected,'tech',datasetnum,fold,seed,".txt",sep = "-"), row.names = F, col.names = F, append = FALSE)
		
		write.table(t(time), file = paste("../results/SVM/time",npi,number_selected,'tech',datasetnum,fold,seed,".txt",sep = "-"), row.names = F, col.names = F, append = FALSE)
	
		# We plot the values of the G function (Only works for uni-dimensional privileged data)
		#pdf(paste("./results/GPC_conf/",i,"_plot_g_function_", npi, ".pdf", sep = ""), width = 9, height = 6)
		#plot_function_g(ret, Xstar_train)
		#dev.off()
		# We plot the Xconf values vs the confidence in the class label by the model for training
                #confX <- predictGPC_privileged_train(ret, Xtrain, Xstar_train) 	
		#pdf(paste("./results/GPC_conf/",i,"_plot_conf_train_vs_Xconf_", npi, ".pdf", sep = ""), width = 7, height = 5)
		#plot_training_confidence_vs_privileged_values(confX, Xstar_train)
		#dev.off()
}
