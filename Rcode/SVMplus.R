#SVM+ using QP
library(quadprog)
library(SparseM)
library(e1071)

#
# This functions computes the kernel between X 
#

computeKmm <- function(X, gamma) {

	m <- nrow(X)
	Q <- matrix(apply(X^2, 1, sum), m, m)
	distance <- Q + t(Q) - 2 * X %*% t(X)
	exp(- gamma * distance)
}

computeKnm <- function(Xnew, X, gamma) {

	n <- nrow(Xnew)
	m <- nrow(X)
	Q <- matrix(apply(Xnew^2, 1, sum), n, m)
	Qbar <- matrix(apply(X^2, 1, sum), n, m, byrow = T)
	distance <- Qbar + Q - 2 * Xnew %*% t(X)
	exp(-gamma * distance)
}

#
# This function computes the solution to the SVM+ problem
#

solve_svm_plus_QP <- function(X, Xstar, Y, C, Cstar, gamma, gammaStar) {
	
	# We prepare the data to call the quadratic solver

	K <- computeKmm(X, gamma)	
	Kstar <- computeKmm(Xstar, gammaStar)
	n <- nrow(X)

	Dmat <- matrix(0, 2 * n, 2 * n)

	Dmat[ 1 : n, 1 : n ] <- K * (Y %*% t(Y))
	Dmat[ (1 : n) + n, (1 : n) + n ] <- Kstar / Cstar

	diag(Dmat) <- diag(Dmat) + 1e-10

	dvec <- c(rep(1, n), rep(0, n))

	b0 <- 0
	b0 <- c(b0, 0)
	b0 <- c(b0, rep(0, n))
	b0 <- c(b0, rep(-C, n))

	meq <- 2 

	A <- matrix(0, 2 + 2 * n, 2 * n)
	
	A[ 1, ] <- c(rep(0, n), rep(1, n))
	A[ 2, ] <- c(Y, rep(0, n))
	A[ 3 : (n + 3 - 1), ] <- cbind(diag(n), matrix(0, n, n))
	A[ (3 + n) : (3 + 2 * n - 1), ] <- cbind(-diag(n), diag(n))

	fail <- FALSE
	result <- NULL
	
	tryCatch(result <- solve.QP(Dmat, dvec, t(A), b0, meq) , error = function(x) fail <<- TRUE)
	
	if (fail == TRUE) {
		result <- list()
		result$solution <- rep(0, 2 * n)
	}

	result$solution[ is.na(result$solution) ] <- 0

	alphas <- result$solution[ 1 : n ]

	betas <- result$solution[ (1 : n) + n ] - alphas + C
	upsilon <- result$solution[ (1 : n) + n ]

	alphas[ abs(alphas) < 1e-5 ] <- 0
	betas[ abs(betas) < 1e-5 ] <- 0

	# We compute the bias as explained in Fast Optimization Algorithms for Solving SVM+ by D.Pechyony and V.Vapnik

	Fi <- c((Y * alphas) %*% K)
	fi <- c(upsilon %*% Kstar)

	sel_pos <- which(alphas > 0 & Y == 1)
	sel_neg <- which(alphas > 0 & Y == -1)

	if (length(sel_pos) > 0) {
		s_pos <- mean((1 - fi / Cstar - Fi)[ sel_pos ])
	} else {
		s_pos <- 0
	}

	if (length(sel_neg) > 0) {
		s_neg <- mean((-1 + fi / Cstar - Fi)[ sel_neg ])
	} else {
		s_neg <- 0
	}

	rho <- - (s_pos + s_neg) / 2

	list(X = X, Y = Y, alphas = alphas * Y, rho = rho, betas = betas, gamma = gamma)
}


##
# This function computes the predictions for new data instances 
#

predict_svm_plus_QP <- function(model, Xnew) {

	if (length(which(model$alphas != 0)) > 0) {
		Knm <- computeKnm(Xnew, model$X[ which(model$alphas != 0) ,, drop = FALSE ], model$gamma)
		Fi <- Knm %*% as.vector(model$alphas[ model$alphas != 0 ]) - model$rho
	} else {
		Fi <- -rep(model$rho, nrow(Xnew))
	}

	list(classes = sign(Fi), decision = Fi)
}



# The performance of SVM+ trained with QP
svm_plus_performance_QP <- function(X, Xstar, Y, C, Cstar, gamma, gammaStar, Xtest, Ytest) {

	errors <- array(0, dim = c(length(gamma), length(C), length(Cstar), length(gammaStar)))

	# We call the method for training 

	for (l in 1 : length(gamma))
		for (k in 1 : length(C))
			for (i in 1 : length(Cstar))
				for (j in 1 : length(gammaStar)) {

				        model <- solve_svm_plus_QP(X, Xstar, Y, C[ k ], Cstar[ i ], gamma[ l ], gammaStar[ j ])
					classes <- predict_svm_plus_QP(model, Xtest)$classes

					errors[ l, k, i, j ] <- mean(classes != Ytest)

					cat(".")
	}
	cat("\n")

	errors	
}


