##
# Function which implements the EP algorithm for the GP classifiar
# based on the FITC approximation.
#
# @param	X	n x d matrix with the data points.
# @param	Xbar	m x d matrix with the pseudo inputs.
# @param	sigma	scalar with the log-amplitude of the GP.
# @param	sigma0	scalar with the log-noise level in the GP.
# @param	l	d-dimensional vector with the log-lengthscales.
# @param	Y	n-dimensional vector with the class labels.
# @param	m0	scalar with the mean of the prior Gaussian Process
# @param 	start	Initial values for the approximate factors.
#
# @return	ret	A list with the following elements:
#
#			a		The posterior approximation.
#			f1Hat		The approximation to the first factor.
#			f2Har		The approximation to the second factor.
#			logZ		The approximation to the log evidence.
#			gradientLogZ	The gradient of the log evidence.
#
# @author Jose Miguel Hernandez-Lobato
#

epGPCInternal <- function(X, Xbar, sigma, sigma0, l, Y, m0 = 0, start = NULL) {

	# We initialize the structure with the problem information

	gFITCinfo <- initGFITCinfo(X, Xbar, sigma, sigma0, l, m0)

	# We initialize the approximate factor to be uniform

	f1Hat <- list(eta1 = rep(0, gFITCinfo$n), eta2 = rep(0, gFITCinfo$n))

	# We check for an initial solution
	
	if (!is.null(start)) 
		f1Hat <- start$f1Hat

	# Main loop of EP

	i <- 1
	damping <- .5
	convergence <- FALSE
	while (!convergence && i < 1e3) {

		f1HatOld <- f1Hat

		# We get the marginals of the posterior approximation

		ret <- computeTitledDistribution(gFITCinfo, f1Hat$eta2, f1Hat$eta1)

		meanMarginals <- ret$mNew
		varMarginals <- ret$vNew

		# We refine the first approximate factor in parallel

		vOld <- (varMarginals^-1 - f1Hat$eta2)^-1
		mOld <- vOld * (meanMarginals  / varMarginals - f1Hat$eta1)

		logZ <- pnorm(Y * mOld / sqrt(vOld), log.p = T)
		ratio <- exp(-logZ + dnorm(mOld / sqrt(vOld), log = T))
		alpha <- ratio * Y / sqrt(vOld)
		beta <- -ratio * (Y * mOld / sqrt(vOld) + ratio) / vOld

		eta2HatNew <- -beta / (1 + beta * vOld)
		eta1HatNew <- (alpha - mOld * beta) / (1 + beta * vOld)

#		eta1HatNew[ eta2HatNew < 0 ] <- f1Hat$eta1[ eta2HatNew < 0 ]
#		eta2HatNew[ eta2HatNew < 0 ] <- f1Hat$eta2[ eta2HatNew < 0 ]
#		eta2HatNew[ eta2HatNew < 0 ] <- 1e-2

		f1Hat$eta2 <- damping * eta2HatNew + (1 - damping) * f1Hat$eta2
		f1Hat$eta1 <- damping * eta1HatNew + (1 - damping) * f1Hat$eta1

		# We check for convergence

		change <- max(abs(f1HatOld$eta1 - f1Hat$eta1))
		change <- max(change, abs(f1HatOld$eta2 - f1Hat$eta2))

		if (change < 1e-4)
			convergence <- T

#		cat("Iteration",  i, change, "\n")

		# Annealed damping scheme

		damping <- damping * 0.99

		i <- i + 1
	}

	# We compute the evidence and its gradient

	logZ <- computeEvidence(f1Hat, gFITCinfo, Y)
	gradientLogZ <- computeDerivativesEvidence(gFITCinfo, f1Hat$eta2, f1Hat$eta1)

	# We are done!

	list(f1Hat = f1Hat, logZ = logZ, gradientLogZ = gradientLogZ, gFITCinfo = gFITCinfo)
}

###
# Function which computes the EP approximation of the log evidence.
#
# @param	f1Hat		The approximation for the first factor.
# @param	gFITCinfo	The list with the problem information.
# @param	Y		The class labels.
#
# @return	logZ		The log evidence.
#

computeEvidence <- function(f1Hat, gFITCinfo, Y) {

	ret <- computeTitledDistribution(gFITCinfo, f1Hat$eta2, f1Hat$eta1)

	meanMarginals <- ret$mNew
	varMarginals <- ret$vNew

	vOld <- (varMarginals^-1 - f1Hat$eta2)^-1
	mOld <- vOld * (meanMarginals  / varMarginals - f1Hat$eta1)

#	logZ <- sum(pnorm(Y * mOld / sqrt(vOld + 1), log.p = T)) -
#		sum(dnorm(0, mean = (mOld - f1Hat$eta1 * f1Hat$eta2^-1), sd = 
#		sqrt(vOld + f1Hat$eta2^-1), log = TRUE))

	logZ <- pnorm(Y * mOld / sqrt(vOld), log.p = T)

	logZ <- sum(logZ) - sum(0.5 * log(varMarginals) - 0.5 * log(2 * pi) - 0.5 * log(vOld) - 
		0.5 * mOld^2 / vOld - 0.5 * f1Hat$eta1^2 / f1Hat$eta2 + 0.5 * meanMarginals^2 / varMarginals) #- sum(0.5 * log(f1HatF$eta2))

	logZret <- logZ + getFITCevidence(gFITCinfo, f1Hat$eta2, f1Hat$eta1)
	
	logZret
}

###
# Function which computes the probability of class 1 on new data.
#
# @param	ret	The list returned by epGPCExternal
# @param	Xtest	The n x d matrix with the new data points.
#
# @return	pOne	The probability of class 1 on the new data.
#

predictGPC <- function(ret, Xtest) {

	retInternal <- predictFITC(ret$gFITCinfo, ret$f1Hat$eta2, ret$f1Hat$eta1, Xtest)
	pnorm(retInternal$m / sqrt(retInternal$v))
}

##
# Function which adjust the kernel parameters by gradient descent.
#
# @param	X	n x d matrix with the data points.
# @param	Y	n-dimensional vector with the class labels.
# @param	m	The number of pseudo-inputs to use.
#
# @return	ret	A list with the following elements:
#
#			a		The posterior approximation.
#			f1Hat		The approximation to the first factor.
#			f2Har		The approximation to the second factor.
#			logZ		The approximation to the log evidence.
#			gradientLogZ	The gradient of the log evidence.
#
#			optimize_flags: booleans corresponding to sigma, sigma0, l and pesudoinputs
#

epGPCExternal <- function(X, Y, m, m0 = 0, sigma = 1, sigma0 = 1, l = 1, optimize_flags = rep(TRUE, 4)) {

	# We initialize the hyper-parameters

	m0 <- m0
	sigma <- log(sigma)
	sigma0 <- log(sigma0)
	l <- rep(log(l), ncol(X)) #-log((apply(X, 2, max) - apply(X, 2, min)) / 2)
#	Xbar <- X[ sample(1 : nrow(X), m), ]
	Xbar <- X[ 1 : m, ]

	# We initialize the gradient optimization process

	ret <- epGPCInternal(X, Xbar, sigma, sigma0, l, Y, m0)
	best <- ret

	cat(0, "New evidence:", ret$logZ, "\n")

	eps <- 0.01
	convergence <- F
	iteration <- 1
	
	while (!convergence && iteration < 100) {

		if (optimize_flags[ 1 ])
			sigma <- sigma + eps * ret$gradientLogZ$dLogZdSigma
		if (optimize_flags[ 2 ])
			sigma0 <- sigma0 + eps * ret$gradientLogZ$dLogZdSigma0
		if (optimize_flags[ 3 ])
			l <- l + eps * sum(ret$gradientLogZ$dLogZdl)
		if (optimize_flags[ 4 ])
			Xbar <- Xbar + eps * ret$gradientLogZ$dLogZdXbar

#		m0 <- m0 + eps * ret$gradientLogZ$dLogZdm0

		# We train the model using the previous solution as the starting point. If that fails (the starting point
		# is unfeasible we start from scracht)


		tryCatch( retNew <<- epGPCInternal(X, Xbar, sigma, sigma0, l, Y, m0, ret)
                , error = function(x)  {
			tryCatch( retNew <<- epGPCInternal(X, Xbar, sigma, sigma0, l, Y, m0)
                	, error = function(x)  convergence <<- TRUE)
			}
                )


		if (is.nan(retNew$logZ))
			return(best)

		if (abs(retNew$logZ - ret$logZ) < 1e-4)
			convergence <- T

		if (retNew$logZ < ret$logZ)
			eps <- eps * 0.5
		else
			eps <- eps * 1.1

		cat(iteration, "New evidence:", retNew$logZ, "eps:", eps, "Change:", abs(retNew$logZ - ret$logZ), l[ 1 ], "\n")

		ret <- retNew

		if (ret$logZ > best$logZ)
			best <- ret

		iteration <- iteration + 1
	}

	best
}


