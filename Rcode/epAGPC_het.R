##
# Function which implements the EP algorithm for the GP classifier
# based on the FITC approximation. Two gaussian process are considered.
# A first GP for the decision boundaries and a second for the noise.
#
# @param	X	n x d matrix with the data points.
# @param	Xstar	n x d matrix with the data points, priviledged information.
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
# @author Daniel Hernandez-Lobato
#

epGPCInternal <- function(X, Xstar, XbarF, XbarG, sigmaF, sigmaG, sigma0F, sigma0G, lF, lG, Y, m0F = 0, 
	m0G = 0, start = NULL) {

	# We initialize the structure with the problem information

	gFITCinfoF <- initGFITCinfo(X, XbarF, sigmaF, sigma0F, lF, m0F)
	gFITCinfoG <- initGFITCinfo(Xstar, XbarG, sigmaG, sigma0G, lG, m0G)

	# We initialize the approximate factor to be uniform

	f1HatF <- list(eta1 = rep(0, gFITCinfoF$n), eta2 = rep(0, gFITCinfoF$n))
	f1HatG <- list(eta1 = rep(0, gFITCinfoG$n), eta2 = rep(0, gFITCinfoG$n))

	# We check for an initial solution
	
	if (!is.null(start))  {
		f1HatF <- start$f1HatF
		f1HatG <- start$f1HatG
	}

	a <- list(gFITCinfoG = gFITCinfoG, gFITCinfoF = gFITCinfoF, f1HatG = f1HatG, f1HatF = f1HatF)
	aOld <- a

	# Main loop of EP

	i <- 1
	damping <- .5
	convergence <- FALSE
	while (!convergence && i < 1e3) {

		update_correct <- FALSE
		damping_inner <- damping

		while(update_correct != TRUE) {

			error <- FALSE

			tryCatch(aNew <<- process_likelihood_factors(a, Y, damping_inner), error = function(x) error <<- TRUE)

			if (error == FALSE) {
				update_correct <- TRUE
			} else {

				if (i == 1)
					stop("Error in first iteration!")

				a <- aOld
				damping_inner <- damping_inner * 0.5
			}
		}

		aOld <- a
		a <- aNew

		# We get the marginals of the posterior approximation
		# We check for convergence

		change <- max(abs(aOld$f1HatF$eta1 - a$f1HatF$eta1))
		change <- max(change, abs(aOld$f1HatF$eta2 - a$f1HatF$eta2))
		change <- max(change, abs(aOld$f1HatG$eta1 - a$f1HatG$eta1))
		change <- max(change, abs(aOld$f1HatG$eta2 - a$f1HatG$eta2))

		if (change < 1e-3)
			convergence <- T

		cat("\tIteration",  i, change, "\n")

		# Annealed damping scheme

		damping <- damping * 0.99

		i <- i + 1
	}

	# We compute the evidence and its gradient

	logZ <- computeEvidence(a$f1HatF, a$f1HatG, a$gFITCinfoF, a$gFITCinfoG, Y)
	gradientLogZF <- computeDerivativesEvidence(a$gFITCinfoF, a$f1HatF$eta2, a$f1HatF$eta1)
	gradientLogZG <- computeDerivativesEvidence(a$gFITCinfoG, a$f1HatG$eta2, a$f1HatG$eta1)

	# We are done!

	list(f1HatF = a$f1HatF, f1HatG = a$f1HatG, logZ = logZ, gradientLogZF = gradientLogZF, gradientLogZG = gradientLogZG, 
		gFITCinfoF = a$gFITCinfoF, gFITCinfoG = a$gFITCinfoG)
}

###
# Function which processes the likelihood
#
# @param	a 		The approximation 
# @param	Y		The class labels.
# @param	damping		The damping factor 
#
# @return	a 		The updated approximation
#
#

process_likelihood_factors <- function(a, Y, damping) {

	# We get the marginals of the posterior approximation

	retF <- computeTitledDistribution(a$gFITCinfoF, a$f1HatF$eta2, a$f1HatF$eta1)
	retG <- computeTitledDistribution(a$gFITCinfoG, a$f1HatG$eta2, a$f1HatG$eta1)

	meanMarginalsF <- retF$mNew
	varMarginalsF <- retF$vNew

	meanMarginalsG <- retG$mNew
	varMarginalsG <- retG$vNew

	# We compute an old distribution

	vOldF <- (varMarginalsF^-1 - a$f1HatF$eta2)^-1
	mOldF <- vOldF * (meanMarginalsF / varMarginalsF - a$f1HatF$eta1)

	vOldG <- (varMarginalsG^-1 - a$f1HatG$eta2)^-1
	mOldG <- vOldG * (meanMarginalsG / varMarginalsG - a$f1HatG$eta1)

	if (any(vOldF < 0)  || any(varMarginalsF < 0) || any(vOldG < 0) || any(varMarginalsG < 0)) {
		stop("Negative variances!")
	}

	# We refine the first approximate factor in parallel

	Z <- computeZ(Y, mOldF, vOldF, mOldG, vOldG)
	logZ <- log(Z)

	alphaF <- computedZdmOldF(Y, mOldF, vOldF, mOldG, vOldG) / Z
	betaF <- computed2ZdmOldF2(Y, mOldF, vOldF, mOldG, vOldG) / Z - alphaF^2

	eta2HatNewF <- -betaF / (1 + betaF * vOldF)
	eta1HatNewF <- (alphaF - mOldF * betaF) / (1 + betaF * vOldF)

	a$f1HatF$eta2 <- damping * eta2HatNewF + (1 - damping) * a$f1HatF$eta2
	a$f1HatF$eta1 <- damping * eta1HatNewF + (1 - damping) * a$f1HatF$eta1

	# We refine the second approximate factor in parallel

	alphaG <- computedZdmOldG(Y, mOldF, vOldF, mOldG, vOldG) / Z
	betaG <- computed2ZdmOldG2(Y, mOldF, vOldF, mOldG, vOldG) / Z - alphaG^2

	eta2HatNewG <- -betaG / (1 + betaG * vOldG)
	eta1HatNewG <- (alphaG - mOldG * betaG) / (1 + betaG * vOldG)

	a$f1HatG$eta2 <- damping * eta2HatNewG + (1 - damping) * a$f1HatG$eta2
	a$f1HatG$eta1 <- damping * eta1HatNewG + (1 - damping) * a$f1HatG$eta1

	a
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
#

computeEvidence <- function(f1HatF, f1HatG, gFITCinfoF, gFITCinfoG, Y) {

	retF <- computeTitledDistribution(gFITCinfoF, f1HatF$eta2, f1HatF$eta1)
	retG <- computeTitledDistribution(gFITCinfoG, f1HatG$eta2, f1HatG$eta1)

	meanMarginalsF <- retF$mNew
	varMarginalsF <- retF$vNew

	meanMarginalsG <- retG$mNew
	varMarginalsG <- retG$vNew

	vOldF <- (varMarginalsF^-1 - f1HatF$eta2)^-1
	mOldF <- vOldF * (meanMarginalsF  / varMarginalsF - f1HatF$eta1)

	vOldG <- (varMarginalsG^-1 - f1HatG$eta2)^-1
	mOldG <- vOldG * (meanMarginalsG  / varMarginalsG - f1HatG$eta1)

	Z <- computeZ(Y, mOldF, vOldF, mOldG, vOldG)

	logZ <- sum(log(Z)) - sum(0.5 * log(varMarginalsF) - 0.5 * log(2 * pi) - 0.5 * log(vOldF) - 
		0.5 * mOldF^2 / vOldF - 0.5 * f1HatF$eta1^2 / f1HatF$eta2 + 0.5 * meanMarginalsF^2 / varMarginalsF) #- sum(0.5 * log(f1HatF$eta2))

	logZ <- logZ - sum(0.5 * log(varMarginalsG) - 0.5 * log(2 * pi) - 0.5 * log(vOldG) - 
		0.5 * mOldG^2 / vOldG - 0.5 * f1HatG$eta1^2 / f1HatG$eta2 + 0.5 * meanMarginalsG^2 / varMarginalsG) #- sum(0.5 * log(f1HatG$eta2))

	logZret <- logZ + getFITCevidence(gFITCinfoF, f1HatF$eta2, f1HatF$eta1) + getFITCevidence(gFITCinfoG, f1HatG$eta2, f1HatG$eta1)

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
#

#predictGPC <- function(ret, Xtest) {
#
#	retInternalF <- predictFITC(ret$gFITCinfoF, ret$f1HatF$eta2, ret$f1HatF$eta1, Xtest)
#	retInternalG <- predictFITC(ret$gFITCinfoG, ret$f1HatG$eta2, ret$f1HatG$eta1, Xtest)
#
#	pnorm(retInternalF$m / sqrt(retInternalF$v)) * pnorm(-retInternalG$m/ sqrt(retInternalG$v + 1)) + 
#		0.5 * pnorm(retInternalG$m / sqrt(retInternalG$v+ 1))
#}

predictGPC <- function(ret, Xtest) {

	retInternalF <- predictFITC(ret$gFITCinfoF, ret$f1HatF$eta2, ret$f1HatF$eta1, Xtest)
	probNoise <- predictGPCNoise(ret, Xtest)

	pnorm(retInternalF$m / sqrt(retInternalF$v)) 
}

predictGPCNoiseStar <- function(ret, Xtest) {
	retInternalG <- predictFITC(ret$gFITCinfoG, ret$f1HatG$eta2, ret$f1HatG$eta1, Xtest)
	pnorm(retInternalG$m / sqrt(retInternalG$v)) / 2
}

predictGPCNoise <- function(ret, Xtest) {
	retInternalG <- pnorm(ret$gFITCinfoG$m0 / sqrt(ret$gFITCinfoG$sigma + ret$gFITCinfoG$sigma0))
}

predictGPCDecision <- function(ret, Xtest) {
	retInternalF <- predictFITC(ret$gFITCinfoF, ret$f1HatF$eta2, ret$f1HatF$eta1, Xtest)
	pnorm(retInternalF$m / sqrt(retInternalF$v))
}

##
# Function which adjust the kernel parameters by gradient descent.
#
# @param	X	n x d matrix with the data points.
# @param	Xstar	n x d matrix with the data points in the priviledged data.
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
#			optimize_flags_F booleans corresponding to sigmaF, sigma0F, lF ans pesudinputs
#
#

epGPCExternal <- function(X, Xstar, Y, m, sigmaF = 0, sigma0F = 0, lF = 0, 
	sigmaG = 0, sigma0G = 0, lG = 0, optimize_flags_F = rep(TRUE, 4), optimize_flags_G = rep(TRUE, 5)) {

	# We initialize the hyper-parameters

	m0F <- 0
	m0G <- -1
	sigmaF <- log(sigmaF)
	sigmaG <- log(sigmaG)
	sigma0F <- log(sigma0F)
	sigma0G <- log(sigma0G)
	lF <- rep(log(lF), ncol(X)) #mean(-log((apply(X, 2, max) - apply(X, 2, min)) / 2))
	lG <- rep(log(lG), ncol(Xstar)) #mean(-log((apply(Xstar, 2, max) - apply(Xstar, 2, min)) / 2))
#	XbarF <- X[ sample(1 : nrow(X), m), ]
#	XbarG <- Xstar[ sample(1 : nrow(Xstar), m), ]
	XbarF <- matrix(c(X[ 1 : m, ]), m, ncol(X))
	XbarG <- matrix(Xstar[ 1 : m, ], m, ncol(Xstar))

	# We initialize the gradient optimization process

	ret <- epGPCInternal(X, Xstar, XbarF, XbarG, sigmaF, sigmaG, sigma0F, sigma0G, lF, lG, Y, m0F, m0G) 
	best <- ret

	cat(0, "New evidence:", ret$logZ, "\n")

	eps <- 0.01
	convergence <- F
	iteration <- 1

	while (!convergence && iteration < 100) {

		if (optimize_flags_F[ 1 ])
			sigmaF <- sigmaF + eps * ret$gradientLogZF$dLogZdSigma

		if (optimize_flags_F[ 2 ])
			sigma0F <- sigma0F + eps * ret$gradientLogZF$dLogZdSigma0
		
		if (optimize_flags_F[ 3 ])
			lF <- lF + eps * sum(ret$gradientLogZF$dLogZdl)

		if (optimize_flags_F[ 4 ])
			XbarF <- XbarF + eps * ret$gradientLogZF$dLogZdXbar

#		m0F <- m0F + eps * ret$gradientLogZF$dLogZdm0

		if (optimize_flags_G[ 1 ])
			sigmaG <- sigmaG + eps * ret$gradientLogZG$dLogZdSigma

		if (optimize_flags_G[ 2 ])
			sigma0G <- sigma0G + eps * ret$gradientLogZG$dLogZdSigma0

		if (optimize_flags_G[ 3 ])
			lG <- lG + eps * sum(ret$gradientLogZG$dLogZdl)
		
		if (optimize_flags_G[ 4 ])
			XbarG <- XbarG + eps * ret$gradientLogZG$dLogZdXbar

		if (optimize_flags_G[ 5 ])
			m0G <- m0G + eps * ret$gradientLogZG$dLogZdm0

		# We train the model using the previous solution as the starting point. If that fails (the starting point
		# is unfeasible we start from scracht)

		cat("\t", "sF:", sigmaF, "s0F:", sigma0F, "lF:", lF[ 1 ], "sG:", sigmaG, "s0G:", sigma0G, "lG:", lG[ 1 ], "m0G:", m0G, "\n")

		vars <- list(sigmaF = sigmaF, sigma0F = sigma0F, sigmaG = sigmaG, sigma0G = sigma0G, m0G = m0G, m0F = m0F, ret = ret,
				X = X, Xstar = Xstar, XbarF = XbarF, XbarG = XbarG, lF = lF, lG = lG, Y = Y)

		tryCatch(
			retNew <<- epGPCInternal(X, Xstar, XbarF, XbarG, sigmaF, sigmaG, sigma0F, sigma0G, lF, lG, Y, m0F, m0G, ret)
                , error = function(x) 
			tryCatch(
                	retNew <<- epGPCInternal(X, Xstar, XbarF, XbarG, sigmaF, sigmaG, sigma0F, sigma0G, lF, lG, Y, m0F, m0G)  
                	, error = function(x) 
			convergence <<- TRUE
                	)
                )

		if (is.nan(retNew$logZ))
			return(best)

		if (abs(retNew$logZ - ret$logZ) < 1e-4)
			convergence <- T

		if (retNew$logZ < ret$logZ)
			eps <- eps * 0.5
		else
			eps <- eps * 1.05

		cat(iteration, "New evidence:", retNew$logZ, "eps:", eps, "Change:", retNew$logZ - ret$logZ, "\n")

		ret <- retNew

		if (ret$logZ > best$logZ)
			best <- ret

		iteration <- iteration + 1
	}

	best
}

##
# We compute the required quantities for a heteroscedastic GP
##

# This computes the normalization value Z of the cavity times the exact factor

computeZ <- function(Y, mOldF, vOldF, mOldG, vOldG) {

	# We marginalize f

	marginal_value <- function(Y, mOldF, vOldF, g) {
		pnorm(Y * mOldF / sqrt(exp(g) + vOldF))
	}

	m <- cbind(Y, mOldF, vOldF, mOldG, vOldG)

	as.vector(apply(m, 1, function(x) {
		Y <- x[ 1 ]
		mOldF <- x[ 2 ]
		vOldF <- x[ 3 ]
		mOldG <- x[ 4 ]
		vOldG <- x[ 5 ]
		integrate(function(g) dnorm(g, mOldG, sqrt(vOldG)) * marginal_value(Y, mOldF, vOldF, g), -Inf, +Inf, rel.tol = 1e-10)$value
	}))
}

# This computes the derivate of the normalization value Z of the cavity times the exact factor with respect to 
# the mean of the cavity for f

computedZdmOldF <- function(Y, mOldF, vOldF, mOldG, vOldG) {

	# First derivative

	dmarginal_valuedmOldF <- function(Y, mOldF, vOldF, g) {
		dnorm(Y * mOldF / sqrt(exp(g) + vOldF)) * Y / sqrt(exp(g) + vOldF)
	}

	m <- cbind(Y, mOldF, vOldF, mOldG, vOldG)

	as.vector(apply(m, 1, function(x) {
		Y <- x[ 1 ]
		mOldF <- x[ 2 ]
		vOldF <- x[ 3 ]
		mOldG <- x[ 4 ]
		vOldG <- x[ 5 ]
		integrate(function(g) dnorm(g, mOldG, sqrt(vOldG)) * dmarginal_valuedmOldF(Y, mOldF, vOldF, g), -Inf, +Inf, rel.tol = 1e-10)$value
	}))
}

# This computes the second derivate of the normalization value Z of the cavity times the exact factor with respect to 
# the mean of the cavity for f

computed2ZdmOldF2 <- function(Y, mOldF, vOldF, mOldG, vOldG) {

	# Second derivative

	d2marginal_valuedmOldF2 <- function(Y, mOldF, vOldF, g) {
		dnorm(Y * mOldF / sqrt(exp(g) + vOldF)) / (exp(g) + vOldF)^(3 / 2) * (- mOldF) * Y
	}

	m <- cbind(Y, mOldF, vOldF, mOldG, vOldG)

	as.vector(apply(m, 1, function(x) {
		Y <- x[ 1 ]
		mOldF <- x[ 2 ]
		vOldF <- x[ 3 ]
		mOldG <- x[ 4 ]
		vOldG <- x[ 5 ]
		integrate(function(g) dnorm(g, mOldG, sqrt(vOldG)) * d2marginal_valuedmOldF2(Y, mOldF, vOldF, g), -Inf, +Inf, rel.tol = 1e-10)$value
	}))
}

# This computes the derivate of the normalization value Z of the cavity times the exact factor with respect to 
# the mean of the cavity for g

computedZdmOldG <- function(Y, mOldF, vOldF, mOldG, vOldG) {

	marginal_value <- function(Y, mOldF, vOldF, g) {
		pnorm(Y * mOldF / sqrt(exp(g) + vOldF))
	}

	m <- cbind(Y, mOldF, vOldF, mOldG, vOldG)

	dgaussiandmean <- function(x, mean, var) {
		dnorm(x, mean, sqrt(var)) * (x - mean) / var
	}

	as.vector(apply(m, 1, function(x) {
		Y <- x[ 1 ]
		mOldF <- x[ 2 ]
		vOldF <- x[ 3 ]
		mOldG <- x[ 4 ]
		vOldG <- x[ 5 ]
		integrate(function(g) dgaussiandmean(g, mOldG, vOldG) * marginal_value(Y, mOldF, vOldF, g), -Inf, +Inf, rel.tol = 1e-10)$value
	}))
}

# This computes the second derivate of the normalization value Z of the cavity times the exact factor with respect to 
# the mean of the cavity for g

computed2ZdmOldG2 <- function(Y, mOldF, vOldF, mOldG, vOldG) {

	marginal_value <- function(Y, mOldF, vOldF, g) {
		pnorm(Y * mOldF / sqrt(exp(g) + vOldF))
	}

	m <- cbind(Y, mOldF, vOldF, mOldG, vOldG)

	d2gaussiandmean2 <- function(x, mean, var) {
		dnorm(x, mean, sqrt(var)) * (x - mean)^2 / var^2 + 
		dnorm(x, mean, sqrt(var)) * - 1 / var
	}

	as.vector(apply(m, 1, function(x) {
		Y <- x[ 1 ]
		mOldF <- x[ 2 ]
		vOldF <- x[ 3 ]
		mOldG <- x[ 4 ]
		vOldG <- x[ 5 ]
		integrate(function(g) d2gaussiandmean2(g, mOldG, vOldG) * marginal_value(Y, mOldF, vOldF, g), -Inf, +Inf, rel.tol = 1e-10)$value
	}))

}



