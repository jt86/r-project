set.seed(0)

for (i in 1 : 10) {

	data <- list()
	print(i)	#i-th train/test split
	
	c1 <- read.table("c1.txt")$V1[ 1 ]
	c2 <- read.table("c2.txt")$V1[ 1 ]
	
	x1 <- as.matrix(read.table(paste("../../data_SUNAttributes",i,"/", c1,".cnn", sep = "")))	
	x2 <- as.matrix(read.table(paste("../../data_SUNAttributes",i,"/", c2,".cnn", sep = ""))) 	
	x <- rbind(x1, x2)
	
	xstar1 <- as.matrix(read.table(paste("../../data_SUNAttributes",i,"/", c1,".conf", sep = "")))
	xstar2 <- as.matrix(read.table(paste("../../data_SUNAttributes",i,"/", c2,".conf", sep = "")))
	xstar <- rbind(xstar1, xstar2)

	labels <- c(rep(1, nrow(x1)), rep(-1, nrow(x2)))

	n_train <- 200
	n_test <- 50

	sel_tr1 <- sample(1 : nrow(x1), n_train)					#200 pos samples 
	sel_tr2 <- sample(1 : nrow(x2), n_train) + nrow(x1)				#200 neg samples 
	
	sel_test1 <- sample((1 : nrow(x1))[ -c(sel_tr1) ], n_test)			#50 pos
	sel_test2 <- sample((1 : nrow(x2))[ -c(sel_tr2 - nrow(x1)) ], n_test) + nrow(x1)#50 neg

	sel_train <- c(sel_tr1, sel_tr2)
	sel_test <- c(sel_test1, sel_test2)

	data$x <- x
	data$x_star <- xstar
	data$itrain <- sel_train
	data$itest <- sel_test
	data$y_train <- labels[  sel_train ]
	data$y_test <- labels[  sel_test ]
	
	s <- sample(1 : length(data$itrain))

	data$itrain <- data$itrain[s]
	data$y_train <- data$y_train[s]

	save(data, file = paste(i,"data.dat", sep = ""))
	data <- NULL
}

