load(paste("/home/j/jt/jt306/Documents/CVPR2016_Rcode/data/1data.dat",sep = ""))
Xtrain <- data$x[ data$itrain, ]
# Ytrain <- as.vector(c(data$y_train))
# Xtest <- data$x[ data$itest, ]
# Ytest <- as.vector(c(data$y_test))
# 
# 
typeof(data$x)
typeof(Xtrain)
dim(Xtrain)

dataset <- (read.table(paste("/Volumes/LocalDataHD/j/jt/jt306/Documents/CVPR2016_Rcode/saved_datasets/tech",datasetnum,"data",sep = "")))
dims=dim(dataset)
is.character(dataset)
is.numeric(dataset)
typeof(dataset)
dims
dataset2 <- matrix(as.double(unlist(dataset)),dims[1],dims[2])
typeof(dataset)
dim(dataset)

# typeof(data$x_star)
# typeof(data$itrain)
# typeof(data$itest)
# typeof(data$y_train)
# typeof(data$y_test)
# 
# 
# typeof(Xtrain)
# typeof(Ytrain)
# typeof(Xtest)
# typeof(Ytest)
# 
# 
# 
# dim(Xtrain)
# dim(Ytrain)
# dim(Xtest)
# dim(Ytest)

# (seq(-3, 3, length = 5))
# print(commandArgs())



