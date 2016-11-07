library(RcppCNPy)
load("/Volumes/LocalDataHD/j/jt/jt306/Documents/CVPR2016_Rcode/data/1data.dat")
print((data$itrain))

dataset <- (read.table("/Volumes/LocalDataHD/j/jt/jt306/Documents/CVPR2016_Rcode/saved_datasets/tech39data"))
labels <-(as.double(unlist(read.table("/Volumes/LocalDataHD/j/jt/jt306/Documents/CVPR2016_Rcode/saved_datasets/tech39labels"))))
dim(labels)
Ytrain <- labels[train_indices]
Ytest <- labels[test_indices]
dim(Ytrain)
length(Ytrain)
length(Ytest)

labels
dim(dataset)
train_indices <- as.integer(read.table('/home/j/jt/jt306/Documents/CVPR2016_Rcode/saved_indices/top500-tech39-3-4-train_instances_indices')[[1]])+1
test_indices <- as.integer(read.table('/home/j/jt/jt306/Documents/CVPR2016_Rcode/saved_indices/top500-tech39-3-4-test_instances_indices')[[1]])+1
selected_indices <- as.integer(read.table('/home/j/jt/jt306/Documents/CVPR2016_Rcode/saved_indices/top500-tech39-3-4-selected_feat_indices')[[1]])+1
priv_indices <- as.integer(read.table('/home/j/jt/jt306/Documents/CVPR2016_Rcode/saved_indices/top500-tech39-3-4-unselected_feat_indices')[[1]])+1
 

dims <- dim(dataset)
dataset =matrix(dataset,dims[1],dims[2])
dim(dataset)
typeof(dataset)
dims

dataset[1,1:130]

length(train_indices)

dim(train_data)
dim(test_data)

print (selected_indices)

length(selected_indices)
length(priv_indices)

train_normal <- dataset[train_indices,selected_indices]
dim(train_normal)
test_normal <- dataset[test_indices,selected_indices]
dim(test_normal)
train_priv <- dataset[train_indices,priv_indices]
dim(train_priv)
test_priv <- dataset[test_indices,selected_indices]
dim(test_priv)


 