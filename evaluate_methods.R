
#This R script evaluate all methods compared in the manuscript 
#This data corresponds to one of 10 repeats when training 
#the attribute classifier 'railing', test scenario A.


commandArgs <- function() 1
source("simulateGPC.R")
source("simulateGPC_conf.R")
source("simulateGPC_plus.R")
source("simulateSVM.R")
source("simulateSVMplusQP.R")
