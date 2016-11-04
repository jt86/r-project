This folder contains the R code of all the methods described in the main 
manuscript. To simulate them simply type

R --no-save < evaluate_methods.R

on a Linux prompt. This will evaluate the performance of each method
using example data from the SUNAttribute experiment. 
This data corresponds to one of 10 repeats when training the attribute classifier 'railing'. 
The result of each method is stored in a different sub-folder, inside the folder results.
For this case, we got: 

Method 	  Error   Accuracy
GPC       0.38      62%
GPCconf   0.30      70%
GPC+	  0.28      72%
SVM	  0.31      69%
SVM+	  0.31      69%

Running time is stated in the Table 3 of the main paper.



