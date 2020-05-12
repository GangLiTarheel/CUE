
## Set the working directory
#setwd(“/YOUR_DIRECTORY/CUE”)
# Replace YOUR_DIRECTORY with your directory where you downloaded CUE.

	root_dir = getwd()


## load sample dataset
	load("PTSD/Sample_Dataset.RData")
	X<-sample_data                            # Input HM450
	m<-dim(X)[2] # number of samples          # Number of smaples

## load required datasets 
	load("PTSD/Annotations.RData")                # load Annotation data
	load("Data/PTSD.neighbors.RData")             # load All neighbors
	load("Data/Probes.RData")                     # load All probes' names

## load required datasets for PTSD imputation
	load("Data/PTSD_CpG_Best_method_list.RData")  # load the best method list 
	load("PTSD/RF/best_RF.RData")                 # load neighbors for RF
	load("PTSD/XGB/best_XGB.RData")               # load neighbors for XGB
	load("PTSD/KNN/best_KNN.RData")               # load neighbors for KNN
	load("PTSD/TCR/best_TCR.RData")               # load neighbors for PFR

## load all the required packages
	library(randomForest)
	library(xgboost)
	library(xgboost)
	source("R/refund_lib.R")

## load the CUE imputation function
	source("R/impute.R")

## Check the input
	CUE_check(X)                        # check if the inuput X match the requirements

## Imputation (Single CPU might takes 4-5 days to imputation for a middle sized methylation dataset (n=~100 samples).)
	m.imputed<-CUE.impute(X,m,"PTSD")

## Save the imputed probes as y_impute.RData
	setwd(root_dir)
	save(m.imputed,file="y_impute.RData")
