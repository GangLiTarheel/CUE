load("Data/PTSD_CpG_Best_method_list.RData")
load("PTSD/Annotations.RData")
#load("PTSD/Sample_Dataset.RData")
load("Data/PTSD.neighbors.RData")
load("Data/Probes.RData")

# line 64 need to change the directory

#X<-sample_data
#  the input dataset should have row as probes, columns as samples.
#m<-dim(X)[2] # number of samples

root_dir = getwd()

# Check the input
CUE_check <- function(X){
    # check if NA
    if (sum(is.nan(X))==0){
        print("Completeness checked! The input HM450 data don't have any missing NAs.")
    }else{
        stop("Error: the input HM450 data contains missingness! Please complete the dataset first!")
    }
    
    # check probe.list
    if (sum(probe.list %in% rownames(X)) == 248421){
        print("Probe list checked! The input HM450 data contain all the required probes.")
    }else{
        stop("Error: the input HM450 data does not contain all the required probes.! Please complete the dataset first!")
    }
    
    m<-dim(X)[2]
    cat("The input data contain all 248,421 required probes of",m,"samples.\n")
}

RF.impute<-function(i){
    cov.probe<-neighbors_RF[[i]]
    cov.select <- X[cov.probe,]
    #temp=c(1,1,1)
    rf<-readRDS(paste("PTSD/RF/model_",CpG_list_RF[i],".rds",sep=""))
    #y_RF[1,]<-predict(rf,t(X))
    return(predict(rf,t(X[cov.probe,])))
}

XGB.impute<-function(i){
    cov.probe<-neighbors_XGB[[i]]
    cov.select <- X[cov.probe,]
    bst<-xgb.load(paste("PTSD/XGB/model_",CpG_list_XGBoost[i],".model",sep=""))
    #y_RF[1,]<-predict(rf,t(X))
    return(predict(bst,t(X[cov.probe,])))
}

KNN.impute<-function(i){
    cov.probe<-neighbors_KNN[[i]]
    cov.select <- X[cov.probe,]
    bst<-xgb.load(paste("PTSD/KNN/model_",CpG_list_KNN[i],".model",sep=""))
    #y_RF[1,]<-predict(rf,t(X))
    return(predict(bst,t(X[cov.probe,])))
}

TCR.impute<-function(i){
    j <- paste(annotation.EPIC[CpG_list_PFR[i],"Relation_to_UCSC_CpG_Island"])
    if (j == "") {
        j <- "NA"
    }
    cov.probe<-neighbors_PFR[[i]]
    cov.select <- X[cov.probe,]
    X.test <- create_predictors(t(cov.select), t(test.funcs[[j]]))
    Y.test <- X.test %*% TCR_model_best[[i]]#fit$coefs  
    return(2 ^ Y.test / (2 ^ Y.test + 1))
}

load("PTSD/RF/best_RF.RData")
load("PTSD/XGB/best_XGB.RData")
load("PTSD/KNN/best_KNN.RData")
source("R/refund_lib.R")
load("PTSD/TCR/best_TCR.RData")

library(randomForest)
library(xgboost)

# impute
CUE.impute <- function(X=X,m=m,tissue="PTSD"){
    ## RF
    y_RF<-sapply(1:length(CpG_list_RF),RF.impute)
    colnames(y_RF)<-CpG_list_RF
    rownames(y_RF)<-colnames(X)
    
    ## XGBoost
    y_XGB<-sapply(1:length(CpG_list_XGBoost),XGB.impute)
    rownames(y_XGB)=colnames(X)
    colnames(y_XGB)<-CpG_list_XGBoost
    
    ## KNN
    y_KNN<-sapply(1:length(CpG_list_KNN),KNN.impute)
    rownames(y_KNN)=colnames(X)
    colnames(y_KNN)<-CpG_list_KNN
    
    ## TCR
    X=log2(X/(1-X))
    # Create the density 
    # Build the functional predictor for each group
    # setwd("PTSD/CUE/TCR/")
    isl_group <- c("", "Island", "N_Shelf", "N_Shore", "S_Shore", "S_Shelf")
    test.dens <- list()
    test.funcs <- list()
    for (l in isl_group){
        if(l==""){j<-"NA"}else{j<-l}
        test.dens[[j]] <- X[ (annotation.450K[, "Relation_to_UCSC_CpG_Island"] == l),]
        test.funcs[[j]]<-apply(test.dens[[j]], 2, function(x){density(x,from=-13.38757,to=13.38757)$y})
    }
    y_PFR<-sapply(1:length(CpG_list_PFR),TCR.impute)
    rownames(y_PFR)=colnames(X)
    colnames(y_PFR)<-CpG_list_PFR
    
    m.impute<-cbind(y_PFR,y_RF,y_XGB,y_KNN)
    return(m.impute)
}

X=sample_data
CUE_check(X)
m.imputed<-CUE.impute(X,m,"PTSD")

setwd(root_dir)
save(m.imputed,file="y_impute.RData")