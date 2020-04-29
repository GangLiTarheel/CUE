load("Data/PTSD_CpG_Best_method_list.RData")
load("Annotations.RData")
load("Sample_Dataset.RData")
load("Data/PTSD.neighbors.RData")

# line 64 need to change the directory

X<-sample_data
m<-dim(X)[2] # number of samples

CUE.impute(X=X,m=m,tissue="PTSD"){
    ## RF
    library(randomForest)
    load("RF/best_RF.RData")

    RF.impute<-function(i){
    cov.probe<-neighbors_RF[[i]]
    cov.select <- X[cov.probe,]
    #temp=c(1,1,1)
    rf<-readRDS(paste("RF/model_",CpG_list_RF[i],".rds",sep=""))
    #y_RF[1,]<-predict(rf,t(X))
    return(predict(rf,t(X[cov.probe,])))
    }
    y_RF<-sapply(1:length(CpG_list_RF),RF.impute)
    colnames(y_RF)<-CpG_list_RF
    rownames(y_RF)<-colnames(X)

    ## XGBoost
    library(xgboost)
    load("XGB/best_XGB.RData")
    XGB.impute<-function(i){
    cov.probe<-neighbors_XGB[[i]]
    cov.select <- X[cov.probe,]
    bst<-xgb.load(paste("XGB/model_",CpG_list_XGBoost[i],".model",sep=""))
    #y_RF[1,]<-predict(rf,t(X))
    return(predict(bst,t(X[cov.probe,])))
    }

    y_XGB<-sapply(1:length(CpG_list_XGBoost),XGB.impute)
    rownames(y_XGB)=colnames(X)
    colnames(y_XGB)<-CpG_list_XGBoost

    ## KNN
    library(xgboost)
    load("XGB/best_XGB.RData")
    load("KNN/best_KNN.RData")
    KNN.impute<-function(i){
    cov.probe<-neighbors_KNN[[i]]
    cov.select <- X[cov.probe,]
    bst<-xgb.load(paste("KNN/model_",CpG_list_KNN[i],".model",sep=""))
    #y_RF[1,]<-predict(rf,t(X))
    return(predict(bst,t(X[cov.probe,])))
    }
    y_KNN<-sapply(1:length(CpG_list_KNN),KNN.impute)
    rownames(y_KNN)=colnames(X)
    colnames(y_KNN)<-CpG_list_KNN

    ## TCR
    X=log2(X/(1-X))

    source("/proj/yunligrp/users/gangli/methylation/EPIC_450K_VA/11_VarSel/refund_lib.R")
    load("TCR/best_TCR.RData")

    # Create the density 
    # Build the functional predictor for each group
    setwd("CUE/TCR/")
    isl_group <- c("", "Island", "N_Shelf", "N_Shore", "S_Shore", "S_Shelf")
    test.dens <- list()
    test.funcs <- list()

    for (l in isl_group) {
    if (l == "") {
        j <- "NA"
    } else {
        j <- l
    }
    # max(X)==13.28757
    test.dens[[j]] <- X[annotation.450K[, "Relation_to_UCSC_CpG_Island"] == l,]
    test.funcs[[j]] <- apply(test.dens[[j]], 2, function(x) {density(x,from=-13.38757,to=13.38757)$y})
    # train.dens[[j]] <- X[annotation.450K[, "Relation_to_UCSC_CpG_Island"] == l,]
    # train.funcs[[j]] <- apply(train.dens[[j]], 2, function(x) {density(x,from=-13.38757,to=13.38757)$y})   
    }
    #filename=paste("funcs.Rdata",sep = "")
    save(test.funcs,file = "funcs.Rdata")

    i=1
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

    y_PFR<-sapply(1:length(CpG_list_PFR),TCR.impute)
    rownames(y_PFR)=colnames(X)
    colnames(y_PFR)<-CpG_list_PFR

    m.impute<-cbind(y_PFR,y_RF,y_XGB,y_KNN)
    return(m.impute)
}

m.imputed<-CUE.impute(data,probe.list)

#setwd("")
save(m.imputed,file="y_impute.RData")