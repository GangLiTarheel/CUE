# !/usr/bin/env Rscript
# 04/23/2020
#
#' @title CUE
#' 
#' @description This function CUE.impute is designed to perform CpG impUtation Ensemble (CUE), which leverages multiple classical statistical and modern machine learning methods, to impute from the Illumina HumanMethylation450 (HM450) BeadChip to the Illumina HumanMethylationEPIC (HM850) BeadChip.
#' It takes as input the beta methylation matrix for the probse of the Illumina HumanMethylation450 (HM450) BeadChip.
#' It outputs the methylation matrix of for imputed probes of the Illumina HumanMethylationEPIC (HM850) BeadChip.
#' @usage CUE.impute(HM450)
#' @param X is matrix of XXXK HM450 probes.
#'
#' @return SMNNcorrect returns the HM850 imputed probes
#' @author Gang Li <franklee@live.unc.edu>, Laura Raffield, Yun Li <yunli@med.unc.edu>
#' @examples 
#' # Load the example data data_CUE
#' data("data_CUE")
#' 
#' # Perform imputation using CUE.impute
#' HM850.impute <- CUE.impute(HM450)
#' @import randomForest, xgboost, class
#' @export
CUE.impute <- function(HM450){
  load("/proj/yunligrp/users/gangli/methylation/Hudson/10_user/ELGAN_CpG_Best_method_list.RData")
    setwd("/proj/yunligrp/users/gangli/methylation/Hudson/10_user/CUE/")
    load("Annotations.RData")
    load("Sample_Dataset.RData")
    load("ELGAN.neighbors.RData")

    X<-sample_data
    m<-dim(X)[2] # number of samples


    ## RF
    library(randomForest)
    load("RF/best_RF.RData")
    # y_RF<-matrix(0,length(CpG_list_RF),m)
    # rownames(y_RF)=CpG_list_RF
    # colnames(y_RF)=colnames(X)

    # i=1
    # cov.probe<-neighbors_RF[[1]]
    # cov.select <- X[cov.probe,]
    # #temp=c(1,1,1)
    # rf<-readRDS(paste("RF/model_",CpG_list_RF[i],".rds",sep=""))
    # y_RF[1,]<-predict(rf,t(X))
    # #Meth<-data.frame(cbind(temp,t(cov.select)))
    # #colnames(Meth)[1]<-probe[i]

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
    # y_XGB<-matrix(0,length(CpG_list_XGBoost),m)
    # rownames(y_XGB)=CpG_list_XGBoost
    # colnames(y_XGB)=colnames(X)

    i=1
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
    # library("class") 
    # load("KNN/best_KNN.RData")
    # y_KNN<-matrix(0,length(CpG_list_KNN),m)
    # rownames(y_KNN)=CpG_list_KNN
    # colnames(y_KNN)=colnames(X)

    # i=1
    # cov.probe<-neighbors_KNN[[1]]
    # cov.select <- X[cov.probe,]
    # temp=c(1,1,1)
    # Meth<-data.frame(cbind(temp,t(cov.select)))
    # colnames(Meth)[1]<-probe[i]
    # knn=knn(train = Meth[,2:dim(cov.select)[1]],
    #     test = Meth[,2:dim(cov.select)[1]],
    #     cl = Meth[,1],
    #     k=50,
    #     prob=T)

    ## LOGIT 
    load("LOGIT/best_LOGIT.RData")

    LOGIT.impute<-function(i){
        cov.probe<-neighbors_LOGIT[[i]]
        cov.select <- X[cov.probe,]
        #temp=c(1,1,1)
        logit<-readRDS(paste("LOGIT/model_",CpG_list_Logistic[i],".rds",sep=""))
        Meth<-data.frame(t(cov.select))
        #y_RF[1,]<-predict(rf,t(X))
        return(predict(logit,Meth))
    }
    y_LOGIT<-sapply(1:length(CpG_list_Logistic),LOGIT.impute)
    rownames(y_LOGIT)=colnames(X)
    colnames(y_LOGIT)<-CpG_list_Logistic

    ## TCR
    X=log2(X/(1-X))

    source("/proj/yunligrp/users/gangli/methylation/EPIC_450K_VA/11_VarSel/refund_lib.R")
    load("TCR/best_TCR.RData")

    # Create the density 
    # Build the functional predictor for each group
    setwd("/proj/yunligrp/users/gangli/methylation/Hudson/10_user/CUE/TCR/")
    isl_group <- c("", "Island", "N_Shelf", "N_Shore", "S_Shore", "S_Shelf")
    # train.dens <- list()
    # train.funcs <- list()
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

    setwd("/proj/yunligrp/users/gangli/methylation/Hudson/10_user/CUE/")
    save(y_PFR,y_RF,y_XGB,y_KNN,Y_LOGIT,file="y_impute.RData")

    return (HM850.impute)
}
