#  the input dataset X should have row as probes, columns as samples.

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
    cat("Check input: Done!!!\n")
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

TCR.input<-function(X){
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
    newlist <- list(X,test.funcs)
    return(newlist)
}

TCR.impute<-function(i){
    j <- paste(annotation.EPIC[CpG_list_PFR[i],"Relation_to_UCSC_CpG_Island"])
    if (j == "") {
        j <- "NA"
    }
    cov.probe<-neighbors_PFR[[i]]
    cov.select <- X_logit[cov.probe,] #logit transformation
    X.test <- create_predictors(t(cov.select), t(test.funcs[[j]]))
    Y.test <- X.test %*% TCR_model_best[[i]]#fit$coefs  
    return(2 ^ Y.test / (2 ^ Y.test + 1)) # inverse logit tranformation
}

# impute
CUE.impute <- function(X=X,m=m,tissue="PTSD"){
    ## RF
    y_RF<-sapply(1:length(CpG_list_RF),RF.impute)
    colnames(y_RF)<-CpG_list_RF
    rownames(y_RF)<-colnames(X)
    
    cat("Imputation with random forests: Done!!!\n")

    ## XGBoost
    y_XGB<-sapply(1:length(CpG_list_XGBoost),XGB.impute)
    rownames(y_XGB)=colnames(X)
    colnames(y_XGB)<-CpG_list_XGBoost
    
    cat("Imputation with XGBoosting: Done!!!\n")

    ## KNN
    y_KNN<-sapply(1:length(CpG_list_KNN),KNN.impute)
    rownames(y_KNN)=colnames(X)
    colnames(y_KNN)<-CpG_list_KNN
    
    cat("Imputation with KNN: Done!!!\n")

    y_PFR<-sapply(1:length(CpG_list_PFR),TCR.impute)
    rownames(y_PFR)=colnames(X)
    colnames(y_PFR)<-CpG_list_PFR
    
    cat("Imputation with penalized functional regression: Done!!!\n")
    
    cat("Imputation with logistic regression: Omitted for PTSD!!!\n")

    m.impute<-cbind(y_PFR,y_RF,y_XGB,y_KNN)
    return(m.impute)
}

