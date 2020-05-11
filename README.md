# CUE
### CpG impUtation Ensemble for DNA Methylation Levels Across the HumanMethylation450 (HM450) BeadChip and HumanMethylation EPIC (HM850) BeadChip Platforms

DNA methylation at CpG dinucleotides is one of the most extensively studied epigenetic marks. 
With technological advancements, geneticists can profile DNA methylation with multiple reliable approaches. 
However, different profiling platforms can differ substantially in the CpGs they assess, consequently hindering integrated analysis across platforms. 
Here, we present CpG impUtation Ensemble (CUE), which leverages multiple classical statistical and modern machine learning methods, to impute from the Illumina HumanMethylation450 (HM450) BeadChip to the Illumina HumanMethylationEPIC (HM850) BeadChip. 

CUE is maintained by Gang Li [franklee@live.unc.edu] and Laura Raffield [laura_raffield@unc.edu].

## News and Updates
* Version 0.0.9 released

## Brief introduction
From this study, we provide two sets of imputation models: one for whole blood and the other for placenta. 
Investigators can therefore complete their own imputation of placental or whole blood HM850 CpG sites using their own HM450 data, without access to their own reference panel. 
Our method is also easily useable for imputation in other tissues, provided the user can provide a reference dataset assayed on both HM450 and HM850.

## Installation 

CUE package can be directly installed as following:
```
git clone https://github.com/GangLiTarheel/CUE.git
cd CUE
```

<!-- ```{r installation}
install.packages("devtools")

devtools::install_github("")
``` -->
<!-- 

## Set up the library
```{r init, message=TRUE}
library("CUE")
``` -->
### Required library for imputation.
Please install the following version of R and R packages before performing imputation.



* R version 3.6.0 (2019-04-26)

* class_7.3-15        
* xgboost_0.82.1      
* randomForest_4.6-14
* parallel

Note: to install an R packages, simply run install.packages("Name_of_package"). For example,
```
install.packages("xgboost")
```

Note: you can also install the version of xgboost we used from the source code:
https://cran.r-project.org/src/contrib/Archive/xgboost/xgboost_0.82.1.tar.gz


```{r init, message=TRUE}
library(randomForest) # RF
library(xgboost) # XGBoost
library("class") # KNN
source("R/refund_lib.R")

library(parallel) # Use this package if you need parallel computing to accelerate computation speed

```

## Download the pretrained imputation models 

Note: The compressed pre-trained models need around ~ 58 GB (ELGAN ~36 GB; and PTSD ~21 GB) memory of storage.


Please save the above pretrained models (tar.gz files) in the same directory of your CUE packages.

To extract the pre-trained models:

Make two directories for two sets of pre-trained models:
```
pwd # this should be your root directory of CUE
mkdir PTSD
mkdir ELGAN
```

Download the pre-trained models amd move the PTSD and ELGAN pretrained models to the corresponding folder.
### Whole blood (PTSD): 
ftp://yunlianon:anon@rc-ns-ftp.its.unc.edu/CUE/PTSD_model.tar.gz
### Placenta (ELGAN): 
ftp://yunlianon:anon@rc-ns-ftp.its.unc.edu/CUE/ELGAN_model.tar.gz


###  then decompress.

```
cd PTSD
tar -xf PTSD_model.tar.gz
cd ../

cd ELGAN
tar -xf ELGAN_model.tar.gz
cd ../
```

## CUE imputation
Here we show the imputation for a toy dataset with three samples, 248K HM450K probes. 
The reason we use only 248K probes is that when we trained models, we only retain the 248,421 HM450 CpG sites (sites overlapping between ELGAN and PTSD, and without missingness in our samples) to train since we don't want the incompleteness or the pre-imputation on HM450 affects our evaluations on HM850 imputation.
The full list of 248K HM450K probes can be found in Probes.RData.
###  The input dataset X should have row as probes, columns as samples.
```{r perform imputation}
setwd(“/YOUR_DIRECTORY/CUE”)
# Replace YOUR_DIRECTORY with your directory where you downloaded CUE.

sample_data<-load("PTSD/Sample_Dataset.RData")
X<-sample_data
m<-dim(X)[2] # number of samples

source("impute.R")

# the output will be saved as y_impute.RData
#save(m.imputed,file="y_impute.RData")

```

<!-- #m.imputed<-CUE.impute(X=X,m=m,tissue="PTSD") -->

Note: we impute all 339K HM850 specific probes which had complete data in our reference whole blood and placenta datasets. Users of CUE must use the following quality control steps to retain the well imputed probes only for use in subsequent analysis (such as epigenome wide association studies).

## Quality Control
We provdies two sets of QC+ probes list for two datsets with the following thresholds:
QC for ELGAN: RMSE<0.1 and Accruracy > 90%; 
QC for PTSD: RMSE<0.05 and Accruracy > 95%.

```{r subset}
QC_probes<-read.csv("Placenta (ELGAN)QC_probe_list.csv") # for Placenta
#QC_probes<-read.csv("Whole Blood (PTSD)QC_probe_list") # for whole blood
load("y_impute.RData")
m.QC <- m.imputed[,paste(QC_probes)]

save(m.QC, file="m.imputed.QC.RData")
```

(Optional): We also provide an shiny app, CUE_QC.R, an visualization tool to select the QC+ probe list.
```{r CUE_QC}
runApp('CUE_QC.R')
```
The shiny app will save the list of well imputed probes. Users can use the previous code to subset the output probes. 


```{r output from CUE}
## Output : DNA methylation matrix
head(m.QC[,1:10])
dim(m.QC)
```

## Citation
CUE

