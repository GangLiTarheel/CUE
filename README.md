# CUE
### CpG impUtation Ensemble for DNA Methylation Levels Across the HumanMethylation450 (HM450) BeadChip and HumanMethylation EPIC (HM850) BeadChip Platforms

DNA methylation at CpG dinucleotides is one of the most extensively studied epigenetic marks. With technological advancements, geneticists can profile DNA methylation with multiple reliable approaches. 
However, different profiling platforms can differ substantially in the CpGs they assess, consequently hindering integrated analysis across platforms. 
Here, we present CpG impUtation Ensemble (CUE), which leverages multiple classical statistical and modern machine learning methods, to impute from the Illumina HumanMethylation450 (HM450) BeadChip to the Illumina HumanMethylationEPIC (HM850) BeadChip. 

CUE is maintained by Gang Li [franklee@live.unc.edu] and Laura Raffield [laura_raffield@unc.edu].

## News and Updates
* Version 0.0.1 released

## Brief introduction
From this study, we provide two sets of imputation models: one for whole blood and the other for placenta. 
Investigators can therefore complete their own imputation of placental or whole blood HM850 CpG sites using their own HM450 data, without access to their own reference panel. 
Our method is also easily useable for imputation in other tissues, provided the user can provide a reference dataset assayed on both HM450 and HM850.

## Installation （Under construction）

CUE package can be directly installed from GitHub with:
```{r installation}
install.packages("devtools")

devtools::install_github("")
```


## Set up the library
```{r init, message=TRUE}
library("CUE")
```
# Required library for imputation.
Please install the following version of R and R packages before performing imputation.

R version 3.6.0 (2019-04-26)

class_7.3-15        
xgboost_0.82.1      
randomForest_4.6-14
parallel

```{r init, message=TRUE}
library(randomForest) # RF
library(xgboost) # XGBoost
library("class") # KNN
source("R/refund_lib.R")

library(parallel) # Use this package if you need parallel computing to accelerate computation speed

```

## Download the pretrained imputation models 

Note: The compressed pre-trained models need around ~ 58 GB (ELGAN ~36 GB; and PTSD ~21 GB) memory of storage.

# Whole blood (PTSD): 
ftp://yunlianon:anon@rc-ns-ftp.its.unc.edu/CUE/PTSD_model.tar.gz
# Placenta (ELGAN): 
ftp://yunlianon:anon@rc-ns-ftp.its.unc.edu/CUE/ELGAN_model.tar.gz

Please save the above pretrained models (tar.gz files) in the same directory of your CUE packages.

To extract the pre-trained models:

pwd # this should be your root directory of CUE
mkdir PTSD
cd PTSD
# copy the pre PTSD pretrained model to this folder and then decompress.
tar -xf PTSD_model.tar.gz
cd ../
mkdir ELGAN
cd ELGAN
# copy the pre ELGAN pretrained model to this folder and then decompress.
tar -xf ELGAN_model.tar.gz
cd ../


## CUE imputation
Here we show the imputation for a toy dataset with three samples, 248K HM450K probes. 
The reason we use only 248K probes is that when we trained models, we only retain the 248,421 HM450 CpG sites (sites overlapping between ELGAN and PTSD, and without missingness in our samples) to train since we don't want the incompleteness or the pre-imputation on HM450 affects our evaluations on HM850 imputation.
The full list of 248K HM450K probes can be found in Probes.RData.
#  The input dataset X should have row as probes, columns as samples.
```{r perform imputation}
sample_data<-load("sample_Data.RData")
X<-sample_data
m<-dim(X)[2] # number of samples
#m.imputed<-CUE.impute(X=X,m=m,tissue="PTSD")
source("impute.R")

save(m.imputed,file="y_impute.RData")

```

Note: we impute all 339K HM850 specific probes which had complete data in our reference whole blood and placenta datasets. Users of CUE must use the following quality control steps to retain the well imputed probes only for use in subsequent analysis (such as epigenome wide association studies).

## Quality Control
Using the shiny app: CUE_QC.R to select the QC+ probe list.
```{r CUE_QC}
runApp(appDir = CUE_QC)
```
The shiny app will save a list of well imputed probes. Users can use the following code to save this list. 

```{r subset}
#load("csv from QC app")
load("y_impute.RData")
m.QC <- m.imputed[,paste(QC_probes)]

save(m.QC, file="m.imputed.QC.RData")
```


```{r output from CUE}
## Output : DNA methylation matrix
head(m.QC[,1:10])
dim(m.QC)
```

## Citation
CUE

