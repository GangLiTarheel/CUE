# CUE
### CpG impUtation Ensemble for DNA Methylation Levels Across the HumanMethylation450 (HM450) BeadChip and HumanMethylation EPIC (HM850) BeadChip Platforms

DNA methylation at CpG dinucleotides is one of the most extensively studied epigenetic marks. With technological advancements, geneticists can profile DNA methylation with multiple reliable approaches. 
However, different profiling platforms can differ substantially in the density of and actual measurements for the CpGs they assess, consequently hindering integrated analysis across platforms. 
Here, we present CpG impUtation Ensemble (CUE), which leverages multiple classical statistical and modern machine learning methods, to impute from the Illumina HumanMethylation450 (HM450) BeadChip to the Illumina HumanMethylationEPIC (HM850) BeadChip. 

CUE is maintained by Gang Li [franklee@live.unc.edu] and Laura Raffield.

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
# Required library for imputaiton.
Please install those R packages before imputations.

```{r init, message=TRUE}
library(randomForest) # RF
library(xgboost) # XGBoost
library("class") # KNN
source("R/refund_lib.R")

library(parallel) # Use this package if you need parallel computing to accelerate computation speed

```

## Download the pretrained imputaiton models 

Note: The compressed pre-trained models need around ~ 58 GB (ELGAN ~36 GB; and PTSD ~21 GB) memory of storage.

Whole blood: ftp://yunlianon:anon@rc-ns-ftp.its.unc.edu/CUE/PTSD_model.tar.gz
ELGAN: ftp://yunlianon:anon@rc-ns-ftp.its.unc.edu/CUE/ELGAN_model.tar.gz

Please save the above pretrained models (tar.gz files) in the same directory of your CUE packages.

To extract the pre-trained models:

tar -xf PTSD_model.tar.gz
tar -xf ELGAN_model.tar.gz

## CUE imputation
Here we show the imputation for a toy dataset with three samples, 248K HM450K probes. 
The reason we use only 248K probes is that when we trained models, we only retain the 248,421 HM450 CpG sites (sites overlapping between ELGAN and PTSD, and without missingness in our samples) to train since we don't want the incompleteness or the pre-imputation on HM450 affects our evaluations on HM850 imputation.
The full list of 248K HM450K probes can be found in Probes.RData.
```{r perform imputation}
sample_data<-load("sample_Data.RData")
X<-sample_data
m<-dim(X)[2] # number of samples
#m.imputed<-CUE.impute(X=X,m=m,tissue="PTSD")
source("impute.R")

save(m.imputed,file="y_impute.RData")

```

Note: we impute all 339K probes, and one must used the following quality control to retain the well imputed model.

## Qualtiy Control
Using the shiny app: CUE_QC.R to select the QC+ probe list.
```{r CUE_QC}
runApp(appDir = CUE_QC)
```
App would save a list of well imputed probes. Users can uses the following code

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

