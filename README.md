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
Please install the following R packages before performing imputation.

```{r init, message=TRUE}
library(randomForest) # RF
library(xgboost) # XGBoost
library("class") # KNN
source("R/refund_lib.R")

library(parallel) # Use this package if you need parallel computing to accelerate computation speed

```

## Download the pretrained imputation models 


ftp -i rc-ns-ftp.its.unc.edu
  user: yunlianon
  secret (upon request)

Whole blood: CUE/PTSD_model.tar.gz
ELGAN: 


## CUE imputation
Here we show how CUE performs imputation for a toy dataset with three samples and 248K HM450K probes.
```{r perform imputation}
sample_data<-load("sample_Data.RData")
#m.imputed<-CUE(data,probe.list)
source("impute.R")

```

Note: we impute all 339K HM850 specific probes which had complete data in our reference whole blood and placenta datasets. Users of CUE must use the following quality control steps to retain the well imputed probes only for use in subsequent analysis (such as epigenome wide association studies).

## Quality Control
Using the shiny app: CUE_QC.R to select the QC+ probe list.
```{r CUE_QC}
runApp(appDir = CUE_QC)
```
The shiny app will save a list of well imputed probes. Users can use the following code to save this list. 

```{r subset}
load("csv from QC app")
load("All probes from CUE impute.")
m.impute <- m.impute[,paste(QC_probes)]

save ("RData")
```


```{r output from CUE}
## Output : DNA methylation matrix
m.imputed[1:10,1:10]
```

## Citation
CUE

