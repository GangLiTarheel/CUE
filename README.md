# CUE
### CpG impUtation Ensemble for DNA Methylation Levels Across the HumanMethylation450 (HM450) BeadChip and HumanMethylation EPIC (HM850) BeadChip Platforms

DNA methylation at CpG dinucleotides is one of the most extensively studied epigenetic marks. With technological advancements, geneticists can profile DNA methylation with multiple reliable approaches. 
However, different profiling platforms can differ substantially in the density of and actual measurements for the CpGs they assess, consequently hindering integrated analysis across platforms. 
Here, we present CpG impUtation Ensemble (CUE), which leverages multiple classical statistical and modern machine learning methods, to impute from the Illumina HumanMethylation450 (HM450) BeadChip to the Illumina HumanMethylationEPIC (HM850) BeadChip. 

CUE is maintained by Gang Li [franklee@live.unc.edu].

## News and Updates
* Version 0.0.1 released

## Brief introduction
From this study, we provide two sets of imputation models: one for whole blood and the other for placenta. 
Investigators can therefore complete their own imputation of placental or whole blood HM850 CpG sites using their own HM450 data, without access to their own reference panel. 
Our method is also easily useable for imputation in other tissues, provided the user can provide a reference dataset assayed on both HM450 and HM850.

## Installation

CUE package can be directly installed from GitHub with:
```{r installation}
install.packages("devtools")

devtools::install_github("")
```


## Set up the library
```{r init, message=TRUE}
library("CUE")
```


