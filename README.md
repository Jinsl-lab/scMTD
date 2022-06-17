# scMTD
A R package.
## scMTD: a statistical multidimensional imputation method for single-cell RNA-seq data leveraging transcriptome dynamic information


## Installation
You can install scMTD from GitHub with:

``` r
install.packages("devtools")         
library(devtools)           
install_github("Jinsl-lab/scMTD")
```

## Quick start
The input data of scMTD is a gene expression matrix, columns and rows represent cells and genes respectively. 

``` r
library("scMTD")
data(data)
imputed_data<-scMTD(data,do.nor=TRUE,do.log=TRUE,cores=2)
