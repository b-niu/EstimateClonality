# A modified package of EstimateClonality

EstimateClonality is an R package which uses read count and copy number information to temporally and clonally dissect SNVs, written by Nicholas McGranahan (nicholas.mcgranahan@cancer.org.uk).

https://bitbucket.org/nmcgranahan/pancancerclonality/src/master/

This version fixes some issues in the original code and makes a little optimization.

## Major changes

1. The original R package is not compatible with sequenza version 3.0. Now the earlyORlate() function is compatible.
2. Added support for FACETS and CNVKIT (non ASCAT data).

## Minor changes

1. Part of the code is simplified.
2. 

## Usage written by the original author:

```
# The following an is an example script to assess the clonality of one TCGA sample, 'TCGA-BT-A42C'.

rm(list=ls())

setwd("~/Downloads/McGranahan_data/Estimate_Clonality_Package/")

install.packages("EstimateClonality_1.0.tar.gz",repos=NULL,type='source')

library("EstimateClonality")

clonality.estimation(mutation.table.loc="BLCA.mutation.table.txt" ,seg.mat.loc="tcga.blca.seg.hg19.rdata" ,data.type='TCGA_BLCA' ,TCGA.barcode="TCGA-BT-A42C" ,ANALYSIS.DIR="example2/")
```
