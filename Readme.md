# A modified package of EstimateClonality

EstimateClonality is an R package which uses read count and copy number information to temporally and clonally dissect SNVs, written by Nicholas McGranahan (nicholas.mcgranahan@cancer.org.uk).

This version fixes some issues in the original code and makes a little optimization.

## Major changes

1. The original R package is not compatible with sequenza version 3.0. Now the earlyORlate() function is compatible.
2. Added support for FACETS and CNVKIT (non ASCAT data).

## Minor changes

1. Part of the code is simplified.
2. 
