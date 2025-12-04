<div align="center">
 
  <img width="120" height="180" alt="image" src="https://github.com/user-attachments/assets/529166f8-c03d-4542-872a-748b64449580" /> 
  
  # Hypatia
  
  **A statistical framework for comparative isoform profiling across cell populations.**

 To get started, please visit [Hypatia's vignette](https://gaolabtools.github.io/Hypatia/vignettes/Hypatia.html). 

</div>


## About

Hypatia (hi-pay-shuh) is a computational toolkit for the investigation of population-specific isoforms from long-read single-cell RNA-sequencing data, featuring three modes of differential analyses:
1) Isoform usage: Identifies differential isoform usage shifts (DIUs).
2) Isoform diversity: Classifies genes according to isoform species heterogeneity and identifies discrepances (DIVs).
3) Isoform expression: Detects differentially expressed isoforms (DEIs).

## Installation

Hypatia is an R library available through Github.

```r
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

devtools::install_github("gaolabtools/Hypatia")
```

Requires R version 4.1 or higher.

## Citation
If you found Hypatia useful in your research, please cite us:
```
[Citation]
```
