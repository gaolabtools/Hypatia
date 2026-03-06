<div align="center">
 
  <img width="120" height="180" alt="image" src="https://github.com/user-attachments/assets/529166f8-c03d-4542-872a-748b64449580" /> 
  
  # Hypatia
  
  **A statistical framework for comparative isoform profiling across cell populations.**

  [![R version](https://img.shields.io/badge/R-4.1-276DC3?style=for-the-badge&logo=r&logoColor=white)](https://www.r-project.org/)
  ![Build](https://img.shields.io/badge/build-passing-brightgreen?style=for-the-badge)
  [![Download](https://img.shields.io/badge/Download-blue?style=for-the-badge&logo=github&logoColor=white)](https://github.com/gaolabtools/Hypatia/releases)
  [![stars](https://img.shields.io/github/stars/gaolabtools/Hypatia?style=for-the-badge&logo=github&color=FFD700)](https://github.com/gaolabtools/Hypatia/stargazers)

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

## Usage

For detailed documentation and tutorial, please visit [Hypatia's vignette](https://gaolabtools.github.io/Hypatia/vignettes/Hypatia.html). 

## Citation
If you found Hypatia useful in your research, please cite our [preprint](https://www.biorxiv.org/content/10.64898/2026.01.13.699341v2):
```
Pan, T., Shiau, C. K., Lu, L., Wang, C., Wang, M., He, Y., Bhimaraj, A., Brat, D., Huse, J.,
Jessica Li, J., & Gao, R. (2026). Hypatia: Comparative Isoform Profiling Across Cell Populations
 from Long-Read Single-Cell Transcriptomes. bioRxiv : the preprint server for biology,
 2026.01.13.699341. https://doi.org/10.64898/2026.01.13.699341```
