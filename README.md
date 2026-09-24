<div align="center">
 
  <img width="120" height="180" alt="Hypatia logo" src="https://github.com/user-attachments/assets/529166f8-c03d-4542-872a-748b64449580" />
  
  # Hypatia
  
  **A statistical framework for comparative isoform profiling across cell populations.**

  [![R ≥ 4.1](https://img.shields.io/badge/R-%E2%89%A54.1-276DC3?style=for-the-badge&logo=r&logoColor=white)](https://www.r-project.org/)
  [![Download](https://img.shields.io/badge/Download-blue?style=for-the-badge&logo=github&logoColor=white)](https://github.com/gaolabtools/Hypatia/releases)
  [![Documentation](https://img.shields.io/badge/Documentation-376B6D?style=for-the-badge)](https://gaolabtools.github.io/Hypatia/vignettes/Hypatia.html)
  [![stars](https://img.shields.io/github/stars/gaolabtools/Hypatia?style=for-the-badge&logo=github&color=FFD700)](https://github.com/gaolabtools/Hypatia/stargazers)

 To get started, please visit [Hypatia's vignette](https://gaolabtools.github.io/Hypatia/vignettes/Hypatia.html). 

</div>


## About

Hypatia (hy-pay-shuh) is a computational toolkit for the investigation of population-specific isoforms from long-read single-cell RNA-sequencing data, featuring three modes of differential analyses:
1) Isoform usage: Measures isoform *shifts* as differential isoform usage (DIU) events.
2) Isoform diversity: Measures isoform *complexity* as differential isoform diversity (DIV) events.
3) Isoform expression: Measures isoform *abundance* as differentially expressed isoforms (DEI) events.

## Installation

Hypatia is an R package available through GitHub and requires R 4.1 or later.

```r
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

devtools::install_github("gaolabtools/Hypatia")
```

<details>
<summary>Troubleshooting: GitHub authentication</summary>

A "Bad credentials" error may indicate an expired or invalid GitHub personal access token (PAT). If the invalid token is stored in your Git credential manager, update it with `gitcreds::gitcreds_set()`, or remove it with `gitcreds::gitcreds_delete()` if it is no longer needed. Then retry the installation.

</details>

## Usage

For detailed documentation and tutorial, please visit [Hypatia's vignette](https://gaolabtools.github.io/Hypatia/vignettes/Hypatia.html). 

## Citation

If you found Hypatia useful in your research, please cite our [preprint](https://www.biorxiv.org/content/10.64898/2026.01.13.699341v2):

Pan, T., Shiau, C. K., Lu, L., Wang, C., Wang, M., He, Y., Bhimaraj, A., Brat, D., Huse, J., Li, J., & Gao, R. (2026). *Hypatia: Comparative Isoform Profiling Across Cell Populations from Long-Read Single-Cell Transcriptomes*. bioRxiv. [doi:10.64898/2026.01.13.699341](https://doi.org/10.64898/2026.01.13.699341).

## License

See [LICENSE](LICENSE) for the terms governing use of Hypatia, including non-commercial research and educational use. For commercial-use inquiries, contact [Dr. Ruli Gao](mailto:ruli.gao@northwestern.edu).

## Support

For questions, bug reports, or feature requests, open a [GitHub issue](https://github.com/gaolabtools/Hypatia/issues). For bug reports, include a minimal reproducible example, the error message, and your `sessionInfo()` output. Do not include private data or access tokens.
