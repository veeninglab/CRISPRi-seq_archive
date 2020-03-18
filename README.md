> Part of: Liu, X. _et al._ (2020) [title, journal]

# CRISPRi-seq off-target evaluation alpha
- the file **function_specificity_score.R** finds all potential off-target sites and gives them a specificity score
- the file **offTarget_effect_scoring.R** calls on **function_specificity_score.R** for the _Streptococcus pneumoniae_ D39V genome and our CRISPRi-seq whole-genome sgRNA library

Requires R packages [CRISPRseek](https://bioconductor.org/packages/release/bioc/html/CRISPRseek.html) & [biomartr](https://cran.r-project.org/package=biomartr).

# References
- Drost HG & Paszkowski J (2017) Biomartr: Genomic data retrieval with R. Bioinformatics 33: 1216–1217  
- R Core Team (2018) R: A language and environment for statistical computing. Available at: https://www.r-project.org/
- Zhu LJ, Holmes BR, Aronin N & Brodsky MH (2014) CRISPRseek: A Bioconductor package to identify target-specific guide RNAs for CRISPR-Cas9 genome-editing systems. PLoS One 9
