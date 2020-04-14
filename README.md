> Part of: Liu, X. _et al._ (2020) [title, journal]

# CRISPRi-seq target evaluation alpha
- The file **function_sgRNAefficiency.R** contains a function to find and score all binding sites of a given list of sgRNAs in a given genome.
- The file **library-efficiency.R** contains a script that performs the complete pipeline of target identification and library efficiency analysis. Outputs an Excel table. Calls **function_sgRNAefficiency.R**. 
- The file **library-efficiency-pneumococcal-strains.pdf** contains an analysis of the results of **library-efficiency.R** for seven penumococcal genomes, given the CRISPRi library of this publication. Source code of the pdf: **library-efficiency-pneumococcal-strains.Rmd**.

Requires R packages [CRISPRseek](https://bioconductor.org/packages/release/bioc/html/CRISPRseek.html) & [biomartr](https://cran.r-project.org/package=biomartr).

# References
- Drost HG & Paszkowski J (2017) Biomartr: Genomic data retrieval with R. Bioinformatics 33: 1216–1217  
- R Core Team (2018) R: A language and environment for statistical computing. Available at: https://www.r-project.org/
- Zhu LJ, Holmes BR, Aronin N & Brodsky MH (2014) CRISPRseek: A Bioconductor package to identify target-specific guide RNAs for CRISPR-Cas9 genome-editing systems. PLoS One 9
