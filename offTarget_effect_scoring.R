#### off-target detection, scoring and reporting for given sgRNA's ####
# Vincent de Bakker©
####

#### SETTINGS; ONLY USER-DEFINED SETTINGS IN SCRIPT ####

# mandatory settings
wd <- "D:/Projects/CRISPRi/sgRNA_design/"
sgRNA_file <- "/sgRNA_seq_all-lib-1499.xlsx" 
max_mismatch <- 8

# optional settings
PAM <- "NGG" 
name_by <- "locus_tag"
accession_nr <- "GCA_003003495.1" # D39V gb accession nr
db <- "genbank"
first_run <- FALSE # set to FALSE after first run to save time
path_ncbi_downloads <- wd # where the genome will be saved, defaults to wd
path_sgRNA_file <- wd # where the sgRNA_file is saved, defaults to wd


#### PRELIMINARIES ####

# load packages (and automatically installs if necessary)
#biomartrURL <- "https://cran.r-project.org/src/contrib/Archive/biomartr/biomartr_0.8.0.tar.gz" # archived
required_packages_CRAN <- c("BiocManager", "readxl", "writexl", "dplyr", "digest")
required_packages_BioC <- c("Biostrings", "CRISPRseek", "biomartr")
for(i in seq.int(required_packages_CRAN)){
  if(!requireNamespace(required_packages_CRAN[i], quietly = TRUE)){install.packages(required_packages_CRAN[i])}
}
for(i in seq.int(required_packages_BioC)){
  if(!requireNamespace(required_packages_BioC[i], quietly = TRUE)){BiocManager::install(required_packages_BioC[i])}
}
#if(!requireNamespace("biomartr", quietly = TRUE)){install.packages(biomartrURL, repos = NULL, type = "source")}
library(Biostrings)
library(biomartr)
library(CRISPRseek)
library(readxl)
library(writexl)

# load user-defined specificity function from file
source(paste0(wd, "function_specificity_score3.R"))

# set working directory
setwd(wd)


#### PREPARE INPUT ####

# Read genome into R
if(first_run){
  genome_path <- getGenome(db = db, 
                           organism = accession_nr, 
                           reference = FALSE)
} 
genome_path <- paste0(path_ncbi_downloads, "_ncbi_downloads/genomes/", accession_nr,"_genomic_", db, ".fna.gz")
genome <- read_genome(genome_path)

# read sgRNA file
sgRNAs <- read_xlsx(paste0(path_sgRNA_file, sgRNA_file))
colnames(sgRNAs) <- c("tag", "sgRNA")
sgRNAs$sgRNA <- apply(sgRNAs, 1, function(x){paste0(toupper(substr(x[2], 1, 20)), PAM)})


#### RUN ANALYSIS ####

# get scores
start <- Sys.time()
SH <- specificity3(sgRNAs = sgRNAs, 
                   return.offTargets = TRUE, 
                   seqs = unlist(genome), seqname = accession_nr, 
                   outfile = paste0(wd, accession_nr), 
                   max.mismatch = max_mismatch, 
                   PAM = PAM, allowed.mismatch.PAM = 1, PAM.pattern = paste0(PAM, "$"))
end <- Sys.time()
end - start
# all 1499 sgRNA library Xue against D39V genome, return.offTargets=TRUE
# Time difference of 1.081435 hours


#### EXPLORE RESULTS IN R ####

# now SH is a list with two elements
length(SH)
# namely 
names(SH)
# so the first element contains the specificity scores with the original sgRNA names
SH$specificity
# that we can put in a data frame for convinience:
df <- data.frame("tag" = sgRNAs$tag, 
                 "sgRNA" = substr(sgRNAs$sgRNA, 1, 20), # here I remove the PAM from the sequence
                 "specificity" = SH$specificity[match(names(SH$specificity), sgRNAs$tag)])
# where I made sure the names of the scores match the correct original sequences
df # check if you like 
# and the second the full off-target hits data frame
summary(SH$offTargetHits)
head(SH$offTargetHits)
# which contains: 
# - for every position whether there is a mismatch (1) or not (0)
#   between the sgRNA and its target
# - the position of the hit on the chromosome: 
#   strand, start and end  bp corresponding to pneumobrowse / given sequence file
# - the original sgRNA name and sequence + attached PAM as in sgRNA_file
# - the number of mismatches (sum of columns 1:20)
# - some obscure score the deeper package originally uses, but we don't

# check e.g. how many times each nr of mismatches occured per sgRNA
tapply(SH$offTargetHits$n.mismatch, SH$offTargetHits$name, table)
## note that in all cases you will get at least 1 target w/ 0 mismatches: this is the on-target
# or show per sgRNA the hits with the least mismatches
min.misms <- lapply(split(SH$offTargetHits, SH$offTargetHits$name), function(x){
  subset(x, x$n.mismatch == min(x$n.mismatch))
})
min.misms$sgRNA0899
## indeed multiple targets w/ 0 mismatches for this sgRNA, hence its 0 specificity score

# you can also view the table in RStudio by running 
#View(SH$offTargetHits)

# or view the table in Excel by writing it to a file, shown below


#### WRITING RESULTS TO FILE ####

## to save the specificity scores, you can run e.g.:
# CSV
write.csv(df, file = paste0(wd, accession_nr, "_specificity_maxmismatch", max_mismatch, ".csv"), row.names = FALSE)
# Excel
write_xlsx(df, path = paste0(wd, accession_nr, "_specificity_maxmismatch", max_mismatch, ".xlsx"))

## to save the off-target hits table, run e.g.:
# CSV
write.csv(SH$offTargetHits, file = paste0(wd, accession_nr, "_offTargetHits_maxmismatch", max_mismatch, ".csv"), row.names = FALSE)
# Excel
write_xlsx(SH$offTargetHits, path = paste0(wd, accession_nr, "_offTargetHits_maxmismatch", max_mismatch, ".xlsx"))