#### Target-site detection and reporting for given sgRNAs in any genome ####
# Vincent de Bakker
# Veening Lab, DMF, FBM, University of Lausanne, Switzerland
# vincent.debakker@unil.ch
# 27 March 2020
####

# Genbank accession numbers:
#   GCA_003003495.1   S. pneumoniae D39V 
#   GCA_000006885.1   S. pneumoniae TIGR4 
#   GCA_000007045.1   S. pneumoniae R6
#   GCA_000019265.1   S. pneumoniae H19A
#   GCA_000019025.1   S. pneumoniae T19F
#   GCA_002813955.1   S. pneumoniae 11A
#   GCA_000019825.1   S. pneumoniae G54

# SETTINGS ####
## User-defined settings ####
sgRNA_file <- "~/sgRNA_seq_all-lib-1499.xlsx"                   # path to file with sgRNAs
accession_nr <- "GCA_000006885.1"                               # Genbank assembly accession of genome to analyze
path_ncbi_downloads <- "~"                                      # directory with / to save _ncbi_downloads folder
file_function_sgRNAefficiency <- "~/function_sgRNAefficiency.R" # path to function file
out_path <- "~"                                                 # directory to store output in

# optional
PAM <- "NGG"
db <- "genbank"   # other options: refseq, ensembl
max_mismatch <- 8 # more takes much longer; not worth it

## Loading packages ####
required_packages_CRAN <- c("BiocManager", "readxl", "writexl")
required_packages_BioC <- c("Biostrings", "CRISPRseek", "biomartr")
for(i in seq.int(required_packages_CRAN)){
  if(!requireNamespace(required_packages_CRAN[i], quietly = TRUE)){install.packages(required_packages_CRAN[i])}
}
for(i in seq.int(required_packages_BioC)){
  if(!requireNamespace(required_packages_BioC[i], quietly = TRUE)){BiocManager::install(required_packages_BioC[i])}
}
lapply(required_packages_CRAN[-1], library, character.only = TRUE)
lapply(required_packages_BioC, library, character.only = TRUE)

# load sgRNAefficiency function
source(file_function_sgRNAefficiency)


# PREPARE DATA ####

## Read sgRNAs ####
sgRNAs <- read_xlsx(sgRNA_file)
colnames(sgRNAs) <- c("tag", "sgRNA")
sgRNAs$sgRNA <- apply(sgRNAs, 1, function(x){paste0(toupper(substr(x[2], 1, 20)), PAM)})

# format for CRISPRseek::searchHits
gRNAs <- DNAStringSet(unlist(sgRNAs[, 2]), use.names = FALSE)
names(gRNAs) <- unlist(sgRNAs[, 1])

## Read genome, download if needed ####
if(file.exists(paste0(path_ncbi_downloads, "_ncbi_downloads/genomes/", accession_nr,"_genomic_", db, ".fna.gz"))){
  genome_path <- paste0(path_ncbi_downloads, "_ncbi_downloads/genomes/", accession_nr,"_genomic_", db, ".fna.gz")
} else{
  genome_path <- getGenome(db = db, 
                           organism = accession_nr, 
                           reference = FALSE)
}
genome <- read_genome(genome_path)

## Read GFF, download if needed ####
if(file.exists(paste0(path_ncbi_downloads, "_ncbi_downloads/annotation/", accession_nr, "_genomic_", db, ".gff.gz"))){
  gffPath <- paste0(path_ncbi_downloads, "_ncbi_downloads/annotation/", accession_nr, "_genomic_", db, ".gff.gz")
} else{
  gffPath <- getGFF(db = db, 
                    organism = accession_nr, 
                    reference = FALSE, 
                    path = paste0(path_ncbi_downloads, "_ncbi_downloads/annotation/"))
}
GFF <- read_gff(gffPath)

# extract all features with a locus tag
genes <- GFF[unlist(lapply(GFF$attribute, grepl, pattern = "locus_tag")), ]

# replace split features with same attributes by first with total range for start and end
#   e.g. TIGR4 SP_0755 is annotated with "joined feature span" 
#       (https://www.ncbi.nlm.nih.gov/genbank/genomesubmit_annotation/)
if(any(duplicated(genes$attribute))){
  # if any, then per unique attribute...
  genes <- do.call(rbind, lapply(split(genes, genes$attribute), function(x){
    # ...find total range of start and end
    tmp_range <- range(x[, c("start", "end")])
    # return only first row of attribute...
    res <- x[1, ]
    # ...with total chromosomal location span
    res$start <- tmp_range[1]
    res$end <- tmp_range[2]
    res
  }))
}

# find duplicates of same feature
genes_tags <- unlist(lapply(genes$attribute, function(x){
  # add ";" at end in case locus_tag is last attribute (regular expression requires ending character)
  sub(paste0(".*?", "locus_tag", "=(.*?);.*"), "\\1", paste0(x, ";"))
}))
# remove duplicates
genes <- genes[!duplicated(genes_tags), ]


# RUN ANALYSIS ####
## for evaluated pneumococcal genomes, took ~1 to 1.5 hours
effic <- sgRNAefficiency(sgRNAs = gRNAs, genes = genes, 
                         reprAct = TRUE, dist2SC = TRUE, 
                         name_by = "locus_tag", 
                         penalties = "qi.mean.per.region", 
                         seqs = unlist(genome), seqname = accession_nr, 
                         outfile = paste0(out_path, accession_nr), 
                         max.mismatch = max_mismatch, 
                         PAM = PAM, allowed.mismatch.PAM = 1, PAM.pattern = paste0(PAM, "$"))

## Save results ####
#write.csv(effic, paste0(out_path, "lib-1499_", accession_nr, ".csv"))
write_xlsx(effic, paste0(out_path, "lib-1499_", accession_nr, ".xlsx"))