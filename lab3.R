# Lab3 - Genome Assemby and Annotation Basics

# Install BiocManager if you don't have it yet
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Install all the requested packages together
BiocManager::install(c("IRanges", "GenomicRanges", "Rsamtools", "GenomicAlignments", "Gviz"))

# load the following package
library("IRanges")
library("GenomicRanges")
library("Rsamtools")
library("GenomicAlignments")
library("Gviz")

# Download the reference genome from refseq
# wget -nc -P data/ https://ftp.ncbi.nlm.nih.gov/genomes/refseq/fungi/Saccharomyces_cerevisiae/latest_assembly_versions/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.fna.gz

# Download the genome annotation file from refseq
# wget -nc -P data/ https://ftp.ncbi.nlm.nih.gov/genomes/refseq/fungi/Saccharomyces_cerevisiae/latest_assembly_versions/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.gff.gz

# Download the cds containing fasta file from refseq
# wget -nc -P data/ https://ftp.ncbi.nlm.nih.gov/genomes/refseq/fungi/Saccharomyces_cerevisiae/latest_assembly_versions/GCF_000146045.2_R64/GCF_000146045.2_R64_cds_from_genomic.fna.gz

# 1. Uncompress and Read the raw NCBI file as simple text
con <- gzfile("data/GCF_000146045.2_R64_genomic.fna.gz", "rt")
ygen <- readLines(con)
close(con)  

# 2. Define the new UCSC-style headers (chrI through chrXVI + chrM)
new_headers <- c(paste0(">chr", as.roman(1:16)), ">chrM")

# 3. Find and Replace
header_lines <- grep(">", ygen)

# Look at the header for chromosome 1 before replacing
print(ygen[header_lines[1]])

# Safety Check: Yeast has 16 nuclear chromosomes + 1 mitochondrial
if(length(header_lines) == 17) {
  ygen[header_lines] <- new_headers
  print("After replacement:")
  print(ygen[header_lines[1]])
  writeLines(ygen, "data/S288C_UCSC.fna")
  message("SUCCESS: Reference renamed to UCSC standards.")
} else {
  stop("Error: Chromosome count mismatch! Check your input file.")
}

# Ensure the output directory exists
# mkdir -p outputs/flye_output

# wget -nc -P data/ ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR153/006/ERR1539006/ERR1539006.fastq.gz

# Run the Flye assembler, note the expected genome size
# flye --nano-raw data/ERR1539006.fastq.gz --out-dir outputs/flye_output --threads 8 --genome-size 12m
