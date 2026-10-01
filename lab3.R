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
library("Biostrings")

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

# Do in terminal

# Ensure the output directory exists
## mkdir -p outputs/flye_output

# wget -nc -P data/ ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR153/006/ERR1539006/ERR1539006.fastq.gz

# Run the Flye assembler, note the expected genome size
# flye --nano-raw data/ERR1539006.fastq.gz --out-dir outputs/flye_output --threads 8 --genome-size 12m


# install medaka
## conda create -n medaka_final -c conda-forge -c bioconda -c nanoporetech medaka python=3.10 -y

# Activate your environment
## conda activate medaka_final

# Ensure directories exist
## mkdir -p outputs/medaka_round1 outputs/medaka_round2

# Round 1
# medaka_consensus -i data/ERR1539006.fastq.gz -d outputs/flye_output/assembly.fasta -o outputs/medaka_round1 -t 8

# Round 2
# medaka_consensus -i data/ERR1539006.fastq.gz -d outputs/medaka_round1/consensus.fasta -o outputs/medaka_round2 -t 8


# Do in R
# Visualizing Alignment
# BAM file path
bam_file <- "outputs/alignment.sorted.bam"

# Open connection to BAM file
bam <- BamFile(bam_file)

# Read all alignments
alns <- readGAlignments(bam)

# Convert to GRanges object
gr_alns <- granges(alns)

# Find the start and end boundaries of all mapped contigs combined
range(gr_alns)

# Peek at the individual contig boundaries
ranges(gr_alns)

# create track showing alignment
aln_track <- AnnotationTrack(gr_alns, name = "Contigs", genome = "sacCer3", chromosome = "chrI")

# Add genome axis for scale
axis_track <- GenomeAxisTrack()

# Plot
plotTracks(list(axis_track, aln_track),
           from = min(start(gr_alns)),
           to = max(end(gr_alns)),
           main = "Contig Alignments to Chromosome 1")


# code to find # of contigs aligning to chromosome 1
width(gr_alns)
max(width(gr_alns))

# Read contigs
readContigs <- readDNAStringSet("outputs/medaka_round2/consensus.fasta")

# Find longest contig
longest <- which.max(width(readContigs))

# Extract longest contig 
longest_contig <- readContigs[longest]

# write to new fasta file
writeXStringSet(longest_contig, "outputs/longest_contig.fasta")


## Do this is in R terminal 

# gunzip -k data/GCF_000146045.2_R64_cds_from_genomic.fna.gz
# minimap2 -ax splice -uf --secondary=no outputs/longest_contig.fasta data/GCF_000146045.2_R64_cds_from_genomic.fna > outputs/cds_aln.sam

## code chunk that converts your cds_aln.sam file to a sorted and indexed BAM file for visualization
# samtools view -bS outputs/cds_aln.sam > outputs/cds_alignment.bam
# samtools sort outputs/cds_alignment.bam -o outputs/cds_aln.sorted.bam
# samtools index outputs/cds_aln.sorted.bam

## Using your sorted BAM file, create Gviz tracks that represent the aligned CDSs

# BAM file path
bam_file1 <- "outputs/cds_aln.sorted.bam"

# Open connection to BAM file
bam1 <- BamFile(bam_file1)

# Read all alignments
alns1 <- readGAlignments(bam1)

# Convert to GRanges object
gr_alns1 <- granges(alns1)

# Find the start and end boundaries of all mapped contigs combined
range(gr_alns1)

# set to false
options(ucscChromosomeNames=FALSE)

# Create tracks showing alignment
aln_track1 <- AlignmentsTrack("outputs/cds_aln.sorted.bam", isPaired = FALSE)

# Add genome axis for scale
axis_track1 <- GenomeAxisTrack()

# Visualize the alignment
# Plot
plotTracks(list(axis_track1, aln_track1),
           from = min(start(gr_alns1)),
           to = max(end(gr_alns1)),
           chromosome = "contig_301",
           main = "CDSs Alignments to Chromosome 1")

# reset to true
options(ucscChromosomeNames=TRUE)

