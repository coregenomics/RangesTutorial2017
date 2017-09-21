## ----Vignette setup, echo = FALSE----------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")

## ----Tutorial installation, eval = FALSE---------------------------------
#  # Change to vignettes directory to additional files.
#  setwd("vignettes/")
#  # Opens the vignette in your editor window.
#  file.edit("granges-intro.Rmd")
#  
#  # Install the packages we need to run the tutorial.
#  #
#  # Update bioconductor installer.
#  source("https://bioconductor.org/biocLite.R")
#  # Devtools makes it convenient to work with packages.
#  if (!requireNamespace("devtools", quietly = TRUE))
#    install.packages("devtools")
#  # If you already have devtools, update if not newest.
#  devtools::install_cran("devtools")
#  # Install everything required for this tutorial.
#  devtools::install(repos = biocinstallRepos(), dependencies = TRUE)

## ----R BED file raw read-------------------------------------------------
file_bed <- system.file(package = "rtracklayer", "tests", "test.bed")
df <- read.delim(file_bed, header = FALSE, skip = 2)
df <- df[, 1:6]
colnames(df) <- c("chr", "start", "end", "name", "score", "strand")
df

## ----R BED file import---------------------------------------------------
suppressPackageStartupMessages(library(rtracklayer))

gr <- import(file_bed)
gr

GRanges("chr11:100-1000")
# What happens when you run this?
try(GRanges("chr11:1000-100"))

## ----Equivalent R workflow-----------------------------------------------
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(GenomicAlignments)
  library(rtracklayer)   # import
  library(OrganismDbi)   # makeTxDbFromGFF, genes
})

file_annotations <- system.file(package = "Rsamtools",
                                "extdata", "example.gtf.gz")
file_reads <- system.file(package = "Rsamtools",
                          "extdata", "ex1.bam")

# Anntations ----------------------------------------------------------

# Less optimal way of extracting genes from annotations.
annox <- import(file_annotations)
unique(mcols(annox)$type)
genes <- annox[mcols(annox)$type == "gene"]

# Better way.
txdb <- makeTxDbFromGFF(file_annotations, organism = "Homo sapiens")
txdb
genes <- genes(txdb)
genes

# Read counts ---------------------------------------------------------
reads <- readGAlignments(file_reads)
# Cobber example data.  Never do this :)
levels <- c("chr1", "chr2")
reads@seqinfo@seqnames <- levels
reads@seqnames@values <- factor(levels)

# Find reads in genes.
counts <- countOverlaps(genes, reads)
counts

# Smooth reads using windows.
genes_windows <- slidingWindows(genes, width = 250, step = 50)
counts_windowed <- relist(
  countOverlaps(unlist(genes_windows), reads),
  genes_windows)
names(counts_windowed) <- names(genes)
counts_windowed

## ----Magic numbers, eval = FALSE-----------------------------------------
#  genes_trimmed  # Good - uses underscores to separate words
#  genesTrimmed   # Bad - uses mixedCase
#  genes.trimmed  # Bad - uses dots

## ----Reverse verb convention, eval = FALSE-------------------------------
#  promoters
#  promoters_windowed  # Good - variant at end
#  win_promoters       # Bad - "win" is not clear
#  windows             # Bad - relation not clear

## ----Test reformatting, eval = FALSE-------------------------------------
#  # *** Reflow the long comment line below with Ctrl + Shift + / ***
#  #
#  # Read multiple BAM alignment files and simplify the GAlignments into a GRangesList object.
#  read_bams<-function( file_names ){
#  as(List(lapply( file_names, GenomicAlignments::readGAlignments )),
#  "GRangesList")
#  }
#  
#  # *** Select the above code and reindent with Ctrl + I ***
#  # *** Select the above code and reformat with Ctrl + Shift + A ***

