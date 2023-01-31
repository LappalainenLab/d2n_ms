##### D2N - Get dosage genes TSS #####
# Use K562 5' scRNA-seq 10X data to get the TSS coordinates of the 4 dosage genes
# (1) Get the exonic coordinates for genes of interest 
# (2) Get mapped reads in these exon regions from ECCITE-seq BAM file
# (3) Per gene, calculate the read coverage per sliding window and return all putative TSS (local maximas)

# JDE, November 2021
# Last modified: January 2023


## Libraries
library(data.table)
library(Rsamtools)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(biomaRt)
library(rtracklayer)
library(BiocParallel)

## Dir
setwd("/gpfs/commons/groups/lappalainen_lab/jdomingo/manuscripts/d2n/d2n_ms/01-library_design/")
folder_name <- basename(getwd())
data_dir = file.path("../../data", folder_name)
processed_data_dir = file.path("../../processed_data", folder_name)


## Functions
# Function to infer TSS sites for one gene (inspired by https://github.com/argschwind/TAPseq)
TSS_site_gene <- function(gene, coverage, TSS_upstream, wdsize, by, extend_upstream, perc_threshold) {
  
  # process annotation and coverage data -----------------------------------------------------------
  
  # get chromosome and strand of gene
  chr <- as.character(unique(seqnames(gene)))
  strand <- as.character(unique(strand(gene)))
  
  # merge exons of gene
  exons <- reduce(gene)
  
  # find putative TSS site if gene is on the '-' strand
  
  if (strand == "-") {
    
    # extend last exon downstream because TSS on minus strand, so plus the end of "last exon"
    ######()here!!!!!!!end of exo
    #extend_upstreamdownstream when estimating TSS sites (default: 0). A reasonable value (e.g. 100-200 bp)
    end(exons[length(exons)]) <- end(exons[length(exons)]) + extend_upstream
    
    # exon ends and starts in spliced (concatenated) transcript
    exon_ends <- cumsum(width(exons))
    exon_starts <- c(1, exon_ends[-length(exon_ends)] + 1)
    tx_exons <- IRanges(exon_starts, exon_ends)
    
    # calculate coverage in sliding windows along gene ---------------------------------------------
    
    # coverage of gene, exons contains the genomic coordinates
    gene_cvrg <- Views(coverage[[chr]], ranges(exons))
    
    # coverage of spliced transcript (concatenated genomic coverage)
    tx_cvrg <- unlist(gene_cvrg)
    
    # define windows based on strand of gene, so that the 3' end is for sure, default wdsize=200, by=1
    # covered by a window, default wdsize=200
    wd_ends <- seq(from = sum(width(exons)), to = wdsize, by = -by)#diff, from higher # to lower #
    wd_starts <- wd_ends - wdsize + 1
    wds <- IRanges(wd_starts, wd_ends)
    
    # calculate coverage in each window, kinda smoothed 
    wd_cvrg <- sum(Views(tx_cvrg, wds))
    
    # define putative TSS site based on local maxima (peaks) -------------------------------------
    
    # calculate local maxima in smoothed coverage
    peaks <- local_maxima(wd_cvrg)
    peaks_cvrg <- wd_cvrg[peaks]
    
    # retain peaks which coverage is >= perc_threshold, default perc_threshold=0.9
    min_cvrg <- as.numeric(stats::quantile(wd_cvrg, probs = perc_threshold))
    peaks <- peaks[peaks_cvrg >= min_cvrg]
    
    # get windows for these peaks
    peak_wds <- wds[peaks]
    
    # get center of peaks
    peak_centers <- round(start(peak_wds) + wdsize / 2 - 1)
    
    # define putative TSS site： TSS_downstream (numeric) How far downstream of a peak in coverage are TSS sites
    #'   expected? Somewhat depends on input DNA fragment size. (default: 100)
    TSS_sites <- peak_centers + TSS_upstream
    
    # if gene is on '-' strand
  } else if (strand == "+") {
    
    # Define the starts and ends of the gene's exons (elongating first exon for upstream bp)
    # this is different by strand
    start(exons[1]) <- start(exons[1]) - extend_upstream # add upstream start of exon to calculate coverage 
    # exon ends and starts in spliced (concatenated) transcript
    exon_ends <- cumsum(width(exons))
    exon_starts <- c(1, exon_ends[-length(exon_ends)] + 1)
    tx_exons <- IRanges(exon_starts, exon_ends)
    
    # calculate coverage in sliding windows along gene ---------------------------------------------
    
    gene_cvrg <- Views(coverage[[chr]], ranges(exons))
    tx_cvrg <- unlist(gene_cvrg)
    
    wd_starts <- seq(from = 1, to = (sum(width(exons)) - wdsize + 1), by = by) # from lower to higher
    wd_ends <- wd_starts + wdsize - 1
    wds <-IRanges(wd_starts, wd_ends)
    
    wd_cvrg <- sum(Views(tx_cvrg, wds))
    
    # define putative TSS site based on local maxima (peaks) -------------------------------------
    
    peaks <- local_maxima(wd_cvrg)
    peaks_cvrg <- wd_cvrg[peaks]
    
    min_cvrg <- as.numeric(stats::quantile(wd_cvrg, probs = perc_threshold))
    peaks <- peaks[peaks_cvrg >= min_cvrg]
    
    peak_wds <- wds[peaks]
    
    peak_centers <- round(start(peak_wds) + wdsize / 2 - 1)
    
    TSS_sites <- peak_centers - TSS_upstream + 1
    
  } else {
    
    stop("Strand not '+' or '-' for at least 1 gene!", call. = FALSE)
    
  }
  
  # translate TSS sites to genomic coordinates
  if (length(TSS_sites) > 0) {
    
    # transform to IRanges object
    TSS_sites <- IRanges(start = TSS_sites, end = TSS_sites)
    
    # add metadata
    cov_bp <- wd_cvrg[peaks] / wdsize # coverage per bp of peaks
    mcols(TSS_sites) <- data.frame("name" = "TSS_site",
                                   "score" = cov_bp)
    
    # sort and remove possible duplicate TSS sites due to rounding of coordinates
    TSS_sites <- unique(sort(TSS_sites))
    
    # get exon overlapping with TSS sites
    TSS_overlaps <- findOverlaps(query = tx_exons, subject = TSS_sites)
    TSS_exons <- S4Vectors::queryHits(TSS_overlaps)
    
    # only retain TSS sites that overlap any exons
    TSS_sites <- TSS_sites[S4Vectors::subjectHits(TSS_overlaps)]
    
    # calculate distance of TSS site to the start of the overlapping exon
    pos_diff <- start(TSS_sites) - start(tx_exons[TSS_exons])
    
    # translate TSS site coordinates to genomic coordinates
    TSS_start <- start(exons[TSS_exons]) + pos_diff
    
    # create GRanges object with genomic coordinates of TSS site
    TSS_coords <- GRanges(seqnames = chr, strand = strand,
                          ranges = IRanges(TSS_start, TSS_start),
                          mcols(TSS_sites))
    
  } else {
    
    # if no TSS sites are found, return empty GRanges object
    TSS_coords <- GRanges()
    
  }
  
  # return output
  return(TSS_coords)
 
}

## FUN: find local maxima in an ordered data series (after post in:
#  https://stackoverflow.com/questions/6836409/finding-local-maxima-and-minima)
local_maxima <- function(x) {
  # Use -Inf instead if x is numeric (non-integer)
  y <- diff(c(-1, x)) > 0L
  # y <- diff(c(-.Machine$integer.max, x)) > 0L
  rle(y)$lengths
  y <- cumsum(rle(y)$lengths)
  y <- y[seq.int(1L, length(y), 2L)]
  if (x[[1]] == x[[2]]) {
    y <- y[-1]
  }
  # return vector with indices of local maxima
  return(y)
}


# Initial variables
input_gene_list =  c("GFI1B", "NFE2", "MYB", "TET2") # Dosage genes



## (1) Get the exons genomic coordinates for genes of interest 

# Convert SYMBOL to Enterez IDs
ensembl <- biomaRt::useEnsembl(
  biomart = "ensembl", 
  dataset = "hsapiens_gene_ensembl",
  host = "https://www.ensembl.org"
)
symbols_entrezid <- as.data.table(select(org.Hs.eg.db, 
                                         keys = input_gene_list, 
                                         columns = c("ENTREZID", "SYMBOL"), 
                                         keytype = "SYMBOL"))

# Get genomic coordinates of all known genes
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# Keep only the coordinates of exons, splitting into different GRanges for each individual gene (gene names are Entrez IDs)
exons_by_gene <- exonsBy(txdb, by = "gene")

# Subset only for the genes of interest
glist_exon <- exons_by_gene[symbols_entrezid$ENTREZID]

# Use only chromosome levels that are in use
seqlevels(glist_exon) <- seqlevelsInUse(glist_exon) # Avoid unusual chromosome levels so that is compatible with BAM file header




# (2) Get mapped reads in these exon regions from ECCITE-seq BAM file

# CellRanger BAM file 
bam_file <- file.path(data_dir, "K562_WT_scRNAseq5cap_possorted_genome_bam.bam")

# ScanBamParam creates a parameter object influencing what fields and which records are imported from a BAM file
gene_range <- ScanBamParam(which=unlist(glist_exon))

# Load reads for only those existing granges in the grange list
reads <- GenomicAlignments::readGAlignments(bam_file, param=gene_range)

# Calculate coverage per base
read_cvrg <- GenomicAlignments::coverage(reads)



# (3) Per gene, calculate the read coverage per sliding window and return all putative TSS (local maximas)

# Register backend for parallelization
register(MulticoreParam(workers = 4))

# Set variable values to infer TSS sites
TSS_upstream = 100
wdsize = 100
by = 1
extend_upstream = 0
perc_threshold = 0.8
parallel = TRUE
min_cvrg = 0


# Find TSS iterating over each individual gene
if (parallel == TRUE) {
  
  TSS_sites <- BiocParallel::bplapply(X = glist_exon, FUN = TSS_site_gene, coverage = read_cvrg,
                                      TSS_upstream = TSS_upstream, wdsize = wdsize, by = by,
                                      extend_upstream = extend_upstream, perc_threshold = perc_threshold)
} else {
  
  TSS_sites <- lapply(X = glist_exon, FUN = TSS_site_gene, coverage = read_cvrg,
                      TSS_upstream = TSS_upstream, wdsize = wdsize, by = by,
                      extend_upstream = extend_upstream, perc_threshold = perc_threshold)
}


# Transform output to GRanges object
TSS_sites <- unlist(GRangesList(TSS_sites))

# filter TSS sites according to min_cvrg
TSS_sites <- TSS_sites[TSS_sites$score >= min_cvrg]

# sort according to position
TSS_sites <- sort(TSS_sites)

# Add gene names to TSS
TSS_sites$gene_enterez <- names(TSS_sites)
TSS_sites$gene_name <- sapply(TSS_sites$gene_enterez, function(x) {symbols_entrezid[ENTREZID == x, SYMBOL]})



# Export TSS bed and exons gtf to visualize with IGV
export(TSS_sites, con = "processed_data/titration_genes_TSS_sites.bed", format = "bed")


# (3) Given stand, find the first TSS of the transcript
CollapseTSS <- function(gr) {
  if (unique(strand(gr)) == "+") {
    tss <- gr[which(min(start(gr)) == start(gr))]
  } else if (unique(strand(gr)) == "-") {
    tss <- gr[which(max(start(gr)) == start(gr))]
  }
}

# Concatenate all TSS into a single GRanges object
Unique_TSS_sites <- unlist(as(BiocParallel::bplapply(X = split(TSS_sites, f=TSS_sites$gene_name), FUN = CollapseTSS), "GRangesList"))
Unique_TSS_sites$name <- paste0(Unique_TSS_sites$gene_name, "_TSS")

# Export TSS bed 
export(Unique_TSS_sites, con = file.path(processed_data_dir, "D2N_TSS_dosage_genes.bed"), format = "bed")







