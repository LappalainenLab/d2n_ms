##### D2N - Design tiling sgRNA library #####
# (0) Create non-overlapping genomic regions where to design sgRNAs
#     - Load identified TSS from 10X BAM files
#     - Define fist region as 50 bp around the TSS
#     - Define upstream and downstream regions where to design guides
# (1) Get all possible guide sequences from those defined regions
# (2) Filter out guides with:
#     - Repeats 
#     - U6 terminators 
#     - Duplicated spacers 
#     - EcoRI restriction sites (full and star activity) 
# (3) Generate bed file with Hg19 coordinates to intersect with K562 phased genome vcf (bedtools) to discard guides that fall in K562 specific SNVs or SNPs
# (4) Generate data table for running FlashFry
# (5) Select best FlashFry scoring guide within each region
# (6) Make sure distance between guides is adequate
# (7) Overlap 'weismanTSS' guides and substitute close or overlapping ones for new extra tiling guides


# JDE, January 2022
# Last modified: February 2023

## Libs
library(data.table)
library(GenomicRanges)
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg38)

## Dirs
setwd("/gpfs/commons/groups/lappalainen_lab/jdomingo/manuscripts/d2n/d2n_ms/01-library_design/")
folder_name <- basename(getwd())
data_dir = file.path("../../data", folder_name)
processed_data_dir = file.path("../../processed_data", folder_name)
plots_dir = file.path("../../plots", folder_name)


## Functions
DefineSingleUpDownGenomicRegions <- function(gene, chr, strand, tss, d_tss, d2, d3, num_up1, num_up2, num_down, PAM_pattern = "NGG") {
  
  # Positive strand
  if (strand == "+") {
    
    # upstream regions 
    s_u = tss-d_tss-1-(d2*num_up1)-(d3*num_up2) # start
    e_u = tss-d_tss-1 # end
    
    # downstream regions
    s_d = tss+d_tss+1
    e_d = tss+d_tss+1+(d2*num_down)
    
    starts = c(s_u, tss-d_tss, s_d)
    ends = c(e_u, tss+d_tss, e_d) 
    
    region = c("up", "tss", "down")
    
    
  } else if (strand == "-") {
    # downstream regions 
    s_d = tss-d_tss-1-(d2*num_down) # start
    e_d = tss-d_tss-1 # ends
    
    # upstream regions
    s_u = tss+d_tss+1
    e_u = tss+d_tss+1+(d2*num_up1)+(d3*num_up2)
    
    starts = c(s_d, tss-d_tss, s_u)
    ends = c(e_d, tss+d_tss, e_u) 
    
    region = c("down", "tss", "up")
    
  } else {
    return(NULL)
  }
  
  # Add features in a single table
  dt <- data.table(gene = gene,
                   chr = chr,
                   strand = strand,
                   start = starts,
                   end = ends,
                   region = region,
                   id = paste0(gene, "_", region))
  
  # Count if there are PAM seqs in each region
  dt[, n_PAMs := vcountPattern(DNAString(PAM_pattern), getSeq(Hsapiens, names = chr, start = start, end = end),fixed = FALSE) +
       vcountPattern(reverseComplement(DNAString(PAM_pattern)), getSeq(Hsapiens, names = chr, start = start, end = end), fixed = FALSE)]
  
  return(dt)
}



## (0) Create non-overlapping genomic regions where to design sgRNAs

# Load data
TSS_dt <- fread(input = file.path(processed_data_dir, "D2N_TSS_dosage_genes.bed"), col.names = c("chr", "TSS_start", "TSS_end", "gene_name", "cvrg_score", "strand"))

# Define variables
# diameter bp dist from TSS to map gRNA in the TSS
d1 = 25
# distance in bp to bin the upstream or downstream regions that are closer to the TSS
d2 = 100
# distance in bp to bin the distal part of the promoter upstream of the TSS
d3 = 200
# num of closest upstream TSS promoter bins
p_up1 = 5
# num of distal upstream TSS promoter bins
p_up2 = 2
# num of closest downstream TSS promoter bins
p_dw = 5


# Define single upstream, TSS and downstream regions per dosage gene
GR_dt <- do.call("rbind", apply(TSS_dt, 1, function(x){
  
  DefineSingleUpDownGenomicRegions(gene = gsub("_TSS", "", x["gene_name"]), 
                                   chr = x["chr"], 
                                   strand = x["strand"], 
                                   tss = as.numeric(x["TSS_start"]), 
                                   d_tss = d1, 
                                   d2 = d2, 
                                   d3 = d3, 
                                   num_up1 = p_up1, 
                                   num_up2 = p_up2,
                                   num_down = p_dw)
}))

# Define sub regions to make several upstream and downstream guides
fwd = c(rep(d3, p_up2), rep(d2, p_up1), d1*2+1, rep(d2, p_dw-1), d2+1)
names(fwd) = c(rep("up", p_up2), rep("up", p_up1), "tss", rep("down", p_dw))
rev = rev(fwd)

sub_GR_l <- list()
for (g in GR_dt[, gene]){
  dt <- GR_dt[gene == g,]
  if (unique(dt$strand == "+")) {
    sub_GR_l[[g]] <- data.table(gene = unique(dt$gene),
                                chr = unique(dt$chr),
                                strand = unique(dt$strand),
                                start = c(min(dt$start), min(dt$start)+1 + cumsum(fwd))[1:length(fwd)],
                                end = c(min(dt$start) + cumsum(fwd)),
                                region = names(fwd))
  } else {
    sub_GR_l[[g]] <- data.table(gene = unique(dt$gene),
                                chr = unique(dt$chr),
                                strand = unique(dt$strand),
                                start = c(min(dt$start), min(dt$start) + cumsum(rev))[1:length(rev)],
                                end = c(min(dt$start)-1 + cumsum(rev)[1:length(rev)-1], max(dt$end)),
                                region = names(rev))
    
  }
}
sub_GR_dt <- rbindlist(sub_GR_l)
sub_GR_dt[, region_id := paste0(gene, "_", region, "_", 1:.N), .(gene, region)]

