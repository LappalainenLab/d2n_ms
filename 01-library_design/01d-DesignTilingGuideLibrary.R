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
library(tidyverse)

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

# Function to generate all possible sgRNAs with matching PAM sequence within a sequence region
GetAllGuidesInRegion <- function(regionToSearch, PAMpattern) {
  
  GetGuidesFromMatch_TopStrand <- function(matchPattern_Output){
    PAMends = end(matchPattern_Output); indsToKeep = which(PAMends >= 23)
    if (length(indsToKeep) == 0) {return(NULL)}
    PAMends_Keep = PAMends[indsToKeep]
    return(data.frame("GUIDE_START" = (PAMends_Keep - 22),
                      "GUIDE_END" = (PAMends_Keep - 3),
                      "GUIDE_STRAND" = "+",
                      "GUIDE_SEQ" = sapply(PAMends_Keep, function(x){return(matchPattern_Output %>% subject %>% as.character %>% substr(start=(x-22),stop=(x)))}),
                      "PAM" = as.data.frame(matchPattern_Output[indsToKeep,])$seq)) 
  }
  
  GetGuidesFromMatch_BottomStrand <- function(matchPattern_Output){
    PAMends = end(matchPattern_Output)
    indsToKeep = which(PAMends <= ((matchPattern_Output %>% subject %>% length)-20))
    if (length(indsToKeep) == 0) {return(NULL)}; PAMends_Keep = PAMends[indsToKeep]
    return(data.frame("GUIDE_START" = (PAMends_Keep + 20),
                      "GUIDE_END" = (PAMends_Keep + 1),
                      "GUIDE_STRAND" = "-",
                      "GUIDE_SEQ" = sapply(PAMends_Keep, function(x){return(matchPattern_Output %>% subject %>% as.character %>% substr(start=(x-2),stop=(x+20)) %>% DNAString %>% reverseComplement %>% as.character)}),
                      "PAM" = as.data.frame(reverseComplement(matchPattern_Output[indsToKeep,]))$seq))
  }
  
  
  temp <- bind_rows(GetGuidesFromMatch_TopStrand(matchPattern(PAMpattern, regionToSearch, fixed=F)),
                    GetGuidesFromMatch_BottomStrand(matchPattern(reverseComplement(PAMpattern), regionToSearch, fixed=F))) %>%
    mutate(Design_PAM = as.character(PAMpattern))
  
  return(temp)
}




### (0) Create non-overlapping genomic regions where to design sgRNAs

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
GRs_dt <- rbindlist(sub_GR_l)
GRs_dt[, region_id := paste0(gene, "_", region, "_", 1:.N), .(gene, region)]




### (1) Get all possible guides within the defined regions


# Fuse all regions to get a singe sequence region per gene
GR_fused_dt <- GRs_dt[, list("start" = min(start), "end" =  max(end)), .(gene, chr, strand)]

# Get all guides within each region per gene
GDT.1 <- GR_fused_dt[, GetAllGuidesInRegion(regionToSearch = getSeq(Hsapiens, names=chr, start=start, end=end), PAMpattern = DNAString("NGG")), gene]




### (2) Filter guides based on different criteria
spacerVec = GDT.1$GUIDE_SEQ

# A) Remove repeats
indsToRemove <- NULL; indsToRemove = c(indsToRemove, grep("AAAAA",spacerVec)); indsToRemove = c(indsToRemove, grep("CCCCC",spacerVec)); indsToRemove = c(indsToRemove, grep("GGGGG",spacerVec))
indsRepeats = unique(indsToRemove); numWithRepeats = length(indsRepeats)

# B) Remove U6 terminators
indsU6 = grep("TTTT",spacerVec); numWithU6Term = length(indsU6)

# C) Num duplicate spacers
indsDups = which(duplicated(spacerVec)); numDups = length(indsDups)

# D) EcoRI restriction sites (full and star activity) 
EcorRI_main = "GAATTC"
EcorRI_star1 = "TAATTC"
EcorRI_star2 = "CAATTC"

cnst_5 = toupper("caccG")
cnst_3 = toupper("gttta")
spacerVec_re <- paste0(cnst_5, substr(spacerVec, 1, 20), cnst_3)

indsREsites <- NULL; indsREsites = c(indsREsites, grep(EcorRI_main, spacerVec_re)); indsREsites = c(indsREsites, grep(EcorRI_star1, spacerVec_re)); indsREsites = c(indsREsites, grep(EcorRI_star2, spacerVec_re))

# All to remove
uniqueIndsToRemove = sort(unique(c(indsRepeats, indsU6, indsDups, indsREsites)))
numToRemove = length(uniqueIndsToRemove)
numUnfiltered = length(spacerVec)

GDT.2 <- GDT.1[-uniqueIndsToRemove,]





### (3) Discard guides that fall in unique genetic variation of K562 (SNVs and SNPs)

# Get actual genomic coordinates of the guides
GDT.3 <- merge.data.table(GDT.2, setnames(GR_fused_dt[, .(gene, chr, start)], new=c("gene", "chr", "regions_start")), by="gene")
GDT.3[, c("coord_start", "coord_end") := list(regions_start+GUIDE_START-1, regions_start+GUIDE_END-1)] # add actual genomic coordinates
GDT.3[,  coord_pam_end := ifelse(GUIDE_STRAND == "+", coord_end+3, coord_end-3)]

# Generate a GRanges object to liftOver the hg38 coordinates to hg19 (vcf of K562 genome assembly is hg19)
gr_GDT.3 <- makeGRangesFromDataFrame(
  as.data.frame(rbind(GDT.3[GUIDE_STRAND == "+", .(chr, coord_start, coord_pam_end, GUIDE_SEQ)], 
                      setnames(GDT.3[GUIDE_STRAND == "-", .(chr, coord_pam_end, coord_start, GUIDE_SEQ)], c("chr", "coord_start", "coord_pam_end", "GUIDE_SEQ")))),
  keep.extra.columns = T,
  ignore.strand=T,
  start.field = "coord_start",
  end.field = "coord_pam_end")

# The chain file for hg38 to hg19 transformation
path = system.file(package="liftOver", "extdata", "hg38ToHg19.over.chain")
ch = import.chain(path)


# liftOver coordinates to Hg19
seqlevelsStyle(gr_GDT.3) = "UCSC"  # necessary
gr_GDT.3_hg19 = liftOver(gr_GDT.3, ch)
gr_GDT.3_hg19 = unlist(gr_GDT.3_hg19)
genome(gr_GDT.3_hg19) = "hg19"

# Generate bed file with guides coordinates
fwrite(as.data.table(gr_GDT.3_hg19)[, .(seqnames, start, end, GUIDE_SEQ)],
       file = file.path(processed_data_dir, "bedtools/D2N_GDT3_hg19.bed"),
       quote = F, row.names = F, col.names = F, sep="\t")

# Get guides to retain 
GDT.3_hg19_K562 <- fread(file.path(processed_data_dir, "bedtools/D2N_GDT3_hg19_intK562.bed"), 
                         col.names = c("chr", "start", "end", "GUIDE_SEQ"))

# Keep guides that K562 does not differe from the human ref
GDT.4 <- GDT.3[GUIDE_SEQ %in% GDT.3_hg19_K562$GUIDE_SEQ, ]




## (4) Generate table to run FlashFry
# These are the flanking random sequences added
# TTCGTACAAA
# TGTCCGCACT
GDT.5 <- GDT.4
GDT.5[, FlasFryID := paste0(gene, "-", chr, "-", coord_start, "-", coord_end)]
GDT.5[, FlashFrySEQ := paste0("TTCGTACAAA", GUIDE_SEQ, "TGTCCGCACT")]

fwrite(GDT.5, file = "processed_data/guides/D2N_titration_ToScore.txt", quote = F, row.names = F, sep="\t")
GDT.5_flashfry_scored <- fread(file = "processed_data/guides/D2N_titration_Scored.txt")

## Remove from the data table the weismanTSS guides
weismanTSS_scored_dt <- fread(input = "processed_data/guides/D2N_weismanTSS_Scored.txt")
distalCRE_scored_dt <- fread(input = "processed_data/guides/D2N_distalCRE_Scored_Selected.txt")

GDT.5_flashfry_scored[, is_weismanTSS := guide %in% weismanTSS_scored_dt$guide]
GDT.5_flashfry_scored[, is_distalCRE := guide %in% distalCRE_scored_dt$GUIDE_SEQ]

sum(GDT.5_flashfry_scored[, is_distalCRE])
sum(GDT.5_flashfry_scored[, is_weismanTSS])












