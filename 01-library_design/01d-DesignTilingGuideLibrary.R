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
                      "PAM" = (matchPattern_Output[indsToKeep,] %>% as.data.frame %>% as.matrix)) %>% mutate_if(is.factor, as.character) %>% mutate(PAM=x) %>% select(-x))
    }
  
  GetGuidesFromMatch_BottomStrand <- function(matchPattern_Output){
    PAMends = end(matchPattern_Output)
    indsToKeep = which(PAMends <= ((matchPattern_Output %>% subject %>% length)-20))
    if (length(indsToKeep) == 0) {return(NULL)}; PAMends_Keep = PAMends[indsToKeep]
    return(data.frame("GUIDE_START" = (PAMends_Keep + 20),
                      "GUIDE_END" = (PAMends_Keep + 1),
                      "GUIDE_STRAND" = "-",
                      "GUIDE_SEQ" = sapply(PAMends_Keep, function(x){return(matchPattern_Output %>% subject %>% as.character %>% substr(start=(x-2),stop=(x+20)) %>% DNAString %>% reverseComplement %>% as.character)}),
                      "PAM" = (matchPattern_Output[indsToKeep,] %>% reverseComplement %>% as.data.frame %>% as.matrix)) %>% mutate_if(is.factor, as.character) %>% mutate(PAM=x) %>% select(-x))
    }
  
  temp <- bind_rows(GetGuidesFromMatch_TopStrand(vmatchPattern(PAMpattern, regionToSearch, fixed=F)),
                    GetGuidesFromMatch_BottomStrand(vmatchPattern(reverseComplement(PAMpattern), regionToSearch, fixed=F))) %>%
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



## (2) Filter guides based on different criterias
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
GDT.3 <- merge.data.table(GDT.2, setnames(GR_dt[, .(gene, chr, start)], new=c("gene", "chr", "regions_start")), by="gene")
GDT.3[, c("coord_start", "coord_end") := list(regions_start+GUIDE_START-1, regions_start+GUIDE_END-1)] # add actual genomic coordinates
GDT.3[,  coord_pam_end := ifelse(GUIDE_STRAND == "+", coord_end+3, coord_end-3)]

# Generate a GRanges object to liftover the hg38 cooordinates to hg19 (vcf of K562 genome assembly is hg19)
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
# fwrite(as.data.table(gr_GDT.3_hg19)[, .(seqnames, start, end, GUIDE_SEQ)], 
#        file = "processed_data/bedtools/D2N_GDT3_hg19.bed",
#        quote = F, row.names = F, col.names = F, sep="\t")


# Get guides to retain 
GDT.3_hg19_K562 <- fread("processed_data/bedtools/D2N_GDT3_hg19_intK562.bed", col.names = c("chr", "start", "end", "GUIDE_SEQ"))

# Keep guides that K562 does not differe from the human ref
GDT.4 <- GDT.3[GUIDE_SEQ %in% GDT.3_hg19_K562$GUIDE_SEQ, ]


## (4) Generate table to run FlashFry
# These are the flanking random sequences added
# TTCGTACAAA
# TGTCCGCACT
GDT.5 <- GDT.4
GDT.5[, FlasFryID := paste0(gene, "-", chr, "-", coord_start, "-", coord_end)]
GDT.5[, FlashFrySEQ := paste0("TTCGTACAAA", GUIDE_SEQ, "TGTCCGCACT")]

# fwrite(GDT.5, file = "processed_data/guides/D2N_titration_ToScore.txt", quote = F, row.names = F, sep="\t")
GDT.5_flashfry_scored <- fread(file = "processed_data/guides/D2N_titration_Scored.txt")

## Remove from the data table the weismanTSS guides
weismanTSS_scored_dt <- fread(input = "processed_data/guides/D2N_weismanTSS_Scored.txt")
distalCRE_scored_dt <- fread(input = "processed_data/guides/D2N_distalCRE_Scored_Selected.txt")

GDT.5_flashfry_scored[, is_weismanTSS := guide %in% weismanTSS_scored_dt$guide]
GDT.5_flashfry_scored[, is_distalCRE := guide %in% distalCRE_scored_dt$GUIDE_SEQ]

sum(GDT.5_flashfry_scored[, is_distalCRE])
sum(GDT.5_flashfry_scored[, is_weismanTSS])


## (5) Select best scoring guide within each subregion 

GDT.6 <- merge(GDT.5_flashfry_scored[is_weismanTSS == F & is_distalCRE == F, .(ID, Doench2014OnTarget, DoenchCFD_specificityscore, Hsu2013, otCount)], GDT.5, by.x = "ID", by.y = "FlasFryID")
GDT.6[, guide_no_pam := substr(GUIDE_SEQ, 1, 20)]

ggplot(GDT.6, aes(x=Doench2014OnTarget, fill=gene)) + geom_histogram(bins = 50) 
ggplot(GDT.6, aes(x=Hsu2013, fill=gene)) + geom_histogram(bins = 50) 
ggplot(GDT.6, aes(x=DoenchCFD_specificityscore, fill=gene)) + geom_histogram(bins = 50) 
ggplot(GDT.6, aes(x=otCount, fill=gene)) + geom_histogram(bins = 50) 


# Classify each guide into each genomic subregion
subGR_dt <- fread("processed_data/20220104_titration_genes_guides_regions.txt")
# some guides will fall in the 'margins' between regions so their label becomes 'NA' for now
GDT.6[, region_id := subGR_dt[min(coord_start, coord_end) > subGR_dt$start & max(coord_start, coord_end) <= subGR_dt$end & gene == gene, region_id], by="GUIDE_SEQ"]


# For each gene subregion, select the best scoring guide by ranking the scores and suming the ranks within each region
titration_selguides_l <- list()
for (n in 1:nrow(subGR_dt)) {
  dt <- GDT.6[region_id == subGR_dt[n, region_id],] # Higher scores indicate higher predicted on-target guide activity.
  dt[, rank_Hsu2013 := frank(-Hsu2013)] # Scores range from 0 to 100; higher scores indicate lower predicted off-target activity. 
  dt[, rank_otCount := frank(otCount)] # A higher score indicates lower predicted off-target activity
  dt[, rank_Doench2014OnTarget := frank(-Doench2014OnTarget)]
  dt[, sum_ranks := sum(rank_Doench2014OnTarget, rank_Hsu2013, rank_otCount), GUIDE_SEQ]
  dt_ord <- dt[base::order(sum_ranks),]
  titration_selguides_l[[n]] <- dt_ord[1,]
}
titration_selguides_dt <- rbindlist(titration_selguides_l)
GDT.6[, selected := ifelse(GUIDE_SEQ %in% titration_selguides_dt$GUIDE_SEQ, TRUE, FALSE)]

# Make sure each region has one selected guide, if not pick manually one that fall in between regions
GDT.6[selected == T, .N, region_id] # GFI1B GFI1B_up_4 is missing
GDT.6[45:50, .(ID, Doench2014OnTarget, Hsu2013, DoenchCFD_specificityscore, otCount, region_id, selected)]
GDT.6[ID == "GFI1B-chr9-132978502-132978483", c("region_id", "selected") := list("GFI1B_up_4", TRUE)]


ggpairs(GDT.6[!is.na(region_id), .(gene, selected, Doench2014OnTarget, Hsu2013, DoenchCFD_specificityscore, otCount)], columns = 3:6, aes(colour=selected, alpha=0.7))

#### (6) Double chech reasonable distances between guides


GDT.7 <- GDT.6[selected == T, ]
GDT.7[, unstrand_coord_start := min(coord_start, coord_end), GUIDE_SEQ]

dist_dt <- cbind(setnames(GDT.7[base::order(unstrand_coord_start), unstrand_coord_start, .(region_id)][1:(nrow(GDT.7)-1)], new = c("region_id_g1", "unstrand_coord_start_g1")),
                 setnames(GDT.7[base::order(unstrand_coord_start), unstrand_coord_start, .(region_id)][2:nrow(GDT.7)], new = c("region_id_g2", "unstrand_coord_start_g2")))
dist_dt[, inter_region_id := paste0(region_id_g1, "-", region_id_g2)]
dist_dt[, dist_bp := unstrand_coord_start_g2 - unstrand_coord_start_g1, ]
dist_dt[, gene := tstrsplit(inter_region_id, "_", keep = 1)]

ggplot(dist_dt[dist_bp<500,], aes(x=dist_bp, fill=gene)) + geom_histogram(bins = 30)


##### (7) Overlap weismanTSS guides and subsitute close or overlaping ones for new extra titration guides


# Plot guides distributed in their genomic locations and show predicted TSS position
TSS_dt <- fread(input = file.path("processed_data/titration_genes_TSS_unique_site.bed"),
                col.names = c("chr", "TSS_start", "TSS_end", "gene_tss", "cvrg_score", "strand"))
TSS_dt[, gene := tstrsplit(gene_tss, "_", keep = 1)]
weismanTSS_scored_dt[, gene := tstrsplit(ID, "-", keep = 1)]
weismanTSS_scored_dt[, c("coord_start", "coord_end") := tstrsplit(ID, "-", keep = c(3,4))]
weismanTSS_scored_dt[, c("coord_start", "coord_end") := list(as.numeric(coord_start), as.numeric(coord_end))]
weismanTSS_scored_dt[, GUIDE_STRAND := ifelse(coord_start < coord_end, "+", "-")]


plot_dt <- merge(rbind(cbind(data.table::melt(GDT.7[, .(ID, gene, coord_start, coord_end, GUIDE_STRAND)], id.vars = c("ID" ,"gene", "GUIDE_STRAND"), value.name = "coord", variable.name = "variable"), data.table(category="titration")),
                       cbind(data.table::melt(weismanTSS_scored_dt[, .(ID, gene, coord_start, coord_end, GUIDE_STRAND)], id.vars = c("ID" ,"gene", "GUIDE_STRAND"), value.name = "coord", variable.name = "variable"), data.table(category="weismanTSS"))),
                 TSS_dt[, .(TSS_start, gene)], by="gene")
ggplot(plot_dt, aes(x = coord, y=category, color=GUIDE_STRAND, group=ID)) +
  geom_line(data = plot_dt[GUIDE_STRAND == "+",], arrow = arrow(length=unit(0.15,"cm"), ends="last", type = "closed")) +
  geom_line(data = plot_dt[GUIDE_STRAND == "-",], arrow = arrow(length=unit(0.15,"cm"), ends="first", type = "closed")) +
  geom_vline(data = filter(plot_dt, gene=="GFI1B"), aes(xintercept = TSS_start), lty = 2) +
  geom_vline(data = filter(plot_dt, gene=="NFE2"), aes(xintercept = TSS_start), lty = 2) +
  geom_vline(data = filter(plot_dt, gene=="MYB"), aes(xintercept = TSS_start), lty = 2) +
  geom_vline(data = filter(plot_dt, gene=="TET2"), aes(xintercept = TSS_start), lty = 2) +
  facet_wrap(. ~ gene, scales = "free", nrow = 4) +
  scale_color_manual("strand", values = c("red","#1b98e0")) +
  theme(legend.position = "bottom")

GDT.7[gene == "GFI1B", .(ID, region_id, Doench2014OnTarget, Hsu2013, otCount)]
weismanTSS_scored_dt[gene == "GFI1B", .(ID, Doench2014OnTarget, Hsu2013, otCount)]

GDT.7[gene == "MYB", .(ID, region_id, Doench2014OnTarget, Hsu2013, otCount)]
weismanTSS_scored_dt[gene == "MYB", .(ID, Doench2014OnTarget, Hsu2013, otCount)]

GDT.7[gene == "NFE2", .(ID, region_id, Doench2014OnTarget, Hsu2013, otCount)]
weismanTSS_scored_dt[gene == "NFE2", .(ID, Doench2014OnTarget, Hsu2013, otCount)]

GDT.7[gene == "TET2", .(ID, region_id, Doench2014OnTarget, Hsu2013, otCount)]
weismanTSS_scored_dt[gene == "TET2", .(ID, Doench2014OnTarget, Hsu2013, otCount)]


## Flag those guides that overlap with the WeismanTSS ones and have worst scores
rid_remove = c("GFI1B_up_7", "GFI1B_tss_1", "NFE2_down_4", "MYB_tss_1", "TET2_down_5")
rid_subst = c("GFI1B_up_0", "GFI1B_down_6", "NFE2_down_0", "MYB_down_6", "TET2_down_7")
GDT.7[, substitute_for_weismanTSS := ifelse(region_id %in% rid_remove, T, F) ]

## Generate new seq regions, design guides and score to choose the new best ones
extraGR_dt <- data.table("gene_rid_subs" = rid_subst,
                         start = c(min(subGR_dt[gene == "GFI1B", start])-200,
                                   max(subGR_dt[gene == "GFI1B", end])+1, 
                                   min(subGR_dt[gene == "NFE2", start])-200,
                                   max(subGR_dt[gene == "MYB", end])+1,
                                   max(weismanTSS_scored_dt[gene == "TET2", coord_end])+30
                         ),
                         end = c(min(subGR_dt[gene == "GFI1B", start])-1, 
                                 max(subGR_dt[gene == "GFI1B", end])+200,
                                 min(subGR_dt[gene == "NFE2", start])-1,
                                 max(subGR_dt[gene == "MYB", end])+200,
                                 max(weismanTSS_scored_dt[gene == "TET2", coord_end])+230
                         )
)
extraGR_dt[, gene := tstrsplit(gene_rid_subs, "_", keep = 1)]
extraGR_dt <- merge.data.table(extraGR_dt, GR_dt[, .(gene, chr, strand)], by = "gene")

# Generate all possible guides
GDT.1b <- extraGR_dt[, GetAllGuidesInRegion(regionToSearch = getSeq(BSgenome.Hsapiens.UCSC.hg38::Hsapiens, names=chr, start=start, end=end), PAMpattern = DNAString("NGG")), gene_rid_subs]

# Filter guides based on different criterias
spacerVec = GDT.1b$GUIDE_SEQ

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

GDT.2b <- GDT.1b[-uniqueIndsToRemove,]

## (3) Discard guides that fall in unique genetic variation of K562 (SNVs and SNPs)

# Get actual genomic coordinates of the guides
GDT.3b <- merge.data.table(GDT.2b, setnames(extraGR_dt[, .(gene_rid_subs, chr, start)], new=c("gene_rid_subs", "chr", "regions_start")), by="gene_rid_subs")
GDT.3b[, c("coord_start", "coord_end") := list(regions_start+GUIDE_START-1, regions_start+GUIDE_END-1)] # add actual genomic coordinates
GDT.3b[,  coord_pam_end := ifelse(GUIDE_STRAND == "+", coord_end+3, coord_end-3)]

# Generate a GRanges object to liftover the hg38 cooordinates to hg19 (vcf of K562 genome assembly is hg19)
gr_GDT.3b <- makeGRangesFromDataFrame(
  as.data.frame(rbind(GDT.3b[GUIDE_STRAND == "+", .(chr, coord_start, coord_pam_end, GUIDE_SEQ)], 
                      setnames(GDT.3b[GUIDE_STRAND == "-", .(chr, coord_pam_end, coord_start, GUIDE_SEQ)], c("chr", "coord_start", "coord_pam_end", "GUIDE_SEQ")))),
  keep.extra.columns = T,
  ignore.strand=T,
  start.field = "coord_start",
  end.field = "coord_pam_end")

# liftOver coordinates to Hg19
seqlevelsStyle(gr_GDT.3b) = "UCSC"  # necessary
gr_GDT.3b_hg19 = liftOver(gr_GDT.3b, ch)
gr_GDT.3b_hg19 = unlist(gr_GDT.3b_hg19)
genome(gr_GDT.3b_hg19) = "hg19"

# Generate bed file with guides coordinates
fwrite(as.data.table(gr_GDT.3b_hg19)[, .(seqnames, start, end, GUIDE_SEQ)],
       file = "processed_data/bedtools/D2N_GDT3b_hg19.bed",
       quote = F, row.names = F, col.names = F, sep="\t")

# Get guides to retain 
GDT.3b_hg19_K562 <- fread("processed_data/bedtools/D2N_GDT3b_hg19_intK562.bed", col.names = c("chr", "start", "end", "GUIDE_SEQ"))

# Keep guides that K562 does not differe from the human ref
GDT.4b <- GDT.3b[GUIDE_SEQ %in% GDT.3b_hg19_K562$GUIDE_SEQ, ]


## (4) Generate table to run FlashFry
# These are the flanking random sequences added
# TTCGTACAAA
# TGTCCGCACT
GDT.5b <- GDT.4b
GDT.5b[, FlasFryID := paste0(gsub("_", "\\.", gene_rid_subs), "-", chr, "-", coord_start, "-", coord_end)]
GDT.5b[, FlashFrySEQ := paste0("TTCGTACAAA", GUIDE_SEQ, "TGTCCGCACT")]

fwrite(GDT.5b, file = "processed_data/guides/D2N_titrationExtra_ToScore.txt", quote = F, row.names = F, sep="\t")
GDT.5b_flashfry_scored <- fread(file = "processed_data/guides/D2N_titrationExtra_Scored.txt")
GDT.5b_flashfry_scored[, is_distalCRE := guide %in% distalCRE_scored_dt$guide]

GDT.6b <- merge(GDT.5b_flashfry_scored[is_distalCRE == F, .(ID, Doench2014OnTarget, DoenchCFD_specificityscore, Hsu2013, otCount)], GDT.5b, by.x = "ID", by.y = "FlasFryID")
GDT.6b[, guide_no_pam := substr(GUIDE_SEQ, 1, 20)]
GDT.6b[, gene := tstrsplit(ID, "\\.", keep = 1)]
setnames(GDT.6b, old="gene_rid_subs", new="region_id")

ggplot(GDT.6b, aes(x=Doench2014OnTarget, fill=region_id)) + geom_histogram(bins = 50) 
ggplot(GDT.6b, aes(x=Hsu2013, fill=region_id)) + geom_histogram(bins = 50) 
ggplot(GDT.6b, aes(x=DoenchCFD_specificityscore, fill=region_id)) + geom_histogram(bins = 50) 
ggplot(GDT.6b, aes(x=otCount, fill=region_id)) + geom_histogram(bins = 50) 


# For each gene extra subregion, select the best scoring guide by ranking the scores and suming the ranks within each region (try to avoid more off-target than on-target activity now)
titration_selguides_l <- list()
rid_vec <- unique(GDT.6b$region_id)
for (n in 1:length(rid_vec)) {
  print(n)
  dt <- GDT.6b[region_id == rid_vec[n],] # Higher scores indicate higher predicted on-target guide activity.
  dt[, rank_Hsu2013 := frank(-Hsu2013)] # Scores range from 0 to 100; higher scores indicate lower predicted off-target activity. 
  dt[, rank_otCount := frank(otCount)] # A higher score indicates lower predicted off-target activity
  dt[, rank_Doench2014OnTarget := frank(-Doench2014OnTarget)]
  dt[, sum_ranks := sum(0.2*rank_Doench2014OnTarget, 0.4*rank_Hsu2013, 0.4*rank_otCount), GUIDE_SEQ]
  dt_ord <- dt[base::order(sum_ranks),]
  titration_selguides_l[[n]] <- dt_ord[1,]
}
titration_selguides_dt <- rbindlist(titration_selguides_l)
GDT.6b[, selected := ifelse(GUIDE_SEQ %in% titration_selguides_dt$GUIDE_SEQ, TRUE, FALSE)]

## Merge both initial titration and extra titration guides in a single data.table
GDT.8 <- rbind(GDT.6[!(GUIDE_SEQ %in% GDT.7[substitute_for_weismanTSS == T, GUIDE_SEQ]), ], GDT.6b)

p1 <- ggpairs(GDT.8[!is.na(region_id), .(gene, selected, Doench2014OnTarget, Hsu2013, DoenchCFD_specificityscore, otCount)], columns = 3:6, aes(colour=selected, alpha=0.7))
ggsave(p1, filename = "plots/20220118_titration_metrics_selectedVsnon.pdf", height = 7, width = 8)

#### ¡¡ FINAL TITRATION GUIDE SELECTION !! ####
GDT.9 <- GDT.8[selected == T,]

GDT.9[, unstrand_coord_start := min(coord_start, coord_end), GUIDE_SEQ]

dist_dt <- cbind(setnames(GDT.9[base::order(unstrand_coord_start), unstrand_coord_start, .(region_id)][1:(nrow(GDT.9)-1)], new = c("region_id_g1", "unstrand_coord_start_g1")),
                 setnames(GDT.9[base::order(unstrand_coord_start), unstrand_coord_start, .(region_id)][2:nrow(GDT.9)], new = c("region_id_g2", "unstrand_coord_start_g2")))
dist_dt[, inter_region_id := paste0(region_id_g1, "-", region_id_g2)]
dist_dt[, dist_bp := unstrand_coord_start_g2 - unstrand_coord_start_g1, ]
dist_dt[, gene := tstrsplit(inter_region_id, "_", keep = 1)]

p2 <- ggplot(dist_dt[dist_bp<500,], aes(x=dist_bp, fill=gene)) + geom_histogram(bins = 30) +
  theme(legend.position = "bottom", legend.direction = "vertical")
ggsave(p2, filename = "plots/20220118_titration_distances_dist.pdf", width = 4, height = 5)



plot_dt <- merge(rbind(cbind(data.table::melt(GDT.9[, .(ID, gene, coord_start, coord_end, GUIDE_STRAND)], id.vars = c("ID" ,"gene", "GUIDE_STRAND"), value.name = "coord", variable.name = "variable"), data.table(category="titration")),
                       cbind(data.table::melt(weismanTSS_scored_dt[, .(ID, gene, coord_start, coord_end, GUIDE_STRAND)], id.vars = c("ID" ,"gene", "GUIDE_STRAND"), value.name = "coord", variable.name = "variable"), data.table(category="weismanTSS"))),
                 TSS_dt[, .(TSS_start, gene)], by="gene")
p3 <- ggplot(plot_dt, aes(x = coord, y=category, color=GUIDE_STRAND, group=ID)) +
  geom_line(data = plot_dt[GUIDE_STRAND == "+",], arrow = arrow(length=unit(0.15,"cm"), ends="last", type = "closed")) +
  geom_line(data = plot_dt[GUIDE_STRAND == "-",], arrow = arrow(length=unit(0.15,"cm"), ends="first", type = "closed")) +
  geom_vline(data = filter(plot_dt, gene=="GFI1B"), aes(xintercept = TSS_start), lty = 2) +
  geom_vline(data = filter(plot_dt, gene=="NFE2"), aes(xintercept = TSS_start), lty = 2) +
  geom_vline(data = filter(plot_dt, gene=="MYB"), aes(xintercept = TSS_start), lty = 2) +
  geom_vline(data = filter(plot_dt, gene=="TET2"), aes(xintercept = TSS_start), lty = 2) +
  facet_wrap(. ~ gene, scales = "free", nrow = 4) +
  scale_color_manual("strand", values = c("red","#1b98e0")) +
  theme(legend.position = "bottom")
p3
ggsave(p3, filename = "plots/20220118_titration_genomic_location.pdf", width = 7.5, height = 5)


# Select WT guide for designing the attenuated ones

AttWT_dt <- rbind(setnames(weismanTSS_scored_dt[ID == "GFI1B-chr9-132978817-132978798", .(ID, guide, GUIDE_STRAND, Doench2014OnTarget, otCount)], new = c("ID", "WT_GUIDE_SEQ", "GUIDE_STRAND", "Doench2014OnTarget_WT", "otCount_WT")),
                  setnames(GDT.9[ID == "NFE2-chr12-54300900-54300919", .(ID, GUIDE_SEQ, GUIDE_STRAND, Doench2014OnTarget, otCount)], new = c("ID", "WT_GUIDE_SEQ", "GUIDE_STRAND", "Doench2014OnTarget_WT", "otCount_WT")),
                  setnames(weismanTSS_scored_dt[ID == "MYB-chr6-135181383-135181402", .(ID, guide, GUIDE_STRAND, Doench2014OnTarget, otCount)], new = c("ID", "WT_GUIDE_SEQ", "GUIDE_STRAND", "Doench2014OnTarget_WT", "otCount_WT")),
                  setnames(GDT.9[ID == "TET2-chr4-105146774-105146755", .(ID, GUIDE_SEQ, GUIDE_STRAND, Doench2014OnTarget, otCount)], new = c("ID", "WT_GUIDE_SEQ", "GUIDE_STRAND", "Doench2014OnTarget_WT", "otCount_WT")))


fwrite(AttWT_dt, file = "processed_data/guides/D2N_WT_attenuated.txt", row.names = F, quote = F)


##### Export table with final selection to collapse with the additional guides of the library #####
# Add general category
GDT.9[, guide_category := "titration"]
GDT.9[, GUIDE_SEQ_NO_PAM := substr(GUIDE_SEQ, 1, 20), GUIDE_SEQ]
GDT.9[, id := paste0(gene, "_", chr, ":", coord_start, "-", coord_end, "_", GUIDE_STRAND), GUIDE_SEQ]


# Write file of selected guides
fwrite(GDT.9[, .(id, coord_start, coord_end, GUIDE_SEQ, GUIDE_STRAND, GUIDE_SEQ_NO_PAM, PAM, Doench2014OnTarget, Hsu2013, otCount, guide_category)],
       file = "processed_data/guides/D2N_titration_Scored_Selected.txt",
       quote = F, row.names = F)











