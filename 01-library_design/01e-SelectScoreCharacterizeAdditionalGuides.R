##### Manually select, score and include genomic features of additional sgRNAs in the d2n library #####
# 1) Select guides from Reploge et al. BiorXiv 2021 (Weisman's lab large-scale Pertub-seq), score and add genomic features
# 2) Include distal CRE guides: John's GFI1B and NFE2 enhancers, score and add genomic features 
# 3) Generate guides for MYB enhancers and select two best scoring
# 4) Include attenuated sgRNAs for each gene and score after missmatches introduced
# 5) Select 5 NT guides
# 6) Collapse all guides with dosage-titration ones into single table

# JDE, January 2022

## Libraries
library(data.table)
library(rtracklayer)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)
library(seqinr)
library(Biostrings)
library(ggplot2); theme_set(theme_bw())
library(GGally)



## Directories
setwd("/gpfs/commons/groups/lappalainen_lab/jdomingo/projects/004-dosage_network/001-library_design")
data_dir = "/gpfs/commons/groups/lappalainen_lab/jdomingo/data/"





##### (1) CRISPRi sgRNAs from Reploge et al. 2021 - Global category = weismanTSS #####

# Load data
Reploge_dt <- do.call("rbind", lapply(list.files(file.path(data_dir, "Replogle_Cell_2022")), function(ff){
  if(grepl(pattern = "^sgRNAs_.*csv$", ff) == T) {
    fread(input = file.path(data_dir, "Replogle_Cell_2022", ff))
  }
}))

# Check if dosage titration genes were targeted in this study
target_genes = c("GFI1B", "NFE2", "MYB", "TET2")
chr = c("chr9", "chr12", "chr6", "chr4")
names(chr) <- target_genes
Sel_dt <- rbind(setnames(Reploge_dt[gene %in% target_genes, .(sgID_A, `targeting sequence A`)], new = c("sgID", "GUIDE_SEQ")),
                setnames(Reploge_dt[gene %in% target_genes, .(sgID_B, `targeting sequence B`)], new = c("sgID", "GUIDE_SEQ")))
Sel_dt[, c("gene", "strand_Reploge", "coord") := tstrsplit(sgID, "_", fixed = TRUE)]
Sel_dt[, coord := as.numeric(gsub("\\.23-P1P2", "", coord))]
Sel_dt[, coord_end := coord +1]
Sel_dt[, chr := chr[gene]]
# The strand of the guides is inverted, correct that to get the same guide start and end coordinates
Sel_dt[, GUIDE_STRAND := ifelse(strand_Reploge == "+", "-", "+")]


# Get coordinates and liftOver to get the Hg38 ones
gr_Sel_hg19 <- makeGRangesFromDataFrame(
  as.data.frame(Sel_dt[, .(chr, coord, coord_end, GUIDE_SEQ)]),
  keep.extra.columns = T,
  ignore.strand=T,
  start.field = "coord",
  end.field = "coord_end")

# The chain file for hg19 to hg38 transformation
ch = import.chain("/gpfs/commons/groups/lappalainen_lab/jdomingo/data/genomes/Hsapiens/hg19ToHg38.over.chain")

# liftOver coordinates to Hg38
seqlevelsStyle(gr_Sel_hg19) = "UCSC"  # necessary
gr_Sel_hg38 = liftOver(gr_Sel_hg19, ch)
gr_Sel_hg38 = unlist(gr_Sel_hg38)
genome(gr_Sel_hg38) = "hg38"

# Create a dataset with the new Hg38 coordinates
Sel_dt_hg38 <- merge(Sel_dt[, .(sgID, GUIDE_SEQ, gene, GUIDE_STRAND, chr)],
                     as.data.table(gr_Sel_hg38)[, .(start, GUIDE_SEQ)], 
                     by="GUIDE_SEQ")
# Correct the start and end of the coordinates from the single position given in the supplementary table
setnames(Sel_dt_hg38, old="start", new="coord_Reploge")
Sel_dt_hg38[, coord_start := ifelse(GUIDE_STRAND == "+", as.integer(coord_Reploge-21), as.integer(coord_Reploge +23))]
Sel_dt_hg38[, coord_end := ifelse(GUIDE_STRAND == "+", as.integer(coord_Reploge-2), as.integer(coord_Reploge+4))]

# Get the exact genomic sequence of the guide in Hg38
Sel_dt_hg38[GUIDE_STRAND == "+", genomic_GUIDE_PAM_SEQ := as.character(getSeq(BSgenome.Hsapiens.UCSC.hg38::Hsapiens, names= chr, start = coord_start, end = coord_end +3))]
Sel_dt_hg38[GUIDE_STRAND == "-", genomic_GUIDE_PAM_SEQ := as.character(reverseComplement(getSeq(BSgenome.Hsapiens.UCSC.hg38::Hsapiens, names= chr, start = coord_end -3, end = coord_start)))]

Sel_dt_hg38[, genomic_GUIDE_SEQ := substr(genomic_GUIDE_PAM_SEQ, 1, 20)]
Sel_dt_hg38[, PAM := substr(genomic_GUIDE_PAM_SEQ, 21, 23)]

# Find if there are mismatches in the designed guide sequence and annotate them
Sel_dt_hg38[, MM_NUM := sum(unlist(strsplit(GUIDE_SEQ, "")) != unlist(strsplit(genomic_GUIDE_SEQ, ""))), GUIDE_SEQ]
Sel_dt_hg38[, MM_POS := ifelse(MM_NUM <= 1,
                               ifelse(MM_NUM == 0, 0, which(unlist(strsplit(GUIDE_SEQ, "")) != unlist(strsplit(genomic_GUIDE_SEQ, "")))),
                               paste(which(unlist(strsplit(GUIDE_SEQ, "")) != unlist(strsplit(genomic_GUIDE_SEQ, ""))), collapse = ",")), GUIDE_SEQ]
Sel_dt_hg38[, MM_WT_NT := ifelse(MM_NUM <= 1,
                               ifelse(MM_NUM == 0, "-", unlist(strsplit(genomic_GUIDE_SEQ, ""))[MM_POS]),
                               paste(unlist(strsplit(genomic_GUIDE_SEQ, ""))[as.numeric(unlist(strsplit(MM_POS, ",")))], collapse = ",")), GUIDE_SEQ]
Sel_dt_hg38[, MM_MUT_NT := ifelse(MM_NUM <= 1,
                                 ifelse(MM_NUM == 0, "-", unlist(strsplit(GUIDE_SEQ, ""))[MM_POS]),
                                 paste(unlist(strsplit(GUIDE_SEQ, ""))[as.numeric(unlist(strsplit(MM_POS, ",")))], collapse = ",")), GUIDE_SEQ]

# Write data.table to run FlashFry and score guides  
DT.1 <- Sel_dt_hg38[, .(gene, GUIDE_STRAND, GUIDE_SEQ, PAM, chr, coord_start, coord_end)]
DT.1[, GUIDE_SEQ := paste0(GUIDE_SEQ, PAM), GUIDE_SEQ]
DT.1[, GUIDE_START := 1]
DT.1[, GUIDE_END := 23]
DT.1[, Design_PAM := "NGG"]
DT.1[, regions_start := 1]
DT.1[, coord_pam_end := ifelse(GUIDE_STRAND == "+", coord_end + 3, coord_end - 3)]
DT.1[, FlasFryID := paste0(gene, "-", chr, "-", coord_start, "-", coord_end)]
DT.1[, FlashFrySEQ := paste0("TTCGTACAAA", GUIDE_SEQ, "TGTCCGCACT")]

# fwrite(DT.1[, .(gene, GUIDE_START, GUIDE_END, GUIDE_STRAND, GUIDE_SEQ, PAM, Design_PAM, chr, regions_start, coord_start, coord_end, coord_pam_end, FlasFryID, FlashFrySEQ)], 
       # file = "processed_data/guides/D2N_weismanTSS_ToScore.txt", quote = F, row.names = F, sep="\t")

# Add features to collapse with the rest of the guides
DT.1[, id := paste0(gene, "_", chr, ":", coord_start, "-", coord_end, "_", GUIDE_STRAND), GUIDE_SEQ]
DT.1[, GUIDE_SEQ_NO_PAM := substr(GUIDE_SEQ, 1, 20), GUIDE_SEQ]
DT.1[, guide_category := "weismanTSS"]
DT.1_scored <- merge.data.table(DT.1, fread("processed_data/guides/D2N_weismanTSS_Scored.txt")[, .(ID, Doench2014OnTarget, Hsu2013, otCount)], by.x = "FlasFryID", by.y = "ID")


# Write file of selected guides
fwrite(DT.1_scored[, .(id, coord_start, coord_end, GUIDE_SEQ, GUIDE_STRAND, GUIDE_SEQ_NO_PAM, PAM, Doench2014OnTarget, Hsu2013, otCount, guide_category)],
       file = "processed_data/guides/D2N_weismanTSS_Scored_Selected.txt",
       quote = F, row.names = F)









##### (2-3) DISTAL CRE guides: SSv1 GFI1B & NFE2 enhancers + Literature MYB enhancers - Global category = distalCRE #####

# Load data and double check guides and coordinates are correct

# FI1B & NFE2 (SSv1) MYB enhancers guides (Li et al. 2021) 
# Load results from UCSC BLAT, merge with ids from fasta file
dd1 <- as.character(Biostrings::readDNAStringSet("processed_data/guides/D2N_distalCRE_guides.fasta"))

distalCRE_dt1 <- merge.data.table(fread("processed_data/guides/D2N_distalCRE_guides.csv"), data.table(id = names(dd1), guide_seq = dd1), by ="id")
distalCRE_dt1[, coord_start := ifelse(GUIDE_STRAND == "+", ucsc_coord_start, ucsc_coord_end)]
distalCRE_dt1[, coord_end := ifelse(GUIDE_STRAND == "+", ucsc_coord_end, ucsc_coord_start)]

# Load MYB additional enhancers guides (Xie et al. 2020)
distalCRE_dt2 <- fread("processed_data/guides/MYB_enh_guides_Li2021.txt", col.names = c("long_id", "guide_seq", "rnd", "enh_coord"))
distalCRE_dt2[,  c("short_id", "guide_coord", "GUIDE_STRAND") := tstrsplit(long_id, "|", fixed = TRUE)]
distalCRE_dt2[, id := paste0(gsub(">MYB-enh", "MYB_distalCRE_Xie2020-Enh", short_id), "-G", seq(1, .N)), by="enh_coord"]
distalCRE_dt2[, c("chr", "guide_coord") := tstrsplit(guide_coord, ":", fixed = TRUE)]

# Coordinates of the guides are in hg19 assembly
distalCRE_dt2[,  c("hg19_coord_start", "hg19_coord_end") := tstrsplit(guide_coord, "-", fixed = TRUE)]

# Get coordinates and liftOver to get the Hg38 ones
gr_distalCRE_dt2_hg19 <- makeGRangesFromDataFrame(
  as.data.frame(distalCRE_dt2[, .(chr, hg19_coord_start, hg19_coord_end, guide_seq)]),
  keep.extra.columns = T,
  ignore.strand=T,
  start.field = "hg19_coord_start",
  end.field = "hg19_coord_end")

# liftOver coordinates to Hg38
seqlevelsStyle(gr_distalCRE_dt2_hg19) = "UCSC"  # necessary
gr_distalCRE_dt2_hg38 = liftOver(gr_distalCRE_dt2_hg19, ch)
gr_distalCRE_dt2_hg38 = unlist(gr_distalCRE_dt2_hg38)
genome(gr_distalCRE_dt2_hg38) = "hg38"

# Put coordinates following same structure as with other guide data.tables
distalCRE_dt3 <- merge(distalCRE_dt2, as.data.table(gr_distalCRE_dt2_hg38)[, .(guide_seq, start, end)], by="guide_seq")
distalCRE_dt3[, coord_start := ifelse(GUIDE_STRAND == "+", start, end)]
distalCRE_dt3[, coord_end := ifelse(GUIDE_STRAND == "+", end, start)]

# Merge the two data.tables into single distalCRE data.table
DT.2 <- rbind(distalCRE_dt1[, .(id, chr, coord_start, coord_end, GUIDE_STRAND, guide_seq)],
              distalCRE_dt3[, .(id, chr, coord_start, coord_end, GUIDE_STRAND, guide_seq)])
DT.2[GUIDE_STRAND == "+", genomic_GUIDE_PAM_SEQ := as.character(getSeq(BSgenome.Hsapiens.UCSC.hg38::Hsapiens, names= chr, start = coord_start, end = coord_end +3))]
DT.2[GUIDE_STRAND == "-", genomic_GUIDE_PAM_SEQ := as.character(reverseComplement(getSeq(BSgenome.Hsapiens.UCSC.hg38::Hsapiens, names= chr, start = coord_end -3, end = coord_start)))]
DT.2[, PAM := substr(genomic_GUIDE_PAM_SEQ, 21, 23)]

# Guides from Xie et al are missing a nucleotide, add at begining or end depending on the strand
DT.2[grep("Xie2020", id), guide_seq := ifelse(GUIDE_STRAND == "+", paste0(substr(genomic_GUIDE_PAM_SEQ, 1, 1), guide_seq), paste0(guide_seq, substr(genomic_GUIDE_PAM_SEQ, 20, 20)))]

# Add remaining features to run FlashFry
DT.2[, gene := tstrsplit(id, "_", keep = 1), guide_seq]
DT.2[, GUIDE_SEQ := paste0(guide_seq, PAM), guide_seq]
DT.2[, GUIDE_START := 1]
DT.2[, GUIDE_END := 23]
DT.2[, Design_PAM := "NGG"]
DT.2[, regions_start := 1]
DT.2[, coord_pam_end := ifelse(GUIDE_STRAND == "+", coord_end + 3, coord_end - 3)]
DT.2[, FlasFryID := paste0(gene, "-", chr, "-", coord_start, "-", coord_end)]
DT.2[, FlashFrySEQ := paste0("TTCGTACAAA", GUIDE_SEQ, "TGTCCGCACT")]


# fwrite(DT.2[, .(gene, GUIDE_START, GUIDE_END, GUIDE_STRAND, GUIDE_SEQ, PAM, Design_PAM, chr, regions_start, coord_start, coord_end, coord_pam_end, FlasFryID, FlashFrySEQ)], 
#        file = "processed_data/guides/D2N_distalCRE_ToScore.txt", quote = F, row.names = F, sep="\t")

## Select MYB guides based on the resulted scored guides

# Add features
DT.3 <- merge(fread(input = "processed_data/guides/D2N_distalCRE_Scored.txt")[, .(ID, Doench2014OnTarget, DoenchCFD_specificityscore, Hsu2013, otCount)], 
              DT.2, by.x = "ID", by.y = "FlasFryID")
DT.3[, gene := tstrsplit(id, "_", keep = 1)]
DT.3[, enh_id := gsub("-G[0-9]+$", "", id)]
DT.3[, guide_no_pam := substr(GUIDE_SEQ, 1, 20)]

# Flag those with restriction enzyme sites
# EcoRI restriction sites (full and star activity) 
EcorRI_main = "GAATTC"
EcorRI_star1 = "TAATTC"
EcorRI_star2 = "CAATTC"
# Constant upstream and downstream regions in the pGC003 to where the guide is cloned
cnst_5 = toupper("caccG")
cnst_3 = toupper("gttta")

DT.3[, plasmid_guide_seq := paste0(cnst_5, guide_no_pam, cnst_3)]
DT.3[, EcorRI_site := grepl(EcorRI_main, plasmid_guide_seq) | grepl(EcorRI_star1, plasmid_guide_seq) | grepl(EcorRI_star2, plasmid_guide_seq) ]
DT.3 <- DT.3[EcorRI_site == F,]

# Give ascending rank to each of the scoring metrics and aggregate into a rank sum statistic to chose best guide
DT.3[, rank_Doench2014OnTarget := frank(-Doench2014OnTarget)] # Higher scores indicate higher predicted on-target guide activity.
DT.3[, rank_Hsu2013 := frank(-Hsu2013)] # Scores range from 0 to 100; higher scores indicate lower predicted off-target activity. 
DT.3[, rank_DoenchCFD_specificityscore := frank(-DoenchCFD_specificityscore)] # A higher score indicates lower predicted off-target activity
DT.3[, rank_otCount := frank(otCount)] # A higher score indicates lower predicted off-target activity

# Sum ranks
DT.3[, sum_ranks := sum(rank_Doench2014OnTarget, rank_Hsu2013, rank_otCount), GUIDE_SEQ]

NumGuidesEnh_dt <- data.table(enh_id = unique(DT.3[, enh_id]))
NumGuidesEnh_dt[, num_guides := ifelse(grepl("MYB", enh_id), 1, 2)]

# Select best guide that will have the best scoring sum of ranks
distalCRE_selguides_l <- list()
for (n in 1:nrow(NumGuidesEnh_dt)) {
  dt <- DT.3[enh_id == NumGuidesEnh_dt[n, enh_id],]
  ord <- frankv(dt[, sum_ranks])
  distalCRE_selguides_l[[n]] <- dt[ord,][1:NumGuidesEnh_dt[n, num_guides]]
}
distalCRE_selguides_dt <- rbindlist(distalCRE_selguides_l)

DT.3[, selected := ifelse(GUIDE_SEQ %in% distalCRE_selguides_dt$GUIDE_SEQ, TRUE, FALSE)]


ggplot(DT.3, aes(x=Doench2014OnTarget, fill=gene)) + geom_histogram(bins = 50) 
ggplot(DT.3, aes(x=Hsu2013, fill=gene)) + geom_histogram(bins = 50) 
ggplot(DT.3, aes(x=DoenchCFD_specificityscore, fill=gene)) + geom_histogram(bins = 50) 
ggplot(DT.3, aes(x=otCount, fill=gene)) + geom_histogram(bins = 50) 

ggpairs(DT.3[, .(gene, selected, Doench2014OnTarget, Hsu2013, DoenchCFD_specificityscore, otCount)], columns = 3:6)


# Add features to the selected guides and write table to collapse with others
DT.3_selected <- DT.3[selected == T, ]
DT.3_selected[, guide_category := "distalCRE"]
DT.3_selected[, id_literature := id]
DT.3_selected[, id := paste0(gene, "_", chr, ":", coord_start, "-", coord_end, "_", GUIDE_STRAND), GUIDE_SEQ]
DT.3_selected[, GUIDE_SEQ_NO_PAM := substr(GUIDE_SEQ, 1, 20), GUIDE_SEQ]

# Write final table
fwrite(DT.3_selected[, .(id, coord_start, coord_end, GUIDE_SEQ, GUIDE_STRAND, GUIDE_SEQ_NO_PAM, PAM, Doench2014OnTarget, Hsu2013, otCount, guide_category)], 
       file = "processed_data/guides/D2N_distalCRE_Scored_Selected.txt", quote = F, row.names = F)











##### (4) Select NT sgRNAs #####
# Load guides used in SSv1, list of 12 NTC guides
SSv1_gRNA_list <- fread(input = file.path(data_dir, "GWAS-CRISPRi_phaseI", "sgRNA_list.txt"), fill = T, 
                        col.names = c("ID", "GUIDE_STRAND", "GUIDE_SEQ_NO_PAM", "PAM", "SNP", "ADT.1", "ADT.2") )

# Discard those with EcoRI restriction sites
SSv1_NTC <- SSv1_gRNA_list[grepl("NTC-", ID), ]
SSv1_NTC[, plasmid_guide_seq := paste0(cnst_5, GUIDE_SEQ_NO_PAM, cnst_3)]
SSv1_NTC[, EcorRI_site := grepl(EcorRI_main, plasmid_guide_seq) | grepl(EcorRI_star1, plasmid_guide_seq) | grepl(EcorRI_star2, plasmid_guide_seq) ]
SSv1_NTC <- SSv1_NTC[EcorRI_site == F,]

# Select 5 random NT guides
set.seed(1234)
NTC_id_vec <- sample(SSv1_NTC[, ID], 5)

DT.4 <- SSv1_NTC[ID %in% NTC_id_vec, .(ID, GUIDE_STRAND, GUIDE_SEQ_NO_PAM)]
DT.4[, c("coord_start", "coord_end") := list(NA, NA), GUIDE_SEQ_NO_PAM]
DT.4[, c("GUIDE_SEQ", "PAM") := list(NA, NA)]
DT.4[, c("Doench2014OnTarget", "Hsu2013", "otCount") := list(NA, NA, NA)]
DT.4[, guide_category := "NTC"]
DT.4[, id := paste0("NTC-", 1:.N)]


# Write final table
fwrite(DT.4[, .(id, coord_start, coord_end, GUIDE_SEQ, GUIDE_STRAND, GUIDE_SEQ_NO_PAM, PAM, Doench2014OnTarget, Hsu2013, otCount, guide_category)], 
       file = "processed_data/guides/D2N_NTC_Selected.txt", quote = F, row.names = F)





##### (5) Select attenuated gRNAs #####
AttWT_dt <- fread(input = "processed_data/guides/D2N_WT_attenuated.txt")
AttWT_dt[, WT_GUIDE_SEQ_NO_PAM := substr(WT_GUIDE_SEQ, 1, 20)]

# Function to generate all possible single mutants given a range of positions and the guide sequence
GetAllSingleMutants <- function(guide_seq, pos_range=2:20) {
  actg = c("A", "C", "T", "G")
  return(do.call("c", lapply(pos_range, function(x){
    guide_nts =unlist(strsplit(guide_seq, ""))
    muts = actg[!(actg %in% guide_nts[x])]
    paste0(substr(guide_seq, 0, x-1), muts, substr(guide_seq, x+1, nchar(guide_seq)))
  })))
}

# Double check if mutations make sense
DT.5 <- setnames(AttWT_dt[, GetAllSingleMutants(WT_GUIDE_SEQ_NO_PAM), .(ID, WT_GUIDE_SEQ, GUIDE_STRAND, Doench2014OnTarget_WT,otCount_WT ,WT_GUIDE_SEQ_NO_PAM)], old = "V1", new="GUIDE_SEQ_NO_PAM")
DT.5[, MM_NUM := sum(unlist(strsplit(GUIDE_SEQ_NO_PAM, "")) != unlist(strsplit(WT_GUIDE_SEQ_NO_PAM, ""))), GUIDE_SEQ_NO_PAM]
DT.5[, MM_POS := which(unlist(strsplit(GUIDE_SEQ_NO_PAM, "")) != unlist(strsplit(WT_GUIDE_SEQ_NO_PAM, ""))), GUIDE_SEQ_NO_PAM]
DT.5[, MM_WT_NT := unlist(strsplit(WT_GUIDE_SEQ_NO_PAM, ""))[MM_POS], GUIDE_SEQ_NO_PAM]
DT.5[, MM_MUT_NT := unlist(strsplit(GUIDE_SEQ_NO_PAM, ""))[MM_POS], GUIDE_SEQ_NO_PAM]
DT.5[, MM_MUT_ID := paste0(MM_WT_NT, MM_POS, MM_MUT_NT), GUIDE_SEQ_NO_PAM]

# Discard guides EcorRI restriction sites
DT.5[, plasmid_guide_seq := paste0(cnst_5, GUIDE_SEQ_NO_PAM, cnst_3)]
DT.5[, EcorRI_site := grepl(EcorRI_main, plasmid_guide_seq) | grepl(EcorRI_star1, plasmid_guide_seq) | grepl(EcorRI_star2, plasmid_guide_seq)]
DT.5 <- DT.5[EcorRI_site == F, ]

# Write table to run FlashFry on each of these single mutants guides
DT.5[, gene := tstrsplit(ID, "-", keep = 1), GUIDE_SEQ_NO_PAM]
DT.5[, GUIDE_SEQ := paste0(GUIDE_SEQ_NO_PAM, substr(WT_GUIDE_SEQ, 21, 23)), GUIDE_SEQ_NO_PAM]
DT.5[, GUIDE_START := 1]
DT.5[, GUIDE_END := 20]
DT.5[, Design_PAM := "NGG"]
DT.5[, regions_start := 1]
DT.5[, coord_pam_end := 23]
DT.5[, FlasFryID := paste0(ID, "-", MM_MUT_ID)]
DT.5[, FlashFrySEQ := paste0("TTCGTACAAA", GUIDE_SEQ, "TGTCCGCACT")]
DT.5[, PAM := substr(WT_GUIDE_SEQ, 21, 23), GUIDE_SEQ_NO_PAM]
DT.5[, chr := tstrsplit(ID, "-", keep = 2), GUIDE_SEQ_NO_PAM]
DT.5[, coord_start := as.numeric(tstrsplit(ID, "-", keep = 3)), GUIDE_SEQ_NO_PAM]
DT.5[, coord_end := as.numeric(tstrsplit(ID, "-", keep = 4)), GUIDE_SEQ_NO_PAM]

# fwrite(DT.5[, .(gene, GUIDE_START, GUIDE_END, GUIDE_STRAND, GUIDE_SEQ, PAM, Design_PAM, chr, regions_start, coord_start, coord_end, coord_pam_end, FlasFryID, FlashFrySEQ)], 
#               file = "processed_data/guides/D2N_attenuated_ToScore.txt", quote = F, row.names = F, sep="\t")

## Chose the best attenuated guides depending on the position of the mutation and the difference in the guide scores vs the WT
DT.6 <- merge(fread(input = "processed_data/guides/D2N_attenuated_Scored.txt")[, .(ID, Doench2014OnTarget, otCount, Hsu2013)], 
                      DT.5, by.x = "ID", by.y = "FlasFryID")

# Check the fold-change on the on-target and off-target activity after the missmatch
DT.6[, fc_Doench2014OnTarget := Doench2014OnTarget/Doench2014OnTarget_WT-1, GUIDE_SEQ]
DT.6[, fc_otCount := otCount/otCount_WT-1, GUIDE_SEQ]

ggplot(DT.6, aes(x= fc_otCount, y=fc_Doench2014OnTarget, color=gene)) + geom_point() +
  ggtitle(paste0("r = ", round(cor(DT.6$fc_Doench2014OnTarget, DT.6$fc_otCount), 2)))
ggplot(DT.6, aes(x= otCount, y=Hsu2013, color=gene)) + geom_point() +
  ggtitle(paste0("r = ", round(cor(DT.6$fc_Doench2014OnTarget, DT.6$fc_otCount), 2)))

# Categorize each attenuated gRNA by region
DT.6[, region_id := ifelse(MM_POS >= 18, "MM_18:20", ifelse(MM_POS >= 14, "MM_14:17", ifelse(MM_POS >= 10, "MM_10:13", ifelse(MM_POS >= 6, "MM_6:9", "MM_2:5") ) ) ), GUIDE_SEQ_NO_PAM]
DT.6[, gene_region_id := paste0(gene, "_", region_id), GUIDE_SEQ ]

ggplot(DT.6, aes(x= fc_otCount, y=fc_Doench2014OnTarget, color=region_id)) + geom_point() +
  ggtitle(paste0("r = ", round(cor(DT.6$fc_Doench2014OnTarget, DT.6$fc_otCount), 2))) +
  facet_wrap(. ~ gene, scales = "free") +
  geom_vline(xintercept = 0, lty=2) +
  geom_hline(yintercept = 0, lty=2)


# Select that guide with bigest decrease in on-target activity but mantaining a maximum 0.5 FC on off-targets
rid_vec <- unique(DT.6$gene_region_id)
att_guides_l <- list()
for (n in 1:length(rid_vec)) {
  dt <- DT.6[gene_region_id == rid_vec[n] & fc_Doench2014OnTarget < 0 & abs(fc_otCount) <= 0.5,]
  dt_ord <- dt[base::order(fc_Doench2014OnTarget),]
  att_guides_l[[n]] <- dt_ord[1,]
}
att_selguides_dt <- rbindlist(att_guides_l)
DT.6[, selected := GUIDE_SEQ %in% att_selguides_dt$GUIDE_SEQ]


DT.6_selected <- DT.6[selected == T, ]

ggplot(DT.6_selected, aes(x=fc_Doench2014OnTarget, fill=gene)) +
  geom_histogram()

# Add features and write final table to collapse with all other guides
DT.6_selected[, id := paste0(gene, "_", chr, ":", coord_start, "-", coord_end, "_", GUIDE_STRAND), GUIDE_SEQ]
DT.6_selected[, guide_category := "attenuated"]



# Write final table
fwrite(DT.6_selected[, .(id, coord_start, coord_end, GUIDE_SEQ, GUIDE_STRAND, GUIDE_SEQ_NO_PAM, PAM, Doench2014OnTarget, Hsu2013, otCount, guide_category)], 
       file = "processed_data/guides/D2N_attenuated_Scored_Selected.txt", quote = F, row.names = F)


DT.6_selected[gene == "TET2", .(id, fc_Doench2014OnTarget, MM_MUT_ID)]









