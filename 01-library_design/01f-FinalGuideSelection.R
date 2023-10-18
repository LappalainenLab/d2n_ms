##### D2N - Library design - Join all designed guides and modify sequence to order #####
# 1) Join all categories of designed guides into a single table
# 2) Identify mismatches (double check that attenuated mismatches are the same as designed)
# 3) Remove, add or substitute any additional guide
# 4) Write final table, include NTC guides, and export the IDT-format ordering table

# JDE, January 2022

## Libraries
library(data.table)
library(rtracklayer)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)
library(Biostrings)
library(ggplot2); theme_set(theme_bw())
library(reshape2)


setwd("/gpfs/commons/groups/lappalainen_lab/jdomingo/projects/004-dosage_network/001-library_design/")


##### (1) Load and cbind all selected guide tables
GDT.1 <- rbind(fread(input = "processed_data/guides/D2N_titration_Scored_Selected.txt"),
               fread(input = "processed_data/guides/D2N_weismanTSS_Scored_Selected.txt"),
               fread(input = "processed_data/guides/D2N_distalCRE_Scored_Selected.txt"),
               fread(input = "processed_data/guides/D2N_attenuated_Scored_Selected.txt"))
GDT.1[, gene := tstrsplit(id, "_", keep = 1), GUIDE_SEQ]
GDT.1[, chr := gsub(".*_", "", tstrsplit(id, ":", keep = 1)), GUIDE_SEQ]

##### (2) Add mismatches features
GDT.1[GUIDE_STRAND == "+", genomic_GUIDE_PAM_SEQ := as.character(getSeq(BSgenome.Hsapiens.UCSC.hg38::Hsapiens, names= chr, start = coord_start, end = coord_end +3))]
GDT.1[GUIDE_STRAND == "-", genomic_GUIDE_PAM_SEQ := as.character(reverseComplement(getSeq(BSgenome.Hsapiens.UCSC.hg38::Hsapiens, names= chr, start = coord_end -3, end = coord_start)))]

GDT.1[, MM_NUM := sum(unlist(strsplit(GUIDE_SEQ, "")) != unlist(strsplit(genomic_GUIDE_PAM_SEQ, ""))), GUIDE_SEQ]
GDT.1[, MM_POS := ifelse(MM_NUM <= 1,
                         ifelse(MM_NUM == 0, 0, which(unlist(strsplit(GUIDE_SEQ, "")) != unlist(strsplit(genomic_GUIDE_PAM_SEQ, "")))),
                         paste(which(unlist(strsplit(GUIDE_SEQ, "")) != unlist(strsplit(genomic_GUIDE_PAM_SEQ, ""))), collapse = ",")), GUIDE_SEQ]
GDT.1[, MM_WT_NT := ifelse(MM_NUM <= 1,
                           ifelse(MM_NUM == 0, "-", unlist(strsplit(genomic_GUIDE_PAM_SEQ, ""))[MM_POS]),
                           paste(unlist(strsplit(genomic_GUIDE_PAM_SEQ, ""))[as.numeric(unlist(strsplit(MM_POS, ",")))], collapse = ",")), GUIDE_SEQ]
GDT.1[, MM_MUT_NT := ifelse(MM_NUM <= 1,
                            ifelse(MM_NUM == 0, "-", unlist(strsplit(GUIDE_SEQ, ""))[MM_POS]),
                            paste(unlist(strsplit(GUIDE_SEQ, ""))[as.numeric(unlist(strsplit(MM_POS, ",")))], collapse = ",")), GUIDE_SEQ]
GDT.1[, MM_ID := paste0(MM_WT_NT, MM_POS, MM_MUT_NT), GUIDE_SEQ]
GDT.1[MM_ID == "-0-", MM_ID := "WT", GUIDE_SEQ]
GDT.1[, ID := paste0(id, "_", MM_ID), GUIDE_SEQ]

##### (3) Substitute any guide?
# Substitute "TET2_chr4:105146774-105146755_-_A3C for the WT version of the WeismanTSS mismatched guide

new_guide <- GDT.1[MM_NUM > 0 & guide_category == "weismanTSS", ]
new_guide[, GUIDE_SEQ := genomic_GUIDE_PAM_SEQ]
new_guide[, GUIDE_SEQ_NO_PAM := substr(GUIDE_SEQ, 1, 20)]
new_guide[, MM_NUM := sum(unlist(strsplit(GUIDE_SEQ, "")) != unlist(strsplit(genomic_GUIDE_PAM_SEQ, ""))), GUIDE_SEQ]
new_guide[, MM_POS := ifelse(MM_NUM <= 1,
                         ifelse(MM_NUM == 0, 0, which(unlist(strsplit(GUIDE_SEQ, "")) != unlist(strsplit(genomic_GUIDE_PAM_SEQ, "")))),
                         paste(which(unlist(strsplit(GUIDE_SEQ, "")) != unlist(strsplit(genomic_GUIDE_PAM_SEQ, ""))), collapse = ",")), GUIDE_SEQ]
new_guide[, MM_WT_NT := ifelse(MM_NUM <= 1,
                           ifelse(MM_NUM == 0, "-", unlist(strsplit(genomic_GUIDE_PAM_SEQ, ""))[MM_POS]),
                           paste(unlist(strsplit(genomic_GUIDE_PAM_SEQ, ""))[as.numeric(unlist(strsplit(MM_POS, ",")))], collapse = ",")), GUIDE_SEQ]
new_guide[, MM_MUT_NT := ifelse(MM_NUM <= 1,
                            ifelse(MM_NUM == 0, "-", unlist(strsplit(GUIDE_SEQ, ""))[MM_POS]),
                            paste(unlist(strsplit(GUIDE_SEQ, ""))[as.numeric(unlist(strsplit(MM_POS, ",")))], collapse = ",")), GUIDE_SEQ]
new_guide[, MM_ID := paste0(MM_WT_NT, MM_POS, MM_MUT_NT), GUIDE_SEQ]
new_guide[MM_ID == "-0-", MM_ID := "WT", GUIDE_SEQ]
new_guide[, ID := paste0(id, "_", MM_ID), GUIDE_SEQ]


GDT.2 <- rbind(GDT.1[!(ID %in% "TET2_chr4:105146774-105146755_-_A3C"), ], new_guide)
GDT.2[guide_category == "weismanTSS" & MM_NUM > 0, guide_category := "attenuated",]


##### (4) Write final table, include NTC guides, and export the IDT-format ordering table

# Write final sgRNA table
fwrite(GDT.2, file = "processed_data/guides/D2N_ALLguides_features.txt", quote = F, row.names = F)


# Create IDT table including NTC guides

GDT.3 <- rbind( GDT.2[base::order(ID), .(ID, GUIDE_SEQ_NO_PAM)], 
                setnames(fread("processed_data/guides/D2N_NTC_Selected.txt")[, .(id, GUIDE_SEQ_NO_PAM)], old="id", new="ID") )
# Add Gibson OH regions
Gibson_OH_5 = "GGAAAGGACGAAACACCG"
Gibson_OH_3 = "GTTTAAGAGCTATGCTGGAAAC"

GDT.3[, Sequence := paste0(Gibson_OH_5, GUIDE_SEQ_NO_PAM, Gibson_OH_3), GUIDE_SEQ_NO_PAM]
GDT.3[, "Well Position" := c(paste0("A", 1:12), paste0("B", 1:12), paste0("C", 1:12), paste0("D", 1:12), paste0("E", 1:12), paste0("F", 1:12), paste0("G", 1:12) , paste0("H", 1:12))   ]
GDT.3[, Name := ID, GUIDE_SEQ_NO_PAM]

# Write csv for IDT order
# fwrite(GDT.3[, .(`Well Position`, Name, Sequence)], row.names = F, quote = F, file = "processed_data/guides/D2N_IDT_order.csv")
fwrite(GDT.3[, .(ID, GUIDE_SEQ_NO_PAM)], row.names = F, quote = F, file = "../005-sgRNA_library_distribution/processed_data/d2n_guides.txt", sep = "\t")

length(unique(GDT.3$Sequence))







##### (5)  PLOTS ####

plot_dt <- as.data.table(melt(table(GDT.2[, .(gene, guide_category)]), id.vars = "gene", variable.name = "guide_category", value.name = "count"))
plot_dt[, guide_category := factor(guide_category, levels = rev(c("attenuated", "distalCRE", "titration", "weismanTSS")))]
p1 <- ggplot(plot_dt, aes(x=gene, y=count)) +
  geom_bar(stat = "identity", aes(fill = guide_category)) +
  geom_text(position = "stack", aes(label = count, hjust = 0.5))
p1
ggsave(p1, filename = "plots/20220118_guide_category_barplot.pdf", width = 4.5, height = 6)


plot_dt <- rbind( cbind( setnames(GDT.2[, .(gene, guide_category, Hsu2013)], old= "Hsu2013", new = "value"), data.table(metric = "Hsu2013") ),
                  cbind( setnames(GDT.2[, .(gene, guide_category,otCount)], old= "otCount", new = "value"), data.table(metric = "otCount") ),
                  cbind( setnames(GDT.2[, .(gene, guide_category, Doench2014OnTarget)], old= "Doench2014OnTarget", new = "value"), data.table(metric = "Doench2014OnTarget") ))
plot_dt[, guide_category := factor(guide_category, levels = rev(c("attenuated", "distalCRE", "titration", "weismanTSS")))]
p2 <- ggplot(plot_dt, aes(x=value, fill=guide_category)) + geom_histogram(bins = 50) + 
  theme(axis.title.x = element_blank()) + theme(legend.position = "bottom") +
  facet_grid(gene ~ metric, scales = "free") 
  
p2
ggsave(p2, filename = "plots/20220118_allguides_properties_dist.pdf", width = 8, height = 5)


GDT.2[guide_category == "attenuated", .(gene, MM_ID)]

