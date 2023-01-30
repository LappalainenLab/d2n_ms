##### D2N - Targeted-seq gene selection #####
# 10X Target Gene Expression custom probes panel are designed based on genome annotations (no cell line unique transcript expression)
# (1) Check the read coverage of all probes of these using WT K562 5' scRNA-seq
#     - For those probes in exon-exon junctions, recalculate depth of the exon regions only
# (2) Filter out probes given coverage thresholds
#     - All probes for non-expressed LHX3 gene 
#     - Probes with median cvrg >= 1 for lowly expressed genes
#     - Probes with median cvrg >= 3 for remaining genes


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
library(spatstat)
library(ggplot2)
library(spatstat)
library(doParallel); registerDoParallel(detectCores()/2)
library(foreach)



## Dir
setwd("/gpfs/commons/groups/lappalainen_lab/jdomingo/manuscripts/d2n/d2n_ms/01-library_design/")
folder_name <- basename(getwd())
data_dir = file.path("../../data", folder_name)
processed_data_dir = file.path("../../processed_data", folder_name)


# Import probes bed file and tranform it into GRangesList
glist_baits <- import(file.path(processed_data_dir, "D2N_10X_baits.bed"))
# Add gene name to the GRanges objects
glist_baits$gene_name <- sapply(glist_baits$name, function(x){
  unlist(strsplit(x, "\\|"))[2]
})

# Convert them into compressed GRangesList grouped by target gene
glist_baits <- split(glist_baits, f=glist_baits$gene_name)

# List of target gene symbols
target_gene_symbol <- names(glist_baits)

# BAM file path - CellRanger output from WT K562 ECCITE-seq experiment
bam_file = file.path(data_dir, "K562_WT_scRNAseq5cap_possorted_genome_bam.bam")

# Get genomic coordinates of all known genes
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
# Keep only the coordinates of exons, spliting into different GRanges for each individual gene (gene names are Entrez IDs)
exons_by_gene <- exonsBy(txdb, by = "gene")

# Get Entrez IDs for each target gene (manually change CCDC173 to CFAP210 when retrieving Entrez IDs, since external gene name is not recognized)
symbols_entrezid <- as.data.table(select(org.Hs.eg.db, 
                                         keys = gsub("CCDC173", "CFAP210", target_gene_symbol), 
                                         columns = c("ENTREZID", "SYMBOL"), 
                                         keytype = "SYMBOL"))
symbols_entrezid[, SYMBOL := gsub("CFAP210", "CCDC173", SYMBOL)]



### (1) Check the read coverage of all probes of these using WT K562 5' scRNA-seq
# For those probes that cover a single exon, use probe coordinates to calculate read coverage
# For those probes that span multiple exons, find median coverage of each exon and combine into a single coverage value by doing a weighted median

# Function to find median read coverage and number of reads that overlap in each probe - run in parallel for each unique gene
MedianCvrgPerProbe <- function(gr, gene_enterz_dt, bam_file_path, grl_exons_by_gene){
  
  # Identify gene name and chromosome
  gene_name <- as.character(unique(gr$gene_name))
  message(gene_name)
  gene_enterz <- gene_enterz_dt[SYMBOL == gene_name, ENTREZID]
  chr <- as.character(unique(seqnames(gr)))
  
  # Split probes into two GRangeLists - those that map to a single exon and those that span multiple exons
  gr_se <- subset(gr, sapply(gr$blocks, length) == 1)
  gr_me <- subset(gr, sapply(gr$blocks, length) > 1)
  
  # Split the GRanges by probe id
  gr_se_split <- split(gr_se, f=gr_se$name)
  gr_me_split <- split(gr_me, f=gr_me$name)
  
  # (1) For those single exon probes, get median coverage using the entire gene GRanges (max 120b)
  if (length(gr_se_split) > 0) {
    # ScanBamParam creates a parameter object influencing what fields and which records are imported from a BAM file
    baits_se_range <- ScanBamParam(which=unlist(gr_se_split))
    
    # Get read alignments (with coordinates - compatible with GRanges) for those genomic regions
    se_reads <- GenomicAlignments::readGAlignments(bam_file_path, param=baits_se_range)
    
    # Calculate coverage per base - creates list of Rle objects, one for each chromosome
    # Only the chromosome of that gene should have read coverage info
    se_read_cvrg <- GenomicAlignments::coverage(se_reads)
    
    # Output mean coverage per probe in a data.table format
    dt_se <- data.table(gene_name = gene_name,
                        bait_id = names(gr_se_split),
                        median_cvrg =  aggregate(se_read_cvrg[[chr]], ranges(unlist(gr_se_split)), FUN = median),
                        # read_overlap = countOverlaps(gr_se_split, reads, maxgap = 100),
                        coord_width = unlist(width(gr_se_split)),
                        overlaping_exons = 1)
  } else {
    dt_se <- data.table()
  }
  
  # (2) For probes that fall into different exons, get the weighted median coverage accross exons
  if (length(gr_me_split) > 0) {
    
    # Get exons for that specific target gene
    exons_gr <- grl_exons_by_gene[[gene_enterz]]
    exons_gr <- keepStandardChromosomes(exons_gr) # Removes non-sandard chromosome names that do not match the BAM file header
    
    # Iterate over each probe's exon to 
    cvrg_dt <- do.call("rbind", lapply(gr_me_split, MedianCvrgAcrossExons, exons_gr, bam_file_path, chr))
    
    # Output mean coverage per probe in a data.table format
    dt_me <- data.table(gene_name = gene_name,
                        bait_id = names(gr_me_split),
                        median_cvrg = cvrg_dt$median_cvrg,
                        coord_width = unlist(width(gr_me_split)),
                        overlaping_exons = cvrg_dt$overlaping_exons)
  } else {
    dt_me <- data.table()
  }
  
  dt <- rbindlist(list(dt_se, dt_me))
  return(dt)
}

# pgr = gr_me_split[[1]]

MedianCvrgAcrossExons <- function(pgr, exons_gr, bam_file_path, chr){
  
  # Get the intersect between the probe coordinates and gene exons coordinates
  probe_exon_gr <- GenomicRanges::intersect(pgr, exons_gr)
  
  ex_len <- length(probe_exon_gr)
  
  if (ex_len  == 0) {
    
    mex_cvrg = 0
    
  } else {
    # Define regions to scan the BAM file
    probe_exon_range <- ScanBamParam(which=probe_exon_gr)
    
    # Get read alignments (with coordinates - compatible with GRanges) for those genomic regions
    probe_exon_reads <- GenomicAlignments::readGAlignments(bam_file_path, param=probe_exon_range)
    
    # Only the chromosome of that gene should have read coverage info
    read_cvrg <- GenomicAlignments::coverage(probe_exon_reads)
    
    # Return the weighted median coverage accoss exons (weiths are probe-exon overlap)
    mex_cvrg = weighted.median(x = aggregate(read_cvrg[[chr]], ranges(probe_exon_gr), FUN = median),
                                         w = width(probe_exon_gr))
  }
  
  dt <- data.table(median_cvrg = mex_cvrg,
             overlaping_exons = ex_len)
  return(dt)
  
}


# Rbind into a single data.table

CVRG <- foreach(x = c(1:length(glist_baits)), .combine = rbind) %dopar% {
  MedianCvrgPerProbe(gr = glist_baits[[x]], gene_enterz_dt = symbols_entrezid, bam_file_path = bam_file, grl_exons_by_gene = exons_by_gene)
}
fwrite(file = file.path(processed_data_dir, "D2N_cvrg_per_10X_bait.txt"), x = CVRG, sep="\t", quote = F, row.names = F)



### (2) Filter out probes given coverage thresholds

cvrg_thrs = 3

# Num probes
sum(CVRG$median_cvrg >= cvrg_thrs)
# Total cost
sum(CVRG$median_cvrg >= cvrg_thrs)*0.75

# Add variable to CVRG
CVRG[, pass_cvrg_thrs := ifelse(median_cvrg >= cvrg_thrs, "include", "exclude")]

# Stats per gene
CVRG_stats_gene_melt <- CVRG[, .N, by = c("gene_name", "pass_cvrg_thrs")]
CVRG_stats_gene <- merge.data.table(dcast.data.table(CVRG_stats_gene_melt, gene_name ~ pass_cvrg_thrs, value.var = "N"), 
                                    CVRG[, list(median_gene_cvrg = median(median_cvrg), 
                                                mean_gene_cvrg = mean(median_cvrg)), 
                                         gene_name],
                                    by = "gene_name")
CVRG_stats_gene[is.na(exclude), exclude := 0]
CVRG_stats_gene[is.na(include), include := 0]
CVRG_stats_gene[, pct_excluded := exclude*100/(exclude + include), gene_name]

genes_lower_thrs <- CVRG_stats_gene[pct_excluded >= 50, gene_name]
cvrg_thrs_low = 1
genes_all_probes <- c("KLK1", "LHX3", "CHST3")

CVRG[gene_name %in% genes_lower_thrs, pass_cvrg_thrs := ifelse(median_cvrg >= cvrg_thrs_low, "include", "exclude") ]
CVRG[gene_name %in% genes_all_probes, pass_cvrg_thrs := "include" ]

CVRG_stats_gene_final <- dcast.data.table(CVRG[, .N, by = c("gene_name", "pass_cvrg_thrs")], gene_name ~ pass_cvrg_thrs, value.var = "N")
CVRG_stats_gene_final[is.na(exclude), exclude := 0]
CVRG_stats_gene_final[is.na(include), include := 0]
CVRG_stats_gene_final[, pct_excluded := exclude*100/(exclude + include), gene_name]


# Final Num probes 
sum(CVRG$pass_cvrg_thrs == "include")
# Final Total cost
sum(CVRG$pass_cvrg_thrs == "include")*0.75




############## Plots ################


# Plot distribution of median coverage per probe
p <- ggplot(CVRG, aes(x=median_cvrg)) +
  geom_histogram(bins = length(unique(CVRG$median_cvrg))/10) +
  theme_classic() +
  scale_x_log10() +
  xlab("Median read coverage per probe") + 
  geom_vline(xintercept = c(cvrg_thrs_low, cvrg_thrs), lty=2, col=c("blue", "red") )
ggsave(file.path(plots_dir, "01b_01_10Xbaits_MedianCvrgDist.pdf"), p, width = 6, height = 4)


# Plot number of included and excluded probes per target gene
plot_dt <- melt.data.table(CVRG_stats_gene_final[, .(gene_name, include, exclude)], id.vars = "gene_name", variable.name = "Probe_type", value.name = "N")
plot_dt[, Probe_type := factor(Probe_type, levels = c("exclude", "include"))]
p <- ggplot(plot_dt, aes(x=gene_name, y=N, color=Probe_type)) +
  geom_point() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(legend.position = c(0.2, 0.8)) +
  theme(axis.title.x = element_blank()) +
  ylab("# Probes")
p
ggsave(file.path(plots_dir, "01b_02_10Xbaits_NumIncludedExcludedPerGene.pdf"), p, width = 12, height = 5)



######### PROBES TABLES #########

# Get probes sequences
baits_seqs <- fread(file.path(processed_data_dir, "D2N_10X_baits_seqs.csv"), sep = ",", header = T)

CVRG_seqs <- merge.data.table(CVRG, baits_seqs, by = "bait_id")

baits_seqs_included <- CVRG_seqs[pass_cvrg_thrs == "include", .(gene_name, bait_id, bait_seq)]
fwrite(baits_seqs_included, sep = ",", file = file.path(processed_data_dir, "D2N_10X_baits_seqs_included.csv"), quote = F, row.names = F)


