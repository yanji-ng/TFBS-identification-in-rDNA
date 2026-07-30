library(data.table)
library(dplyr)
library(tidyr)
library(Rsamtools)
library(GenomicAlignments)
library(stringr)
library(ggplot2)

setwd("./_rdna_chip_ubtf")

###################################################
### Collect data
###################################################

ref_df <- as_tibble(fread("./analysis_20260701/haplotype_pos_mef.csv"))
ref_df <- ref_df %>% left_join(
    ref_df %>% filter(!haplotype %in% c("other", "A", "C")) %>% dplyr::rename(haplotype_group = haplotype) %>% 
    dplyr::select(-allele), by = c("pos_ref"))
ref_df <- ref_df %>% mutate(haplotype_group = ifelse(is.na(haplotype_group), "A/C", haplotype_group))
ref_df <- ref_df %>% mutate(pos_loop = ifelse(pos_ref >= 1, pos_ref + 3008, pos_ref + 1 + 3008))

list_samples <- list.files(full = FALSE, path = "./bam_BK000964_3_looped_3008") 
list_samples <- list_samples[!grepl(".bai", list_samples)]
list_samples <- gsub(".dedup.BK000964_3_looped_3008.bam$", "", list_samples)

result_df <- tibble(haplotype_group = as.character(), haplotype = as.character())
for (i in 1:length(list_samples)){
    pick_sample <- list_samples[i]
    print(paste0("Processing sample ", pick_sample, "."))
    path_file <- paste0(
        "./bam_BK000964_3_looped_3008/", 
        pick_sample, 
        ".dedup.BK000964_3_looped_3008.bam")
     
    bam <- BamFile(path_file)
    aln <- readGAlignments(bam, use.names = FALSE, param = ScanBamParam(
        what = c("qname", "mapq", "seq"),
        flag = scanBamFlag(
            isUnmappedQuery = FALSE,
            isSecondaryAlignment = FALSE,
            isDuplicate = FALSE,
            isNotPassingQualityControls = FALSE)))
    aln <- as_tibble(aln)
    aln <- aln %>% dplyr::select(qname, start, end, cigar, seq)

    pos_loop <- ref_df %>% pull(pos_loop) %>% unique()
    aln_mod <- aln %>% rowwise() %>% mutate(
        pos_loop = list(pos_loop[pos_loop >= start & pos_loop <= end])) %>%
        ungroup() %>% unnest(pos_loop)
    
    aln_mod <- aln_mod %>% mutate(seq = as.character(sequenceLayer(
        DNAStringSet(aln_mod$seq), aln_mod$cigar)))
    aln_mod <- aln_mod %>% mutate(seq = str_sub(seq, pos_loop - start + 1, pos_loop - start + 1))
    aln_mod <- aln_mod %>% dplyr::select(qname, seq, pos_loop) %>% left_join(ref_df, by = "pos_loop")
  
    aln_mod <- aln_mod %>% filter(seq == allele) %>% group_by(haplotype_group, haplotype) %>% summarize(num_reads = length(haplotype), .groups = "drop")
    aln_mod <- aln_mod %>% group_by(haplotype_group) %>% summarize(haplotype = haplotype, num_reads = num_reads, freq_reads = round(num_reads/sum(num_reads), 3), .groups = "drop") 
    
    result_df <- result_df %>% full_join(aln_mod %>% dplyr::select(-num_reads) %>% dplyr::rename(!!pick_sample := freq_reads), by = c("haplotype_group", "haplotype"))}
    
result_df <- result_df %>% dplyr::rename(rep1 = SRR8279916, rep2 = SRR8279917, input = SRR8279918) 
result_df <- result_df %>% mutate(rep1_norm = rep1/input, rep2_norm = rep2/input)
  
###################################################
### Figures
###################################################

data_gg <- result_df %>% dplyr::select(-c(haplotype_group, rep1, rep2, input)) %>% dplyr::rename('1' = rep1_norm, '2' = rep2_norm) %>% tidyr::pivot_longer(-haplotype, names_to = "rep", values_to = "freq")
data_gg <- data_gg %>% filter(!haplotype %in% c("other", "A", "C")) %>% mutate(haplotype = gsub("_", "", haplotype), haplotype = sub("^(.).", "\\1", haplotype))

plot_gg <- data_gg %>% ggplot() +
    geom_point(aes(x = factor(haplotype), y = freq, shape = rep, col = haplotype),     
        position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.8), size = 3) +
    scale_shape_manual(values = c("1" = 16, "2" = 15)) +
labs(x = "", y = "Relative Frequency in ChIP-Seq \n UBTF vs Input", shape = "Replicate") +  guides(col = "none") +
theme_bw() + theme( 
    panel.grid.minor = element_blank(), 
    plot.margin = margin(20, 20, 20, 20), 
    plot.title = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(size = 7),
        #axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 7),
    axis.title = element_text(size = 7, face = "bold"),
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 7, face = "bold"),
    legend.key.size = unit(0.5, "lines"),
    strip.text = element_text(size = 7))
ggsave("./analysis_20260701/20260709-panel-01-ubtf-dedup-fran-method.pdf", plot = plot_gg, width = 3, height = 2.8)

###################################################
### Download data
###################################################

path_from="./_rdna_chip_ubtf/analysis_20260701/*.pdf"
path_to="/Users/yzk/Desktop/TFBS-Analysis"
scp -i ~/.ssh/create_msc k21218585@hpc.create.kcl.ac.uk:"${path_from}" "${path_to}"