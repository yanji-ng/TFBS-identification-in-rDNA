
#srun -p cpu --ntasks 1 --nodes 1 --mem 50G --time 0-04:00 --pty /bin/bash -l
#conda activate cookbook_v02_env
#cd /scratch/prj/mmg_holland_sandbox/Mila/_expression_atlas

library(data.table)
library(dplyr)
library(tidyr)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(rtracklayer)

setwd("/scratch/prj/mmg_holland_sandbox/Mila/_expression_atlas")

################################################# 
### Combine all expression profiles
#################################################

path_list <- list.files(path = "./data_raw", full = TRUE)

count_df <- tibble(gene_id = as.character())
counter <- 0
path_tsv <- path_list[1]

for (path_tsv in path_list){
    counter <- counter + 1
    print(counter)

    sample_id <- gsub(".tsv", "", basename(path_tsv))
    tmp_df <- as_tibble(fread(path_tsv))
        n_rows <- tmp_df %>% nrow()
        if (n_rows == 0){next}
    
    n_transcripts <- tmp_df %>% filter(if_any(where(is.character), ~ grepl("ENSMUST", .))) %>% nrow()
    if (n_transcripts != 0){

        tmp_gene_id <- tmp_df %>% dplyr::select(where(is.character) & where(~ any(grepl("ENSMUSG", .))))
        tmp_transcript_id <- tmp_df %>% dplyr::select(where(is.character) & where(~ any(grepl("ENSMUST", .))))
        if (colnames(tmp_gene_id) == colnames(tmp_transcript_id)){
            tmp_gene_id <- tibble(gene_id = gsub(".*(ENSMUSG[^|]+).*", "\\1", tmp_gene_id %>% pull(1)))
            tmp_transcript_id <- tibble(transcript_id = gsub(".*(ENSMUST[^|]+).*", "\\1", tmp_transcript_id %>% pull(1)))}
        if (ncol(tmp_gene_id) == 0){
            transcript2gene <- mapIds(org.Mm.eg.db, 
                keys = tmp_transcript_id %>% pull(1) %>% unlist(), 
                column = "ENSEMBL", keytype = "ENSEMBLTRANS", multiVals = "first")
            tmp_transcript_id <- tibble(gene_id = transcript2gene)}
                
        tmp_count <- tmp_df %>% dplyr::select(where(is.numeric) & (contains("expected_count") | contains("est_counts") | contains("ENCSR")))
            if (ncol(tmp_count) == 0 & ncol(tmp_df) == 2){tmp_count <- tmp_df %>% dplyr::select(where(is.numeric))}
        
        tmp_mod_df <- tmp_gene_id %>% cbind(tmp_count) %>% as_tibble() %>% setNames(c("gene_id", "count"))
            if (ncol(tmp_mod_df) != 2){print("error")}
        tmp_mod_df <- tmp_mod_df %>% group_by(gene_id) %>% summarize(count = sum(count, na.rm = TRUE), .groups = "keep") %>% ungroup()
    }        
        
    if (n_transcripts == 0){
        n_genes <- tmp_df %>% filter(if_any(where(is.character), ~ grepl("ENSMUSG", .))) %>% nrow()
            if (n_genes == 0){next}
        
        tmp_gene_id <- tmp_df %>% dplyr::select(where(is.character) & where(~ any(grepl("ENSMUSG", .))))
        tmp_count <- tmp_df %>% dplyr::select(where(is.numeric) & (contains("count") | contains(sample_id)))
            if (ncol(tmp_count) == 0 & ncol(tmp_df) == 2){tmp_count <- tmp_df %>% dplyr::select(where(is.numeric))}
        
        tmp_mod_df <- tmp_gene_id %>% cbind(tmp_count) %>% as_tibble() %>% setNames(c("gene_id", "count"))
            if (ncol(tmp_mod_df) != 2){print("error")}
    }
                   
    count_df <- count_df %>% full_join(
        tmp_mod_df, by = "gene_id") %>% dplyr::rename(!!sample_id := count)}

#save(count_df, file = "./data_shiny/input.RData")

################################################# 
### Filter profiles
#################################################

load("./data_shiny/input.RData")

count_filter <- count_df %>% arrange(gene_id)

# rm novel transcripts, probably from long-read RNA
count_filter <- count_filter %>% filter(!grepl("ENC", gene_id)) 

# rm spike-ins, as just part of samples have it
count_filter <- count_filter %>% filter(!grepl("ERCC", gene_id) & !(grepl("Spikein", gene_id))) 

# rm tRNAs (all the above can be simply run with this one line)
count_filter <- count_filter %>% filter(grepl("ENSMUSG", gene_id)) # 75,776 genes

# merge genes with different versions 
count_merge <- count_filter %>% mutate(gene_id = gsub("[.].*", "", gene_id)) 
count_merge <- count_merge %>% group_by(gene_id) %>% summarize(across(where(is.numeric), ~ round(sum(.x, na.rm = TRUE))), .groups = "keep") %>% ungroup()

# save under proper names

data_raw <- count_df
data_count <- count_merge
    data_count <- data_count %>% dplyr::select(-ENCFF606SUS.1)
    
# create ref_gene

ref_gene <- mapIds(org.Mm.eg.db,
  keys = data_count %>% pull(gene_id), 
  column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
ref_gene <- tibble(gene_id = names(ref_gene), gene_name = ref_gene)

# create ref_sample

ref_sample <- as_tibble(read.table("encode_transcription_mouse_filtered.tsv", sep = "\t", header = TRUE)) 
ref_sample <- ref_sample %>% filter(Files %in% colnames(data_count)[-1])

#save(data_raw, data_count, ref_gene, ref_sample, file = "./data_shiny/input.RData")

################################################# 
### Normalise: TPM, FPKM, DESeq2 norm
#################################################

# wget https://www.encodeproject.org/files/gencode.vM21.primary_assembly.annotation_UCSC_names/@@download/gencode.vM21.primary_assembly.annotation_UCSC_names.gtf.gz
# untar gencode.vM21.primary_assembly.annotation_UCSC_names.gtf.gz

gtf <- import("gencode.vM21.primary_assembly.annotation_UCSC_names.gtf") # ENCODE gene/transcript gtf
    gtf <- gtf[gtf$type == "exon"]
    gene_mod <- gsub("[.].*", "", gtf$gene_id)  
    gtf <- unlist(reduce(split(gtf, gene_mod)))
    gtf$gene_id <- names(gtf)
    gtf <- as_tibble(gtf) %>% group_by(gene_id) %>% summarize(width = sum(width, na.rm = TRUE), .groups = "keep") %>% ungroup()

ref_gene <- ref_gene %>% left_join(gtf, by = "gene_id")
ref_gene <- ref_gene %>% arrange(gene_id)

data_count <- data_count %>% arrange(gene_id)
data_tpm <- apply(X = data_count %>% dplyr::select(-gene_id), MAR = 2, 
    FUN = function(q, l = ref_gene %>% pull(width)){
        nom <- q/l
        denom <- sum(nom, na.rm = TRUE)
        res <- 10^6 * nom/denom
        return(res)})
data_tpm <- tibble(gene_id = data_count %>% pull(gene_id)) %>% cbind(as_tibble(data_tpm)) %>% as_tibble()
data_tpm <- data_tpm %>% mutate(across(where(is.numeric), ~ tidyr::replace_na(., 0)))

#save(data_raw, data_count, data_tpm, ref_gene, ref_sample, file = "./data_shiny/input.RData")

#avr_tpm <- data_tpm %>% summarize(across(-gene_id, \(x) mean(x[x > 10], na.rm = TRUE)))
#sd_tpm <- data_tpm %>% summarize(across(-gene_id, \(x) sd(x[x > 10], na.rm = TRUE)))










