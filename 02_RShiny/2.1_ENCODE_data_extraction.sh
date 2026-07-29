
library(dplyr)
library(tidyr)

setwd("./_expression_atlas")

ref_df <- as_tibble(read.table("encode_transcription_mouse.tsv", 
    skip = 1, sep = "\t", header = TRUE))

#colnames(ref_df)
#ref_df %>% pull("Biosample.class") %>% table()

# black 6 mouse strain
ref_df <- ref_df %>% filter(
    grepl("C57BL", Simple.Biosample.summary) &
    !grepl("CAST", Simple.Biosample.summary) &
    !grepl("EiJ", Simple.Biosample.summary))
ref_df <- ref_df %>% filter(!Simple.Biosample.summary %in% c(
    "strain Patski M.spretus x C57BL/6J", 
    "strain Bruce4 C57BL/6"))

# untreated
ref_df <- ref_df %>% filter(Biosample.treatment == "") 

# polyA/long-read/total RNA-seq (remove microRNA-seq)
ref_df <- ref_df %>% filter(Assay.name != "microRNA-seq") 

# mm10 genome assembly
ref_df <- ref_df %>% filter(Genome.assembly == "mm10") 

# only cell type/tissue
ref_df <- ref_df %>% filter(Biosample.classification != "cell line") 

# make shorter version (141 samples only)

ref_df <- ref_df %>% 
    dplyr::select(Biosample.accession, Assay.name, Description, Lab, Biosample.classification, Biosample.term.name, Life.stage, Biosample.age, Simple.Biosample.summary, Files) %>% 
    dplyr::rename(Accession = Biosample.accession, Assay = Assay.name, Biosample.class = Biosample.classification, Biosample.term = "Biosample.term.name", Stage = Life.stage, Age = Biosample.age, Strain = Simple.Biosample.summary)

ref_df <- ref_df %>% mutate(Files = strsplit(Files, ","), .groups = "keep") %>% unnest(Files)
ref_df <- ref_df %>% mutate(Files = gsub("/files/", "", Files), Files = gsub("/", "", Files))

link_df <- ref_df %>% dplyr::select(Files) %>% rowwise() %>% mutate(Link = paste0(
    "https://www.encodeproject.org/files/", Files, "/@@download/", Files, ".tsv")) %>% dplyr::select(Link)

write.table(ref_df, "encode_transcription_mouse_filtered.tsv", 
    sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
write.table(link_df, "encode_transcription_mouse_filtered_links.tsv", 
    sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

################################################# Bash

cd /scratch/prj/mmg_holland_sandbox/Mila/_expression_atlas
mkdir -p data_raw

for LINK in $(cat encode_transcription_mouse_filtered_links.tsv); do
    BASE=$(basename ${LINK})

    if ls ./data_raw | grep -q -e ${BASE}; then continue; fi
    wget "${LINK}" -P ./data_raw/ 
done





