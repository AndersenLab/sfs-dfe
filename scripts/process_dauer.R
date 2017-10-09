#!/usr/env/Rscript

library(tidyverse)
try(setwd(dirname(rstudioapi::getActiveDocumentContext()$path)))
setwd(system("git rev-parse --show-toplevel", intern=T))

gene_ids <- readr::read_csv("ftp://ftp.wormbase.org/pub/wormbase/species/c_elegans/annotation/geneIDs/c_elegans.PRJNA13758.WS256.geneIDs.txt.gz",
                            col_names = c("n", "wb_id", "locus", 'ID', 'status'))

df <- readxl::read_excel("data/dauer_genes.xls") %>% 
      dplyr::select(`ID` = 2) %>%
      dplyr::mutate(ID = str_extract(ID, "[^#]+")) %>%
      dplyr::inner_join(gene_ids) 
      
readr::write_tsv(df %>% dplyr::select(wb_id), 'data/dauer_genes.txt', col_names = F)
