try(setwd(dirname(dirname(dirname(rstudioapi::getActiveDocumentContext()$path)))))

gene_ids <- readr::read_csv("ftp://ftp.wormbase.org/pub/wormbase/species/c_elegans/PRJNA13758/annotation/geneIDs/c_elegans.PRJNA13758.current.geneIDs.txt.gz", col_names = c("taxon_id", "wbgene", "locus", "fosmid_id", "status"))

gene_ids <- gene_ids %>% dplyr::filter(status == "Live")

genes <- readxl::read_excel("data/gene_set/dauer_genes.xls", col_names = c("id", "fosmid_id"), skip = 1) %>%
         dplyr::left_join(gene_ids) %>%
         dplyr::filter(!is.na(wbgene)) %>%
         dplyr::select(wbgene)

readr::write_tsv(genes, path = "data/gene_set/dauer_genes.txt", col_names = F)
