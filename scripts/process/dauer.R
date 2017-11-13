try(setwd(dirname(dirname(dirname(rstudioapi::getActiveDocumentContext()$path)))))

genes <- readxl::read_excel("data/gene_set/dauer_genes.xls", col_names = c("id", "fosmid_id"), skip = 1)
