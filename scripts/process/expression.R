try(setwd(dirname(dirname(dirname(rstudioapi::getActiveDocumentContext()$path)))))
library(tidyverse)

expr <- readr::read_tsv("expression/expression.txt.gz", skip = 5) %>%
  dplyr::select(
    "gene" = `### Gene`,
    "expr_count" = `WBls:0000002 count`,
    "expr_median" = `WBls:0000002 median`,
    "expr_mean" = `WBls:0000002 mean`
  ) %>%
  dplyr::filter(expr_mean > 1e-10) %>%
  dplyr::mutate(p_rank = dplyr::percent_rank(expr_mean))

expr %>%
  dplyr::filter(p_rank <= 0.25) %>% 
  dplyr::select(gene) %>%
  readr::write_tsv("expression/q1.expression.tsv", col_names=FALSE)

expr %>%
  dplyr::filter(p_rank > 0.25 & p_rank <= 0.5) %>% 
  dplyr::select(gene) %>%
  readr::write_tsv("expression/q2.expression.tsv", col_names=FALSE)

expr %>%
  dplyr::filter(p_rank > 0.50 & p_rank <= 0.75) %>% 
  dplyr::select(gene) %>%
  readr::write_tsv("expression/q3.expression.tsv", col_names=FALSE)

expr %>%
  dplyr::filter(p_rank > 0.75) %>% 
  dplyr::select(gene) %>%
  readr::write_tsv("expression/q4.expression.tsv", col_names=FALSE)
