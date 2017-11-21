#!/usr/env/Rscript

try(setwd(dirname(rstudioapi::getActiveDocumentContext()$path)))
setwd(system("git rev-parse --show-toplevel", intern=T))


# Life stages
stages = list(
               "WBls:0000002" = "all_stages",
               "WBls:0000003" = 'embryo',
               "WBls:0000032" = 'dauer',
               "WBls:0000024" = 'L1',
               "WBls:0000027" = 'L2',
               "WBls:0000035" = 'L3',
               "WBls:0000038" = 'L4',
               "WBls:0000041" = "adult"
              )

rna_seq = readr::read_tsv("data/expression/RNASeq_FPKM.txt.gz", skip = 5) %>%
          dplyr::rename(`gene` = `### Gene`)

# Fix names
purrr::map2(names(stages),
             stages, 
             function(.x, .y) {
              name_set <- names(rna_seq)[grepl(.x, names(rna_seq))]
              names(rna_seq)[grepl(.x, names(rna_seq))] <<- gsub(.x, .y, name_set)
             }
)

mean_values <- rna_seq %>% 
               dplyr::select(gene, dplyr::contains('mean')) %>%
               tidyr::gather('stage', 'expression', -gene) %>%
               dplyr::filter(!(stage %in% c("Total stages mean", "dauer mean", "all_stages mean"))) %>%
               dplyr::group_by(gene) %>%
               dplyr::mutate(log_10_expr = log10(expression)) %>% 
               dplyr::filter(log_10_expr > -10) %>%
               dplyr::mutate(z = scale(expression),
                             n = n()) %>%
               dplyr::filter(n == 6) %>%
               dplyr::mutate(mean_values, diff = sum(abs(log_10_expr - dplyr::lag(log_10_expr))))
               dplyr::ungroup() %>%
               dplyr::group_by(gene) %>%
               dplyr::filter(n() == 6) %>%
               dplyr::mutate(sd_z = sd(z)) %>%
               dplyr::filter(sd_z < 0.005)
               

ggplot(mean_values %>% dplyr::mutate(l = dplyr::lag(z)) %>% dplyr::summarize(sum_diff = sum(l, na.rm = T))) +
  geom_histogram(aes(x =sum_diff))
ggplot(mean_values) +
  geom_histogram(aes(x = z))
ggplot(mean_values) +
  geom_histogram(aes(x = (sd_z)))+
  xlim(-0.5,0.5)
