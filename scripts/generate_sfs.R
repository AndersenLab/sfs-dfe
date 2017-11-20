try(setwd(dirname(rstudioapi::getActiveDocumentContext()$path)))
setwd(system("git rev-parse --show-toplevel", intern =T))

library(tidyverse)
library(glue)

outgroup = "QX1211"

df <- readr::read_tsv(glue::glue("data/df_outgroup/{outgroup}.tsv.gz"))

# Spectra with restricted regions
spectra <- function(df) {
    df %>%
    dplyr::group_by(ancestral_allele_count, add=TRUE) %>%
    dplyr::summarize(neutral = sum(fold__4), selected = sum(fold__4 == F)) %>%
    dplyr::arrange(ancestral_allele_count)
}

multi_dfe_out <- function(df, fname) {
  writeLines(
    c(
      nrow(df),
      paste(df$selected, collapse =" "),
      paste(df$neutral, collapse =" ")
    ),
    con = glue::glue(fname), 
    sep = "\n"
  )
}

#=========#
# Operons #
#=========#
spectra(df %>% dplyr::filter(operon__operon == T)) %>%
        multi_dfe_out(fname="data/spectra/{outgroup}/operon_T.sfs")
spectra(df %>% dplyr::filter(operon__operon == T)) %>%
        multi_dfe_out(fname="data/spectra/{outgroup}/operon_F.sfs")


#===========#
# Chrom Arm #
#===========#
spectra(df %>% dplyr::filter(aoc__arm == T)) %>%
  multi_dfe_out(fname="data/spectra/{outgroup}/operon_T.sfs")
spectra(df %>% dplyr::filter(aoc__center == T)) %>%
  multi_dfe_out(fname="data/spectra/{outgroup}/operon_F.sfs")



df %>% 
  dplyr::group_by(operon__operon) %>%
  dplyr::group_by(ancestral_allele_count) %>%
  dplyr::select(fold__0, fold__4) %>%
  dplyr::summarize(f0 = sum(fold__0), f4 = sum(fold__4)) -> out

df %>% 
  dplyr::group_by(ancestral_allele_count, operon__operon) %>%
  dplyr::select(ancestral_allele_count, operon__operon, fold__0, fold__4) %>%
  dplyr::summarize(f0 = sum(fold__0), f4 = sum(fold__4)) %>% View()



# Genome-wide spectra
df %>% dplyr::filter()