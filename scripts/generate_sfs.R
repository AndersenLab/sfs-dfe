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

#============#
# Chromosome #
#============#
sapply(c("I", "II", "III", "IV", "V", "X"), function(x) {
  outfile=glue::glue("data/spectra/{outgroup}/chromosome_{x}.sfs")
  df %>% dplyr::filter(chrom == x) %>%
    multi_dfe_out(fname=outfile)
})

#=========#
# Operons #
#=========#
spectra(df %>% dplyr::filter(is__gene == T, operon__operon == T)) %>%
        multi_dfe_out(fname="data/spectra/{outgroup}/operon_T.sfs")
spectra(df %>% dplyr::filter(is__gene == T, operon__operon == F)) %>%
        multi_dfe_out(fname="data/spectra/{outgroup}/operon_F.sfs")

#===========#
# Chrom Arm #
#===========#
spectra(df %>% dplyr::filter(aoc__arm == T)) %>%
  multi_dfe_out(fname="data/spectra/{outgroup}/arm_T.sfs")
spectra(df %>% dplyr::filter(aoc__center == T)) %>%
  multi_dfe_out(fname="data/spectra/{outgroup}/center_T.sfs")

#============#
# Expression #
#============#
sapply(c("q1", "q2", "q3", "q4"), function(x) {
  outfile=glue::glue("data/spectra/{outgroup}/expression_{x}.sfs")
  expr_var = rlang::sym(glue::glue("expression__{x}"))
  df %>% dplyr::filter((!!expr_var) == T) %>%
    spectra() %>%
    multi_dfe_out(fname=outfile)
})

#============#
# Dauer Gene #
#============#
spectra(df %>% dplyr::filter(is__gene == T, dauer__gene == T)) %>%
  multi_dfe_out(fname="data/spectra/{outgroup}/dauer_T.sfs")
spectra(df %>% dplyr::filter(is__gene == T, dauer__gene == F)) %>%
  multi_dfe_out(fname="data/spectra/{outgroup}/dauer_F.sfs")

#=====================#
# Genome-wide classes #
#=====================#

# For these, I examine a single genome-wide, neutral spectrum

neutral_spectrum <- df %>%
  dplyr::group_by(ancestral_allele_count, add=TRUE) %>%
  dplyr::summarize(neutral = sum(fold__4)) 

sapply(names(df)[grepl("biotype", names(df))], function(x) {
  x == rlang::sym(x)
  df %>% dplyr::filter((!!x) == T)
}
