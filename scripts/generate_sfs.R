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
  # An extra 0 is added to the end of each spectra for 
  writeLines(
    c(
      nrow(df),
      paste(c(df$selected, 0), collapse =" "),
      paste(c(df$neutral, 0), collapse =" ")
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
  spectra(df %>% dplyr::filter(chrom == x)) %>%
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

#=========#
# Biotype #
#=========#

sapply(names(df)[grepl("biotype", names(df))], function(x) {
  x = rlang::sym(x)
  df_out = df %>% dplyr::filter((!!x) == T) %>%
    dplyr::group_by(ancestral_allele_count, add=TRUE) %>%
    dplyr::summarize(selected = sum(fold__4 == F)) %>%
    dplyr::full_join(neutral_spectrum) %>%
    dplyr::arrange(ancestral_allele_count) %>%
    dplyr::mutate(selected = ifelse(is.na(selected), 0, selected))
  if(sum(df_out$selected) > 0) {
    print(glue::glue("data/spectra/{outgroup}/{x}.sfs"))
    df_out %>% multi_dfe_out(fname=glue::glue("data/spectra/{outgroup}/{x}.sfs"))
  }
})

#========#
# Effect #
#========#

sapply(names(df)[grepl("effect", names(df))], function(x) {
  x = rlang::sym(x)
  df_out = df %>% dplyr::filter((!!x) == T) %>%
    dplyr::group_by(ancestral_allele_count, add=TRUE) %>%
    dplyr::summarize(selected = sum(fold__4 == F)) %>%
    dplyr::full_join(neutral_spectrum) %>%
    dplyr::arrange(ancestral_allele_count) %>%
    dplyr::mutate(selected = ifelse(is.na(selected), 0, selected))
  if(sum(df_out$selected) > 0) {
    print(glue::glue("data/spectra/{outgroup}/{x}.sfs"))
    df_out %>% multi_dfe_out(fname=glue::glue("data/spectra/{outgroup}/{x}.sfs"))
  }
})
