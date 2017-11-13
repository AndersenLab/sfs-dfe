df <- readr::read_tsv("out.test")


# Spectra with restricted regions

spectra <- function(df, sfs = TRUE) {
    dout <- df %>%
    dplyr::group_by(ancestral_allele_count, add=TRUE) %>%
    dplyr::summarize(neutral = sum(fold__0), selected = n() - sum(fold__0)) %>%
    dplyr::arrange(ancestral_allele_count)
    
    if (sfs == TRUE) {
      selected <- paste0(dout$selected, collapse = " ")
      neutral <- paste0(dout$neutral, collapse = " ")
      dout <- glue::glue("{selected}\n{neutral}")
    }
    dout
}


operon = spectra(df %>% dplyr::filter(operon__operon == T))

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