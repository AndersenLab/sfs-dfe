try(setwd(dirname(rstudioapi::getActiveDocumentContext()$path)))
setwd(system("git rev-parse --show-toplevel", intern =T))

df <- readr::read_tsv("data/df_outgroup/QX1211.tsv.gz")

# Spectra with restricted regions
spectra <- function(df) {
    df %>%
    dplyr::group_by(ancestral_allele_count, add=TRUE) %>%
    dplyr::summarize(neutral = sum(fold__4), selected = sum(fold__4 == F)) %>%
    dplyr::arrange(ancestral_allele_count)
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