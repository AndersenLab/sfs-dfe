#!/usr/bin/env Rscript

library(tidyverse)
library(glue)

multi_dfe_out <- function(df, fname) {
  writeLines(
    c(
      nrow(df) - 1,
      paste(df$Selected, collapse =" "),
      paste(df$Neutral, collapse =" ")
    ),
    con = glue::glue(fname), 
    sep = "\n"
  )
}

# writing function
process_sfs <- function(neut, sel) {
  
  # c_regions <- gsub(".tsv", "", strsplit(neut, split = "_")[[1]][4])
  n_deets <- strsplit(neut, split = "_")[[1]][3]
  s_deets <- gsub(".tsv", "",strsplit(sel, split = "_SELECTED_")[[1]][2])
  
  save_name <- glue::glue("{n_deets}_{s_deets}")
  
  neutral_sites <- data.table::fread(neut, col.names = c("CHROM","POS","AA","GT","COUNT","AA_CLASS")) %>%
    dplyr::group_by(CHROM,POS) %>%
    dplyr::mutate(DERIVED_AF = ifelse(AA_CLASS == "DERIVED", COUNT/sum(COUNT), NA)) %>%
    na.omit() %>%
    dplyr::group_by(DERIVED_AF) %>% 
    dplyr::summarise(Neutral = n())
  
  selected_sites <- data.table::fread(sel, col.names = c("CHROM","POS","AA","GT","COUNT","AA_CLASS")) %>%
    dplyr::group_by(CHROM,POS) %>%
    dplyr::mutate(DERIVED_AF = ifelse(AA_CLASS == "DERIVED", COUNT/sum(COUNT), NA)) %>%
    na.omit() %>%
    dplyr::group_by(DERIVED_AF) %>% 
    dplyr::summarise(Selected = n())
  
  sfs_df <- dplyr::left_join(neutral_sites, selected_sites , by = "DERIVED_AF")
  
  # Replace NAs with 0
  sfs_df[is.na(sfs_df)] <- 0
  
  # save file
  multi_dfe_out(df = sfs_df,
                fname = glue::glue("{save_name}.sfs"))
  
  # plot
  sfs_plot <- sfs_df%>%
    tidyr::gather(SITES, COUNTS, -DERIVED_AF) %>%
    dplyr::group_by(SITES) %>%
    dplyr::mutate(frq = COUNTS/sum(COUNTS)) %>%
    ggplot()+
    aes(x = DERIVED_AF, y = frq, color = SITES)+
    geom_line()+
    theme_bw(15)+
    labs(x = "Derived AF", y = "Frequency")
 
  ggsave(sfs_plot, filename = glue::glue("{save_name}.pdf"), height = 6, width = 8) 
}

# prep file names
neutral <- grep("NEUTRAL", list.files(), value = T)
neutral_arms <- grep("ARMS", neutral, value = T)
neutral_centers <- grep("CENTERS", neutral, value = T)
neutral_genome <- grep("GENOME", neutral, value = T)

selected <- grep("SELECTED", list.files(), value = T)
selected_arms <- grep("ARMS", selected, value = T)
selected_centers <- grep("CENTERS", selected, value = T)
selected_genome <- grep("GENOME", selected, value = T)


for(nt_f in 1:length(neutral_genome)) {
  for(st_f in 1:length(selected_genome)) {
    process_sfs(neutral_genome[nt_f], selected_genome[st_f]) 
    process_sfs(neutral_arms[nt_f], selected_arms[st_f]) 
    process_sfs(neutral_centers[nt_f], selected_centers[st_f]) 
  }
}




