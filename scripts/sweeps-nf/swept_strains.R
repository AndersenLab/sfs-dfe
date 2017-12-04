# load packages required for processing
library(data.table)
library(plyr)
library(tidyverse)
library(igraph)
library(dplyr)
library(ggplot2)

try(setwd(dirname(rstudioapi::getActiveDocumentContext()$path)))

#works, returns an intervals object, output goes to make.graphs
make.breaks <- function(matchset, remove_small, size) {
  match = matchset  
  starts <- match$start
  stops <-  match$stop
  breakpoints <- unique(c(starts,stops))
  breakpoints <- breakpoints[order(breakpoints)]
  rownames(breakpoints)=NULL
  shifted <- breakpoints[-1]
  breakpoints <- breakpoints[-length(breakpoints)]
  startstop <- data.frame(start=breakpoints, stop=shifted)
  if(remove_small == TRUE) {
    print("removing small haplotypes")
    startstop <- startstop[(startstop$stop - startstop$start) > size,]
  }
  rownames(startstop) = NULL  
  return(startstop)
}

make.graphs <- function(intervals, matches, strainset){
  test1 = any(order(intervals$start) != 1:length(intervals$start))
  test2 = any(order(matches$start) != 1:length(matches$start))
  if (test1 | test2){
    stop("intervals and matches must be sorted by start position, fool!")
  }
  starts = intervals$start
  stops = intervals$stop
  match_starts <- matches$start
  match_stops <- matches$stop
  #segment_graphs = list()
  segment_table_full <- c()
  segment_membership_full <- c()
  HH_full <- c()
  search_start = 1
  for (i in 1:length(starts)){
    start_point <- starts[i]
    stop_point <- stops[i]
    matched <- c()
    for (j in search_start:length(match_starts)){
      if (match_starts[j] <= start_point & match_stops[j] >= stop_point){
        matched <- c(matched, j)
        #print (matched)
      }
      else if (match_starts[j] > start_point){
        
        break
      }
    }
    
    segment_graph = graph.edgelist(matrix(c(as.character(matches$strain1[matched]), as.character(matches$strain2[matched])), ncol = 2), directed = F)
    segment_strains = V(segment_graph)$name
    missing_strains = !(strainset %in% segment_strains)
    segment_graph <- add.vertices(segment_graph, nv = sum(missing_strains), name = strainset[missing_strains])
    #segment_graphs[[i]] <- add.vertices(segment_graph, nv = sum(missing_strains), name = strainset[missing_strains])
    
    seg_clust <- clusters(segment_graph)
    samplesize <- length(strainset)
    haplotype_number <- seg_clust$no
    strain_order = match(strainset, V(segment_graph)$name)
    segment_membership = seg_clust$membership[strain_order]
    names(segment_membership) <- strainset
    sizes <- seg_clust$csize
    segstarts <- rep(start_point, length(sizes))
    segstops <- rep(stop_point, length(sizes))
    segment_table <- cbind(segstarts, segstops, sizes/samplesize)
    HH <- sum((seg_clust$csize/samplesize)^2)
    HH_full <- rbind(HH_full, HH)
    segment_membership_full <- rbind(segment_membership_full, segment_membership)
    segment_table_full <- rbind(segment_table_full, segment_table)
  }
  HH_full <- cbind(starts, stops, HH_full)
  segment_membership_full <- cbind(starts, stops, segment_membership_full)
  return(list(HH = HH_full, membership=segment_membership_full, frequencies=segment_table_full))
}

analyze.segs<- function(matches, strains, remove_small=F, remove_size) {
  full_membership <- c()
  full_frequencies <- c()
  full_HH <- c()
  splitmatch <- split(matches, matches$chrom)
  for(i in 1:length(splitmatch)) {
    print(i)
    sorted <- splitmatch[[i]][order(splitmatch[[i]]$start),]
    startstop <- make.breaks(splitmatch[[i]], remove_small, size = remove_size)
    statistics <- make.graphs(startstop, sorted, strains)
    chrom_key <- rep(i, nrow(statistics$frequencies))
    statistics$frequencies <- cbind(chrom_key, statistics$frequencies)
    colnames(statistics$frequencies) <- c("chromosome", "start", "stop", "frequency")
    chrom_key <- rep(i, nrow(startstop))
    statistics$HH <- cbind(chrom_key, statistics$HH)
    colnames(statistics$HH) <- c("chromosome", "start", "stop", "frequency")
    statistics$membership <- cbind(chrom_key, statistics$membership)
    full_membership <- rbind(full_membership, statistics$membership)
    full_frequencies <- rbind(full_frequencies, statistics$frequencies)
    full_HH <- rbind(full_HH, statistics$HH)
    full_frequencies <- data.frame(full_frequencies)
  }
  full_frequencies <- data.frame(full_frequencies)
  full_frequencies$chromosome <- factor(full_frequencies$chromosome, levels=c("1", "2", "3", "4", "5", "6"), labels=c("I", "II", "III", "IV", "V", "X"), ordered=T)
  full_HH <- data.frame(full_HH, row.names=NULL)
  full_HH$chromosome <- factor(full_HH$chromosome, levels=c("1", "2", "3", "4", "5", "6"), labels=c("I", "II", "III", "IV", "V", "X"), ordered=T)
  colnames(full_frequencies) <- c("Chromosome", "Start", "Stop", "Frequency")
  colnames(full_HH) <- c("Chromosome", "Start", "Stop", "Frequency")
  
  return(list(frequencies=full_frequencies, membership=full_membership, HH=full_HH)) 
}

process_haps <- function(haplotype_calls, directory, number_strains) {
  
  chromosome <- factor(haplotype_calls$membership[, "chrom_key"], levels=1:6, 
                       ordered = T, labels = c("I", "II", "III", "IV", "V","X"))
  start <- haplotype_calls$membership[, "starts"]
  stop <- haplotype_calls$membership[, "stops"]
  strains <- colnames(haplotype_calls$membership)[4:ncol(haplotype_calls$membership)]
  strains[strains == "PX174"] <- "RC301"
  
  haplotypes <- data.frame(chromosome, start, stop)
  
  haps_all <- matrix(as.character(haplotype_calls$membership[, 4:ncol(haplotype_calls$membership)]), ncol=length(4:ncol(haplotype_calls$membership)))
  colnames(haps_all) <- strains
  
  
  haps_noSingle <- t(apply(haps_all, 1, function(hrow){hrow[table(hrow)[hrow] == 1 ] <- NA; return(hrow)}))
  colnames(haps_noSingle) <- strains
  
  haps <- haps_noSingle
  ## rename haplotypes for some consistency
  
  # rank by similarity to LSJ1
  # similarity <- sort(colSums((haps== haps[, "LSJ1"])), dec = T)
  # sorted_strains <- names(similarity)
  
  #rank by total haplotype sharing
  similarity <- laply(strains, function(S){sum( (haplotypes$stop - haplotypes$start) * (haps == haps[,S]), na.rm = T )})
  sorted_strains <- strains[order(similarity, decreasing = T)]
  
  # replace numbers with strain names 
  for (strain in sorted_strains){
    haps[haps == haps[,strain] & !(haps %in% sorted_strains)] <- strain
  }
  
  strain_alpha = .9 + (1:number_strains %% 2) * .1
  
  haplotypes <- cbind(haplotypes, haps)
  
  colorscheme = hsv(h = 0:(number_strains-1)/number_strains * .85, v = c(.8,1,1), s = c(1,1, .6))
  
  
  # Melt for plotting with ggplot
  
  haplotype_melt <- gather(haplotypes, haplotype, value,-chromosome,-start,-stop)
  
  
  haplotype_melt$value <- factor(haplotype_melt$value, levels = sorted_strains, ordered = T)
  haplotype_melt$sat <- strain_alpha[match(haplotype_melt$value,levels(haplotype_melt$value))]
  #get positions for plotting
  haplotype_melt$plotpoint = match(haplotype_melt$haplotype, sorted_strains)
  
  hap_counts = ddply(haplotype_melt, .variables=.(chromosome, start, stop, value), .fun= nrow)
  
  
  dir.create(directory)
  setwd(directory)
  
  hap_counts <- dplyr::tbl_df(hap_counts)
  chromosome <- dplyr::tbl_df(chromosome)
  haplotype_melt <- dplyr::tbl_df(haplotype_melt)
  sorted_strains <- dplyr::tbl_df(sorted_strains)
  colorscheme <- dplyr::tbl_df(colorscheme)
  
  save(hap_counts, file = "complete_hap_counts.Rda")
  save(chromosome, file = "complete_chromosome_haps.Rda")
  save(haplotype_melt, file = "complete_melted_haps.Rda")
  save(sorted_strains, file = "complete_sorted_strains.Rda")
  save(colorscheme, file = "colorscheme.Rda")
  
  return(list(hap_counts,chromosome,haplotype_melt,sorted_strains,colorscheme))
}

proc_germ <- readr::read_tsv("haplotype_out.ibd", col_names = c("strain1", "ind1", "strain2", "ind2", "chrom", "start", "stop", 'LOD'))
strain <- unique(c(unique(proc_germ$strain1),unique(proc_germ$strain2)))

haps <- analyze.segs(proc_germ, strains = strain)

processed_haps <- process_haps(haps, 
                               directory = getwd(),
                               number_strains = length(strain))
save(processed_haps, filename='processed_haps.Rda')

color_plotpoint <- processed_haps[[5]] %>% 
  dplyr::mutate(cvalue = row_number()) %>%
  dplyr::rename(color = value)

plot_df <- 
  processed_haps[[3]] %>% 
  dplyr::mutate(cvalue = as.integer(value),
                hap_length= stop - start) %>%
  dplyr::group_by() %>%
  dplyr::filter(hap_length > 1E5) %>%
  dplyr::rename(isotype=value) %>%
  dplyr::left_join(color_plotpoint, by = c("cvalue"))

strain_labels <- (plot_df %>%
                    dplyr::select(haplotype, plotpoint) %>%
                    dplyr::distinct() %>%
                    dplyr::arrange(plotpoint) %>%
                    dplyr::select(haplotype))$haplotype

#==========#
# Analysis #
#==========#

# Pull out max haplotypes
chrom_hap_value <- plot_df %>%
  dplyr::group_by(chromosome, haplotype, isotype) %>%
  dplyr::summarize(hap_len = sum(stop - start)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(chromosome, isotype) %>%
  dplyr::mutate(hap_max = max(hap_len),
                shared_w_max = hap_len / hap_max,
                hap_count = n(),
                hap_prevalence = n()/length(strain)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(chromosome) %>%
  dplyr::mutate(max_haplotype_count = max(hap_count),
                max_haplotype = hap_count == max_haplotype_count)

plot_df <- plot_df %>% dplyr::left_join(chrom_hap_value, by=c("chromosome", "haplotype", "isotype"))

# Output distrib of haplotype sharing with max.
ggplot(chrom_hap_value) +
  geom_histogram(aes(x = shared_w_max)) +
  labs(x = "Shared with Max", y = "Count")

ggsave("Shared_w_max.png")

#========================#
# Filter for only swept? #
#========================#

#======#
# Plot #
#======#


plot_df %>% 
  dplyr::filter(chromosome == "I") %>%
  dplyr::arrange(shared_w_max) %>%
  ggplot(.,
         aes(xmin = start/1E6, xmax = stop/1E6, 
             ymin = plotpoint - 0.5, ymax = plotpoint + 0.5,
             fill = max_haplotype)) +
  geom_rect() +
  scale_fill_manual(values = c("Gray", "Red")) + 
  scale_y_continuous(breaks = 1:length(strain),
                     labels = strain_labels) + 
  xlab("Position (Mb)") + 
  theme_bw() +
  facet_grid(.~chromosome, scales="free", space="free")

ggsave(paste("I.png"),
       width = 8,
       height = 16)

plot_df %>%
  dplyr::mutate(val_name = as.character(value)) %>% 
  dplyr::group_by(chromosome, haplotype, value) %>%
  dplyr::mutate(len = sum(stop - start)) %>%
  dplyr::left_join(max_hap, by = c("chromosome", "haplotype", "val_name")) %>% 
  ggplot(.,
         aes(xmin = start/1E6, xmax = stop/1E6, 
             ymin = plotpoint - 0.5, ymax = plotpoint + 0.5,
             fill = max_haplotype)) +
  geom_rect() +
  scale_fill_manual(values = c("Gray", "Red")) + 
  scale_y_continuous(breaks = 1:length(strain),
                     labels = strain_labels) + 
  xlab("Position (Mb)") + 
  theme_bw() +
  facet_grid(.~chromosome, scales="free", space="free")

#===========#
# Heuristic #
#===========#

ggplot(plot_df %>% dplyr::filter(chromosome == "I", stop - start < 5000),
       aes(xmin = start/1E6, xmax = stop/1E6, 
           ymin = plotpoint - 0.5, ymax = plotpoint + 0.5,
           fill = color)) +
  geom_rect() +
  scale_y_continuous(breaks = 1:length(strain),
                     labels = strain_labels) + 
  xlab("Position (Mb)") + 
  theme_bw()

ggsave(paste("I.png"),
       width = 16,
       height = 28)

lapply(unique(processed_haps[[2]]), function(x) {
  ggplot(procssed_haps)
})
