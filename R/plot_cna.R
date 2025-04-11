plot_cna <- function(binned_gr, chrom_sizes_gr,
                     where=NULL,
                     ylim = c(0, 75),
                     annotate_cnv = FALSE,
                     cnv_gr = NULL,
                     cnv_ypos = ylim[2]-5) {
  # note that the keep_sequence_levels depends on binned_gr;
  # keep chr1-22 color 
  require(RColorBrewer)
  
  # filter
  if (!is.null(where)) {
    binned_gr <- binned_gr %>%
      plyranges::filter_by_overlaps(where)
  }
  
  # sanity check
  if (length(binned_gr) < 1) 
    stop(paste0(as.character(where), ' is not within the binned_gr.'))
  
  
  keep_seq_levels <- seqlevels(binned_gr)
  
  # define color palette
  chrom_palette <-
    colorRampPalette(brewer.pal(12, 'Set3'))(length(keep_seq_levels))
  
  # Assign colors to chromosome names
  chrom_colors <- setNames(chrom_palette, keep_seq_levels)
  
  # must have a sanity check on the seq-levels, must agree with keep_seq_levels
  chrom_sizes_gr <- keepSeqlevels(chrom_sizes_gr,
                                  value=keep_seq_levels,
                                  pruning.mode='coarse')
  
  # prepare data.frame for binned depth: get global starts
  df <- as.data.frame(binned_gr) %>%
    dplyr::inner_join(as.data.frame(chrom_sizes_gr), by='seqnames',
                      suffix=c('_bin', '_chrom_sizes')) %>%
    dplyr::mutate(global_start = start_bin + cum_start) 
  
  # plot CNA: labels based on chrom_sizes_gr's seqlevels
  gg <- df %>%
    ggplot(aes(y=score, x=global_start, color=seqnames)) +
    geom_point(size=0.2, alpha=0.8, show.legend = FALSE) +
    scale_x_continuous(
      breaks = chrom_sizes_gr$cum_midpoint,  
      labels = as.character(seqnames(chrom_sizes_gr))
    ) +
    labs(title = "5k Binned HiFiCNV Depth", 
         x = "Chromosome", y = "Depth") +
    theme_void() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1), 
          legend.position = "bottom")  +
    scale_y_continuous(limits = ylim) + 
    scale_color_manual(values = chrom_colors) 
    
  
  if (annotate_cnv) {
    gg <- .annotate_cnv(gg, 
                        where = where,
                        keep_seq_levels = keep_seq_levels, 
                        chrom_sizes_gr = chrom_sizes_gr, 
                        cnv_gr, 
                        cnv_ypos=ylim[2]-5)
  }
  
  return(gg)
}

.annotate_cnv <- function(gg, 
                          where = NULL,
                          keep_seq_levels,
                          chrom_sizes_gr, 
                          cnv_gr, 
                          cnv_ypos=ylim[2]-5) {
  
  # sanity check
  # tidy up cnv_gr
  cnv_gr <- cnv_gr %>% 
    plyranges::filter(!is.na(SVTYPE))
  cnv_gr <- keepSeqlevels(cnv_gr,
                          value=keep_seq_levels,
                          pruning.mode='coarse')
  if (!is.null(where))
    cnv_gr <- cnv_gr %>%
      plyranges::filter_by_overlaps(where)
  
  # sanity check to conform chrom_size_gr
  chrom_sizes_gr <- keepSeqlevels(chrom_sizes_gr,
                                  value=keep_seq_levels,
                                  pruning.mode='coarse')
    
  # convert to data.frame and get global_start and global_end
  cnv_df <- as.data.frame(cnv_gr) %>%
    dplyr::mutate(SVTYPE = factor(SVTYPE)) %>%
    dplyr::left_join(as.data.frame(chrom_sizes_gr), by='seqnames',
                      suffix=c('_cnv', '_chrom_sizes')) %>%
    dplyr::mutate(global_start = start_cnv + cum_start,
                  global_end = end_cnv + cum_start)
        
  gg <- gg +
    geom_rect(data = cnv_df,
              aes(xmin = global_start, xmax = global_end,
                  ymin = cnv_ypos - 5, ymax = cnv_ypos, fill = SVTYPE),
              inherit.aes = FALSE) +
    scale_fill_manual(values = c(
      DEL = "skyblue",
      DUP = "tomato"
    )) +
    theme(
      legend.key.size = unit(0.4, "cm"),       # smaller boxes
      legend.text = element_text(size = 8),    # smaller text
      legend.title = element_text(size = 9),   # smaller title
      legend.spacing.y = unit(0.1, "cm")       # less vertical spacing
    )
  
  return(gg)
}


