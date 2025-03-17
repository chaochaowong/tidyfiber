plot_cna <- function(binned_gr, chrom_sizes_gr,
                     where=NULL) {
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
    colorRampPalette(brewer.pal(12,"Set3"))(length(keep_seq_levels))

  # Assign colors to chromosome names
  chrom_colors <- setNames(chrom_palette, keep_seq_levels)

  # must have a sanity check on the seq-levels
  chrom_sizes_gr <- keepSeqlevels(chrom_sizes_gr,
                                  value=keep_seq_levels,
                                  pruning.mode='coarse')

  df <- as.data.frame(binned_gr) %>%
    dplyr::inner_join(as.data.frame(chrom_sizes_gr), by='seqnames',
                      suffix=c('_bin', '_chrom_sizes')) %>%
    dplyr::mutate(global_start = start_bin + cum_start)

  df %>%
    ggplot(aes(y=score, x=global_start, color=seqnames)) +
    geom_point(size=0.2, alpha=0.8) +
    scale_x_continuous(
      breaks = chrom_sizes_gr$cum_midpoint,
      labels = seqlevels(chrom_sizes_gr)
    ) +
    labs(title = "5k Binned HiFiCNV Depth",
         x = "Chromosome", y = "Depth") +
    theme_void() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          legend.position = "none")  +
    scale_y_continuous(limits = c(0, 25)) +
    scale_color_manual(values = chrom_colors)
}

