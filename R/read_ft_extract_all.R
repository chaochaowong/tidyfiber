# ` @export
read_ft_extract_all <- function(file, col_names = NULL,
                                genome_info = NULL,
                                overlap_ranges = NULL) {
  args <- plyranges:::norm_args_reader(genome_info)
  extra_cols <- c(
    sam_flag = "character",
    HP = "character",
    RG = "character",
    fiber_length = "numeric",
    fiber_sequence = "DNAStringSet",
    ec = "character",
    rq = "character",
    total_AT_bp = "numeric",
    total_m6a_bp = "numeric",
    total_nuc_bp = "numeric",
    total_msp_bp = "numeric",
    total_5mC_bp = "numeric",
    nuc_starts = "character",
    nuc_lengths = "character",
    ref_nuc_starts = "character",
    ref_nuc_lengths = "character",
    msp_starts = "character",
    msp_lengths = "character",
    fire = "character",
    ref_msp_starts = "character",
    ref_msp_lengths = "character",
    m6a = "character",
    ref_m6a = "character",
    m6a_qual = "character",
    `5mC` = "character",
    ref_5mC = "character",
    `5mC_qual` = "character"
  )

  rtracklayer::import.bed(file,
    colnames = col_names,
    extraCols = extra_cols,
    genome = args$genome_info,
    seqinfo = args$seq_info,
    which = overlap_ranges
  )
}
