#` @export
read_ft_extract_all <- function(file) {

  extra_cols <- col_names <- c(sam_flag = 'character', 
                               HP = 'character',
                               RG = 'character',
                               fiber_length = 'numeric', 
                               fiber_sequence = 'DNAStringSet',  
                               ec = 'character',
                               rq = 'character',
                               total_AT_bp = 'numeric',
                               total_m6a_bp = 'numeric',
                               total_nuc_bp = 'numeric',
                               total_msp_bp = 'numeric',
                               total_5mC_bp = 'numeric',
                               nuc_starts = 'character',
                               nuc_lengths = 'character',
                               ref_nuc_starts = 'character',
                               ref_nuc_lengths = 'character', 
                               msp_starts = 'character', 
                               msp_lengths = 'character',
                               fire = 'character',
                               ref_msp_starts = 'character',
                               ref_msp_lengths = 'character', 
                               m6a = 'character', 
                               ref_m6a = 'character', 
                               m6a_qual = 'character', 
                               `5mC` = 'character',
                               ref_5mC = 'character', 
                               `5mC_qual` = 'character')
  rtracklayer::import.bed(file, extraCols = extra_cols)
  
}