vcf_to_granges <- function(vcf) {
  require(plyranges)
  require(VariantAnnotation)

  width <- sapply(info(vcf)$SVLEN, function(x) {
    x[1]
  })

  data.frame(
    chr = as.character(seqnames(rowRanges(vcf))),  # Chromosome
    start = start(rowRanges(vcf)),  # Start position
    end = start(rowRanges(vcf))+ width,
    type = info(vcf)$SVTYPE) %>%
    plyranges::as_grannges(seqnames=chr)
}
