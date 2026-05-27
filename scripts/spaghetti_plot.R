getSpagPlot <- function(CpG, mvals, pheno, anno) {
  # error handling
  stopifnot(CpG %in% rownames(mvals))
  stopifnot(CpG %in% rownames(anno))
  stopifnot("rgName" %in% names(pheno))
  stopifnot(all(c("sample_agedys","maskid") %in% names(pheno)))
  
  # match matrix column names by rgName in pheno
  samp <- intersect(colnames(mvals), pheno$rgName)
  if (length(samp) == 0) 
    stop("No overlap between colnames(matrix) and pheno$rgName")
  
  # build plotting data
  toPlot <- pheno[match(samp, pheno$rgName), , drop = F] # match sample names to pheno
  toPlot$CpG_value <- as.numeric(mvals[CpG, samp]) # match CpG & sample name to matrix
  toPlot$beta_value <- m2beta(toPlot$CpG_value)
  
  # annotation for title
  a <- anno[CpG, , drop = F]
  title <- paste0(
    CpG, " (", a$UCSC_RefGene_Name, ") ",
    a$chr, ": ", prettyNum(a$pos, big.mark = ","))
  
  ggplot(toPlot, 
         aes(x = sample_agedys, y = beta_value, group = maskid)) +
    geom_line(alpha = 0.4) +
    # scale_y_continuous(
    #   breaks = seq(0, 1, by = 0.10),
    #   limits = c(0, 1)) +
    theme_bw() +
    ylab("Beta Value") +
    labs(title = title) +
    theme(plot.title = element_text(size = 12))
}