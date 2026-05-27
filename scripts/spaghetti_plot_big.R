# helper function for plotting predicted trajectories from splines
make_spline_pred_df <- function(
    CpG, fitted_model, value_type = c("beta", "m")) 
{
  value_type <- match.arg(value_type)
  
  fitted_model %>%
    filter(.data$CpG == !!CpG) %>%
    select(CpG, starts_with("pred_age_")) %>%
    pivot_longer(
      cols = starts_with("pred_age_"),
      names_to = "age_col",
      values_to = "pred_m") %>%
    mutate(
      age_yrs = as.numeric(str_remove(age_col, "pred_age_")),
      sample_agedys = age_yrs * 365.25,
      pred_y = if (value_type == "beta") m2beta(pred_m) else pred_m,
      model = "Spline fit")
}


getSpagPlot_big <- function(
    CpG,
    mvals_bm,
    pheno,
    anno,
    probe_ids,
    sample_ids,
    pred_df = NULL,
    value_type = c("beta", "m"),
    age_unit = c("days", "years"),
    y_limits = NULL,
    group_col = NULL,
    facet_by_group = FALSE,
    legend = "bottom") 
{
  value_type <- match.arg(value_type)
  age_unit <- match.arg(age_unit)
  
  stopifnot(is.character(CpG), length(CpG) == 1)
  stopifnot(CpG %in% probe_ids)
  stopifnot(CpG %in% rownames(anno))
  stopifnot("rgName" %in% names(pheno))
  stopifnot(all(c("sample_agedys", "maskid") %in% names(pheno)))
  
  if (!is.null(group_col) && !group_col %in% names(pheno)) {
    stop("group_col not found in pheno: ", group_col)
  }
  
  samp <- intersect(sample_ids, pheno$rgName)
  if (length(samp) == 0) stop("No overlap between sample_ids and pheno$rgName")
  
  probe_idx <- match(CpG, probe_ids)
  samp_idx  <- match(samp, sample_ids)
  
  toPlot <- pheno[match(samp, pheno$rgName), , drop = FALSE]
  toPlot$CpG_value <- as.numeric(mvals_bm[probe_idx, samp_idx])
  
  if (!is.null(group_col)) {
    toPlot$group_plot <- as.factor(toPlot[[group_col]])
  }
  
  if (value_type == "beta") {
    toPlot$yval <- m2beta(toPlot$CpG_value)
    ylab_txt <- "Beta Value"
  } else {
    toPlot$yval <- toPlot$CpG_value
    ylab_txt <- "M-value"
  }
  
  if (age_unit == "years") {
    toPlot$x_age <- toPlot$sample_agedys / 365.25
    xlab_txt <- "Age (years)"
  } else {
    toPlot$x_age <- toPlot$sample_agedys
    xlab_txt <- "Age (days)"
  }
  
  a <- anno[CpG, , drop = FALSE]
  title <- paste0(
    CpG, " (", a$UCSC_RefGene_Name, ") ",
    a$chr, ": ", prettyNum(a$pos, big.mark = ","))
  
  p <- ggplot(toPlot, aes(x = x_age, y = yval, group = maskid))
  
  if (!is.null(group_col)) {
    p <- p +
      geom_line(aes(color = group_plot), alpha = 0.45) +
      #geom_point(aes(color = group_plot), alpha = 0.6, size = 1) +
      labs(color = group_col)
  } else {
    p <- p +
      geom_line(alpha = 0.4) 
      #geom_point(alpha = 0.6, size = 1)
  }
  
  p <- p +
    theme_bw() +
    ylab(ylab_txt) +
    xlab(xlab_txt) +
    labs(title = title) +
    theme(
      plot.title = element_text(size = 12),
      legend.position = legend)
  
  if (!is.null(group_col) && facet_by_group) {
    p <- p + facet_wrap(~ group_plot)
  }
  
  if (!is.null(pred_df)) {
    pred_df <- pred_df %>%
      mutate(
        x_age = if (age_unit == "years") sample_agedys / 365.25 else sample_agedys)
    
    p <- p +
      geom_line(
        data = pred_df,
        aes(x = x_age, y = pred_y, linetype = model),
        linewidth = 0.9,
        color = "firebrick",
        inherit.aes = FALSE)
  }
  
  if (age_unit == "years") {
    p <- p + scale_x_continuous(breaks = seq(0, 6, by = 1))
  }
  
  if (!is.null(y_limits)) {
    p <- p + coord_cartesian(ylim = y_limits)
  }
  
  return(p)
}