# Generate CpG report with Metadata
#
# This function extracts methylation values for a specified CpG site from an
# M-value (or Beta-value) matrix, merges those values with phenotype metadata,
# identifies outlier samples using either an interquartile range (IQR) or
# z-score method, and returns a structured report containing annotation,
# summaries, and sample-level results.
#
# The function is designed for Illumina methylation array data (e.g., EPIC/850K)
# where CpG probes are stored as rownames of the annotation dataset and samples
# are matched between methylation and phenotype data via a shared identifier
# (rgName in this use case).
#
# Parameters:
# **CpG** Character string specifying the CpG probe ID of interest.
# **mvals** Numeric matrix of methylation values (M-values or Beta-values)
#   with CpGs stored either in rows or columns and samples in the opposite
#   dimension.
# **pheno** Data frame containing sample-level phenotype or metadata.
# **anno** Data frame containing CpG annotation with CpG probe IDs stored
#   as rownames.
# **pheno_id** Character string specifying the column in `pheno` that
#   corresponds to sample identifiers in `mvals`. Default is `"rgName"`.
# **outlier_method** Character specifying the outlier detection approach:
#   `"iqr"` for interquartile range rule (default) or `"z"` for z-score rule.
# **z_thresh** Numeric threshold defining the number of standard deviations
#   from the mean required to classify a value as an outlier when using the
#   z-score method. Default is 3.
#
# Returns anamed list containing:
# 
#   $CpG - The CpG probe ID analyzed.
#   $anno - Annotation row corresponding to the CpG probe.
#   $outlier_rule - Details of the applied outlier detection rule.
#   $summary - Summary statistics of CpG DNAm values.
#   $n - Number of samples.
#   $n_outliers - Number of detected outliers.
#   $outliers - Data frame of outliers with phenotype metadata.
#   $lowest - Top 10 lowest CpG values.
#   $highest - Top 10 highest CpG values.
#   $data - Merged sample-level dataset containing DNAm values,
#   phenotype data, and outliers.
# 
#
# The function detects whether CpGs are stored in rows of the DNAm Matrix
# and aligns samples with phenotype metadata using specified identifier column. 
# Missing values are handled using `na.rm = TRUE` during summary and outlier 
# calculations.

cpg_report <- function(
    CpG, # CpG(s) of interest
    mvals, # M-value matrix
    pheno, # phenotype data set
    anno, # annotation dataset - 850k
    pheno_id = "rgName", # column in pheno file that matches matrix columns
    outlier_method = c("iqr","z"), # outlier calculation
    z_thresh = 3, # number of deviations away from mean for z-test 
    DNAm_values = c("M", "beta")  # return M- or beta-values
    ) 
  {
  # enforce
  outlier_method <- match.arg(outlier_method)
  DNAm_values  <- match.arg(DNAm_values)
  
  # M to Beta function
  m_to_beta <- function(M) {
    p <- 2^M
    p / (p + 1)
  }
  
  # error handling 
  if (!CpG %in% rownames(anno)) stop("CpG not in rownames(anno): ", CpG)
  
  cpg_in_rows <- CpG %in% rownames(mvals)
  cpg_in_cols <- CpG %in% colnames(mvals)
  if (!cpg_in_rows && !cpg_in_cols) stop("CpG not found in mvals: ", CpG)
  
  if (!pheno_id %in% names(pheno)) stop("pheno_id column not in pheno: ", pheno_id)
  
  # extract values from M-Matrix 
  vals <- if (cpg_in_rows) mvals[CpG, ] else mvals[, CpG]
  vals <- as.numeric(vals)
  names(vals) <- if (cpg_in_rows) colnames(mvals) else rownames(mvals)
  
  # convert M -> Beta if chosen
  if (DNAm_values == "beta") {
    vals <- m_to_beta(vals)
  }
  
  # match to pheno dataset 
  keep <- intersect(names(vals), pheno[[pheno_id]])
  df <- pheno[match(keep, pheno[[pheno_id]]), , drop = F]
  df$CpG_value <- vals[keep]
  
  # dentify outliers 
  x <- df$CpG_value
  if (outlier_method == "iqr") {
    q1 <- quantile(x, 0.25, na.rm = T)
    q3 <- quantile(x, 0.75, na.rm = T)
    iqr <- q3 - q1
    lo <- q1 - 1.5 * iqr
    hi <- q3 + 1.5 * iqr
    df$is_outlier <- x < lo | x > hi
    df$outlier_score <- NA_real_
    outlier_rule <- list(method="iqr", lo=lo, hi=hi, q1=q1, q3=q3, iqr=iqr)
  } else {
    z <- (x - mean(x, na.rm = T)) / sd(x, na.rm = T)
    df$is_outlier <- abs(z) >= z_thresh
    df$outlier_score <- z
    outlier_rule <- list(method="z", z_thresh=z_thresh)
  }
  
  # annotation  
  a <- anno[CpG, , drop = F]
  
  # summaries 
  df_sorted <- df[order(df$CpG_value), ]
  outliers <- df[df$is_outlier %in% T, , drop = F]
  
  # convert to a large list for organization
  list(
    CpG = CpG,
    DNAm_scale = DNAm_values,
    anno = a,
    outlier_rule = outlier_rule,
    summary = summary(df$CpG_value),
    n = nrow(df),
    n_outliers = nrow(outliers),
    outliers = outliers,
    lowest = head(df_sorted, 10),
    highest = head(df_sorted[rev(seq_len(nrow(df_sorted))), ], 10),
    data = df)
  }
