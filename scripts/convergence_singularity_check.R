
bad_cpgs <- unique(test_res$CpG[!is.na(test_res$error_message)])

cpgs <- bad_cpgs

# 1) Precompute formula ONCE
form_cpg <- update(form, CpG ~ .)

# 2) Keep only columns used by the formula (big speed win)
vars_needed <- setdiff(all.vars(form_cpg), "CpG")
base_dat <- as.data.frame(colData(se.filt)[, vars_needed, drop = FALSE])

# 3) Pull the response matrix ONCE
Y <- assay(se.filt, "Mvalue")[cpgs, , drop = FALSE]

# 4) Chunk work to reduce Snow overhead
chunk_size <- 50
idx <- split(seq_along(cpgs), ceiling(seq_along(cpgs) / chunk_size))

res_list <- unlist(
  bplapply(
    idx,
    function(ii, Y, base_dat, form_cpg) {
      out <- vector("list", length(ii))
      for (k in seq_along(ii)) {
        i <- ii[k]
        dat <- base_dat
        dat$CpG <- as.numeric(Y[i, ])
        
        m <- tryCatch(
          suppressWarnings(
            lmer(form_cpg, data = dat, REML = TRUE)
          ),
          error = function(e) NULL
        )
        
        if (is.null(m)) {
          out[[k]] <- list(is_sing = NA, conv_msg = "lmer error")
        } else {
          msg <- m@optinfo$conv$lme4$messages
          out[[k]] <- list(
            is_sing  = lme4::isSingular(m, tol = 1e-4),
            conv_msg = if (length(msg)) paste(unlist(msg), collapse = " | ") else ""
          )
        }
      }
      out
    },
    Y = Y, base_dat = base_dat, form_cpg = form_cpg,
    BPPARAM = param
  ),
  recursive = FALSE
)

is_sing  <- vapply(res_list, `[[`, logical(1), "is_sing")
conv_msg <- vapply(res_list, `[[`, character(1), "conv_msg")

sing_df <- data.frame(CpG = cpgs, 
                      is_singular = is_sing, 
                      lmer_messages = conv_msg)

table(sing_df$is_singular, useNA = "ifany")

# write to csv
write.csv(sing_df, file = "/Users/ryan_summers/GitHub/teddy-dnam-analysis/results/sing_df.csv")