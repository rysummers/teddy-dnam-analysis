
# ------------------------------------------------------------
# count_grp_obs()
#
# Calculates counts and percentages for one or more categorical
# variables in a dataset. The function supports both:
#
# 1) Separate summaries for each variable
# 2) Interaction (joint combination) summaries across variables
#
# Selected variables are standardized, optionally reshaped to long
# format, and summarized into a tidy table with counts and
# proportions for each category or interaction level.
#
# Key features:
# - Works with data.frame and tibble inputs
# - Handles factor, character, and labelled categorical variables
# - Supports multiple variables simultaneously
# - Allows interaction (cross-classified) group summaries
# - Can compute observation-level or subject-level counts
# - Missing values optionally included and labeled "Missing"
# - Interaction summaries can optionally exclude rows with missing
#   values in any grouping variable
#
# Arguments:
# data        : data frame containing variables of interest
# vars        : character vector of categorical variable names
# include_na  : logical; if TRUE (default), include missing values
# id          : optional subject identifier; if provided, counts
#               represent unique subjects instead of observations
# mode        : "separate" (default) or "interaction"
#               - separate    → counts each variable independently
#               - interaction → counts joint combinations of vars
# sep         : separator used when constructing interaction labels
#
# Returns:
# A tibble containing category or interaction labels with counts
# and proportions. Output column names differ slightly depending
# on observation-level vs subject-level counting.
#
# Examples:
# # Separate summaries
# count_grp_obs(df, vars = c("sex", "race"))
#
# # Subject-level summaries
# count_grp_obs(df, vars = c("sex", "race"), id = participant_id)
#
# # Interaction summaries
# count_grp_obs(df, vars = c("sex", "race"), mode = "interaction")
#
# # Interaction summaries excluding missing values
# count_grp_obs(df, vars = c("sex", "race"), mode = "interaction",
#               include_na = FALSE)
# ------------------------------------------------------------


count_grp_obs <- function(data, vars, include_na = TRUE, id = NULL,
                          mode = c("separate", "interaction"),
                          sep = " × ") {
  stopifnot(is.data.frame(data))
  stopifnot(is.character(vars), length(vars) >= 1)
  
  mode <- match.arg(mode)
  
  # helper: standardize missing + coerce to character
  standardize_cat <- function(df, cols) {
    df %>%
      dplyr::mutate(
        dplyr::across(
          dplyr::all_of(cols),
          ~ dplyr::if_else(is.na(.x), "Missing", as.character(.x))))
  }
  
  # --- MODE 1: separate counts (your original behavior) ---
  if (mode == "separate") {
    dat <- dplyr::as_tibble(data) %>%
      dplyr::select(dplyr::all_of(vars)) %>%
      standardize_cat(vars)
    
    dat_long <- dat %>%
      tidyr::pivot_longer(
        cols = dplyr::everything(),
        names_to = "variable",
        values_to = "category")
    
    if (!include_na) dat_long <- dplyr::filter(dat_long, category != "Missing")
    
    if (is.null(id)) {
      return(
        dat_long %>%
          dplyr::count(variable, category, name = "n") %>%
          dplyr::group_by(variable) %>%
          dplyr::mutate(pct = n / sum(n)) %>%
          dplyr::ungroup())
    } else {
      dat_id <- dplyr::as_tibble(data) %>%
        dplyr::select({{ id }}, dplyr::all_of(vars)) %>%
        standardize_cat(vars) %>%
        tidyr::pivot_longer(
          cols = -{{ id }},
          names_to = "variable",
          values_to = "category"
        )
      
      if (!include_na) dat_id <- dplyr::filter(dat_id, category != "Missing")
      
      return(
        dat_id %>%
          dplyr::distinct({{ id }}, variable, category) %>%
          dplyr::count(variable, category, name = "n_subjects") %>%
          dplyr::group_by(variable) %>%
          dplyr::mutate(pct_subjects = n_subjects / sum(n_subjects)) %>%
          dplyr::ungroup())
    }
  }
  
  # --- MODE 2: interaction (joint combinations across vars) ---
  dat_int <- dplyr::as_tibble(data) %>%
    dplyr::select(dplyr::all_of(vars)) %>%
    standardize_cat(vars)
  
  if (!include_na) {
    # drop rows where ANY of the interaction vars are Missing
    dat_int <- dat_int %>%
      dplyr::filter(dplyr::if_all(dplyr::everything(), ~ .x != "Missing"))
  }
  
  # create a single interaction label like "Male × SiteA × Arm1"
  dat_int <- dat_int %>%
    tidyr::unite("interaction", dplyr::all_of(vars), sep = sep, remove = FALSE)
  
  if (is.null(id)) {
    dat_int %>%
      dplyr::count(interaction, name = "n") %>%
      dplyr::mutate(pct = n / sum(n)) %>%
      dplyr::arrange(dplyr::desc(n))
  } else {
    dplyr::as_tibble(data) %>%
      dplyr::select({{ id }}, dplyr::all_of(vars)) %>%
      standardize_cat(vars) %>%
      { if (!include_na)
        dplyr::filter(., dplyr::if_all(dplyr::all_of(vars), ~ .x != "Missing"))
        else .
      } %>%
      tidyr::unite("interaction", dplyr::all_of(vars), sep = sep, remove = FALSE) %>%
      dplyr::distinct({{ id }}, interaction) %>%
      dplyr::count(interaction, name = "n_subjects") %>%
      dplyr::mutate(pct_subjects = n_subjects / sum(n_subjects)) %>%
      dplyr::arrange(dplyr::desc(n_subjects))
  }
}