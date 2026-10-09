# =============================================================================
# MNAR-CART SENSITIVITY ANALYSIS: 2020 ANES PRESIDENTIAL VOTE
# -----------------------------------------------------------------------------
# Conducts a sensitivity analysis of 2020 ANES presidential vote using
# multiple imputation by chained equations with CART. Covariates are imputed
# under MAR. Vote is imputed with party-stratified MNAR-CART, which applies an
# exponential donor-sampling tilt toward Trump support that depends on each
# missing row's imputed party-identification stratum.
#
# The response indicator A is fully observed, never imputed, and used as a
# predictor for every imputation target. The analysis appends 8,088 blank rows
# with A = 0 to represent unit nonrespondents.
#
# A 27-cell sensitivity grid varies the Republican, Democratic, and
# Independent sensitivity parameters over {0, 0.5, 1.0}. The marginal Trump
# share is compared with an external benchmark of 0.46, and the benchmark-
# closest cells are evaluated using conditional estimates and logistic models.
#
# Expected directory contents:
#   mnar_cart_supplementals/
#     mnar_vote_sensitivity.R
#     anes_timeseries_2020_csv_20220210.csv
#
# Usage (run from any working directory):
#   Rscript mnar_cart_supplementals/mnar_vote_sensitivity.R
# =============================================================================

# Paths and libraries ---------------------------------------------------------
# The data file and output directory are located relative to the
# supplemental-materials folder via here().
suppressMessages({
  library(tidyverse)
  library(rpart)
  library(here)   
})
options(dplyr.summarise.inform = FALSE)

ANES_PATH  <- here("anes_timeseries_2020_csv_20220210.csv")
RESULT_DIR <- here("results")
dir.create(RESULT_DIR, showWarnings = FALSE)   # create results/ if missing


# =============================================================================
# 1. LOAD ANES AND BUILD THE WORKING VARIABLES
# =============================================================================
anes20 <- read_csv(ANES_PATH, show_col_types = FALSE)

vars_of_interest <- c(
  "V200001",   # id
  "V202110x",  # vote choice (raw)
  "V201507x",  # age
  "V201511x",  # education (1-5)
  "V201549x",  # race (1-6)
  "V201600",   # sex (1-2)
  "V201200",   # ideology (1-7)
  "V201617x",  # income
  "V201231x"   # party id (1-7)
)

race_levels <- c("White", "Black", "Hispanic", "Asian/PI",
                 "Native American", "Other/Multiple")
sex_levels  <- c("Male", "Female")

#' Recode ANES missing-value codes as missing values.
#'
#' @param x A vector containing ANES coded values.
#' @return The vector with designated missing-value codes replaced by NA.
na_code <- function(x) {
  x %>% na_if(-9) %>% na_if(-8) %>% na_if(-7) %>% na_if(-6) %>%
    na_if(-5) %>% na_if(-4) %>% na_if(-3) %>% na_if(-2) %>%
    na_if(-1) %>% na_if(99)
}

raw <- anes20 %>% select(all_of(vars_of_interest))

# A is the vote-response indicator: it equals 1 when a valid presidential vote
# was reported. It is calculated from the raw vote code before recoding so it
# remains fully observed after data cleaning.
A_vec <- as.integer(!is.na(raw$V202110x) & raw$V202110x > 0)

clean <- raw %>% mutate(across(everything(), na_code))

anes_raw <- tibble(
  id          = clean$V200001,
  A           = A_vec,
  ap_vote     = factor(if_else(clean$V202110x == 2, "Trump", "Other"),
                       levels = c("Other", "Trump")),
  a_age       = as.numeric(clean$V201507x),
  ap_income   = as.numeric(clean$V201617x),
  a_edu       = as.numeric(clean$V201511x),    # ordinal score 1-5
  ap_ideology = as.numeric(clean$V201200),     # ordinal score 1-7
  ap_partyid  = as.numeric(clean$V201231x),    # ordinal score 1-7
  a_sex       = factor(clean$V201600,   levels = 1:2, labels = sex_levels),
  a_race      = factor(clean$V201549x,  levels = 1:6, labels = race_levels)
)


# =============================================================================
# 2. APPEND UNIT-NONRESPONSE ROWS
# -----------------------------------------------------------------------------
# Appended rows represent unit nonrespondents. All variables except A are
# missing and are imputed; id remains missing.
# =============================================================================
n_add <- 8088

blank_row <- tibble(
  id = NA_real_, A = 0L,
  ap_vote = factor(NA, levels = c("Other", "Trump")),
  a_age = NA_real_, ap_income = NA_real_, a_edu = NA_real_,
  ap_ideology = NA_real_, ap_partyid = NA_real_,
  a_sex = factor(NA, levels = sex_levels),
  a_race = factor(NA, levels = race_levels)
)
blank_rows <- blank_row[rep(1, n_add), ]

dat_aug <- bind_rows(anes_raw, blank_rows) %>% as.data.frame()

cat(sprintf("Augmented data: %d rows | observed vote: %d | missing vote: %d\n",
            nrow(dat_aug), sum(!is.na(dat_aug$ap_vote)),
            sum(is.na(dat_aug$ap_vote))))


# =============================================================================
# 3. VARIABLE ROLES
# =============================================================================
all_vars  <- c("a_age", "ap_income", "a_edu", "ap_ideology", "ap_partyid",
               "a_sex", "a_race", "ap_vote")          # Imputation targets.
cat_vars  <- c("a_sex", "a_race")                     # Variables imputed using classification CART.
covar_order <- c("a_age", "a_sex", "a_race", "a_edu",
                 "ap_income", "ap_ideology", "ap_partyid")
# Party identification is imputed last among the covariates, and vote is
# imputed after all covariates. Thus, the party-identification value used to
# select each missing row's sensitivity parameter has been freshly imputed.

#' Return the predictor set for an imputation target.
#'
#' Each target is predicted by all other analysis variables and the response
#' indicator A. The identifier is not used as a predictor, and A is not an
#' imputation target.
#'
#' @param target The name of the target variable.
#' @return Names of the predictor variables for the target.
predictors_of <- function(target) c(setdiff(all_vars, target), "A")


# =============================================================================
# 4. IMPUTATION HELPERS
# =============================================================================

#' Verify that the predictors required for an imputation step are complete.
#'
#' Same as in the simulation code. Hot-deck initialization completes all missing values before each sweep. A
#' target is re-blanked only immediately before it is re-imputed, so all other
#' predictors are expected to be complete at the time of this check.
#'
#' @param data A data frame used for imputation.
#' @param pred_cols Names of predictor columns.
#' @return The input data invisibly.
assert_predictors_complete <- function(data, pred_cols) {
  bad <- pred_cols[vapply(pred_cols, function(v) anyNA(data[[v]]), logical(1))]
  if (length(bad) > 0)
    stop(sprintf("Predictor(s) unexpectedly contain NA: %s — hotdeck/re-blank invariant broken.",
                 paste(bad, collapse = ", ")))
  invisible(data)
}

#' Return the observed donor values assigned to a tree leaf.
#'
#' Same as in the simulation code. A tree fitted on observed rows must contain
#' observed donors in every leaf; an empty donor pool indicates an error.
#'
#' @param leaf_id Identifier of the leaf of interest.
#' @param leaf_obs Leaf identifiers for observed donor rows.
#' @param obs_vals Observed donor values.
#' @return The observed donor values in the specified leaf.
pool_for_leaf <- function(leaf_id, leaf_obs, obs_vals) {
  in_leaf <- leaf_obs == leaf_id
  if (!any(in_leaf))
    stop(sprintf("Leaf %s has no observed donors — leaf-conditioning broken.",
                 leaf_id))
  obs_vals[in_leaf]
}

#' Extract terminal-node identifiers for rows used to fit a tree.
#'
#' Same as in the simulation code.
#'
#' @param tree A fitted rpart tree.
#' @return Numeric terminal-node identifiers for observed rows.
node_ids_obs <- function(tree) as.numeric(rownames(tree$frame))[tree$where]

#' Determine terminal-node identifiers for new rows.
#'
#' Same as in the simulation code. Node outcome values are temporarily
#' replaced with their identifiers so vector predictions return the leaf
#' reached by each row.
#'
#' @param tree A fitted rpart tree.
#' @param newdata A data frame of rows to assign to leaves.
#' @return Numeric terminal-node identifiers for new rows.
node_ids_new <- function(tree, newdata) {
  t2 <- tree
  t2$frame$yval <- as.numeric(rownames(t2$frame))
  as.numeric(predict(t2, newdata = newdata, type = "vector"))
}

#' Impute a variable using CART donor sampling under MAR.
#'
#' Same as in simulation code. A tree is estimated from rows with observed target values after temporary
#' predictor completion. Each missing value is then imputed by sampling with
#' replacement from observed donors in its assigned terminal node.
#'
#' @param data A data frame containing the target variable.
#' @param target Name of the target variable.
#' @param type Either "anova" for continuous variables or "class" for factors.
#' @return An imputed vector for the target variable.
mar_cart_generic <- function(data, target, type) {
  
  # Identify rows with observed and missing target values.
  obs_idx <- !is.na(data[[target]])
  mis_idx <-  is.na(data[[target]])
  
  # Return the original target when no values are missing.
  if (!any(mis_idx)) return(data[[target]])
  
  # Select eligible predictors.
  preds <- predictors_of(target)
  preds <- preds[vapply(preds, function(v) !all(is.na(data[[v]])), logical(1))]
  
  assert_predictors_complete(data, preds)

  # Define training and prediction data sets.
  obs_data <- data[obs_idx, , drop = FALSE]
  mis_data <- data[mis_idx, , drop = FALSE]
  
  # Specify the CART model.
  form <- as.formula(paste(target, "~", paste(preds, collapse = " + ")))
  ctrl <- rpart.control(cp = 0.001, minsplit = 20, minbucket = 10)
  
  if (type == "anova") {
    # Fit the tree using observed outcomes.
    tree     <- rpart(form, data = obs_data, method = "anova", control = ctrl)
    leaf_obs <- node_ids_obs(tree)
    leaf_mis <- node_ids_new(tree, mis_data)
    obs_y    <- obs_data[[target]]
    
    imp      <- numeric(sum(mis_idx))
    # Sample donor values within each terminal node.
    for (lid in unique(leaf_mis)) {
      idx      <- leaf_mis == lid
      pool     <- pool_for_leaf(lid, leaf_obs, obs_y)
      imp[idx] <- sample(pool, sum(idx), replace = TRUE)
    }
    out <- data[[target]]; out[mis_idx] <- imp
    return(out)
    
  } else {  # Classification tree.
    tree     <- rpart(form, data = obs_data, method = "class", control = ctrl)
    leaf_obs <- node_ids_obs(tree)
    leaf_mis <- node_ids_new(tree, mis_data)
    obs_y    <- as.character(obs_data[[target]])
    imp      <- character(sum(mis_idx))
    for (lid in unique(leaf_mis)) {
      idx      <- leaf_mis == lid
      pool     <- pool_for_leaf(lid, leaf_obs, obs_y)
      imp[idx] <- sample(pool, sum(idx), replace = TRUE)
    }
    out <- data[[target]]
    out[mis_idx] <- factor(imp, levels = levels(data[[target]]))
    return(out)
  }
}

#' Collapse the seven-point party-identification score into three strata.
#'
#' Scores 1--3 define the Democratic stratum, scores 5--7 define the
#' Republican stratum, and all remaining values define the Independent stratum.
#'
#' @param pid_score Numeric party-identification scores.
#' @return A character vector of Democratic, Republican, or Independent strata.
collapse_party_score <- function(pid_score) {
  case_when(
    pid_score %in% c(1, 2, 3) ~ "Democrat",
    pid_score %in% c(5, 6, 7) ~ "Republican",
    TRUE                      ~ "Independent"
  )
}

#' Impute vote using party-stratified MNAR-CART.
#'
#' The procedure follows MAR CART through tree estimation and leaf assignment.
#' Within each leaf, Trump donors receive weight exp(gamma_k), whereas Other
#' donors receive weight 1. The row-specific gamma_k is selected from the
#' imputed party-identification stratum; gamma_k = 0 reproduces MAR donor
#' sampling. The donor tilt is applied to all rows with missing vote values.
#'
#' @param data A data frame containing ap_vote and ap_partyid.
#' @param g_rep MNAR sensitivity parameter for the Republican stratum.
#' @param g_dem MNAR sensitivity parameter for the Democratic stratum.
#' @param g_ind MNAR sensitivity parameter for the Independent stratum.
#' @return An imputed ap_vote factor.
mnar_cart_vote_pid <- function(data, g_rep, g_dem, g_ind) {
  obs_idx <- !is.na(data$ap_vote)
  mis_idx <-  is.na(data$ap_vote)
  if (!any(mis_idx)) return(data$ap_vote)
  
  preds <- predictors_of("ap_vote")
  preds <- preds[vapply(preds, function(v) !all(is.na(data[[v]])), logical(1))]
  
  assert_predictors_complete(data, preds)
  obs_data <- data[obs_idx, , drop = FALSE]
  mis_data <- data[mis_idx, , drop = FALSE]
  
  form <- as.formula(paste("ap_vote ~", paste(preds, collapse = " + ")))
  ctrl <- rpart.control(cp = 0.001, minsplit = 20, minbucket = 10)
  
  tree     <- rpart(form, data = obs_data, method = "class", control = ctrl)
  leaf_obs <- node_ids_obs(tree)
  leaf_mis <- node_ids_new(tree, mis_data)
  obs_vote <- as.character(obs_data$ap_vote)
  
  # Assign a sensitivity parameter from each missing row's party stratum.
  mis_party <- collapse_party_score(mis_data$ap_partyid)
  gamma_vec <- case_when(
    mis_party == "Republican" ~ g_rep,
    mis_party == "Democrat"   ~ g_dem,
    TRUE                      ~ g_ind
  )
  
  imp <- character(sum(mis_idx))
  for (lid in unique(leaf_mis)) {
    in_leaf <- which(leaf_mis == lid)
    pool    <- pool_for_leaf(lid, leaf_obs, obs_vote)
    for (k in in_leaf) { 
      # Apply the row-specific sensitivity parameter.
      g_k <- gamma_vec[k]
      pc  <- pool
      if (g_k == 0) {
        imp[k] <- sample(pc, 1)                     # MAR donor sampling.
      } else {
        # Add an opposite-category donor in homogeneous leaves so the
        # sensitivity parameter can modify donor-sampling probabilities.
        if (length(unique(pc)) == 1) {
          opp <- if (pc[1] == "Trump") "Other" else "Trump"
          pc  <- c(pc, opp)
        }
        # Compute normalized donor-sampling weights.
        w      <- ifelse(pc == "Trump", exp(g_k), 1); w <- w / sum(w)
        imp[k] <- sample(pc, 1, prob = w)
      }
    }
  }
  out <- data$ap_vote
  out[mis_idx] <- factor(imp, levels = c("Other", "Trump"))
  out
}

#' Initialize chained equations using random hot-deck donor sampling.
#'
#' Same as in the simulation code. Every missing target value is replaced by a
#' draw from observed values of the same variable.
#'
#' @param data A data frame containing imputation targets.
#' @param na_flags A named list marking the original missing values.
#' @return A data frame with initialized target values.
hotdeck_init <- function(data, na_flags) {
  curr <- data
  for (v in all_vars) {
    flag <- na_flags[[v]]
    if (any(flag)) {
      donors <- curr[[v]][!is.na(curr[[v]])]
      curr[[v]][flag] <- sample(donors, sum(flag), replace = TRUE)
    }
  }
  curr
}


# =============================================================================
# 5. ESTIMANDS COMPUTED FOR EACH COMPLETED DATA SET
# =============================================================================

#' Estimate marginal Trump support across all rows.
#'
#' Same as in the simulation code.
#'
#' @param data A completed data frame containing ap_vote.
#' @return A one-row tibble with the estimand and within-imputation variance.
est_marginal <- function(data) {
  v <- as.integer(data$ap_vote == "Trump"); n <- length(v); Q <- mean(v)
  tibble(term = "marginal", Q = Q, within_var = Q * (1 - Q) / n)
}

#' Estimate Trump support within levels of a grouping variable.
#'
#' @param data A completed data frame containing ap_vote.
#' @param gvar Name of the grouping variable.
#' @return A tibble of group-specific estimates and within-imputation variances.
est_trump_by <- function(data, gvar) {
  data %>%
    mutate(v = as.integer(ap_vote == "Trump")) %>%
    group_by(level = .data[[gvar]]) %>%
    summarise(n = n(), Q = mean(v), .groups = "drop") %>%
    transmute(term = paste0(gvar, "=", level),
              Q = Q, within_var = Q * (1 - Q) / n)
}

#' Estimate a logistic regression model for Trump support.
#'
#' @param data A completed data frame containing vote and covariates.
#' @return A tibble of coefficient estimates and squared standard errors.
est_logistic <- function(data) {
  d <- data %>%
    transmute(
      y        = as.integer(ap_vote == "Trump"),
      pid      = as.numeric(ap_partyid),
      ideology = as.numeric(ap_ideology),
      edu      = as.numeric(a_edu),
      age      = as.numeric(a_age),
      income   = as.numeric(ap_income),
      sex      = factor(a_sex,  levels = sex_levels),
      race     = factor(a_race, levels = race_levels)
    )
  tryCatch({
    fit <- suppressWarnings(
      glm(y ~ pid + ideology + edu + age + sex + race + income,
          data = d, family = binomial()))
    co  <- summary(fit)$coefficients
    tibble(term = rownames(co), Q = co[, "Estimate"],
           within_var = co[, "Std. Error"]^2)
  }, error = function(e)
    tibble(term = character(), Q = numeric(), within_var = numeric()))
}


# =============================================================================
# 6. RUBIN'S RULES POOLING
# -----------------------------------------------------------------------------
# Q_bar is the mean estimate, U_bar is the mean within-imputation variance,
# B is the between-imputation variance, and Total_var is the pooled variance.
# Confidence intervals use the normal approximation.
# =============================================================================

#' Pool imputation-specific estimates using Rubin's rules.
#'
#' Same as in the simulation code.
#'
#' @param res_list A list of imputation-specific estimate tibbles.
#' @return A tibble of pooled estimates, variances, and confidence intervals.
pool_rubin <- function(res_list) {
  bind_rows(res_list, .id = "imp") %>%
    group_by(term) %>%
    summarise(
      m     = sum(!is.na(Q)),
      Q_bar = mean(Q, na.rm = TRUE),
      U_bar = mean(within_var, na.rm = TRUE),
      B     = if (sum(!is.na(Q)) > 1) var(Q, na.rm = TRUE) else 0,
      .groups = "drop"
    ) %>%
    mutate(Total_var = U_bar + (1 + 1 / pmax(m, 1)) * B,
           SE = sqrt(Total_var),
           ci_lo = Q_bar - 1.96 * SE, ci_hi = Q_bar + 1.96 * SE)
}


# =============================================================================
# 7. GRID-CELL RUNNER
# -----------------------------------------------------------------------------
# Each call evaluates one (g_rep, g_dem, g_ind) cell using m completed data
# sets and maxit chained-equation sweeps. Each sweep imputes covariates under
# MAR-CART and then imputes vote under party-stratified MNAR-CART. Original
# missing values are re-blanked before re-imputation, as in standard MICE.
#
# The marginal Trump share is recorded after each sweep to monitor convergence.
# =============================================================================

#' Run multiple imputation for one sensitivity-parameter grid cell.
#'
#' @param data A data frame containing missing values to impute.
#' @param g_rep MNAR sensitivity parameter for the Republican stratum.
#' @param g_dem MNAR sensitivity parameter for the Democratic stratum.
#' @param g_ind MNAR sensitivity parameter for the Independent stratum.
#' @param m Number of completed data sets to generate.
#' @param maxit Number of chained-equation sweeps.
#' @param base_seed Base random-number seed for the grid cell.
#' @return A named list of pooled estimands and convergence traces.
run_scenario <- function(data, g_rep, g_dem, g_ind, m, maxit, base_seed) {
  na_flags <- setNames(lapply(all_vars, function(v) is.na(data[[v]])), all_vars)
  
  store <- list(marginal = vector("list", m), pid = vector("list", m),
                ideo = vector("list", m), edu = vector("list", m),
                logit = vector("list", m))
  conv_trace <- vector("list", m)
  
  for (i in seq_len(m)) {
    # Use common random numbers so cell differences reflect gamma values.
    set.seed(base_seed + i * 100)
    curr    <- hotdeck_init(data, na_flags)
    trace_i <- numeric(maxit)
    
    for (iter in seq_len(maxit)) {
      for (v in covar_order) {
        if (!any(na_flags[[v]])) next
        curr[[v]][na_flags[[v]]] <- NA
        type <- if (v %in% cat_vars) "class" else "anova"
        curr[[v]] <- mar_cart_generic(curr, v, type)
      }
      if (any(na_flags$ap_vote)) {
        curr$ap_vote[na_flags$ap_vote] <- NA
        curr$ap_vote <- mnar_cart_vote_pid(curr, g_rep, g_dem, g_ind)
      }
      trace_i[iter] <- mean(curr$ap_vote == "Trump")
    }
    
    conv_trace[[i]] <- tibble(imp = i, iter = seq_len(maxit),
                              trump_share = trace_i)
    
    store$marginal[[i]] <- est_marginal(curr)
    store$pid[[i]]      <- est_trump_by(curr, "ap_partyid")
    store$ideo[[i]]     <- est_trump_by(curr, "ap_ideology")
    store$edu[[i]]      <- est_trump_by(curr, "a_edu")
    store$logit[[i]]    <- est_logistic(curr)
  }
  
  out <- lapply(store, pool_rubin)
  out$convergence <- bind_rows(conv_trace)
  out
}


# =============================================================================
# 8. GRID CONFIGURATION
# -----------------------------------------------------------------------------
# The sensitivity grid varies each of g_rep, g_dem, and g_ind over
# {0, 0.5, 1.0}, yielding 27 cells. A positive sensitivity parameter increases
# the odds of Trump donor selection within a missing row's terminal node.
# =============================================================================
target    <- 0.46     # External benchmark for the marginal estimate.
base_seed <- 9000

gamma_vals <- c(0, 0.5, 1.0)
m_run      <- 10
maxit_run  <- 6
grid <- expand_grid(g_rep = gamma_vals, g_dem = gamma_vals, g_ind = gamma_vals)


# =============================================================================
# 9. RUN THE GRID
# -----------------------------------------------------------------------------
# Results are checkpointed after each completed cell so interrupted runs can
# resume. The checkpoint is deleted after the analysis finishes; final results
# remain available as session objects.
# =============================================================================
dir.create(RESULT_DIR, showWarnings = FALSE, recursive = TRUE)
ckpt_file <- file.path(RESULT_DIR, "mnar_grid_ckpt.rds")

cell_results <- vector("list", nrow(grid))
surface_rows <- vector("list", nrow(grid))
k_done       <- 0L
if (file.exists(ckpt_file)) {
  ck <- readRDS(ckpt_file)
  cell_results <- ck$cell_results; surface_rows <- ck$surface_rows
  k_done <- ck$k_done
  cat(sprintf("Resuming grid from checkpoint: %d/%d cells done\n",
              k_done, nrow(grid)))
}

total_start <- proc.time()
for (k in seq_len(nrow(grid))) {
  if (k <= k_done) next
  gr <- grid$g_rep[k]; gd <- grid$g_dem[k]; gi <- grid$g_ind[k]
  cat(sprintf("cell %d/%d | g_rep=%.2f g_dem=%.2f g_ind=%.2f ... ",
              k, nrow(grid), gr, gd, gi))
  
  pooled <- run_scenario(dat_aug, gr, gd, gi, m_run, maxit_run, base_seed)
  cell_results[[k]] <- list(g_rep = gr, g_dem = gd, g_ind = gi,
                            pooled = pooled)
  
  mg <- pooled$marginal
  surface_rows[[k]] <- tibble(g_rep = gr, g_dem = gd, g_ind = gi,
                              Q_bar = mg$Q_bar, SE = mg$SE,
                              ci_lo = mg$ci_lo, ci_hi = mg$ci_hi)
  cat(sprintf("Q_bar=%.4f\n", mg$Q_bar))
  
  saveRDS(list(cell_results = cell_results, surface_rows = surface_rows,
               k_done = k),
          ckpt_file)
}
cat(sprintf("Grid elapsed: %.1f sec\n",
            (proc.time() - total_start)["elapsed"]))


# =============================================================================
# 10. MARGINAL SENSITIVITY SURFACE
# =============================================================================
marginal_surface <- bind_rows(surface_rows) %>%
  mutate(dist_to_benchmark = abs(Q_bar - target))

cat("\n── Marginal sensitivity surface ──\n")
marginal_surface %>%
  arrange(g_rep, g_dem, g_ind) %>%
  mutate(across(c(Q_bar, SE, ci_lo, ci_hi, dist_to_benchmark),
                ~ round(.x, 4))) %>%
  print(n = Inf)


# =============================================================================
# 11. CONVERGENCE SUMMARY
# -----------------------------------------------------------------------------
# For each grid cell, the summary reports the mean marginal Trump share at
# each sweep and the absolute change between the final two sweeps. Small final
# changes indicate stabilization within maxit sweeps. Full traces are retained
# in cell_results[[k]]$pooled$convergence.
# =============================================================================
conv_by_iter <- map_dfr(seq_len(nrow(grid)), function(k) {
  cr <- cell_results[[k]]
  cr$pooled$convergence %>%
    group_by(iter) %>%
    summarise(mean_share = mean(trump_share), .groups = "drop") %>%
    mutate(g_rep = cr$g_rep, g_dem = cr$g_dem, g_ind = cr$g_ind)
})

conv_summary <- conv_by_iter %>%
  group_by(g_rep, g_dem, g_ind) %>%
  arrange(iter, .by_group = TRUE) %>%
  summarise(
    share_first = first(mean_share),
    share_last  = last(mean_share),
    last_drift  = abs(last(mean_share) - nth(mean_share, n() - 1)),
    .groups = "drop"
  )

cat("\n── Convergence: mean marginal by sweep (first/last) and final-sweep drift ──\n")
conv_summary %>%
  mutate(across(c(share_first, share_last, last_drift), ~ round(.x, 4))) %>%
  print(n = Inf)

if (any(conv_summary$last_drift > 0.005))
  cat("NOTE: some cells drifted > 0.005 between the last two sweeps;",
      "consider increasing maxit.\n")


# =============================================================================
# 12. BENCHMARK-CLOSEST CELLS: CONDITIONAL AND LOGISTIC ANALYSIS
# -----------------------------------------------------------------------------
# Conditional estimates and logistic models are evaluated for cells closest to
# the external marginal benchmark of 0.46. Because multiple parameter
# combinations can yield similar marginal estimates, the analysis reports all
# selected cells and summarizes stability across them.
# =============================================================================
TOP_K <- 4    # Number of benchmark-closest cells analyzed in detail.

pid_labels  <- c("Strong Dem", "Not strong Dem", "Ind-Dem", "Independent",
                 "Ind-Rep", "Not strong Rep", "Strong Rep")
ideo_labels <- c("Ext liberal", "Liberal", "Slt liberal", "Moderate",
                 "Slt conservative", "Conservative", "Ext conservative")
edu_labels  <- c("< HS", "HS", "Some college", "Bachelor's", "Graduate")
label_map   <- list(ap_partyid = pid_labels, ap_ideology = ideo_labels,
                    a_edu = edu_labels)

#' Apply readable labels to conditional-estimate output.
#'
#' @param df A pooled conditional-estimate tibble.
#' @param gvar Name of the grouping variable.
#' @return A tibble with readable group labels and rounded estimates.
label_conditional <- function(df, gvar) {
  labs <- label_map[[gvar]]
  df %>%
    mutate(level = as.integer(sub(".*=", "", term)), group = labs[level]) %>%
    arrange(level) %>%
    transmute(group, Q_bar = round(Q_bar, 3),
              ci_lo = round(ci_lo, 3), ci_hi = round(ci_hi, 3))
}

#' Retrieve a stored sensitivity-grid cell by its parameter values.
#'
#' @param gr Republican-stratum sensitivity parameter.
#' @param gd Democratic-stratum sensitivity parameter.
#' @param gi Independent-stratum sensitivity parameter.
#' @return The matching stored cell, or NULL if no match is found.
find_cell <- function(gr, gd, gi) {
  for (c in cell_results)
    if (!is.null(c) && c$g_rep == gr && c$g_dem == gd && c$g_ind == gi)
      return(c)
  NULL
}

bench_cells <- marginal_surface %>%
  arrange(dist_to_benchmark) %>% head(TOP_K)

cat("\n", strrep("=", 70),
    "\nBENCHMARK-CLOSEST CELLS (gamma calibrated to marginal 0.46)\n",
    strrep("=", 70), "\n", sep = "")
print(bench_cells %>%
        mutate(across(c(Q_bar, ci_lo, ci_hi, dist_to_benchmark),
                      ~ round(.x, 4))))

for (r in seq_len(nrow(bench_cells))) {
  gr <- bench_cells$g_rep[r]; gd <- bench_cells$g_dem[r]
  gi <- bench_cells$g_ind[r]
  cell <- find_cell(gr, gd, gi)
  cat(sprintf(
    "\n\n■ Cell  g_rep=%.2f  g_dem=%.2f  g_ind=%.2f  (marginal Q_bar=%.4f)\n",
    gr, gd, gi, cell$pooled$marginal$Q_bar))
  cat("\n  P(Trump | party id):\n")
  print(label_conditional(cell$pooled$pid,  "ap_partyid"))
  cat("\n  P(Trump | ideology):\n")
  print(label_conditional(cell$pooled$ideo, "ap_ideology"))
  cat("\n  P(Trump | education):\n")
  print(label_conditional(cell$pooled$edu,  "a_edu"))
  cat("\n  Logistic (vote ~ pid + ideology + edu + age + sex + race + income):\n")
  print(cell$pooled$logit %>%
          transmute(term, coef = round(Q_bar, 3),
                    ci_lo = round(ci_lo, 3), ci_hi = round(ci_hi, 3)))
}

#' Summarize estimate ranges across benchmark-closest cells.
#'
#' The reported spread is the maximum minus the minimum pooled estimate across
#' the selected cells.
#'
#' @param extract A function that extracts a pooled-estimate tibble from a cell.
#' @return A tibble of minimums, maximums, and spreads by estimand.
stability <- function(extract) {
  lapply(seq_len(nrow(bench_cells)), function(r) {
    cell <- find_cell(bench_cells$g_rep[r], bench_cells$g_dem[r],
                      bench_cells$g_ind[r])
    extract(cell$pooled)
  }) %>%
    bind_rows(.id = "cell") %>%
    group_by(term) %>%
    summarise(min = min(Q_bar), max = max(Q_bar),
              spread = max(Q_bar) - min(Q_bar), .groups = "drop")
}

cat("\n", strrep("=", 70),
    "\nSTABILITY ACROSS BENCHMARK-CLOSEST CELLS (spread = max - min)\n",
    strrep("=", 70), "\n", sep = "")
cat("\n  P(Trump | party id):\n")
print(stability(function(p) p$pid) %>%
        mutate(across(where(is.numeric), ~ round(.x, 3))))
cat("\n  Logistic coefficients:\n")
print(stability(function(p) p$logit) %>%
        mutate(across(where(is.numeric), ~ round(.x, 3))))

# Remove the temporary checkpoint after the analysis completes.
if (file.exists(ckpt_file)) file.remove(ckpt_file)
cat("\nDone. Results are in: marginal_surface, cell_results, bench_cells,",
    "conv_summary.\n")
