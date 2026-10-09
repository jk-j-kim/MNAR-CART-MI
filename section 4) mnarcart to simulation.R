# =============================================================================
# MNAR-CART SIMULATION: MARGINAL, CONDITIONAL, AND REGRESSION ESTIMANDS
# -----------------------------------------------------------------------------
# Evaluates MNAR-CART multiple imputation in a synthetic population with known
# truth. Vote is missing not at random (MNAR), and the seven covariates are
# missing completely at random (MCAR) at a rate of 15% each.
#
# Two imputation strategies are compared:
#   AO (All at Once): the sensitivity parameter is applied at every MICE
#     iteration.
#   GL (Giusti and Little): a MAR chain is run to convergence, followed by one
#     sensitivity-parameter pass for vote.
#
# Imputation conditions on age, income, education, ideology, party
# identification, sex, and race. 

# The analysis estimates the following four estimands.
#   marginal            Trump share across all rows.
#   trump_by_pid        Trump share within each level of party identification.
#   trump_by_ideology   Trump share within each level of ideology.
#   trump_by_edu        Trump share within each level of education.
#   logistic            Coefficients of a logistic regression of the Trump
#                       indicator on party identification, ideology, education,
#                       age, sex, race, and income.
#
# Truth for every estimand is the same quantity computed on the full
# population, using the same function that is applied to each completed data
# set. 
#
# Each data-generation parameter has its own list of imputation parameters, 
# set in GAMMA_IMP_BY_GEN below:
#     gamma_gen = -2.00   ->   gamma_imp in {-2.00, -1.00, 0.00}
#     gamma_gen =  0.00   ->   gamma_imp in {-2.00,  0.00, 2.00}
#     gamma_gen =  2.00   ->   gamma_imp in { 0.00,  1.00, 2.00}
#
#   Every scenario gets a label describing what it represents: correct
#   specification, MAR with no correction, a partial correction in the right
#   direction, or a correction where none was needed. Three of the nine cells
#   are MAR, so MAR appears once under each data-generation parameter.
#
# Trees are estimated using rows with observed target values. Rows with missing
# targets are assigned to leaves using node identifiers obtained from the
# fitted tree.
#
# Hot-deck initialization completes all imputation targets before each chain
# begins. At each step, only the current target is re-blanked; therefore,
# predictors are complete when rpart estimates and applies its trees.
#
# The function pool_for_leaf stops when a leaf contains no observed donors.
# Because trees are fitted on observed rows with minbucket = 10, this condition
# should not occur under the intended procedure.
#
# A single, row-invariant sensitivity parameter is used for vote. Consequently,
# the leaf-level vectorized sampling procedure is distributionally equivalent
# to sampling donors separately for each row.
#
# A is the fully observed vote-response indicator and is never imputed. It may
# be included as a predictor for each target by setting USE_A <- TRUE. Because
# A equals 1 for all observed-vote training rows, it cannot split the vote tree
# but can inform the covariate trees.
#
# Expected directory contents (the folder containing mnarcart_supplementals.Rproj):
#     section 4) mnarcart to simulation.R
#     sim_pop.csv
#
# Usage: open mnarcart_supplementals.Rproj in RStudio (or set the working
# directory to this folder), then either source the script or run, from a terminal
# in this folder:
#   Rscript "section 4) mnarcart to simulation.R" <gamma_gen> <batch>
#   Rscript "section 4) mnarcart to simulation.R" <gamma_gen> <batch> <gamma_imp_csv>
#   Rscript "section 4) mnarcart to simulation.R" <gamma_gen> <batch> <gamma_imp_csv> <reps_per_batch>
#
# The third argument overrides the scenario list for one run, for example
# "0,0.25,0.5". Leave it off (pass "") to use the list defined below.
# The fourth argument overrides how many reps this batch covers (default 250),
# so one can split a g_gen's total reps across more, smaller array tasks.
# =============================================================================

# Paths and libraries ---------------------------------------------------------
# The population data (sim_pop.csv) is created by
# "section 4) simulation synthetic population creation.R" and stored in the
# supplemental-materials folder. Paths are relative to that folder via here().
suppressMessages({
  library(tidyverse)
  library(rpart)
  library(here)  
})
options(dplyr.summarise.inform = FALSE)

POP_PATH <- here("sim_pop.csv")


# Sensitivity-parameter scenarios ---------------------------------------------
# Data-generation parameter, formatted to two decimals. 
# GAMMA_IMP_BY_GEN <- list(
#   "-0.50" = c(-0.50, -0.25,  0.00),
#   "0.00"  = c(-0.50,  0.00,  0.50),
#  "0.50"  = c( 0.00,  0.25,  0.50)
# )

GAMMA_IMP_BY_GEN <- list(
  "-2.00" = c(-2.00, -1.00,  0.00),
  "0.00"  = c(-2.00,  0.00,  2.00),
  "2.00"  = c( 0.00,  1.00,  2.00)
  )


#' Describe what an imputation parameter represents for a given truth.
#'
#' @param g_imp Imputation sensitivity parameter.
#' @param g_gen Data-generation sensitivity parameter.
#' @return A short label used in output and in saved results.
spec_label <- function(g_imp, g_gen) {
  same <- function(a, b) abs(a - b) < 1e-8
  if (same(g_gen, 0) && same(g_imp, 0))  return("MAR, correct")
  if (same(g_imp, g_gen))                return("correct")
  if (same(g_imp, 0))                    return("MAR, no correction")
  if (same(g_gen, 0))                    return("correction not needed")
  if (sign(g_imp) == sign(g_gen) && abs(g_imp) < abs(g_gen))
    return("partial correction")
  "misspecified"
}


# Run configuration -----------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
gamma_imp_values <- NULL

if (length(args) >= 2) {
  # Batch-execution configuration.
  g_gen          <- as.numeric(args[1])
  batch          <- as.integer(args[2])
  reps_per_batch <- if (length(args) >= 4 && nzchar(args[4])) as.integer(args[4]) else 250
  n_sample       <- 5000
  m              <- 10
  maxit          <- 6
  if (length(args) >= 3 && nzchar(args[3])) {
    gamma_imp_values <- as.numeric(strsplit(args[3], ",")[[1]])
  }
} else {
  # Default configuration.
  g_gen          <- 2
  batch          <- 1
  reps_per_batch <- 10
  n_sample       <- 5000
  m              <- 3
  maxit          <- 3
}

if (is.null(gamma_imp_values)) {
  key <- sprintf("%.2f", g_gen)
  if (is.null(GAMMA_IMP_BY_GEN[[key]])) {
    stop(sprintf("No scenario list for gamma_gen = %s. Add one to GAMMA_IMP_BY_GEN or pass a third argument.",
                 key))
  }
  gamma_imp_values <- GAMMA_IMP_BY_GEN[[key]]
}

spec_labels <- vapply(gamma_imp_values, spec_label, character(1), g_gen = g_gen)
USE_A <- FALSE

r_start <- (batch - 1) * reps_per_batch + 1
r_end   <-  batch      * reps_per_batch
CHECKPOINT_EVERY <- 10

cat(sprintf("Config: g_gen=%.2f | batch=%d (reps %d-%d) | n=%d | m=%d | maxit=%d\n",
            g_gen, batch, r_start, r_end, n_sample, m, maxit))
cat("Scenarios:\n")
for (k in seq_along(gamma_imp_values)) {
  cat(sprintf("  g_imp=%+.2f  (%s)\n", gamma_imp_values[k], spec_labels[k]))
}


# =============================================================================
# 1. LOAD POPULATION AND SET VARIABLE TYPES
# =============================================================================
pop_raw <- read_csv(POP_PATH, show_col_types = FALSE)

edu_levels  <- c("Less than HS", "HS diploma", "Some post-high school",
                 "Bachelor's", "Graduate degree")
ideo_levels <- c("Extremely liberal", "Liberal", "Slightly liberal", "Moderate",
                 "Slightly conservative", "Conservative", "Extremely conservative")
pid_levels  <- c("Strong Democrat", "Not very strong Democrat",
                 "Independent-Democrat", "Independent", "Independent-Republican",
                 "Not very strong Republican", "Strong Republican")
race_levels <- c("White", "Black", "Hispanic", "Asian/PI",
                 "Native American", "Other/Multiple")
sex_levels  <- c("Male", "Female")

#' Convert an ordered categorical variable to integer scores.
#'
#' @param x A vector of categorical values.
#' @param lv The ordered category levels.
#' @return An integer vector corresponding to the supplied category levels.
to_score <- function(x, lv) as.integer(factor(x, levels = lv))

pop <- pop_raw %>%
  transmute(
    a_age       = as.numeric(a_age),
    ap_income   = as.numeric(ap_income),
    a_edu       = to_score(a_edu,       edu_levels),    # Scores from 1 to 5.
    ap_ideology = to_score(ap_ideology, ideo_levels),   # Scores from 1 to 7.
    ap_partyid  = to_score(ap_partyid,  pid_levels),    # Scores from 1 to 7.
    a_sex       = factor(a_sex,  levels = sex_levels),
    a_race      = factor(a_race, levels = race_levels),
    ap_vote     = factor(ap_vote, levels = c("Other", "Trump"))
  )

if (anyNA(pop)) {
  cat("WARNING: NAs in population after coercion (check label spellings):\n")
  print(colSums(is.na(pop)))
  pop <- drop_na(pop)
}
pop <- as.data.frame(pop)

cat(sprintf("Population: %d rows | marginal Trump = %.4f\n",
            nrow(pop), mean(pop$ap_vote == "Trump")))

# Variable roles --------------------------------------------------------------
all_vars  <- c("a_age", "ap_income", "a_edu", "ap_ideology", "ap_partyid",
               "a_sex", "a_race", "ap_vote")          # Imputation targets.
cat_vars  <- c("a_sex", "a_race")                     # Variables imputed using classification CART.

# Covariates are imputed first, with party identification last; vote is imputed last.
covar_order <- c("a_age", "a_sex", "a_race", "a_edu",
                 "ap_income", "ap_ideology", "ap_partyid")

#' Return the predictor set for an imputation target.
#'
#' @param target The name of the target variable.
#' @return Names of all other analysis variables and, when enabled, A.
predictors_of <- function(target) {
  c(setdiff(all_vars, target), if (USE_A) "A")
}


# =============================================================================
# 2. IMPUTATION HELPERS
# =============================================================================

#' Draw values from a pool with replacement.
#'
#' Written out rather than calling sample() directly, because sample() on a
#' single numeric value samples from 1:x instead of returning that value.
#'
#' @param pool Values to draw from.
#' @param n Number of draws.
#' @return A vector of n draws.
sample_pool <- function(pool, n) pool[sample.int(length(pool), n, replace = TRUE)]

#' Draw values from a pool with replacement, using supplied weights.
#'
#' @param pool Values to draw from.
#' @param n Number of draws.
#' @param prob Sampling weights, one per pool element.
#' @return A vector of n draws.
sample_pool_w <- function(pool, n, prob) {
  pool[sample.int(length(pool), n, replace = TRUE, prob = prob)]
}

#' Verify that the predictors required for an imputation step are complete.
#'
#' Hot-deck initialization completes all missing values before each sweep. A
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
#' @param tree A fitted rpart tree.
#' @return Numeric terminal-node identifiers for observed rows.
node_ids_obs <- function(tree) as.numeric(rownames(tree$frame))[tree$where]

#' Determine terminal-node identifiers for new rows.
#'
#' Node outcome values are temporarily replaced with their node identifiers so
#' that vector predictions return the leaf reached by each row.
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
#' A tree is estimated from rows with observed target values. Each missing
#' value is then imputed by sampling with replacement from observed donors in
#' its assigned terminal node.
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
      imp[idx] <- sample_pool(pool, sum(idx))
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
      imp[idx] <- sample_pool(pool, sum(idx))
    }
    out <- data[[target]]
    out[mis_idx] <- factor(imp, levels = levels(data[[target]]))
    return(out)
  }
}

#' Impute vote using MNAR-CART with a scalar sensitivity parameter.
#'
#' The procedure follows MAR CART through tree estimation and leaf assignment.
#' Within each leaf, Trump donors receive weight exp(gamma_imp), whereas Other
#' donors receive weight 1. Thus, gamma_imp = 0 reproduces MAR donor sampling.
#'
#' @param data A data frame containing ap_vote.
#' @param gamma_imp Scalar MNAR sensitivity parameter.
#' @return An imputed ap_vote factor.
mnar_cart_vote <- function(data, gamma_imp) {
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
  
  imp <- character(sum(mis_idx))
  for (lid in unique(leaf_mis)) {
    idx  <- leaf_mis == lid
    pool <- pool_for_leaf(lid, leaf_obs, obs_vote)
    n_i  <- sum(idx)
    if (gamma_imp == 0) {
      imp[idx] <- sample_pool(pool, n_i)                # MAR donor sampling.
    } else {
      # Add an opposite-category donor in homogeneous leaves so the sensitivity
      # parameter can modify the donor-sampling distribution.
      if (length(unique(pool)) == 1) {
        opp  <- if (pool[1] == "Trump") "Other" else "Trump"
        pool <- c(pool, opp)
      }
      # Compute normalized donor-sampling weights.
      w        <- ifelse(pool == "Trump", exp(gamma_imp), 1)
      w        <- w / sum(w)
      imp[idx] <- sample_pool_w(pool, n_i, w)
    }
  }
  out <- data$ap_vote
  out[mis_idx] <- factor(imp, levels = c("Other", "Trump"))
  out
}

#' Initialize chained equations using random hot-deck donor sampling.
#'
#' Every missing target value is replaced by a draw from the observed values of
#' the same variable. The fully observed response indicator A is not modified.
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
      curr[[v]][flag] <- sample_pool(donors, sum(flag))
    }
  }
  curr
}


# =============================================================================
# 3. MISSINGNESS INJECTION
# =============================================================================

#' Inject MNAR missingness into vote.
#'
#' Trump voters have higher missingness probabilities when gamma_true is
#' positive, resulting in underrepresentation of Trump voters among observed
#' vote responses.
#'
#' @param data A complete data frame containing ap_vote.
#' @param gamma_true Data-generation MNAR sensitivity parameter.
#' @param target_rate Target marginal vote-missingness rate.
#' @return The data frame with missing values introduced in ap_vote.
inject_mnar_vote <- function(data, gamma_true, target_rate = 0.20) {
  v  <- as.integer(data$ap_vote == "Trump")
  b0 <- uniroot(function(b) mean(plogis(b + gamma_true * v)) - target_rate,
                c(-20, 20))$root
  p_miss <- plogis(b0 + gamma_true * v)
  data$ap_vote[rbinom(nrow(data), 1, p_miss) == 1] <- NA
  data
}

#' Inject vote MNAR missingness and covariate MCAR missingness.
#'
#' The response indicator A is created after vote missingness is injected. It
#' equals 1 for observed vote responses, is never imputed, and may be used as
#' a predictor in the imputation models.
#'
#' @param data A complete data frame.
#' @param gamma_true Data-generation MNAR sensitivity parameter for vote.
#' @param rate_vote Target vote-missingness rate.
#' @param rate_cov Target missingness rate for each covariate.
#' @return The data frame with simulated missingness and response indicator A.
inject_missingness <- function(data, gamma_true,
                               rate_vote = 0.20, rate_cov = 0.15) {
  n    <- nrow(data)
  data <- inject_mnar_vote(data, gamma_true, rate_vote)
  for (v in covar_order) {
    idx <- sample(n, round(rate_cov * n))
    data[[v]][idx] <- NA
  }
  data$A <- as.integer(!is.na(data$ap_vote))
  data
}


# =============================================================================
# 4. ESTIMANDS COMPUTED FOR EACH COMPLETED DATA SET
# -----------------------------------------------------------------------------
# The same functions are applied to the population to obtain truth, so the
# target of each estimand is defined by the code rather than by a separate
# derivation.
# =============================================================================

#' Estimate marginal Trump support across all rows.
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

#' Compute every estimand for one completed data set.
#'
#' @param data A completed data frame.
#' @return A tibble with columns estimand, term, Q, and within_var.
compute_estimands <- function(data) {
  bind_rows(
    est_marginal(data)                  %>% mutate(estimand = "marginal"),
    est_trump_by(data, "ap_partyid")    %>% mutate(estimand = "trump_by_pid"),
    est_trump_by(data, "ap_ideology")   %>% mutate(estimand = "trump_by_ideology"),
    est_trump_by(data, "a_edu")         %>% mutate(estimand = "trump_by_edu"),
    est_logistic(data)                  %>% mutate(estimand = "logistic")
  ) %>%
    select(estimand, term, Q, within_var)
}


# =============================================================================
# 5. RUBIN POOLING
# =============================================================================

#' Pool imputation-specific estimates using Rubin's rules.
#'
#' Pooling is done separately for every estimand and term. A term that is
#' absent from some completed data sets, such as a dropped regression
#' coefficient, is pooled over the data sets where it appears, and m records
#' how many those were.
#'
#' @param res_list A list of imputation-specific estimate tibbles.
#' @return A tibble of pooled estimates, variances, and confidence intervals.
pool_rubin <- function(res_list) {
  bind_rows(res_list, .id = "imp") %>%
    group_by(estimand, term) %>%
    summarise(
      m     = sum(!is.na(Q)),
      Q_bar = mean(Q, na.rm = TRUE),
      U_bar = mean(within_var, na.rm = TRUE),
      B     = if (sum(!is.na(Q)) > 1) var(Q, na.rm = TRUE) else 0,
      .groups = "drop"
    ) %>%
    mutate(
      Total_var = U_bar + (1 + 1 / pmax(m, 1)) * B,
      SE        = sqrt(Total_var),
      ci_lo     = Q_bar - 1.96 * SE,
      ci_hi     = Q_bar + 1.96 * SE
    )
}


# =============================================================================
# 6. ONE MULTIPLE-IMPUTATION RUN (AO OR GL)
# -----------------------------------------------------------------------------
# Each sweep imputes covariates using MAR-CART in covar_order and then imputes
# vote using MNAR-CART. Original missing values are re-blanked immediately
# before re-imputation, as in standard MICE. AO applies the sensitivity
# parameter in every sweep; GL performs MAR sweeps and one final
# sensitivity-parameter pass for vote.
# =============================================================================

#' Perform one multiple-imputation run for a specified strategy.
#'
#' @param data_miss A data frame with simulated missing values.
#' @param gamma_imp Scalar MNAR sensitivity parameter for vote.
#' @param m Number of completed data sets to generate.
#' @param maxit Number of chained-equation sweeps.
#' @param base_seed Base random-number seed for this run.
#' @param variant Imputation strategy: "AO" or "GL".
#' @return A tibble of Rubin-pooled estimates for every estimand and term.
run_mi <- function(data_miss, gamma_imp, m, maxit, base_seed,
                   variant = c("AO", "GL")) {
  variant  <- match.arg(variant)
  na_flags <- setNames(lapply(all_vars, function(v) is.na(data_miss[[v]])),
                       all_vars)
  
  store <- vector("list", m)
  
  for (i in seq_len(m)) {
    set.seed(base_seed + i * 100)
    curr <- hotdeck_init(data_miss, na_flags)
    
    for (iter in seq_len(maxit)) {
      for (v in covar_order) {                       # Impute covariates under MAR.
        if (!any(na_flags[[v]])) next
        curr[[v]][na_flags[[v]]] <- NA
        type <- if (v %in% cat_vars) "class" else "anova"
        curr[[v]] <- mar_cart_generic(curr, v, type)
      }
      if (any(na_flags$ap_vote)) {                    # Impute vote after all covariates.
        curr$ap_vote[na_flags$ap_vote] <- NA
        g_use <- if (variant == "GL") 0 else gamma_imp
        curr$ap_vote <- mnar_cart_vote(curr, g_use)
      }
    }
    
    if (variant == "GL" && any(na_flags$ap_vote)) {    # Apply the final GL sensitivity pass.
      curr$ap_vote[na_flags$ap_vote] <- NA
      curr$ap_vote <- mnar_cart_vote(curr, gamma_imp)
    }
    
    store[[i]] <- compute_estimands(curr)
  }
  
  pool_rubin(store)
}


# =============================================================================
# 7. POPULATION TRUTH
# -----------------------------------------------------------------------------
# Truth is the same set of estimands computed on the full population.
# =============================================================================
population_truth <- compute_estimands(pop) %>%
  select(estimand, term, truth = Q)

cat(sprintf("Truth: %d terms across %d estimands | marginal = %.4f\n",
            nrow(population_truth), n_distinct(population_truth$estimand),
            population_truth$truth[population_truth$term == "marginal"]))

#' Evaluate pooled estimates against known population truth.
#'
#' @param pooled A tibble of pooled estimates and confidence intervals.
#' @param truth_df A tibble of true values by estimand and term.
#' @return The pooled estimates augmented with bias and coverage indicators.
eval_estimand <- function(pooled, truth_df) {
  left_join(pooled, truth_df, by = c("estimand", "term")) %>%
    mutate(bias = Q_bar - truth,
           covered = (ci_lo <= truth) & (truth <= ci_hi))
}


# =============================================================================
# 8. MAIN SIMULATION LOOP WITH CHECKPOINTING
# =============================================================================
res_dir <- file.path(dirname(POP_PATH), "results_estimands")
dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)

K        <- length(gamma_imp_values)
reps_idx <- r_start:r_end
R        <- length(reps_idx)

acc <- list()   # One tibble per replication and scenario, both variants stacked.

ckpt_file <- file.path(res_dir, sprintf("ckpt_gen%.2f_b%d.rds", g_gen, batch))
j_done <- 0L
if (file.exists(ckpt_file)) {
  ck <- readRDS(ckpt_file)
  acc <- ck$acc; j_done <- ck$j_done
  cat(sprintf("  resuming from checkpoint: %d/%d reps done\n", j_done, R))
}

total_start <- proc.time()

for (j in seq_along(reps_idx)) {
  if (j <= j_done) next
  r <- reps_idx[j]
  
  # Sampling and missingness depend only on replication and generation parameter.
  set.seed(r * 1000 + round(1000 * g_gen))
  samp      <- as.data.frame(pop[sample(nrow(pop), n_sample), ])
  data_miss <- inject_missingness(samp, gamma_true = g_gen)
  
  for (k in seq_len(K)) {
    g_imp <- gamma_imp_values[k]
    
    AO <- run_mi(data_miss, g_imp, m, maxit, base_seed = r * 100,         variant = "AO")
    GL <- run_mi(data_miss, g_imp, m, maxit, base_seed = r * 100 + 50000, variant = "GL")
    
    eAO <- eval_estimand(AO, population_truth) %>% mutate(variant = "AO")
    eGL <- eval_estimand(GL, population_truth) %>% mutate(variant = "GL")
    
    acc[[length(acc) + 1]] <- bind_rows(eAO, eGL) %>%
      mutate(rep = r, g_gen = g_gen, g_imp = g_imp, spec = spec_labels[k])
    
    # Progress line: the marginal, plus average absolute bias over the
    # regression coefficients, so a problem in the regression shows up early.
    mAO <- eAO %>% filter(estimand == "marginal")
    mGL <- eGL %>% filter(estimand == "marginal")
    lAO <- mean(abs(eAO$bias[eAO$estimand == "logistic"]), na.rm = TRUE)
    lGL <- mean(abs(eGL$bias[eGL$estimand == "logistic"]), na.rm = TRUE)
    
    cat(sprintf("  rep %d (%d/%d) g_imp=%+.2f (%s) | marginal AO=%.3f (%+.4f) GL=%.3f (%+.4f) | logit mean|bias| AO=%.3f GL=%.3f\n",
                r, j, R, g_imp, spec_labels[k],
                mAO$Q_bar, mAO$bias, mGL$Q_bar, mGL$bias, lAO, lGL))
  }
  
  if (j %% CHECKPOINT_EVERY == 0 || j == R) {
    saveRDS(list(acc = acc, j_done = j, g_gen = g_gen,
                 gamma_imp_values = gamma_imp_values, spec_labels = spec_labels,
                 truth = population_truth, R = R), ckpt_file)
  }
}
cat(sprintf("\nTotal elapsed: %.1f sec\n", (proc.time() - total_start)["elapsed"]))


# =============================================================================
# 9. SUMMARY ACROSS REPLICATIONS
# =============================================================================
all_res <- bind_rows(acc)

#' Summarise replication-level results by estimand and term.
#'
#' @param res The stacked replication-level results.
#' @return A tibble of bias, RMSE, standard deviation, and coverage.
summarise_terms <- function(res) {
  res %>%
    group_by(g_gen, g_imp, spec, variant, estimand, term) %>%
    summarise(
      R_used    = sum(!is.na(Q_bar)),
      truth     = first(truth),
      mean_est  = mean(Q_bar, na.rm = TRUE),
      mean_bias = mean(bias, na.rm = TRUE),
      rmse      = sqrt(mean(bias^2, na.rm = TRUE)),
      sd_est    = sd(Q_bar, na.rm = TRUE),
      coverage  = mean(covered, na.rm = TRUE),
      .groups   = "drop"
    )
}

#' Collapse a term-level summary to one row per estimand family.
#'
#' Absolute bias is averaged across terms, because signed biases within a
#' family can cancel and hide a problem.
#'
#' @param term_summary Output of summarise_terms.
#' @return A tibble with one row per scenario, variant, and estimand.
summarise_families <- function(term_summary) {
  term_summary %>%
    group_by(g_gen, g_imp, spec, variant, estimand) %>%
    summarise(
      n_terms       = n(),
      mean_abs_bias = mean(abs(mean_bias)),
      mean_rmse     = mean(rmse),
      mean_coverage = mean(coverage),
      .groups       = "drop"
    )
}

term_summary   <- summarise_terms(all_res)
family_summary <- summarise_families(term_summary)

fam_order <- c("marginal", "trump_by_pid", "trump_by_ideology",
               "trump_by_edu")

cat(sprintf("\n=== g_gen=%.2f | R=%d | m=%d | maxit=%d ===\n", g_gen, R, m, maxit))

for (k in seq_len(K)) {
  cat(sprintf("\n-- g_imp=%+.2f (%s) --\n", gamma_imp_values[k], spec_labels[k]))
  for (vr in c("AO", "GL")) {
    for (fam in fam_order) {
      row <- family_summary %>%
        filter(abs(g_imp - gamma_imp_values[k]) < 1e-8,
               variant == vr, estimand == fam)
      if (nrow(row) == 0) next
      cat(sprintf("  [%s] %-18s terms=%2d  mean|bias|=%.4f  rmse=%.4f  cov=%.3f\n",
                  vr, fam, row$n_terms, row$mean_abs_bias,
                  row$mean_rmse, row$mean_coverage))
    }
  }
}

#' Print one estimand's results, one row per term, scenario, and variant.
#'
#' @param fam Name of the estimand family.
#' @param heading Text to print above the table.
#' @param cols Columns to show.
#' @return No value. Prints a table.
print_terms <- function(fam, heading,
                        cols = c("term", "g_imp", "spec", "variant", "truth",
                                 "mean_est", "mean_bias", "rmse", "coverage")) {
  cat("\n--", heading, "--\n")
  print(term_summary %>%
          filter(estimand == fam) %>%
          select(all_of(cols)) %>%
          arrange(term, g_imp, variant) %>%
          mutate(across(where(is.numeric), ~ round(.x, 4))),
        n = Inf)
}

cat("\n-- Marginal Trump support, by scenario --\n")
print(term_summary %>%
        filter(estimand == "marginal") %>%
        select(g_imp, spec, variant, truth, mean_est, mean_bias, rmse,
               sd_est, coverage) %>%
        arrange(g_imp, variant) %>%
        mutate(across(where(is.numeric), ~ round(.x, 4))),
      n = Inf)

print_terms("trump_by_pid",      "Trump support within party identification levels")
print_terms("trump_by_ideology", "Trump support within ideology levels")
print_terms("trump_by_edu",      "Trump support within education levels")

print_terms("logistic",          "Logistic regression coefficients")

cat(sprintf("\nTotal elapsed: %.1f sec\n", (proc.time() - total_start)["elapsed"]))

if (file.exists(ckpt_file)) file.remove(ckpt_file)


# =============================================================================
# 10. SAVE FINAL RESULTS
# -----------------------------------------------------------------------------
# The replication-level file is the one to use for plots and for pooling
# batches. The summary files are convenience views of the same numbers.
# =============================================================================
out_rds <- file.path(res_dir, sprintf("reps_gen%.2f_b%d.rds", g_gen, batch))
saveRDS(list(reps = all_res, term_summary = term_summary,
             family_summary = family_summary, truth = population_truth,
             g_gen = g_gen, gamma_imp_values = gamma_imp_values,
             spec_labels = spec_labels, R = R, m = m, maxit = maxit,
             n_sample = n_sample),
        out_rds)

write_csv(term_summary,
          file.path(res_dir, sprintf("term_summary_gen%.2f_b%d.csv", g_gen, batch)))

cat(sprintf("Saved: %s\n", out_rds))

