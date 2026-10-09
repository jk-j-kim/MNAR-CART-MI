# =============================================================================
# BUILDING A SYNTHETIC POPULATION FROM CPS AND ANES 2020
# -----------------------------------------------------------------------------
# The synthetic population combines two sources. The Current Population Survey
# (CPS) supplies the demographic backbone: one record per adult, with age, sex,
# race, and education. The American National Election Studies (ANES) 2020
# supplies everything else. Political variables are absent from the CPS, so
# they are drawn for each CPS person from a model fitted to ANES respondents.
#
# Each political variable is modeled once, on the respondents who answered that
# question, and the fitted tree is then used twice. It fills in the respondents
# who did not answer, so that later models have complete predictors, and it
# generates a value for every CPS record.
#
# Variables are handled in sequence, and each model uses the demographics and
# every variable handled before it:
#
#   income     ~ age + sex + race + education
#   party ID   ~ age + sex + race + education + income
#   ideology   ~ age + sex + race + education + income + party ID
#   vote       ~ age + sex + race + education + income + party ID + ideology
#
# -----------------------------------------------------------------------------
# THE STEPS
# -----------------------------------------------------------------------------
#
#   STEP 1   Demographics
#            Translate age, sex, race, and education from the CPS coding into
#            the ANES coding. Every CPS record carries all four.
#            Set aside the ANES respondents missing any of the four, so that
#            the demographics are complete on both sides from here on.
#
#   STEP 2   Income
#            Fit on   ANES respondents with income observed.
#            Predict  age + sex + race + education.
#            Use it   (a) to fill in ANES respondents missing income
#                     (b) to draw income for every CPS record
#
#   STEP 3   Party identification
#            Fit on   ANES respondents with party ID observed, income now
#                     complete.
#            Predict  age + sex + race + education + income.
#            Use it   (a) to fill in ANES respondents missing party ID
#                     (b) to draw party ID for every CPS record
#
#   STEP 4   Ideology
#            Fit on   ANES respondents with ideology observed, income and
#                     party ID now complete.
#            Predict  age + sex + race + education + income + party ID.
#            Use it   (a) to fill in ANES respondents missing ideology
#                     (b) to draw ideology for every CPS record
#
#   STEP 5   Vote
#            Fit on   ANES respondents who reported a vote, all covariates
#                     now complete.
#            Predict  age + sex + race + education + income + party ID +
#                     ideology.
#            Use it   to draw a vote for every CPS record.
#
#            Vote is the one variable whose survey gaps stay open. Vote
#            nonresponse is the subject of the study and is taken to be
#            missing not at random, so respondents who declined are never
#            filled in and never serve as donors.
#
# -----------------------------------------------------------------------------
#
# Filling in the survey gaps lets each model use the whole sample. Requiring
# every answer to be present would discard a large share of the respondents,
# and would assume that the people who answered everything resemble the people
# who did not.
#
# Each record is assigned to a terminal node, and one observed value is drawn
# at random from the respondents in that node. Drawing rather than predicting
# keeps the spread of each variable realistic instead of collapsing everyone in
# a node onto the same value.
#
# Output: one CSV with the columns the simulation script expects.
# =============================================================================
library(tidyverse)
library(rpart)
library(here)  

# File locations and the random seed ------------------------------------------
# All paths are relative to the mnarcart_supplementals folder.
CPS_PATH  <- here("oct20pub.csv.gz")
ANES_PATH <- here("anes_timeseries_2020_csv_20220210.csv")
OUT_PATH  <- here("sim_pop.csv")

SEED <- 42
set.seed(SEED)

# Tree settings, shared by every model.
TREE_CTRL <- rpart.control(minsplit = 20, minbucket = 7, cp = 0.001)

# =============================================================================
# 1. PREPARE THE ANES VARIABLES
# =============================================================================

anes_raw <- read_csv(ANES_PATH, show_col_types = FALSE)

# ANES variable numbers and what they hold.
anes_vars <- c(
  "V202110x",   # Vote choice, post-election
  "V201507x",   # Age
  "V201511x",   # Education
  "V201549x",   # Race and ethnicity
  "V201600",    # Sex
  "V201200",    # Ideology, seven levels
  "V201617x",   # Income, pre-election
  "V201231x"    # Party identification, seven levels
)

# Vote nonresponse in ANES, recorded before invalid codes become NA. This is
# the nonresponse the simulation study is about, so it is worth reporting.
vote_status <- anes_raw %>%
  transmute(status = case_when(
    V202110x == -9                  ~ "refused",
    V202110x == -8                  ~ "don't know",
    V202110x == -1                  ~ "not asked",
    !is.na(V202110x) & V202110x > 0 ~ "observed",
    TRUE                            ~ "other missing"
  ))

cat("--- ANES vote response status ---\n")
print(table(vote_status$status))

# Negative codes mark refusals, don't-knows, and cases outside the universe.
# Code 99 marks missing values in several ANES variables. All become NA.
anes <- anes_raw %>%
  select(all_of(anes_vars)) %>%
  mutate(across(everything(),
                ~ .x %>% na_if(-9) %>% na_if(-8) %>% na_if(-5) %>%
                  na_if(-4) %>% na_if(-2) %>% na_if(-1) %>% na_if(99)))

# Category labels. The synthetic population uses the same labels throughout.
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

anes <- anes %>%
  transmute(
    # Vote is Trump against everything else, including Biden and third parties.
    vote      = factor(if_else(V202110x == 2, "Trump", "Other"),
                       levels = c("Other", "Trump")),
    age       = as.numeric(V201507x),
    education = factor(V201511x, levels = 1:5, labels = edu_levels,
                       ordered = TRUE),
    race      = factor(V201549x, levels = 1:6, labels = race_levels),
    sex       = factor(V201600,  levels = 1:2, labels = sex_levels),
    ideology  = factor(V201200,  levels = 1:7, labels = ideo_levels,
                       ordered = TRUE),
    income    = as.numeric(V201617x),   # 22-point ordinal scale
    party_id  = factor(V201231x, levels = 1:7, labels = pid_levels,
                       ordered = TRUE)
  ) %>%
  as.data.frame()

cat("\nANES rows:", nrow(anes), "\n")
cat("Missing values per variable:\n")
print(colSums(is.na(anes)))

# The four demographics are the starting point for every model, so they must be
# present. Respondents missing any of them are set aside. Their number is small
# relative to the sample; if it is not, this decision needs revisiting.
n_before <- nrow(anes)
anes <- anes %>% drop_na(age, education, race, sex)
cat(sprintf("\nRespondents set aside for missing demographics: %d of %d (%.1f%%)\n",
            n_before - nrow(anes), n_before,
            100 * (n_before - nrow(anes)) / n_before))
cat("Respondents available:", nrow(anes), "\n")


# =============================================================================
# 2. BUILD THE DEMOGRAPHIC BACKBONE FROM CPS
# -----------------------------------------------------------------------------
# The CPS codes each variable differently from ANES, so each one is translated
# into the ANES categories. The two files must agree on categories, because the
# trees are fitted on ANES and applied to the CPS.
# =============================================================================

cps_raw <- read_csv(CPS_PATH, show_col_types = FALSE)

# ANES surveys voting-age adults, so the CPS is restricted to match.
cps_adults <- cps_raw %>% filter(prtage >= 18)

cat("\nCPS rows:", nrow(cps_raw), "| adults 18 and over:", nrow(cps_adults), "\n")

# CPS records Hispanic ethnicity separately from race. ANES folds it into one
# variable, so Hispanic ethnicity takes priority over the racial category.
hispanic <- if ("pehspnon" %in% names(cps_adults)) {
  cps_adults$pehspnon
} else {
  rep(NA_real_, nrow(cps_adults))
}
cat("Hispanic indicator available in CPS:", !all(is.na(hispanic)), "\n")

pop <- cps_adults %>%
  transmute(
    
    # CPS reports 80 for ages 80 to 84 and 85 for 85 and over. ANES tops out at
    # 80, so both are collapsed to 80.
    a_age = ifelse(prtage >= 80, 80L, as.integer(prtage)),
    
    a_sex = factor(case_when(pesex == 1 ~ "Male",
                             pesex == 2 ~ "Female"),
                   levels = sex_levels),
    
    # CPS race codes: 1 White, 2 Black, 3 Native American, 4 Asian,
    # 5 Hawaiian or Pacific Islander, 6 through 26 multiple races.
    # Hawaiian or Pacific Islander joins the multiple-race group, matching the
    # residual ANES category.
    a_race = factor(case_when(
      hispanic == 1                  ~ "Hispanic",
      ptdtrace == 1                  ~ "White",
      ptdtrace == 2                  ~ "Black",
      ptdtrace == 3                  ~ "Native American",
      ptdtrace == 4                  ~ "Asian/PI",
      ptdtrace >= 5 & ptdtrace <= 26 ~ "Other/Multiple"
    ), levels = race_levels),
    
    # CPS education codes: 31 to 38 no diploma, 39 high school or equivalent,
    # 40 to 42 some college or associate degree, 43 bachelor's,
    # 44 to 46 graduate or professional degree.
    a_edu = factor(case_when(
      peeduca >= 31 & peeduca <= 38 ~ "Less than HS",
      peeduca == 39                 ~ "HS diploma",
      peeduca >= 40 & peeduca <= 42 ~ "Some post-high school",
      peeduca == 43                 ~ "Bachelor's",
      peeduca >= 44 & peeduca <= 46 ~ "Graduate degree"
    ), levels = edu_levels, ordered = TRUE)
  ) %>%
  as.data.frame()

cat("\nMissing demographics before filtering:\n")
print(colSums(is.na(pop)))

pop <- pop %>% drop_na()
cat("Population rows with complete demographics:", nrow(pop), "\n")


# =============================================================================
# 3. CHECK THAT THE TWO SOURCES AGREE
# -----------------------------------------------------------------------------
# The demographic distributions do not have to match. The CPS is a census-based
# sample of all adults and the ANES is a survey sample, so some difference is
# expected. Large gaps would mean a coding mistake in the translation above.
# =============================================================================

#' Print one variable's distribution in both sources.
#'
#' @param anes_vec The ANES vector.
#' @param pop_vec The population vector.
#' @param label Name to print.
#' @return No value. Prints two tables of proportions.
compare_dist <- function(anes_vec, pop_vec, label) {
  cat("\n---", label, "---\n")
  cat("ANES:\n"); print(round(prop.table(table(anes_vec, useNA = "ifany")), 3))
  cat("POP: \n"); print(round(prop.table(table(pop_vec,  useNA = "ifany")), 3))
}

cat("\n\n========== DEMOGRAPHICS: ANES AND POPULATION ==========\n")
compare_dist(anes$sex,       pop$a_sex,  "Sex")
compare_dist(anes$race,      pop$a_race, "Race")
compare_dist(anes$education, pop$a_edu,  "Education")

cat("\n--- Age ---\n")
cat("ANES mean:", round(mean(anes$age), 1), "\n")
cat("POP  mean:", round(mean(pop$a_age), 1), "\n")

# A race and sex combination that appears in the population but never in ANES
# has no respondents to draw from. The trees route such rows through surrogate
# splits, so they still receive a value, but it is worth knowing how many.
anes_combos <- anes %>% distinct(a_race = race, a_sex = sex)

missing_combos <- pop %>%
  distinct(a_race, a_sex) %>%
  anti_join(anes_combos, by = c("a_race", "a_sex"))

if (nrow(missing_combos) == 0) {
  cat("\nEvery race and sex combination in the population appears in ANES.\n")
} else {
  n_affected <- pop %>% semi_join(missing_combos, by = c("a_race", "a_sex")) %>% nrow()
  cat("\nCombinations present in the population but not in ANES:\n")
  print(missing_combos)
  cat("Rows affected:", n_affected,
      sprintf("(%.2f%% of the population)\n", 100 * n_affected / nrow(pop)))
}


# =============================================================================
# 4. FIT A TREE, THEN DRAW DONORS FROM IT
# -----------------------------------------------------------------------------
# One tree per variable, fitted on the respondents who answered the question.
# Filling a survey gap and generating a value for a CPS record are the same
# operation applied to different rows: place the row in the tree, then draw one
# of its neighbours' answers.
# =============================================================================

#' Find the terminal node of each row used to fit a tree.
#'
#' @param tree A fitted rpart tree.
#' @return One node number per training row.
leaf_ids_train <- function(tree) as.integer(rownames(tree$frame))[tree$where]

#' Find the terminal node of each new row.
#'
#' A tree stores one row per node, named by the node number, and a prediction
#' returns that row's fitted value. Replacing the fitted values with the node
#' numbers therefore makes the prediction return the node itself. Rows with a
#' missing or unfamiliar predictor value follow the tree's surrogate splits, so
#' every row reaches a node.
#'
#' @param tree A fitted rpart tree.
#' @param newdata Rows to place in the tree.
#' @return One node number per row of newdata.
leaf_ids_new <- function(tree, newdata) {
  t2 <- tree
  t2$frame$yval <- as.numeric(rownames(t2$frame))
  as.integer(predict(t2, newdata = newdata, type = "vector"))
}

#' Draw values from a pool, with replacement.
#'
#' Indices are drawn rather than values, because sampling directly from a pool
#' holding one number would sample from the sequence up to that number instead.
#'
#' @param pool Values to draw from.
#' @param n Number of draws.
#' @return n drawn values.
sample_pool <- function(pool, n) pool[sample.int(length(pool), n, replace = TRUE)]

#' Fit one donor tree on the respondents who answered a question.
#'
#' The tree is kept together with its training rows, because the training rows
#' are the donor pool.
#'
#' @param data The survey data.
#' @param outcome Name of the outcome column.
#' @param predictors Names of the predictor columns, all complete.
#' @param method Either "class" for categories or "anova" for numbers.
#' @param label Text describing the model, for the progress line.
#' @return A list holding the tree, its training rows, and the variable names.
fit_donor_tree <- function(data, outcome, predictors, method, label = outcome) {
  answered <- !is.na(data[[outcome]])
  train    <- as.data.frame(data[answered, c(predictors, outcome), drop = FALSE])
  form     <- as.formula(paste(outcome, "~", paste(predictors, collapse = " + ")))
  tree     <- rpart(form, data = train, method = method, control = TREE_CTRL)
  
  cat(sprintf("\n  %s\n", label))
  cat(sprintf("    fitted on %d respondents who answered | terminal nodes: %d\n",
              nrow(train), sum(tree$frame$var == "<leaf>")))
  
  list(tree = tree, train = train, outcome = outcome, predictors = predictors)
}

#' Give each row a value drawn from the respondents in its terminal node.
#'
#' @param model A model from fit_donor_tree.
#' @param newdata Rows needing a value, holding the model's predictor columns.
#' @return One drawn value per row of newdata, keeping the outcome's type.
draw_donors <- function(model, newdata) {
  y      <- model$train[[model$outcome]]
  is_fac <- is.factor(y)
  values <- if (is_fac) as.character(y) else as.numeric(y)
  
  id_train <- leaf_ids_train(model$tree)
  id_new   <- leaf_ids_new(model$tree,
                           as.data.frame(newdata[, model$predictors, drop = FALSE]))
  
  out <- if (is_fac) character(nrow(newdata)) else numeric(nrow(newdata))
  
  for (node in unique(id_new)) {
    rows <- id_new == node
    pool <- values[id_train == node]
    out[rows] <- sample_pool(pool, sum(rows))
  }
  
  if (is_fac) factor(out, levels = levels(y), ordered = is.ordered(y)) else out
}

#' Fill the gaps in a survey variable using its own donor tree.
#'
#' @param model A model from fit_donor_tree.
#' @param data The survey data.
#' @return The variable with its missing values filled in.
fill_gaps <- function(model, data) {
  out  <- data[[model$outcome]]
  gaps <- is.na(out)
  if (!any(gaps)) {
    cat("    nothing missing in the survey\n")
    return(out)
  }
  out[gaps] <- draw_donors(model, data[gaps, , drop = FALSE])
  cat(sprintf("    filled in %d missing survey values\n", sum(gaps)))
  out
}


# =============================================================================
# 5. DRAW EACH VARIABLE IN SEQUENCE
# -----------------------------------------------------------------------------
# For each variable: fit one tree on the respondents who answered, use it to
# fill the survey gaps, and use the same tree to generate a value for every CPS
# record. The filled survey variable then becomes a predictor for the next one.
#
# The population columns carry the ANES names while the models run, so that one
# formula applies to both data sets.
# =============================================================================

cat("\n\n========== DRAWING THE POLITICAL VARIABLES ==========\n")

# Distributions among respondents who answered, kept for comparison after the
# gaps are filled.
observed_dist <- list(
  income   = anes$income[!is.na(anes$income)],
  party_id = anes$party_id[!is.na(anes$party_id)],
  ideology = anes$ideology[!is.na(anes$ideology)]
)

stems <- c("age", "education", "race", "sex")

pop_x <- data.frame(age       = pop$a_age,
                    education = pop$a_edu,
                    race      = pop$a_race,
                    sex       = pop$a_sex)

# Income ----------------------------------------------------------------------
model_income  <- fit_donor_tree(anes, "income", stems, "anova",
                                "Income, from demographics")
anes$income   <- as.integer(fill_gaps(model_income, anes))
pop$ap_income <- as.integer(draw_donors(model_income, pop_x))
pop_x$income  <- pop$ap_income

# Party identification --------------------------------------------------------
model_pid      <- fit_donor_tree(anes, "party_id", c(stems, "income"), "class",
                                 "Party identification, from demographics and income")
anes$party_id  <- fill_gaps(model_pid, anes)
pop$ap_partyid <- draw_donors(model_pid, pop_x)
pop_x$party_id <- pop$ap_partyid

# Ideology --------------------------------------------------------------------
model_ideology  <- fit_donor_tree(anes, "ideology", c(stems, "income", "party_id"),
                                  "class",
                                  "Ideology, from demographics, income, party ID")
anes$ideology   <- fill_gaps(model_ideology, anes)
pop$ap_ideology <- draw_donors(model_ideology, pop_x)
pop_x$ideology  <- pop$ap_ideology

# Vote ------------------------------------------------------------------------
# The survey gaps are deliberately left open here. Only the population is drawn.
model_vote  <- fit_donor_tree(anes, "vote",
                              c(stems, "income", "party_id", "ideology"), "class",
                              "Vote, from all covariates")
pop$ap_vote <- draw_donors(model_vote, pop_x)
cat("    survey vote left as it is\n")

cat("\nMissing values in ANES after this step:\n")
print(colSums(is.na(anes)))

# Filling in gaps should not shift a variable's distribution much. A large
# shift means the model is doing more than completing the record.
cat("\n--- Income: answered only, then filled in ---\n")
cat("answered mean:", round(mean(observed_dist$income), 2),
    "| sd:", round(sd(observed_dist$income), 2), "\n")
cat("filled   mean:", round(mean(anes$income), 2),
    "| sd:", round(sd(anes$income), 2), "\n")

cat("\n--- Party identification: answered only, then filled in ---\n")
print(round(rbind(answered = prop.table(table(observed_dist$party_id)),
                  filled   = prop.table(table(anes$party_id))), 3))

cat("\n--- Ideology: answered only, then filled in ---\n")
print(round(rbind(answered = prop.table(table(observed_dist$ideology)),
                  filled   = prop.table(table(anes$ideology))), 3))


# =============================================================================
# 6. CHECK THE POPULATION AGAINST ANES
# -----------------------------------------------------------------------------
# Two questions. Do the variables have realistic distributions, and do they
# relate to each other the way they do in the survey? The second question is
# the important one, because the simulation study estimates associations, not
# just marginal shares.
# =============================================================================

anes_voters <- anes %>% filter(!is.na(vote))

cat("\n\n========== DISTRIBUTIONS: ANES AND POPULATION ==========\n")
compare_dist(anes$party_id,   pop$ap_partyid,  "Party identification")
compare_dist(anes$ideology,   pop$ap_ideology, "Ideology")
compare_dist(anes_voters$vote, pop$ap_vote,    "Vote")

cat("\n--- Income ---\n")
cat("ANES mean:", round(mean(anes$income), 2),
    "| sd:",      round(sd(anes$income),   2), "\n")
cat("POP  mean:", round(mean(pop$ap_income), 2),
    "| sd:",      round(sd(pop$ap_income),   2), "\n")

cat("\n\n========== ASSOCIATIONS: ANES AND POPULATION ==========\n")

#' Rank correlation between two variables, treating categories as scores.
#'
#' @param x First variable.
#' @param y Second variable.
#' @return A single correlation.
score_cor <- function(x, y) {
  suppressWarnings(cor(as.numeric(x), as.numeric(y),
                       method = "spearman", use = "complete.obs"))
}

associations <- tibble(
  pair = c("party ID and ideology", "party ID and income",
           "ideology and income", "party ID and education",
           "ideology and education"),
  anes = c(score_cor(anes$party_id, anes$ideology),
           score_cor(anes$party_id, anes$income),
           score_cor(anes$ideology, anes$income),
           score_cor(anes$party_id, anes$education),
           score_cor(anes$ideology, anes$education)),
  pop  = c(score_cor(pop$ap_partyid,  pop$ap_ideology),
           score_cor(pop$ap_partyid,  pop$ap_income),
           score_cor(pop$ap_ideology, pop$ap_income),
           score_cor(pop$ap_partyid,  pop$a_edu),
           score_cor(pop$ap_ideology, pop$a_edu))
) %>%
  mutate(difference = pop - anes,
         across(where(is.numeric), ~ round(.x, 3)))

print(associations)


cat("\n--- Trump share by party identification ---\n")
print(pop %>%
        group_by(ap_partyid) %>%
        summarise(n = n(), trump = round(mean(ap_vote == "Trump"), 3),
                  .groups = "drop"))

cat("\n--- Trump share by ideology ---\n")
print(pop %>%
        group_by(ap_ideology) %>%
        summarise(n = n(), trump = round(mean(ap_vote == "Trump"), 3),
                  .groups = "drop"))

cat("\n--- Trump share by education ---\n")
print(pop %>%
        group_by(a_edu) %>%
        summarise(n = n(), trump = round(mean(ap_vote == "Trump"), 3),
                  .groups = "drop"))

cat("\nANES respondents who reported a vote, Trump share:",
    round(mean(anes_voters$vote == "Trump"), 4), "\n")
cat("Population, Trump share:                          ",
    round(mean(pop$ap_vote == "Trump"), 4), "\n")


# =============================================================================
# 7. SAVE THE POPULATION
# =============================================================================

sim_pop <- pop %>%
  select(a_age, a_sex, a_race, a_edu, ap_income, ap_partyid, ap_ideology, ap_vote)

stopifnot(!anyNA(sim_pop))

dir.create(dirname(OUT_PATH), showWarnings = FALSE, recursive = TRUE)
write_csv(sim_pop, OUT_PATH)

cat("\nSaved:", OUT_PATH, "|", nrow(sim_pop), "rows and", ncol(sim_pop), "columns\n")