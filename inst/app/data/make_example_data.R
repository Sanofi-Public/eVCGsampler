# Generates the two example files used to try the app out:
#
#   data/example_POOL.csv  - 12 historical control studies, 240 animals
#   data/example_TG.csv    - the treated animals of the current study, 40
#
# Study-11 and Study-12 are deliberately shifted on the outcome variables so
# the outlier study test on the Pre-tests step has something to find, and the
# pool covariates are shifted against the TG so the balancing has work to do.
#
#   source("data/make_example_data.R")

set.seed(20260807)

n_study   <- 12
per_study <- 20

pool <- data.frame(
    animal_id = sprintf("P%03d", seq_len(n_study * per_study)),
    study     = rep(sprintf("Study-%02d", seq_len(n_study)), each = per_study),
    sex       = sample(c("M", "F"), n_study * per_study, replace = TRUE),
    body_weight = rnorm(n_study * per_study, 305, 28),
    age_days    = round(rnorm(n_study * per_study, 78, 7)),
    food_intake = rnorm(n_study * per_study, 24, 3.5),
    stringsAsFactors = FALSE
)

# Outcome / technical variables - not to be balanced, but used by the outlier
# study test. Study-11 is shifted on ALT, Study-12 on both.
pool$ALT <- rnorm(nrow(pool), 42, 6)
pool$AST <- rnorm(nrow(pool), 88, 11)
pool$ALT[pool$study == "Study-11"] <- rnorm(per_study, 71, 6)
pool$ALT[pool$study == "Study-12"] <- rnorm(per_study, 68, 6)
pool$AST[pool$study == "Study-12"] <- rnorm(per_study, 132, 11)

# The treated animals: same covariates, systematically lighter and younger
# than the pool, so the pre-test finds an imbalance.
tg <- data.frame(
    animal_id = sprintf("T%03d", seq_len(40)),
    study     = "Current",
    sex       = rep(c("M", "F"), each = 20),
    body_weight = rnorm(40, 288, 16),
    age_days    = round(rnorm(40, 72, 4)),
    food_intake = rnorm(40, 22, 2.2),
    ALT         = rnorm(40, 44, 6),
    AST         = rnorm(40, 90, 11),
    stringsAsFactors = FALSE
)

num <- vapply(pool, is.numeric, logical(1))
pool[num] <- lapply(pool[num], function(x) round(x, 2))
num <- vapply(tg, is.numeric, logical(1))
tg[num] <- lapply(tg[num], function(x) round(x, 2))

dir.create("data", showWarnings = FALSE)
write.csv(pool, "data/example_POOL.csv", row.names = FALSE)
write.csv(tg,   "data/example_TG.csv",   row.names = FALSE)

# A single-file variant: pool and treated animals already stacked, with the
# group variable in place.
combined <- rbind(
    transform(pool, group = "POOL"),
    transform(tg,   group = "TG")
)
write.csv(combined, "data/example_combined.csv", row.names = FALSE)

message("Wrote data/example_POOL.csv, data/example_TG.csv and data/example_combined.csv")
