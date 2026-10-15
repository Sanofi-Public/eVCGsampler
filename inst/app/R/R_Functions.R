# =====================================================================
# R_Functions.R - shared helpers for the eVCGsampler Shiny app
#
# Everything in here is pure R: no Shiny reactivity, no side effects on
# the session. That keeps the statistics testable outside the app and
# keeps the module servers readable.
# =====================================================================


#' Null-coalescing helper for optional inputs
`%||%` <- function(a, b) if (is.null(a)) b else a

#' Has an action button been pressed yet?
#'
#' A live Shiny session initialises an actionButton to 0, but the input is
#' NULL until the UI reports in - and NULL inside `if (x == 0)` is an error
#' that aborts the whole reactive flush. Always ask through this helper.
clicked <- function(x) !is.null(x) && length(x) == 1 && !is.na(x) && x > 0


# ------ FILE READING ---------------------------------------------------

#' Read an uploaded data file (csv / txt / xlsx / xls)
#'
#' The separator and the decimal mark of delimited files are guessed from
#' the first lines, because study exports arrive both as "," / "." and as
#' ";" / "," depending on the regional settings of the exporting system.
#'
#' @param path     Path to the file on disk (fileInput datapath).
#' @param filename Original file name - only used for the extension.
#' @return A data.frame.
read_upload <- function(path, filename) {
    ext <- tolower(tools::file_ext(filename))

    if (ext %in% c("csv", "txt", "tsv")) {
        first_lines <- readLines(path, n = 10, warn = FALSE)
        n_comma <- sum(vapply(first_lines, function(x) lengths(gregexpr(",",  x, fixed = TRUE)), numeric(1)))
        n_semi  <- sum(vapply(first_lines, function(x) lengths(gregexpr(";",  x, fixed = TRUE)), numeric(1)))
        n_tab   <- sum(vapply(first_lines, function(x) lengths(gregexpr("\t", x, fixed = TRUE)), numeric(1)))

        sep <- if (n_tab > n_comma && n_tab > n_semi) "\t" else if (n_semi > n_comma) ";" else ","
        # With ";" as field separator the decimal mark is almost always ","
        # (German / French Excel export). Everything else uses ".".
        dec <- if (identical(sep, ";")) "," else "."

        out <- utils::read.csv(path, sep = sep, dec = dec, stringsAsFactors = FALSE,
                               check.names = TRUE, na.strings = c("NA", "", ".", "NaN"))
        # A ";"-separated file written with "." decimals still has to work,
        # so retry with "." if the "," guess produced no numeric column. This
        # only catches the case where every column is affected; a file mixing
        # whole-number columns (which parse under either "dec") with decimal
        # ones (which don't) needs the per-column fallback below instead.
        if (identical(dec, ",") && !any(vapply(out, is.numeric, logical(1)))) {
            out <- utils::read.csv(path, sep = sep, dec = ".", stringsAsFactors = FALSE,
                                   check.names = TRUE, na.strings = c("NA", "", ".", "NaN"))
        }
        return(coerce_numeric_like(as.data.frame(out)))
    }

    if (ext == "xlsx") return(coerce_numeric_like(as.data.frame(readxl::read_xlsx(path))))
    if (ext == "xls")  return(coerce_numeric_like(as.data.frame(readxl::read_xls(path))))

    stop("Unsupported file type: .", ext,
         " - please upload a .csv, .txt, .xlsx or .xls file.", call. = FALSE)
}


#' Rescue numeric columns that were read in as text
#'
#' A wrong decimal-mark guess (e.g. a ";"-separated file that turns out to use
#' "." decimals after all - see the retry above, which only catches the case
#' where it breaks *every* column) leaves the affected columns as character,
#' even though the values are perfectly numeric. Lab panels are the classic
#' case: whole-number analytes parse fine under any "dec" guess and mask the
#' problem, while the decimal ones next to them silently become text.
#'
#' Every character column with more than 4 distinct values is tried against
#' both a "." and a "," decimal mark; the version that parses more values
#' wins. The column is only converted if that keeps at least half of its
#' original non-missing values - anything below that is treated as
#' genuinely non-numeric (an ID, a study label, free text) and left as text.
#' Columns with 4 or fewer distinct values are left alone on purpose - those
#' are almost always a category (sex, dose group, ...) that happens to be
#' encoded as numbers, and turning them into a covariate/stratum candidate is
#' exactly what should happen; nothing to fix there.
#'
#' @param data A data.frame, as produced by read.csv() / readxl.
coerce_numeric_like <- function(data) {
    is_chr <- vapply(data, is.character, logical(1))
    if (!any(is_chr)) return(data)

    for (v in names(data)[is_chr]) {
        x <- trimws(data[[v]])
        x[!nzchar(x)] <- NA_character_
        n_before <- sum(!is.na(x))
        if (n_before == 0) next
        if (length(unique(x[!is.na(x)])) <= 4) next

        cand_dot   <- suppressWarnings(as.numeric(x))
        cand_comma <- suppressWarnings(as.numeric(gsub(",", ".", x, fixed = TRUE)))
        n_dot   <- sum(!is.na(cand_dot))
        n_comma <- sum(!is.na(cand_comma))

        best   <- if (n_comma > n_dot) cand_comma else cand_dot
        n_best <- max(n_dot, n_comma)
        if (n_best >= 0.5 * n_before) data[[v]] <- best
    }
    data
}


# ------ VARIABLE CLASSIFICATION ---------------------------------------

#' Number of distinct non-missing values per column
n_unique <- function(data) {
    vapply(data, function(x) length(unique(x[!is.na(x)])), numeric(1))
}

#' Short human-readable type label for a column
var_type <- function(x) {
    if (is.factor(x))                    return("factor")
    if (is.logical(x))                   return("logical")
    if (is.character(x))                 return("character")
    if (is.integer(x))                   return("integer")
    if (is.numeric(x))                   return("numeric")
    class(x)[1]
}

#' Candidate group variables: exactly two distinct non-missing values
binary_vars <- function(data) {
    u <- n_unique(data)
    names(data)[u == 2]
}

#' Candidate stratum variables: few distinct values, and never continuous
#'
#' The upper bound of 8 mirrors eVCGsampler, which refuses more than eight
#' strata (VCG_sampler / energy_test both stop at nlevels > 8).
stratum_vars <- function(data, max_levels = 8) {
    u <- n_unique(data)
    keep <- u >= 2 & u <= max_levels
    names(data)[keep]
}

#' Candidate animal-ID variables
#'
#' An identifier has to be complete (no missing values) and has to single out
#' one animal. Uniqueness is checked *within each group* rather than over the
#' whole data set, because POOL and TG usually come from different files and
#' each may number its animals from 1 - an ID that repeats across the two
#' files still identifies a POOL animal unambiguously, and only POOL animals
#' can end up in the VCG. A globally unique ID passes this test as well.
#'
#' @param data  The combined data frame.
#' @param group Optional name of the POOL / TG variable.
#' @return Column names, in the order they appear - the first one is what the
#'   app pre-selects.
id_vars <- function(data, group = NULL) {
    if (is.null(data) || nrow(data) == 0) return(character(0))
    cand <- names(data)[!vapply(data, anyNA, logical(1))]
    cand <- setdiff(cand, group)
    if (length(cand) == 0) return(character(0))

    g <- if (!is.null(group) && group %in% names(data)) {
        as.factor(data[[group]])
    } else {
        factor(rep(1L, nrow(data)))
    }
    keep <- vapply(cand, function(v) {
        all(tapply(data[[v]], g, function(x) !anyDuplicated(x)))
    }, logical(1))
    cand[keep]
}


#' Move the ID column to the front of a data frame
#'
#' The identifier is the first thing you look for in a result table - it is
#' what lets you pull the sampled animals out of the source data.
id_first <- function(d, id = NULL) {
    if (is.null(d) || is.null(id) || !nzchar(id) || !id %in% names(d)) return(d)
    d[, c(id, setdiff(names(d), id)), drop = FALSE]
}


#' Candidate covariates: anything that can be turned into a number
#'
#' eVCGsampler standardises covariates with robust_scale(), so numeric and
#' low-cardinality categorical variables both work; free text and ID-like
#' columns (one distinct value per row) do not.
covariate_vars <- function(data) {
    u <- n_unique(data)
    ok <- vapply(seq_along(data), function(i) {
        x <- data[[i]]
        is.numeric(x) || is.logical(x) || is.factor(x) ||
            (is.character(x) && u[i] <= 8)
    }, logical(1))
    names(data)[ok & u >= 2]
}


# ------ DATA OVERVIEW --------------------------------------------------

#' Per-variable overview table shown on the Import step
#'
#' @param data   The (possibly combined) data frame.
#' @param origin Named character vector mapping variable -> "both" /
#'   "POOL only" / "TG only" / "single file". Optional.
variable_overview <- function(data, origin = NULL) {
    if (is.null(data) || ncol(data) == 0) return(NULL)

    rows <- lapply(names(data), function(v) {
        x <- data[[v]]
        n_miss <- sum(is.na(x))
        summ <- if (is.numeric(x) && any(!is.na(x))) {
            paste0("median ", signif(stats::median(x, na.rm = TRUE), 4),
                   " [", signif(min(x, na.rm = TRUE), 4), " - ",
                   signif(max(x, na.rm = TRUE), 4), "]")
        } else {
            lv <- unique(as.character(x[!is.na(x)]))
            if (length(lv) == 0) "-" else
                paste0(paste(utils::head(lv, 4), collapse = ", "),
                       if (length(lv) > 4) paste0(", ... (", length(lv), " levels)") else "")
        }
        data.frame(
            Variable = v,
            Type     = var_type(x),
            Distinct = length(unique(x[!is.na(x)])),
            Missing  = paste0(n_miss, " (", round(100 * n_miss / max(1, length(x))), "%)"),
            Summary  = summ,
            stringsAsFactors = FALSE
        )
    })
    out <- do.call(rbind, rows)

    if (!is.null(origin)) {
        out$Source <- ifelse(out$Variable %in% names(origin),
                             origin[out$Variable], "-")
        out <- out[, c("Variable", "Source", "Type", "Distinct", "Missing", "Summary")]
    }
    # Suggest a role so the user sees at a glance which columns are usable.
    bins  <- binary_vars(data)
    covs  <- covariate_vars(data)
    strat <- stratum_vars(data)
    out$Usable_as <- vapply(out$Variable, function(v) {
        roles <- c(if (v %in% bins) "group", if (v %in% covs) "covariate",
                   if (v %in% strat) "stratum")
        if (length(roles) == 0) "not usable" else paste(roles, collapse = " / ")
    }, character(1))
    out
}


#' How the POOL and TG uploads line up, column by column
#'
#' combine_data() silently keeps only the intersection of the column names,
#' so this table makes explicit which variables survive the merge, which are
#' dropped, and where the two files disagree about the type.
combine_check <- function(pool, tg) {
    if (is.null(pool) || is.null(tg)) return(NULL)

    all_names <- union(names(pool), names(tg))
    rows <- lapply(all_names, function(v) {
        in_p <- v %in% names(pool)
        in_t <- v %in% names(tg)
        tp <- if (in_p) var_type(pool[[v]]) else "-"
        tt <- if (in_t) var_type(tg[[v]])   else "-"

        if (in_p && in_t) {
            num_p <- in_p && (is.numeric(pool[[v]]) || is.logical(pool[[v]]))
            num_t <- in_t && (is.numeric(tg[[v]])   || is.logical(tg[[v]]))
            if (identical(tp, tt)) {
                status <- "combined"
                note   <- "Name and type match in both files."
            } else if (num_p == num_t) {
                status <- "combined"
                note   <- "Types differ slightly but both are compatible."
            } else {
                status <- "type mismatch"
                note   <- paste0("POOL is ", tp, " but TG is ", tt,
                                 " - check for text in a numeric column (e.g. \"<0.1\").")
            }
        } else if (in_p) {
            status <- "dropped"
            note   <- "Only in the POOL file - combine_data() keeps the overlap only."
        } else {
            status <- "dropped"
            note   <- "Only in the TG file - combine_data() keeps the overlap only."
        }

        data.frame(Variable = v, POOL = tp, TG = tt,
                   Status = status, Note = note, stringsAsFactors = FALSE)
    })
    out <- do.call(rbind, rows)
    out[order(factor(out$Status, levels = c("type mismatch", "dropped", "combined")),
              out$Variable), ]
}


# ------ FORMULA BUILDING ----------------------------------------------

#' Build the eVCGsampler formula from the user's selection
#'
#' @param group      Name of the binary POOL / TG variable (left-hand side).
#' @param covariates Character vector of covariates to balance.
#' @param stratum    Optional character vector of stratum variables.
#' @param balanced   If TRUE, use "<group>_balanced" as the left-hand side,
#'   which is the variable VCG_sampler() adds for the post-hoc test.
#' @return A formula, or NULL if the selection is incomplete.
build_formula <- function(group, covariates, stratum = NULL, balanced = FALSE) {
    if (is.null(group) || !nzchar(group) || length(covariates) == 0) return(NULL)
    lhs <- if (balanced) paste0(group, "_balanced") else group
    rhs <- paste(covariates, collapse = " + ")
    if (length(stratum) > 0) rhs <- paste0(rhs, " | ", paste(stratum, collapse = " + "))
    stats::as.formula(paste(lhs, "~", rhs), env = globalenv())
}

#' The formula as a plain string, for display and for the report
formula_text <- function(group, covariates, stratum = NULL, balanced = FALSE) {
    f <- build_formula(group, covariates, stratum, balanced)
    if (is.null(f)) return("")
    paste(deparse(f), collapse = " ")
}


#' Parse a comma-separated list of VCG sizes
#'
#' One value balances the whole data set (or is recycled over the strata by
#' VCG_sampler); several values give one size per stratum.
parse_sizes <- function(txt) {
    if (is.null(txt) || !nzchar(trimws(txt))) return(integer(0))
    v <- suppressWarnings(as.integer(trimws(strsplit(txt, ",", fixed = TRUE)[[1]])))
    v[!is.na(v) & v > 0]
}

#' Parse the covariate importance weights
#'
#' eVCGsampler ignores a c_w of length < 2 and rejects values <= 0, so the
#' same rules are applied here before the value ever reaches the sampler.
parse_weights <- function(txt, n_cov) {
    if (is.null(txt) || !nzchar(trimws(txt))) return(NULL)
    v <- suppressWarnings(as.numeric(trimws(strsplit(txt, ",", fixed = TRUE)[[1]])))
    v <- v[!is.na(v) & is.finite(v) & v > 0]
    if (length(v) < 2 || length(v) != n_cov) return(NULL)
    v
}


# ------ SAMPLE SIZE OVERVIEW ------------------------------------------

#' The analysis data set: selected variables only, complete cases only
#'
#' VCG_sampler() and energy_test() both drop incomplete rows internally.
#' Doing it up front means the sample sizes shown in the app are exactly
#' the ones the sampler will work with.
#'
#' The ID rides along as a passthrough column. VCG_sampler() keeps every
#' column it was given (it merges its results back onto the untouched input),
#' and energy_distance() / energy_test() subset to the formula variables, so
#' the identifier reaches the result without taking part in the balancing.
#'
#' @param id Optional name of the animal-ID column.
analysis_data <- function(data, group, covariates, stratum = NULL, id = NULL) {
    if (is.null(data) || is.null(group) || length(covariates) == 0) return(NULL)
    if (!is.null(id) && !nzchar(id)) id <- NULL
    vars <- unique(c(id, group, stratum, covariates))
    vars <- vars[vars %in% names(data)]
    stats::na.omit(data[, vars, drop = FALSE])
}


#' POOL / TG counts, overall and per stratum, next to the requested VCG size
#'
#' @param data   Output of analysis_data().
#' @param group  Group variable name.
#' @param stratum Optional stratum variables.
#' @param sizes  Requested VCG size(s) from parse_sizes().
size_overview <- function(data, group, stratum = NULL, sizes = integer(0)) {
    if (is.null(data) || nrow(data) == 0) return(NULL)

    g <- as.factor(data[[group]])
    lv <- levels(g)
    # VCG_sampler() balances the FIRST level (the pool) towards the second.
    pool_lv <- lv[1]
    tg_lv   <- if (length(lv) > 1) lv[2] else NA_character_

    make_row <- function(idx, label) {
        gg <- g[idx]
        data.frame(Stratum = label,
                   N_POOL  = sum(gg == pool_lv, na.rm = TRUE),
                   N_TG    = if (is.na(tg_lv)) 0L else sum(gg == tg_lv, na.rm = TRUE),
                   stringsAsFactors = FALSE)
    }

    if (length(stratum) > 0 && all(stratum %in% names(data))) {
        strat <- interaction(data[stratum], drop = TRUE)
        out <- do.call(rbind, lapply(levels(strat), function(k)
            make_row(which(strat == k), k)))
    } else {
        out <- make_row(seq_len(nrow(data)), "(all)")
    }

    # Recycle a single requested size over the strata, exactly as the
    # sampler does; leave it blank when the counts do not line up.
    if (length(sizes) == 1L) {
        out$VCG_size <- sizes
    } else if (length(sizes) == nrow(out)) {
        out$VCG_size <- sizes
    } else {
        out$VCG_size <- NA_integer_
    }

    out$Share_of_POOL <- ifelse(is.na(out$VCG_size) | out$N_POOL == 0, NA_real_,
                                round(100 * out$VCG_size / out$N_POOL, 1))
    out$Status <- vapply(seq_len(nrow(out)), function(i) {
        if (is.na(out$VCG_size[i]))                    return("size missing")
        if (out$N_TG[i] == 0)                          return("no TG animals")
        if (out$VCG_size[i] > out$N_POOL[i])           return("size > POOL")
        if (out$VCG_size[i] >= 0.75 * out$N_POOL[i])   return("size >= 75% of POOL")
        if (out$N_POOL[i] < 10)                        return("very small POOL")
        "ok"
    }, character(1))

    attr(out, "levels") <- c(pool = pool_lv, tg = tg_lv)
    out
}

#' Free-text explanation of the size overview, reused in the report
#'
#' @param ov    Output of size_overview().
#' @param n_cov Number of covariates being balanced, for the pool-size /
#'   covariate-count ratio warning. Optional.
size_messages <- function(ov, n_cov = NA_integer_) {
    if (is.null(ov)) return(character(0))
    msg <- character(0)
    if (any(ov$Status == "size missing"))
        msg <- c(msg, paste0("Give either one VCG size for all strata or exactly ",
                             nrow(ov), " comma-separated sizes (one per stratum)."))
    if (any(ov$Status == "no TG animals"))
        msg <- c(msg, "At least one stratum has no TG animals - there is nothing to balance towards there.")
    if (any(ov$Status == "size > POOL"))
        msg <- c(msg, "A requested VCG size is larger than the available POOL in that stratum.")
    if (any(ov$Status == "size >= 75% of POOL"))
        msg <- c(msg, paste0("A requested VCG size takes 75% or more of the POOL. ",
                             "The sampler will run, but with little room to balance - widen the pool if you can."))
    if (any(ov$Status == "very small POOL"))
        msg <- c(msg, "A stratum has fewer than 10 POOL animals; the balancing result will be unstable.")
    # Same rule of thumb eVCGsampler applies internally: fewer than five pool
    # animals per covariate leaves too little room to balance them all.
    if (!is.na(n_cov) && n_cov > 0 && any(ov$N_POOL / n_cov < 5))
        msg <- c(msg, paste0("Fewer than 5 POOL animals per covariate (", n_cov,
                             " covariates). Balancing that many covariates on this pool ",
                             "will be hard - drop the less important ones."))
    if (!is.na(n_cov) && n_cov > 21)
        msg <- c(msg, paste0("You are balancing ", n_cov, " covariates. Every extra ",
                             "covariate makes the balancing harder, not better - select ",
                             "the ones that matter."))
    msg
}


#' Drop the sampler warnings the app already reports itself
#'
#' VCG_sampler() warns about too many covariates, a small POOL per covariate
#' and a VCG size >= 75% of the POOL. size_messages() runs the same checks
#' (per stratum, in the app's wording) on the Setup step and in the report, so
#' those warnings are removed here instead of being shown twice.
clean_sampler_notes <- function(notes) {
    if (length(notes) == 0) return(notes)
    dup <- grepl("POOL size", notes, fixed = TRUE) |
           grepl("You are trying to balance", notes, fixed = TRUE)
    notes[!dup]
}

#' TRUE when the current selection can be handed to VCG_sampler()
size_is_runnable <- function(ov) {
    !is.null(ov) && nrow(ov) > 0 &&
        !any(ov$Status %in% c("size missing", "no TG animals", "size > POOL"))
}


# ------ RESULT SUMMARIES ----------------------------------------------

#' Add the TG / VCG / POOL plotting group to a VCG_sampler() result
#'
#' VCG_sampler() codes its VCG column as NA for the TG animals, 0 for the
#' unselected pool and 1 for the selected pool. "POOL" below is the whole
#' pool (selected and unselected), which is the comparison the summary
#' table is meant to show.
add_plot_group <- function(dat2) {
    tg   <- dat2[is.na(dat2$VCG), , drop = FALSE]
    vcg  <- dat2[!is.na(dat2$VCG) & dat2$VCG == 1, , drop = FALSE]
    pool <- dat2[!is.na(dat2$VCG), , drop = FALSE]
    if (nrow(tg))   tg$plot_g   <- 1L
    if (nrow(vcg))  vcg$plot_g  <- 2L
    if (nrow(pool)) pool$plot_g <- 3L
    out <- rbind(tg, vcg, pool)
    out$plot_g <- factor(out$plot_g, levels = 1:3, labels = c("TG", "VCG", "POOL"))
    out
}


#' The animals that were drawn into the VCG
#'
#' VCG_sampler() codes its VCG column as 1 for the selected POOL animals, 0
#' for the POOL animals it passed over, and NA for the TG animals. Only the
#' 1s are the result; everything else is the material it worked from. Used by
#' the Sampled data table, both download buttons and the report so all four
#' agree on what "the VCG" means.
#'
#' @param dat2 Output of VCG_sampler()[[1]].
vcg_selected <- function(dat2) {
    if (is.null(dat2)) return(dat2)
    cols <- vcg_cols(dat2)
    if (length(cols)) {
        # multiSampler() result: every animal that made it into at least one
        # of the VCGs, with the number of VCGs it was drawn into.
        m <- as.matrix(dat2[, cols, drop = FALSE])
        times <- rowSums(m == 1, na.rm = TRUE)
        out <- dat2[times > 0, , drop = FALSE]
        out$Times_selected <- times[times > 0]
        return(out)
    }
    if (!"VCG" %in% names(dat2)) return(dat2)
    dat2[!is.na(dat2$VCG) & dat2$VCG == 1, , drop = FALSE]
}


# ------ MULTIPLE VCGs -------------------------------------------------

#' The VCG_1 ... VCG_N indicator columns multiSampler() adds
vcg_cols <- function(dat2) {
    cols <- grep("^VCG_[0-9]+$", names(dat2), value = TRUE)
    cols[order(as.integer(sub("VCG_", "", cols, fixed = TRUE)))]
}

#' One VCG out of a multiSampler() result, in VCG_sampler() layout
#'
#' multiSampler() keeps one indicator column per VCG and drops VCG and
#' <group>_balanced. Rebuilding them for VCG number i lets the post-hoc
#' energy_test(), plot_var() and median_mad_table() run on each VCG exactly
#' as they do on a single VCG_sampler() result. <group>_balanced is coded
#' the way VCG_sampler() codes it: TG -> second level, selected -> first
#' level, unselected POOL -> NA.
#'
#' @param dat2  Output of multiSampler()[[1]].
#' @param i     Which VCG.
#' @param group Name of the POOL / TG variable.
vcg_view <- function(dat2, i, group) {
    cols <- vcg_cols(dat2)
    v    <- dat2[[cols[i]]]
    lvls <- levels(as.factor(dat2[[group]]))
    out  <- dat2[, setdiff(names(dat2), cols), drop = FALSE]
    out$VCG <- v
    bal <- ifelse(is.na(v), 1L, ifelse(v == 1, 0L, NA_integer_))
    out[[paste0(group, "_balanced")]] <- factor(bal, levels = 0:1, labels = lvls)
    out
}

#' Average overlap between the VCGs, as multiSampler() reports it
#'
#' Of the POOL animals selected at least once, the mean share of VCGs each
#' one was drawn into, in percent.
vcg_overlap <- function(dat2) {
    cols <- vcg_cols(dat2)
    if (length(cols) < 2) return(NA_real_)
    m <- as.matrix(dat2[, cols, drop = FALSE])
    m <- m[!is.na(m[, 1]), , drop = FALSE]
    freq <- rowSums(m == 1)
    if (!any(freq > 0)) return(NA_real_)
    round(100 * mean(freq[freq > 0]) / length(cols), 1)
}

#' Median (MAD) per covariate for TG, every VCG and POOL
#'
#' @param dat2  Output of multiSampler()[[1]].
#' @param vars  Covariates to describe.
#' @param group Name of the POOL / TG variable.
median_mad_table_multi <- function(dat2, vars, group) {
    cols <- vcg_cols(dat2)
    if (is.null(dat2) || length(cols) == 0 || length(vars) == 0) return(NULL)
    tabs <- lapply(seq_along(cols), function(i) {
        t <- median_mad_table(vcg_view(dat2, i, group), vars)
        if (!is.null(t)) t$Group[t$Group == "VCG"] <- cols[i]
        t
    })
    if (is.null(tabs[[1]])) return(NULL)
    vcg_rows <- do.call(rbind, lapply(tabs, function(t) t[!t$Group %in% c("TG", "POOL"), , drop = FALSE]))
    first <- tabs[[1]]
    out <- rbind(first[first$Group == "TG", , drop = FALSE], vcg_rows,
                 first[first$Group == "POOL", , drop = FALSE])
    rownames(out) <- NULL
    out
}


#' Median (MAD) per covariate for TG, VCG and POOL
#'
#' @param dat2 Output of VCG_sampler()[[1]].
#' @param vars Covariates to describe.
#' @return A data.frame with one row per group; first column is the group.
median_mad_table <- function(dat2, vars) {
    if (is.null(dat2) || length(vars) == 0) return(NULL)
    d <- add_plot_group(dat2)
    vars <- vars[vars %in% names(d)]
    if (length(vars) == 0) return(NULL)

    cell <- function(x, g) {
        m  <- tapply(x, g, function(z) stats::median(as.numeric(z), na.rm = TRUE))
        s  <- tapply(x, g, function(z) stats::mad(as.numeric(z),    na.rm = TRUE))
        paste0(signif(m, 3), " (", signif(s, 3), ")")
    }

    body <- vapply(vars, function(v) cell(d[[v]], d$plot_g), character(nlevels(d$plot_g)))
    body <- matrix(body, nrow = nlevels(d$plot_g),
                   dimnames = list(levels(d$plot_g), vars))

    out <- data.frame(Group = rownames(body), body, check.names = FALSE,
                      stringsAsFactors = FALSE)
    out$N <- as.integer(table(d$plot_g)[out$Group])
    rownames(out) <- NULL
    out
}


#' Turn an energy_test() htest into a one-row data.frame
#'
#' energy_test() returns list(htest, plot) when plot = TRUE and a bare
#' htest otherwise; both are accepted here.
energy_test_table <- function(res, label = "") {
    if (is.null(res)) return(NULL)
    if (is.list(res) && !inherits(res, "htest") && length(res) >= 1 &&
        inherits(res[[1]], "htest")) res <- res[[1]]
    if (!inherits(res, "htest")) return(NULL)

    data.frame(
        Comparison     = label,
        e_distance     = signif(as.numeric(res$estimate), 4),
        critical_value = signif(as.numeric(res$critical.value[2]), 4),
        p_value        = signif(res$p.value, 3),
        Permutations   = res$n.permutations,
        stringsAsFactors = FALSE
    )
}

#' The p-value of an energy_test() result, whatever shape it came in
energy_test_p <- function(res) {
    if (is.null(res)) return(NA_real_)
    if (is.matrix(res)) return(suppressWarnings(min(as.numeric(res[, "p-value"]), na.rm = TRUE)))
    if (is.list(res) && !inherits(res, "htest") && inherits(res[[1]], "htest")) res <- res[[1]]
    if (inherits(res, "htest")) return(as.numeric(res$p.value))
    NA_real_
}

#' The stratum-wise energy_test() matrix as a tidy data.frame
strat_test_table <- function(res) {
    if (!is.matrix(res)) return(NULL)
    out <- as.data.frame(res, stringsAsFactors = FALSE)
    names(out) <- c("Stratum", "e_distance", "critical_value", "p_value")
    out
}


# ------ SMALL UI HELPERS ----------------------------------------------

#' A coloured status pill (see .vcg-pill in www/sass/app.scss)
pill <- function(text, level = c("ok", "warn", "bad", "neutral")) {
    level <- match.arg(level)
    as.character(shiny::span(class = paste0("vcg-pill pill-", level), text))
}

#' Map a size-overview status to a pill level
status_level <- function(status) {
    switch(status,
           "ok"                   = "ok",
           "size >= 75% of POOL"  = "warn",
           "very small POOL"      = "warn",
           "bad")
}

#' One metric tile (see .vcg-tile in www/sass/app.scss)
tile <- function(label, value, sub = NULL, accent = FALSE) {
    shiny::div(
        class = paste("vcg-tile", if (accent) "tile-accent" else ""),
        shiny::div(class = "vcg-tile-label", label),
        shiny::div(class = "vcg-tile-value", value),
        if (!is.null(sub)) shiny::div(class = "vcg-tile-sub", sub)
    )
}

#' A row of metric tiles
tiles <- function(...) shiny::div(class = "vcg-tiles", ...)

#' Verdict wording for a POOL-vs-TG or TG-vs-VCG energy test
#'
#' The test is one-sided on the energy distance, so a small p-value means
#' "the two groups differ more than permutation noise would explain".
#' Before sampling that is expected; after sampling it is a red flag.
test_verdict <- function(p, phase = c("before", "after"), alpha = 0.05) {
    phase <- match.arg(phase)
    if (is.na(p)) return(list(text = "not available", level = "neutral"))
    if (phase == "before") {
        if (p < alpha) list(text = "POOL and TG differ - balancing is worthwhile", level = "warn")
        else           list(text = "no detectable imbalance to start with",        level = "ok")
    } else {
        if (p < alpha) list(text = "VCG still differs from TG - balancing did not succeed", level = "bad")
        else           list(text = "VCG matches TG - balancing successful",                 level = "ok")
    }
}

#' Standard DT table for this app
#'
#' @param scrollx Horizontal scrolling. Off by default: the app's scrollbars
#'   are styled purple and full width, which is distracting under a table
#'   that only overflows by a few pixels. Switch it on for the genuinely
#'   wide tables (raw data previews).
vcg_table <- function(data, dom = "t", page = 25, rownames = FALSE,
                      scrollx = FALSE, ...) {
    DT::datatable(
        data,
        rownames  = rownames,
        selection = "none",
        options   = list(dom = dom, pageLength = page, scrollX = scrollx,
                         columnDefs = list(list(className = "dt-left", targets = "_all"))),
        ...
    )
}
