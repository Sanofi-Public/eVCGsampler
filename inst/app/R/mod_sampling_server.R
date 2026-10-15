sampling_server <- function(
        id,
        setup        # list returned by setup_server()
) { moduleServer(id, function(input, output, session) {

    ns <- session$ns

    # ------ RECAP ---------------------------------------------------------
    output$recap <- renderUI({
        ov <- setup$overview()
        req(ov)
        div(
            class = "mb-3",
            div(class = "vcg-formula mb-3", setup$formula_text()),
            HTML(paste0(
                pill(paste0("VCG n = ",
                            if (all(is.na(ov$VCG_size))) "?" else
                                paste(ov$VCG_size, collapse = " / ")), "neutral"),
                " out of ", paste(ov$N_POOL, collapse = " / "), " POOL animals, ",
                "balanced towards ", paste(ov$N_TG, collapse = " / "), " TG animals."
            ))
        )
    })

    # ------ SAMPLING ------------------------------------------------------
    # One eventReactive holds the whole result: the sampled data, the balance
    # plot and the post-hoc test. That way every output below shows results
    # from the same draw, and nothing recomputes when a tab is switched.
    result <- eventReactive(input$run, {
        ad <- setup$adata()
        req(ad)
        grp   <- setup$group()
        covs  <- setup$covariates()
        strat <- setup$stratum()
        n     <- setup$sizes()
        c_w   <- setup$weights()
        id    <- setup$id()
        req(grp, length(covs) > 0, length(n) > 0)

        f <- build_formula(grp, covs, strat)

        # Several VCGs only exist in random mode - the best subset is unique.
        n_vcg <- if (isTRUE(input$random)) suppressWarnings(as.integer(input$n_vcg)) else 1L
        if (length(n_vcg) != 1 || is.na(n_vcg) || n_vcg < 1) n_vcg <- 1L
        n_vcg <- min(n_vcg, 100L)
        if (n_vcg > 1) {
            return(sample_multi(ad, f, grp, covs, strat, n, c_w, id, n_vcg))
        }

        withProgress(message = "Drawing the VCG", value = 0.1, {
            # The balance plot and the random-subset option both draw random
            # numbers, so the seed is what makes a run reproducible.
            if (is.finite(input$seed)) set.seed(as.integer(input$seed))

            notes <- character(0)
            out <- tryCatch(
                withCallingHandlers(
                    eVCGsampler::VCG_sampler(f, data = ad, n = n, c_w = c_w,
                                             random = isTRUE(input$random), plot = TRUE),
                    warning = function(w) {
                        notes <<- c(notes, conditionMessage(w))
                        invokeRestart("muffleWarning")
                    }
                ),
                error = function(e) e
            )
            if (inherits(out, "error")) return(list(error = conditionMessage(out)))
            setProgress(0.55, message = "Testing the balance")

            # The sampler appends its VCG / e_weights columns at the end, so
            # the ID lands wherever it happened to sit in the input. Put it
            # back in front - it is the first thing you look for in a result.
            dat2  <- id_first(out[[1]], id)
            bplot <- out[[2]]

            # Post-hoc test on the variable VCG_sampler() created: does the
            # selected VCG still differ from the TG?
            f_bal <- build_formula(grp, covs, strat, balanced = TRUE)
            post <- tryCatch(
                eVCGsampler::energy_test(f_bal, data = dat2, R = input$R,
                                         plot = length(strat) == 0),
                error = function(e) e
            )
            setProgress(0.9)

            post_err <- if (inherits(post, "error")) conditionMessage(post) else NULL
            if (!is.null(post_err)) post <- NULL

            list(
                data        = dat2,
                # The drawn VCG on its own. Computed once here so the Sampled
                # data table, both download buttons, the report's Excel export
                # and the report template cannot drift apart on what counts as
                # "selected" - the Rmd renders in its own environment and
                # cannot reach vcg_selected().
                selected    = vcg_selected(dat2),
                cov_data    = dat2,
                multi       = FALSE,
                n_vcg       = 1L,
                plot        = bplot,
                summary     = median_mad_table(dat2, covs),
                post        = post,
                post_p      = if (is.null(post)) NA_real_ else energy_test_p(post),
                post_table  = if (is.null(post)) NULL else
                    if (length(strat) > 0) strat_test_table(post)
                    else energy_test_table(post, "TG vs VCG (after sampling)"),
                post_error  = post_err,
                stratified  = length(strat) > 0,
                covariates  = covs,
                stratum     = strat,
                group       = grp,
                id          = id,
                sizes       = n,
                weights     = c_w,
                random      = isTRUE(input$random),
                seed        = input$seed,
                R           = input$R,
                notes       = clean_sampler_notes(notes),
                formula     = formula_text(grp, covs, strat)
            )
        })
    })

    # Several random VCGs in one go via multiSampler(). Same result list as
    # the single draw, with multi = TRUE: data carries VCG_1 ... VCG_N instead
    # of VCG / e_weights, the balance plot is replaced by the overlap plot(s),
    # and the post-hoc test runs once per VCG.
    sample_multi <- function(ad, f, grp, covs, strat, n, c_w, id, n_vcg) {
        withProgress(message = paste("Drawing", n_vcg, "VCGs"), value = 0.05, {
            if (is.finite(input$seed)) set.seed(as.integer(input$seed))

            notes <- character(0)
            out <- tryCatch(
                withCallingHandlers(
                    eVCGsampler::multiSampler(f, data = ad, n = n, c_w = c_w,
                                              Nsamples = n_vcg, plot = TRUE),
                    warning = function(w) {
                        notes <<- c(notes, conditionMessage(w))
                        invokeRestart("muffleWarning")
                    }
                ),
                error = function(e) e
            )
            if (inherits(out, "error")) return(list(error = conditionMessage(out)))

            dat2 <- id_first(out[[1]], id)
            # Stratified: one overlap plot per stratum, named by stratum.
            overlap <- if (inherits(out[[2]], "ggplot")) list(out[[2]]) else out[[2]]
            if (length(strat) > 0 && "in_stratum" %in% names(dat2))
                names(overlap) <- levels(as.factor(dat2$in_stratum))[seq_along(overlap)]

            # Post-hoc test, once per VCG.
            f_bal <- build_formula(grp, covs, strat, balanced = TRUE)
            post_ps <- rep(NA_real_, n_vcg)
            post_rows <- vector("list", n_vcg)
            post_err <- NULL
            for (i in seq_len(n_vcg)) {
                setProgress(0.1 + 0.8 * (i - 1) / n_vcg,
                            message = paste0("Testing the balance of VCG ", i, " / ", n_vcg))
                post <- tryCatch(
                    eVCGsampler::energy_test(f_bal, data = vcg_view(dat2, i, grp),
                                             R = input$R, plot = FALSE),
                    error = function(e) e
                )
                if (inherits(post, "error")) {
                    post_err <- conditionMessage(post)
                    next
                }
                post_ps[i] <- energy_test_p(post)
                tab <- if (length(strat) > 0) strat_test_table(post)
                       else energy_test_table(post, "TG vs VCG")
                if (!is.null(tab)) post_rows[[i]] <- cbind(VCG = paste0("VCG_", i), tab)
            }
            post_table <- do.call(rbind, post_rows)
            setProgress(0.95)

            list(
                data        = dat2,
                selected    = vcg_selected(dat2),
                # Covariate plots in the report show the first VCG.
                cov_data    = vcg_view(dat2, 1, grp),
                multi       = TRUE,
                n_vcg       = n_vcg,
                overlap     = overlap,
                overlap_pct = vcg_overlap(dat2),
                plot        = NULL,
                summary     = median_mad_table_multi(dat2, covs, grp),
                post        = NULL,
                post_p      = if (any(!is.na(post_ps))) min(post_ps, na.rm = TRUE) else NA_real_,
                post_ps     = post_ps,
                post_table  = post_table,
                post_error  = post_err,
                stratified  = length(strat) > 0,
                covariates  = covs,
                stratum     = strat,
                group       = grp,
                id          = id,
                sizes       = n,
                weights     = c_w,
                random      = TRUE,
                seed        = input$seed,
                R           = input$R,
                notes       = clean_sampler_notes(notes),
                formula     = formula_text(grp, covs, strat)
            )
        })
    }

    # ------ BALANCE TAB ---------------------------------------------------
    output$balance_plot <- renderPlot({
        res <- result()
        req(res, is.null(res$error))
        res$plot
    })

    output$dt_summary <- DT::renderDataTable({
        res <- result()
        req(res, is.null(res$error), res$summary)
        vcg_table(res$summary, dom = if (isTRUE(res$multi)) "tp" else "t", page = 25)
    })

    output$post_plot <- renderPlot({
        res <- result()
        req(res, is.null(res$error), !res$stratified, !is.null(res$post))
        res$post[[2]]
    })

    output$dt_post <- DT::renderDataTable({
        res <- result()
        req(res, is.null(res$error), !is.null(res$post_table))
        vcg_table(res$post_table, dom = if (isTRUE(res$multi)) "tp" else "t", page = 25)
    })

    output$balance_body <- renderUI({
        if (!clicked(input$run)) {
            return(div(class = "vcg-note note-info",
                       "Press \"Draw the VCG\" to sample a balanced control group."))
        }
        res <- result()
        req(res)
        if (!is.null(res$error)) {
            return(div(class = "vcg-note note-bad",
                       paste0("Sampling failed: ", res$error)))
        }
        if (isTRUE(res$multi)) return(balance_body_multi(res))
        v <- test_verdict(res$post_p, "after")
        n_vcg <- sum(!is.na(res$data$VCG) & res$data$VCG == 1)

        tagList(
            tiles(
                tile("VCG selected", n_vcg, "animals", accent = TRUE),
                tile("Post-hoc p-value",
                     if (is.na(res$post_p)) "-" else signif(res$post_p, 3),
                     if (res$stratified) "smallest across strata" else "TG vs VCG"),
                tile("Mode", if (res$random) "random subset" else "best subset",
                     paste0("seed ", res$seed))
            ),
            div(class = "mt-3", HTML(paste0(
                pill(v$text, v$level),
                " A large p-value here is the goal: it means the selected VCG is no ",
                "longer distinguishable from the TG on the balanced covariates."
            ))),
            if (length(res$notes)) div(class = "vcg-note mt-3",
                                       paste(unique(res$notes), collapse = " ")),
            h5("Balance target plot", class = "mt-4 vcg-section-title"),
            p(class = "vcg-section-hint",
              paste("Each arrow starts at the POOL-vs-TG difference and points at the",
                    "VCG-vs-TG difference, in median (x) and MAD (y). The purple ellipse",
                    "is the 95% range of differences seen under random permutation, so an",
                    "arrow tip inside it means the remaining difference is ordinary noise.")),
            plotOutput(ns("balance_plot"), height = "520px"),
            h5("Median (MAD) per group", class = "mt-4 vcg-section-title"),
            p(class = "vcg-section-hint",
              "POOL is the whole pool, VCG the selected subset. The VCG row should sit closer to TG than POOL does."),
            DT::dataTableOutput(ns("dt_summary")),
            h5("Post-hoc energy test - TG vs VCG", class = "mt-4 vcg-section-title"),
            if (!is.null(res$post_error))
                div(class = "vcg-note note-bad", paste0("The post-hoc test failed: ", res$post_error))
            else tagList(
                DT::dataTableOutput(ns("dt_post")),
                if (!res$stratified)
                    div(class = "mt-3", plotOutput(ns("post_plot"), height = "420px"))
            )
        )
    })

    # Several VCGs: overlap plot(s) instead of the balance target plot, and
    # the post-hoc test as one row per VCG.
    observeEvent(result(), {
        res <- result()
        if (is.null(res) || !is.null(res$error) || !isTRUE(res$multi)) return()
        lapply(seq_along(res$overlap), function(k) {
            output[[paste0("overlap_plot_", k)]] <- renderPlot(res$overlap[[k]])
        })
    })

    balance_body_multi <- function(res) {
        n_ok   <- sum(!is.na(res$post_ps) & res$post_ps >= 0.05)
        first  <- res$data[[vcg_cols(res$data)[1]]]
        n_each <- sum(first == 1, na.rm = TRUE)
        # Overlap plots draw one bar per POOL animal, so they grow with the pool.
        pool   <- res$data[!is.na(first), , drop = FALSE]
        n_pool <- if ("in_stratum" %in% names(pool)) max(table(pool$in_stratum)) else nrow(pool)
        h <- max(420, min(1600, 8 * n_pool))
        tagList(
            tiles(
                tile("VCGs drawn", res$n_vcg, paste(n_each, "animals each"), accent = TRUE),
                tile("Balanced VCGs", paste0(n_ok, " / ", res$n_vcg),
                     "post-hoc p >= 0.05"),
                tile("Smallest p-value",
                     if (is.na(res$post_p)) "-" else signif(res$post_p, 3),
                     "across all VCGs"),
                tile("Overlap", if (is.na(res$overlap_pct)) "-" else paste0(res$overlap_pct, "%"),
                     paste0("seed ", res$seed))
            ),
            div(class = "vcg-note mt-3",
                paste("Several random VCGs drawn with multiSampler(). Use them as a set,",
                      "e.g. for a proof-of-concept study - they are not meant for picking",
                      "the best-looking VCG afterwards. For a single VCG, set the number of",
                      "VCGs to 1.")),
            if (length(res$notes)) div(class = "vcg-note mt-3",
                                       paste(unique(res$notes), collapse = " ")),
            h5("VCG overlap", class = "mt-4 vcg-section-title"),
            p(class = "vcg-section-hint",
              paste("How often each POOL animal was drawn across the VCGs. The overlap is",
                    "the mean share of VCGs an animal is in, among the animals drawn at",
                    "least once - a high value means the VCGs are largely the same animals.")),
            lapply(seq_along(res$overlap), function(k) tagList(
                if (!is.null(names(res$overlap)))
                    h6(paste("Stratum", names(res$overlap)[k]), class = "mt-3"),
                plotOutput(ns(paste0("overlap_plot_", k)), height = paste0(h, "px"))
            )),
            h5("Median (MAD) per group", class = "mt-4 vcg-section-title"),
            p(class = "vcg-section-hint",
              "One row per VCG. Every VCG row should sit closer to TG than POOL does."),
            DT::dataTableOutput(ns("dt_summary")),
            h5("Post-hoc energy test - TG vs each VCG", class = "mt-4 vcg-section-title"),
            if (!is.null(res$post_error))
                div(class = "vcg-note note-bad", paste0("The post-hoc test failed: ", res$post_error)),
            if (!is.null(res$post_table)) DT::dataTableOutput(ns("dt_post"))
        )
    }

    # ------ COVARIATES TAB ------------------------------------------------
    # One plot_var() panel per balanced covariate. plot_var() stacks the
    # strata itself when the data carries in_stratum, so the plot gets taller
    # the more strata there are.
    plot_height <- reactive({
        res <- result()
        if (is.null(res) || !is.null(res$error)) return(320)
        nl <- if ("in_stratum" %in% names(res$data)) nlevels(as.factor(res$data$in_stratum)) else 1
        rows <- if (nl <= 3) nl else ceiling(nl / 2)
        320 * max(1, rows)
    })

    observeEvent(result(), {
        res <- result()
        if (is.null(res) || !is.null(res$error)) return()
        lapply(seq_along(res$covariates), function(i) {
            output[[paste0("cov_plot_", i)]] <- renderPlot({
                eVCGsampler::plot_var(cov_data(), what = res$covariates[i], group = "VCG")
            })
        })
    })

    # The data behind the covariate plots: the single VCG, or the VCG picked
    # in the selector when several were drawn.
    cov_data <- reactive({
        res <- result()
        req(res, is.null(res$error))
        if (!isTRUE(res$multi)) return(res$data)
        i <- suppressWarnings(as.integer(sub("VCG_", "", input$cov_vcg %||% "VCG_1", fixed = TRUE)))
        if (is.na(i) || i < 1 || i > res$n_vcg) i <- 1L
        vcg_view(res$data, i, res$group)
    })

    output$covariate_body <- renderUI({
        if (!clicked(input$run)) {
            return(div(class = "vcg-note note-info",
                       "The per-covariate distributions appear here once the VCG has been drawn."))
        }
        res <- result()
        req(res)
        if (!is.null(res$error)) return(NULL)
        h <- plot_height()
        tagList(
            p(class = "vcg-section-hint",
              "TG, VCG and POOL side by side for every balanced covariate. The VCG should track the TG, not the POOL."),
            if (isTRUE(res$multi))
                selectInput(ns("cov_vcg"), "Show VCG",
                            choices = vcg_cols(res$data), selected = "VCG_1"),
            lapply(seq_along(res$covariates), function(i) {
                div(
                    class = "mb-4",
                    h5(res$covariates[i], class = "vcg-section-title"),
                    plotOutput(ns(paste0("cov_plot_", i)), height = paste0(h, "px"))
                )
            })
        )
    })

    # ------ DATA TAB ------------------------------------------------------
    # Only the animals that made it into the VCG: VCG == 1. The unselected
    # POOL animals (VCG == 0) and the TG animals (VCG is NA) are what the
    # sampler worked with, not what it produced, so they would only get in
    # the way here. Both downloads still contain the whole data set.
    vcg_rows <- reactive({
        res <- result()
        req(res, is.null(res$error))
        res$selected
    })

    output$dt_result <- DT::renderDataTable({
        vcg_table(vcg_rows(), dom = "tp", page = 15, scrollx = TRUE)
    })

    output$data_body <- renderUI({
        if (!clicked(input$run)) {
            return(div(class = "vcg-note note-info",
                       "The sampled data set appears here once the VCG has been drawn."))
        }
        res <- result()
        req(res)
        if (!is.null(res$error)) return(NULL)
        if (isTRUE(res$multi)) return(tagList(
            p(class = "vcg-section-hint",
              HTML(paste0(
                  "The <b>", nrow(vcg_rows()), " animals drawn into at least one of the ",
                  res$n_vcg, " VCGs</b>. <code>VCG_1</code> ... <code>VCG_", res$n_vcg,
                  "</code> mark membership in each VCG (<code>1</code> = selected), and",
                  " <code>Times_selected</code> counts the VCGs an animal is in.",
                  if (!is.null(res$id))
                      paste0(" <code>", res$id, "</code> in the first column identifies",
                             " each one in your source data.")
                  else
                      paste0(" No ID variable was selected on the Setup step, so these rows",
                             " cannot be traced back to individual animals."),
                  " The downloads on the left contain exactly these rows."
              ))),
            DT::dataTableOutput(ns("dt_result"))
        ))
        tagList(
            p(class = "vcg-section-hint",
              HTML(paste0(
                  "The <b>", nrow(vcg_rows()), " animals selected into the VCG</b> - the",
                  " rows where <code>VCG</code> is <code>1</code>.",
                  if (!is.null(res$id))
                      paste0(" <code>", res$id, "</code> in the first column identifies",
                             " each one in your source data.")
                  else
                      paste0(" No ID variable was selected on the Setup step, so these rows",
                             " cannot be traced back to individual animals."),
                  " <code>e_weights</code> holds the energy weight each one received, and",
                  " <code>", res$group, "_balanced</code> is the group variable used for the",
                  " post-hoc test. The downloads on the left contain exactly these rows."
              ))),
            DT::dataTableOutput(ns("dt_result"))
        )
    })

    # ------ DOWNLOADS -----------------------------------------------------
    observeEvent(result(), {
        res <- result()
        ok <- !is.null(res) && is.null(res$error)
        if (ok) {
            shinyjs::enable("download_csv")
            shinyjs::enable("download_xlsx")
            shinyjs::enable("continue")
        } else {
            shinyjs::disable("download_csv")
            shinyjs::disable("download_xlsx")
            shinyjs::disable("continue")
        }
    })

    # Both exports carry the drawn VCG only - the same rows the Sampled data
    # table shows. The unselected POOL animals and the TG animals are inputs
    # to the sampling, not its result, and the user already has them in their
    # own source files.
    output$download_csv <- downloadHandler(
        filename = function() paste0("VCG_selected_", gsub("-", "_", Sys.Date(), fixed = TRUE), ".csv"),
        content = function(file) {
            utils::write.table(vcg_rows(), file, row.names = FALSE, sep = ";", dec = ".", na = "")
        }
    )

    output$download_xlsx <- downloadHandler(
        filename = function() paste0("VCG_selected_", gsub("-", "_", Sys.Date(), fixed = TRUE), ".xlsx"),
        content = function(file) {
            res <- result()
            sheets <- list(VCG_selected = vcg_rows(), Median_MAD = res$summary)
            if (isTRUE(res$multi) && !is.null(res$post_table)) sheets$Posthoc_test <- res$post_table
            openxlsx::write.xlsx(sheets, file)
        }
    )

    # ------ RETURNS -------------------------------------------------------
    list(
        continue = reactive(input$continue),
        result   = reactive(if (!clicked(input$run)) NULL else result())
    )

})}
