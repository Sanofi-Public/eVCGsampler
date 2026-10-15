pretest_server <- function(
        id,
        data,        # reactive: the full imported data set
        setup        # list returned by setup_server()
) { moduleServer(id, function(input, output, session) {

    ns <- session$ns

    # ------ SHOW THE RESULTS OF WHAT WAS JUST RUN -------------------------
    # The two tests live on separate tabs but share one settings sidebar, so
    # pressing a Run button while the other tab is open would compute results
    # the user never sees. Follow the button.
    observeEvent(input$run_energy, {
        updateTabsetPanel(session, "result_tabs", selected = "energy")
    })
    observeEvent(input$run_outlier, {
        updateTabsetPanel(session, "result_tabs", selected = "outlier")
    })

    # ------ ENERGY TEST (POOL vs TG) --------------------------------------
    # Run on demand: a permutation test with a few thousand replicates is far
    # too slow to re-run on every keystroke in the Setup step.
    energy_res <- eventReactive(input$run_energy, {
        ad <- setup$adata()
        req(ad)
        grp   <- setup$group()
        covs  <- setup$covariates()
        strat <- setup$stratum()
        alpha <- input$alpha
        R     <- input$R
        req(grp, length(covs) > 0, is.finite(alpha), is.finite(R))

        f_flat <- build_formula(grp, covs)                 # never stratified
        f_full <- build_formula(grp, covs, strat)

        withProgress(message = "Running the permutation test", value = 0.1, {
            # energy_distance() refuses "|", so the overall distance is always
            # computed on the flat formula, with the per-stratum distances
            # added separately below.
            dist_all <- tryCatch(eVCGsampler::energy_distance(f_flat, data = ad),
                                 error = function(e) NA_real_)
            setProgress(0.25)

            if (length(strat) > 0) {
                # With strata energy_test() returns a character matrix and
                # cannot draw the permutation histogram.
                tab <- tryCatch(
                    eVCGsampler::energy_test(f_full, data = ad, alpha = alpha,
                                             R = R, plot = FALSE),
                    error = function(e) e
                )
                setProgress(0.9)
                if (inherits(tab, "error")) {
                    return(list(error = conditionMessage(tab)))
                }
                strat_int <- interaction(ad[strat], drop = TRUE)
                per_strat <- vapply(levels(strat_int), function(k) {
                    tryCatch(eVCGsampler::energy_distance(f_flat,
                                                          data = ad[strat_int == k, , drop = FALSE]),
                             error = function(e) NA_real_)
                }, numeric(1))
                return(list(
                    stratified   = TRUE,
                    distance     = dist_all,
                    table        = strat_test_table(tab),
                    per_stratum  = data.frame(Stratum = names(per_strat),
                                              e_distance = signif(per_strat, 4),
                                              stringsAsFactors = FALSE),
                    p_value      = energy_test_p(tab),
                    alpha        = alpha,
                    R            = R,
                    formula      = formula_text(grp, covs, strat)
                ))
            }

            res <- tryCatch(
                eVCGsampler::energy_test(f_flat, data = ad, alpha = alpha,
                                         R = R, plot = TRUE),
                error = function(e) e
            )
            setProgress(0.9)
            if (inherits(res, "error")) return(list(error = conditionMessage(res)))
            list(
                stratified = FALSE,
                distance   = dist_all,
                htest      = res[[1]],
                plot       = res[[2]],
                table      = energy_test_table(res, "POOL vs TG (before sampling)"),
                p_value    = energy_test_p(res),
                alpha      = alpha,
                R          = R,
                formula    = formula_text(grp, covs)
            )
        })
    })

    output$energy_formula <- renderUI({
        div(class = "vcg-formula mb-3", setup$formula_text())
    })

    output$energy_tiles <- renderUI({
        res <- energy_res()
        req(res)
        if (!is.null(res$error)) return(NULL)
        tiles(
            tile("Energy distance", signif(res$distance, 4),
                 "0 = identical distributions", accent = TRUE),
            tile("p-value",
                 if (is.na(res$p_value)) "-" else signif(res$p_value, 3),
                 if (isTRUE(res$stratified)) "smallest across strata" else "one-sided permutation test"),
            tile("Permutations", res$R, "per test"),
            tile("Alpha", res$alpha, "significance level")
        )
    })

    output$energy_verdict <- renderUI({
        res <- energy_res()
        req(res)
        if (!is.null(res$error)) {
            return(div(class = "vcg-note note-bad mt-3",
                       paste0("The energy test failed: ", res$error)))
        }
        v <- test_verdict(res$p_value, "before", res$alpha)
        div(
            class = "mt-3",
            HTML(paste0(
                pill(v$text, v$level),
                " A small p-value means the POOL and the TG differ more than random ",
                "variation would explain. That is the situation VCG sampling is for; ",
                "if there is no imbalance to begin with, a plain random draw is enough."
            ))
        )
    })

    output$energy_plot <- renderPlot({
        res <- energy_res()
        req(res, isFALSE(res$stratified))
        res$plot
    })

    output$dt_energy <- DT::renderDataTable({
        res <- energy_res()
        req(res, isTRUE(res$stratified))
        vcg_table(res$table, dom = "t")
    })

    output$dt_energy_strat <- DT::renderDataTable({
        res <- energy_res()
        req(res, isTRUE(res$stratified))
        vcg_table(res$per_stratum, dom = "t")
    })

    output$energy_body <- renderUI({
        if (!clicked(input$run_energy)) {
            return(div(class = "vcg-note note-info",
                       "Press \"Run energy test\" to compare the POOL with the TG on the selected covariates."))
        }
        res <- energy_res()
        req(res)
        if (!is.null(res$error)) return(NULL)

        if (isTRUE(res$stratified)) {
            tagList(
                h5("Test per stratum", class = "mt-4 vcg-section-title"),
                p(class = "vcg-section-hint",
                  "Balancing runs inside each stratum, so the test does too. The",
                  "permutation histogram is only available without stratification."),
                DT::dataTableOutput(ns("dt_energy")),
                h5("Energy distance per stratum", class = "mt-4 vcg-section-title"),
                DT::dataTableOutput(ns("dt_energy_strat"))
            )
        } else {
            tagList(
                h5("Permutation distribution", class = "mt-4 vcg-section-title"),
                p(class = "vcg-section-hint",
                  "Purple bars are the energy distances obtained by randomly reshuffling",
                  "the group labels. Orange marks the critical value, cyan the observed",
                  "distance."),
                plotOutput(ns("energy_plot"), height = "420px")
            )
        }
    })

    # ------ OUTLIER STUDIES ------------------------------------------------
    # Candidate grouping variables: anything that looks like a study or batch
    # label. find_outliers() compares every group against every other, so a
    # column with hundreds of levels is not useful here.
    output$study_ui <- renderUI({
        d <- tryCatch(data(), error = function(e) NULL)
        req(d)
        u  <- n_unique(d)
        ch <- names(d)[u >= 2 & u <= 50]
        ch <- setdiff(ch, c(group_indicator, setup$group()))
        selectInput(
            inputId = ns("study"),
            label   = "Study / batch variable",
            choices = ch
        )
    })

    output$outvars_ui <- renderUI({
        d <- tryCatch(data(), error = function(e) NULL)
        req(d)
        ch <- names(d)[vapply(d, is.numeric, logical(1))]
        ch <- setdiff(ch, c(group_indicator, setup$group(), input$study))
        selectizeInput(
            inputId  = ns("outvars"),
            label    = "Outcome or technical variables",
            choices  = ch,
            multiple = TRUE,
            options  = multi_select_options
        )
    })

    output$poolonly_ui <- renderUI({
        req(setup$group())
        checkboxInput(
            inputId = ns("pool_only"),
            label   = "POOL animals only",
            value   = TRUE
        )
    })

    # Balancing covariates do not belong in this test - warn rather than
    # silently filter, because the user may have a reason.
    outlier_warn <- reactive({
        ov <- intersect(input$outvars %||% character(0), setup$covariates() %||% character(0))
        if (length(ov) == 0) return(NULL)
        paste0("You selected ", paste(ov, collapse = ", "),
               " - these are balancing covariates and are handled by the sampling. ",
               "This test is meant for outcome or technical variables.")
    })

    outlier_res <- eventReactive(input$run_outlier, {
        d <- tryCatch(data(), error = function(e) NULL)
        req(d, input$study, length(input$outvars) > 0)

        dd <- d
        if (isTRUE(input$pool_only) && !is.null(setup$group()) && setup$group() %in% names(dd)) {
            lv <- levels(as.factor(dd[[setup$group()]]))
            dd <- dd[as.character(dd[[setup$group()]]) == lv[1], , drop = FALSE]
        }
        f <- stats::as.formula(paste(input$study, "~", paste(input$outvars, collapse = " + ")),
                               env = globalenv())

        withProgress(message = "Comparing the studies", value = 0.2, {
            notes <- character(0)
            res <- tryCatch(
                withCallingHandlers(
                    eVCGsampler::find_outliers(f, data = dd, cutoff = input$cutoff,
                                               R = input$R_out, plot = TRUE),
                    warning = function(w) {
                        notes <<- c(notes, conditionMessage(w))
                        invokeRestart("muffleWarning")
                    }
                ),
                error = function(e) e
            )
            setProgress(0.95)
            if (inherits(res, "error")) return(list(error = conditionMessage(res)))
            res$notes    <- notes
            res$n_used   <- nrow(stats::na.omit(dd[, c(input$study, input$outvars), drop = FALSE]))
            res$study    <- input$study
            res$vars     <- input$outvars
            res$cutoff   <- input$cutoff
            res$pool_only <- isTRUE(input$pool_only)
            res
        })
    })

    output$dt_outlier <- DT::renderDataTable({
        res <- outlier_res()
        req(res, is.null(res$error))
        s <- res$summary
        show <- data.frame(
            Study           = s$group,
            Median_distance = signif(s$median_distance, 4),
            Verdict         = vapply(s$is_outlier, function(o)
                pill(if (isTRUE(o)) "outlier" else "consistent",
                     if (isTRUE(o)) "bad" else "ok"), character(1)),
            stringsAsFactors = FALSE
        )
        names(show) <- c("Study", "Median energy distance", "Verdict")
        vcg_table(show, dom = "t", page = 25, escape = FALSE)
    })

    output$outlier_heatmap <- renderPlot({
        res <- outlier_res()
        req(res, is.null(res$error))
        res$heatmap
    })

    output$outlier_barplot <- renderPlot({
        res <- outlier_res()
        req(res, is.null(res$error))
        res$barplot
    })

    output$outlier_body <- renderUI({
        if (!clicked(input$run_outlier)) {
            return(tagList(
                div(class = "vcg-note note-info",
                    paste("Pick a study variable and one or more outcome variables, then press",
                          "\"Run outlier test\". Studies whose outcome distribution sits far from",
                          "the rest are flagged - consider excluding them from the pool before sampling.")),
                if (!is.null(outlier_warn())) div(class = "vcg-note mt-3", outlier_warn())
            ))
        }
        res <- outlier_res()
        req(res)
        if (!is.null(res$error)) {
            return(div(class = "vcg-note note-bad",
                       paste0("The outlier test failed: ", res$error)))
        }
        n_out <- sum(res$summary$is_outlier)
        tagList(
            div(
                HTML(paste0(
                    pill(paste0(n_out, if (n_out == 1) " outlier study" else " outlier studies"),
                         if (n_out == 0) "ok" else "warn"),
                    " ", nrow(res$summary), " studies compared on ",
                    paste(res$vars, collapse = ", "),
                    " (", format(res$n_used, big.mark = " "), " animals",
                    if (isTRUE(res$pool_only)) ", POOL only" else "", ")."
                ))
            ),
            if (!is.null(outlier_warn())) div(class = "vcg-note mt-3", outlier_warn()),
            if (length(res$notes)) div(class = "vcg-note mt-3",
                                       paste(unique(res$notes), collapse = " ")),
            div(class = "mt-3", DT::dataTableOutput(ns("dt_outlier"))),
            h5("Median distance to the other studies", class = "mt-4 vcg-section-title"),
            p(class = "vcg-section-hint",
              paste0("The dashed line is the ", round(100 * res$cutoff),
                     "th percentile of the permuted distances - bars beyond it are flagged.")),
            plotOutput(ns("outlier_barplot"), height = "360px"),
            h5("Pairwise energy distances", class = "mt-4 vcg-section-title"),
            plotOutput(ns("outlier_heatmap"), height = "460px")
        )
    })

    # ------ RETURNS -------------------------------------------------------
    list(
        continue   = reactive(input$continue),
        energy     = reactive(if (!clicked(input$run_energy))  NULL else energy_res()),
        outlier    = reactive(if (!clicked(input$run_outlier)) NULL else outlier_res())
    )

})}
