function(input, output, session) {

    # ------ HELPERS -----------------------------------------------------------
    shinyhelper::observe_helpers(withMathJax = TRUE)

    # ------ NAVIGATION --------------------------------------------------------
    # The workflow lives in a hidden tabset (see ui.R). The breadcrumb in the
    # header is the only way through it, and each step unlocks the next one.
    observeEvent(input$tabcontent, {
        shinyjs::removeClass(selector = "#menu a", class = "active")
        shinyjs::addClass(id = glue::glue("link_{input$tabcontent}"), class = "active")

        if (input$tabcontent == "home") {
            shinyjs::hide("header")
        } else {
            shinyjs::show("header")
        }
    })

    # ------ MODULES -----------------------------------------------------------

    # ------ * Home module -----------------------------------------------------
    home_module <- home_server(
        id = "home"
    )

    # ------ * Import module (step 1) ------------------------------------------
    import_module <- import_server(
        id = "import"
    )
    # Setup becomes reachable as soon as the uploaded data could be balanced.
    observeEvent(import_module$ready(), {
        if (isTRUE(import_module$ready())) shinyjs::enable("link_setup")
        else                               shinyjs::disable("link_setup")
    })

    # ------ * Setup module (step 2) -------------------------------------------
    setup_module <- setup_server(
        id   = "setup",
        data = import_module$data
    )
    # A complete, runnable selection unlocks the pre-tests and the sampling.
    # The pre-tests are informative rather than mandatory, so the Sampling step
    # is opened at the same time instead of behind them.
    observeEvent(setup_module$ready(), {
        if (isTRUE(setup_module$ready())) {
            shinyjs::enable("link_pretest")
            shinyjs::enable("link_sampling")
        } else {
            shinyjs::disable("link_pretest")
            shinyjs::disable("link_sampling")
        }
    })

    # ------ * Pre-test module (step 3) ----------------------------------------
    pretest_module <- pretest_server(
        id    = "pretest",
        data  = import_module$data,
        setup = setup_module
    )

    # ------ * Sampling module (step 4) ----------------------------------------
    sampling_module <- sampling_server(
        id    = "sampling",
        setup = setup_module
    )
    # The report is built around the sampling result, so it only opens once a
    # VCG has actually been drawn.
    observeEvent(sampling_module$result(), {
        res <- sampling_module$result()
        if (!is.null(res) && is.null(res$error)) shinyjs::enable("link_report")
        else                                    shinyjs::disable("link_report")
    }, ignoreNULL = FALSE)

    # ------ * Report module (step 5) ------------------------------------------
    report_module <- report_server(
        id       = "report",
        import   = import_module,
        setup    = setup_module,
        pretest  = pretest_module,
        sampling = sampling_module
    )

}
