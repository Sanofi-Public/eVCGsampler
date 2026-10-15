report_server <- function(
        id,
        import,      # list returned by import_server()
        setup,       # list returned by setup_server()
        pretest,     # list returned by pretest_server()
        sampling     # list returned by sampling_server()
) { moduleServer(id, function(input, output, session) {

    ns <- session$ns

    # The sampling result is what the report is about; without it there is
    # nothing to document.
    sample_res <- reactive({
        res <- sampling$result()
        if (is.null(res) || !is.null(res$error)) return(NULL)
        res
    })

    # ------ CONTENT OVERVIEW ----------------------------------------------
    output$content_table <- renderUI({
        res  <- sample_res()
        ener <- pretest$energy()
        outl <- pretest$outlier()

        row <- function(label, ok, detail) {
            div(
                class = "d-flex justify-content-between align-items-center py-2 border-bottom",
                div(
                    div(class = "fw-bold", label),
                    div(class = "vcg-section-hint mb-0", detail)
                ),
                div(HTML(pill(if (ok) "included" else "not run", if (ok) "ok" else "neutral")))
            )
        }

        tagList(
            row("Data import", !is.null(import$data()),
                paste0(if (isTRUE(import$combined()))
                    "POOL and TG uploaded separately and combined"
                    else "single uploaded file",
                    if (!is.null(import$data())) paste0(" - ", nrow(import$data()), " rows") else "")),
            row("Setup", isTRUE(setup$ready()), setup$formula_text()),
            row("Energy test (POOL vs TG)", !is.null(ener) && is.null(ener$error),
                "permutation test on the covariates before sampling"),
            row("Outlier studies", isTRUE(input$inc_outlier) && !is.null(outl) && is.null(outl$error),
                "energy distance between the source studies"),
            row("VCG sampling", !is.null(res),
                "balance plot, median (MAD) table and post-hoc test"),
            row("Covariate plots", isTRUE(input$inc_covplots) && !is.null(res),
                "TG / VCG / POOL distribution per covariate"),
            row("Sampled animals", isTRUE(input$inc_data) && !is.null(res),
                "the rows selected into the VCG"),
            if (is.null(res)) div(
                class = "vcg-note note-info mt-4",
                "Draw a VCG on the previous step first - the report is built around the sampling result."
            )
        )
    })

    observe({
        if (is.null(sample_res())) {
            shinyjs::disable("download_html")
            shinyjs::disable("download_xlsx")
        } else {
            shinyjs::enable("download_html")
            shinyjs::enable("download_xlsx")
        }
    })

    # ------ PARAMETERS FOR THE RMD ----------------------------------------
    report_params <- reactive({
        res  <- sample_res()
        req(res)
        ener <- pretest$energy()
        outl <- pretest$outlier()
        if (!isTRUE(input$inc_outlier)) outl <- NULL
        if (!is.null(ener) && !is.null(ener$error)) ener <- NULL
        if (!is.null(outl) && !is.null(outl$error)) outl <- NULL

        ov <- setup$overview()

        list(
            meta = list(
                app_version = app_version,
                date        = format(Sys.Date(), "%Y-%m-%d"),
                study       = input$study,
                author      = input$author,
                comment     = input$comment,
                pool_name   = import$pool_name(),
                tg_name     = import$tg_name(),
                n_pool_file = import$n_pool(),
                n_tg_file   = import$n_tg(),
                combined    = isTRUE(import$combined()),
                pkg_version = as.character(utils::packageVersion("eVCGsampler")),
                # Absolute path: the Rmd is rendered from a temp copy, so a
                # relative image reference would not resolve.
                logo        = normalizePath("www/img/logo_vcg.png", mustWork = FALSE)
            ),
            import = list(
                variables = variable_overview(import$data()),
                check     = import$check(),
                n_rows    = nrow(import$data()),
                n_vars    = ncol(import$data())
            ),
            setup = list(
                formula    = setup$formula_text(),
                group      = setup$group(),
                id         = setup$id(),
                covariates = setup$covariates(),
                stratum    = setup$stratum(),
                sizes      = setup$sizes(),
                weights    = setup$weights(),
                overview   = ov,
                messages   = size_messages(ov, length(setup$covariates())),
                n_complete = nrow(setup$adata())
            ),
            energy   = ener,
            outlier  = outl,
            sampling = res,
            opts = list(
                covplots = isTRUE(input$inc_covplots),
                data     = isTRUE(input$inc_data)
            )
        )
    })

    # ------ DOWNLOADS -----------------------------------------------------
    output$download_html <- downloadHandler(
        filename = function() {
            paste0("VCGsampler_report_", gsub("-", "_", Sys.Date(), fixed = TRUE), ".html")
        },
        content = function(file) {
            withProgress(message = "Rendering the report", value = 0.2, {
                # Rendered from a copy in the session temp folder so a
                # multi-user deployment never writes into the app directory.
                src <- normalizePath("www/rmd/report.Rmd")
                tmp <- file.path(tempdir(), "report.Rmd")
                file.copy(src, tmp, overwrite = TRUE)
                file.copy(normalizePath("www/rmd/report.css"),
                          file.path(tempdir(), "report.css"), overwrite = TRUE)
                setProgress(0.4)
                rmarkdown::render(
                    tmp,
                    output_file   = file,
                    params        = report_params(),
                    envir         = new.env(parent = globalenv()),
                    intermediates_dir = tempdir(),
                    quiet         = TRUE
                )
                setProgress(1)
            })
        }
    )

    # The complete data set: the drawn VCG, the POOL animals that were passed
    # over and the TG animals, all with their energy weights. The Sampling
    # step already exports the selection on its own, so duplicating that here
    # would be redundant - this is the "give me everything" export.
    output$download_xlsx <- downloadHandler(
        filename = function() {
            paste0("VCG_all_data_", gsub("-", "_", Sys.Date(), fixed = TRUE), ".xlsx")
        },
        content = function(file) {
            res <- sample_res()
            req(res)
            openxlsx::write.xlsx(
                list(All_data   = res$data,
                     Median_MAD = res$summary,
                     Sizes      = setup$overview()),
                file
            )
        }
    )

    # ------ RETURNS -------------------------------------------------------
    list(
        restart = reactive(input$restart)
    )

})}
