fluidPage(
    # ------ BOOSTRAP THEME ----------------------------------------------------
    theme = bslib::bs_theme(
        version = 5,
        primary = "#7A00E6",
        secondary = '#23004C',
        light = "#f4f2f6",
        dark = "#BBBAC6",
        "btn-border-radius" = "200px",
        "disabled-background" = "#333",
        "h1-font-size" = "3.15rem",
        "btn-padding-y" = "14px",
        "btn-padding-x" = "20px",
        "btn-line-height" = "1",
        'h6-font-size' = "12px",
    ),
    # Shiny + SCSS
    shinyjs::useShinyjs(),
    includeCSS("www/css/app.min.css"),
    includeScript("www/js/bindEnter.js"),
    includeScript("www/js/hack.js"),
    shinyjs::hidden(shinyhelper::helper("")),
    # ------ HEADER -----------------------------------------------------------
    # Hidden on the welcome screen, which shows the large logo instead.
    shinyjs::hidden(
        tags$header(
            id = "header",
            fluidRow(
                column(
                    width = 12,
                    fixedPage(
                        fluidRow(
                            class = "align-items-center",
                            column(
                                width = 2,
                                id = "menu-logo",
                                div(
                                    img(
                                        src = "img/logo_vcg.png",
                                        class = "logo logo-badge",
                                        alt = "VCG sampler",
                                        title = "Back to the welcome screen",
                                        onclick = "$('#tabcontent').data('shiny-input-binding').setValue($('#tabcontent'), 'home');"
                                    )
                                )
                            ),
                            # ------ * Workflow breadcrumb ------------------------
                            column(
                                id = "menu",
                                width = 10,
                                div(
                                    actionLink(
                                        inputId = "link_import",
                                        label = "Data",
                                        class = "fs-5",
                                        title = "Step 1 - upload the POOL and TG data",
                                        onclick = "$('#tabcontent').data('shiny-input-binding').setValue($('#tabcontent'), 'import');"
                                    ),
                                    span(" > "),
                                    shinyjs::disabled(
                                        actionLink(
                                            inputId = "link_setup",
                                            label = "Setup",
                                            class = "fs-5",
                                            title = "Step 2 - group, covariates, strata and VCG size",
                                            onclick = "$('#tabcontent').data('shiny-input-binding').setValue($('#tabcontent'), 'setup');"
                                        )
                                    ),
                                    span(" > "),
                                    shinyjs::disabled(
                                        actionLink(
                                            inputId = "link_pretest",
                                            label = "Pre-tests",
                                            class = "fs-5",
                                            title = "Step 3 - energy distance, permutation test and outlier studies",
                                            onclick = "$('#tabcontent').data('shiny-input-binding').setValue($('#tabcontent'), 'pretest');"
                                        )
                                    ),
                                    span(" > "),
                                    shinyjs::disabled(
                                        actionLink(
                                            inputId = "link_sampling",
                                            label = "Sampling",
                                            class = "fs-5",
                                            title = "Step 4 - draw the VCG and check the balance",
                                            onclick = "$('#tabcontent').data('shiny-input-binding').setValue($('#tabcontent'), 'sampling');"
                                        )
                                    ),
                                    span(" > "),
                                    shinyjs::disabled(
                                        actionLink(
                                            inputId = "link_report",
                                            label = "Report",
                                            class = "fs-5",
                                            title = "Step 5 - download the HTML report and the sampled data",
                                            onclick = "$('#tabcontent').data('shiny-input-binding').setValue($('#tabcontent'), 'report');"
                                        )
                                    )
                                )
                            )
                        )
                    )
                )
            ),
            # ------ * Help (top right) -----------------------------------
            div(
                id = "help",
                title = "Open the help guide",
                shinyhelper::helper(
                    shiny_tag = p("?"),
                    icon = "question-circle",
                    size = "l",
                    type = "markdown",
                    title = "Help",
                    content = "help",
                    buttonLabel = "Dismiss"
                )
            )
        )
    ),
    # ------ CONTENT -----------------------------------------------------------
    fluidRow(
        class = "bg-light full-screen d-flex flex-column justify-content-center",
        id = "content",
        tabsetPanel(
            id = "tabcontent",
            type = "hidden",
            tabPanelBody(
                value = "home",
                home_ui("home")
            ),
            tabPanelBody(
                value = "import",
                import_ui("import")
            ),
            tabPanelBody(
                value = "setup",
                setup_ui("setup")
            ),
            tabPanelBody(
                value = "pretest",
                pretest_ui("pretest")
            ),
            tabPanelBody(
                value = "sampling",
                sampling_ui("sampling")
            ),
            tabPanelBody(
                value = "report",
                report_ui("report")
            )
        )
    )
)
