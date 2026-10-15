pretest_ui <- function(id) {

    ns <- NS(id)

    fluidRow(
        class = "full-screen-home border-top border-2 border-dark",
        column(
            width = 12,
            fixedPage(
                fluidRow(
                    class = "mt-4",
                    # ------ SETTINGS -------------------------------------------
                    column(
                        width = 4,
                        div(
                            class = "card h-100",
                            div(
                                class = "card-container vcg-side",
                                h4("3. Pre-tests", class = "vcg-section-title"),
                                p(
                                    class = "vcg-section-hint",
                                    "Two questions before sampling: is there an imbalance worth",
                                    "correcting, and does any source study stand out from the rest?"
                                ),
                                h6("Energy test - POOL vs TG", class = "fw-bold text-uppercase mt-3"),
                                numericInput(
                                    inputId = ns("alpha"),
                                    label   = "Significance level (alpha)",
                                    value   = 0.05, min = 0.001, max = 0.2, step = 0.005
                                ),
                                numericInput(
                                    inputId = ns("R"),
                                    label   = "Permutations",
                                    value   = 1000, min = 100, max = 10000, step = 100
                                ),
                                actionButton(
                                    inputId = ns("run_energy"),
                                    label   = "Run energy test",
                                    class   = "btn-primary"
                                ),
                                tags$hr(class = "my-4"),
                                h6("Outlier studies", class = "fw-bold text-uppercase"),
                                p(
                                    class = "vcg-section-hint",
                                    "Meant for outcome or technical variables, not for the covariates",
                                    "you are balancing - those are handled by the sampling itself."
                                ),
                                uiOutput(ns("study_ui")),
                                uiOutput(ns("outvars_ui")),
                                uiOutput(ns("poolonly_ui")),
                                numericInput(
                                    inputId = ns("cutoff"),
                                    label   = "Cut-off percentile",
                                    value   = 0.99, min = 0.5, max = 0.999, step = 0.01
                                ),
                                numericInput(
                                    inputId = ns("R_out"),
                                    label   = "Permutations",
                                    value   = 500, min = 100, max = 2000, step = 50
                                ),
                                actionButton(
                                    inputId = ns("run_outlier"),
                                    label   = "Run outlier test",
                                    class   = "btn-primary"
                                )
                            )
                        )
                    ),
                    # ------ RESULTS --------------------------------------------
                    column(
                        width = 8,
                        div(
                            class = "card",
                            div(
                                class = "card-container",
                                tabsetPanel(
                                    id = ns("result_tabs"),
                                    tabPanel(
                                        title = "Energy test",
                                        value = "energy",
                                        div(
                                            class = "mt-3",
                                            uiOutput(ns("energy_formula")),
                                            uiOutput(ns("energy_tiles")),
                                            uiOutput(ns("energy_verdict")),
                                            shinycssloaders::withSpinner(
                                                uiOutput(ns("energy_body"))
                                            )
                                        )
                                    ),
                                    tabPanel(
                                        title = "Outlier studies",
                                        value = "outlier",
                                        div(
                                            class = "mt-3",
                                            shinycssloaders::withSpinner(
                                                uiOutput(ns("outlier_body"))
                                            )
                                        )
                                    )
                                )
                            )
                        )
                    )
                ),
                p(
                    class = "text-end mt-4 mb-5",
                    actionButton(
                        inputId = ns("continue"),
                        label = "Continue",
                        class = "btn-primary text-uppercase",
                        onclick = "$('#tabcontent').data('shiny-input-binding').setValue($('#tabcontent'), 'sampling');"
                    )
                )
            )
        )
    )
}
