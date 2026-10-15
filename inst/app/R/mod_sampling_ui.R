sampling_ui <- function(id) {

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
                                h4("4. VCG sampling", class = "vcg-section-title"),
                                p(
                                    class = "vcg-section-hint",
                                    "The sampler selects whole animals from the POOL so that the",
                                    "selection sits as close to the TG as the pool allows."
                                ),
                                uiOutput(ns("recap")),
                                checkboxInput(
                                    inputId = ns("random"),
                                    label   = "Random subset sampling",
                                    value   = FALSE
                                ),
                                p(
                                    class = "vcg-section-hint",
                                    "Off: the best possible subset. On: animals are drawn with a",
                                    "probability derived from the energy weights, which adds randomness",
                                    "at the cost of balance quality."
                                ),
                                conditionalPanel(
                                    condition = "input.random",
                                    ns = ns,
                                    numericInput(
                                        inputId = ns("n_vcg"),
                                        label   = "Number of VCGs",
                                        value   = 1, min = 1, max = 100, step = 1
                                    ),
                                    p(
                                        class = "vcg-section-hint",
                                        "More than one draws several independent random VCGs",
                                        "(multiSampler), e.g. for a proof-of-concept study. Not",
                                        "meant for picking the nicest VCG afterwards."
                                    )
                                ),
                                numericInput(
                                    inputId = ns("seed"),
                                    label   = "Random seed",
                                    value   = 1234, min = 1, step = 1
                                ),
                                numericInput(
                                    inputId = ns("R"),
                                    label   = "Permutations for the post-hoc test",
                                    value   = 1000, min = 100, max = 10000, step = 100
                                ),
                                actionButton(
                                    inputId = ns("run"),
                                    label   = "Draw the VCG",
                                    class   = "btn-primary"
                                ),
                                tags$hr(class = "my-4"),
                                shinyjs::disabled(
                                    downloadButton(
                                        outputId = ns("download_csv"),
                                        label    = "Download VCG (CSV)",
                                        class    = "btn-primary"
                                    )
                                ),
                                div(
                                    class = "mt-2",
                                    shinyjs::disabled(
                                        downloadButton(
                                            outputId = ns("download_xlsx"),
                                            label    = "Download VCG (Excel)",
                                            class    = "btn-primary"
                                        )
                                    )
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
                                        title = "Balance",
                                        value = "balance",
                                        div(
                                            class = "mt-3",
                                            shinycssloaders::withSpinner(
                                                uiOutput(ns("balance_body"))
                                            )
                                        )
                                    ),
                                    tabPanel(
                                        title = "Covariates",
                                        value = "covariates",
                                        div(
                                            class = "mt-3",
                                            uiOutput(ns("covariate_body"))
                                        )
                                    ),
                                    tabPanel(
                                        title = "Sampled data",
                                        value = "data",
                                        div(
                                            class = "mt-3",
                                            uiOutput(ns("data_body"))
                                        )
                                    )
                                )
                            )
                        )
                    )
                ),
                p(
                    class = "text-end mt-4 mb-5",
                    shinyjs::disabled(
                        actionButton(
                            inputId = ns("continue"),
                            label = "Continue",
                            class = "btn-primary text-uppercase",
                            onclick = "$('#tabcontent').data('shiny-input-binding').setValue($('#tabcontent'), 'report');"
                        )
                    )
                )
            )
        )
    )
}
