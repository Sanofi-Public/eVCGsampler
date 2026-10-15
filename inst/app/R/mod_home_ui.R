home_ui <- function(id) {

    ns <- NS(id)

    fluidRow(
        class = "full-screen-home",
        column(
            width = 12,
            class = "full-screen-home",
            fixedPage(
                fixedRow(
                    column(
                        width = 5,
                        div(
                            class = "home-logo-wrap",
                            img(
                                src = "img/logo_vcg.png",
                                class = "home-logo",
                                alt = "eVCGsampler"
                            )
                        )
                    ),
                    column(
                        width = 7,
                        div(
                            h1(
                                "VCG sampler",
                                class = "mb-4 fs-1 fw-bold"
                            ),
                            div(
                                class = "mb-4",
                                p(
                                    "Draw a Virtual Control Group from your historical control pool that is",
                                    "balanced to the treated animals, using energy-distance based covariate",
                                    "balancing. Whole, unweighted animals are selected, so no special",
                                    "downstream analysis is needed."
                                ),
                                p(
                                    class = "text-muted",
                                    "Balancing is a final refinement step - it does not replace careful",
                                    "pre-filtering of the historical control data, and it does not address",
                                    "technical heterogeneity between studies. Your data is processed",
                                    "temporarily in this session and is not uploaded, stored or shared."
                                )
                            ),
                            # ------ WORKFLOW OVERVIEW -------------------------
                            div(
                                class = "card p-4 mb-4",
                                h3("The five steps", class = "fs-4 fw-bold mb-3"),
                                tags$ol(
                                    class = "mb-0",
                                    tags$li(tags$b("Data"), " - upload the POOL and the treated (TG) animals, in one file or two."),
                                    tags$li(tags$b("Setup"), " - pick the group variable, the covariates to balance, optional strata and the VCG size."),
                                    tags$li(tags$b("Pre-tests"), " - is there an imbalance worth correcting, and are any source studies outliers?"),
                                    tags$li(tags$b("Sampling"), " - draw the VCG and check the balance covariate by covariate."),
                                    tags$li(tags$b("Report"), " - download an HTML report and the sampled data set.")
                                )
                            ),
                            p(
                                class = "text-end mb-0",
                                actionButton(
                                    inputId = ns("start"),
                                    label = "Start sampling",
                                    class = "btn-primary text-uppercase",
                                    onclick = "$('#tabcontent').data('shiny-input-binding').setValue($('#tabcontent'), 'import');"
                                )
                            )
                        )
                    )
                )
            )
        )
    )
}
