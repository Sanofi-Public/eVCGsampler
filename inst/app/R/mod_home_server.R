home_server <- function(
        id
) { moduleServer(id, function(input, output, session) {

    # The welcome screen is static; the "Start sampling" button switches the
    # tab client-side (see the onclick in home_ui). The observer is only kept
    # so the parent server can react to the start of a workflow if needed.
    list(
        start = reactive(input$start)
    )

})}
