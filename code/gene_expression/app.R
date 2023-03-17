library(shiny)
library(ggplot2)

matrix_data <- as.data.frame(matrix(rnorm(30), ncol = 3))
rownames(matrix_data) <- paste0("R",seq_along(1:nrow(matrix_data)))

ui <- fluidPage(
  titlePanel("Matrix Visualization"),

  sidebarLayout(
    sidebarPanel(
      selectInput("row_name", "Select Row:", rownames(matrix_data))
      ),

    mainPanel(
      plotOutput("scatterplot")
    )
  )
)

server <- function(input, output) {
  data <- reactive({
    row_num <- input$row_name
    plot_data <- matrix_data[row_num, , drop=F] %>%
      gather()
  })

  output$scatterplot <- renderPlot({
    ggplot(data(), aes(x = key, y = value)) +
      geom_point()
  })
}

shinyApp(ui, server)
