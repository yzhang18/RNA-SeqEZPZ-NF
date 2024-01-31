library(shiny)
ui <- pageWithSidebar(
 headerPanel("Click the button"),
 sidebarPanel(
  sliderInput("obs", "Number of observations:",
              min = 0, max = 1000, value = 500),
  actionButton("goButton", "Go!")
 ),
 mainPanel(
  plotOutput("distPlot")
 )
)
server <- function(input, output) {
 react.obs <- reactive({
  input$obs
 })
 output$distPlot <- renderPlot({
  
  # Take a dependency on input$goButton
  input$goButton
  
  # Use isolate() to avoid dependency on input$obs
  obs <- isolate(react.obs())
  dist <- rnorm(obs)
  hist(dist)
 })
}


shinyApp(ui, server)