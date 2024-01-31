library(shiny)

ui <- fluidPage(
 conditionalPanel(
 condition = "output.plot != null",
  plotOutput("plot")
 )
)

server <- function(input, output) {
 # Server logic goes here
 
 output$plot <- renderPlot({
 
 })
}

shinyApp(ui, server)
