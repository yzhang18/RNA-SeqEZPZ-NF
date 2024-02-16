library(shiny)

ui <- fluidPage(
 plotOutput("plot"),
 downloadButton("download_plot", "Download Plot as PDF")
)

server <- function(input, output) {
 output$plot <- renderPlot({
  # Your plot code here
  plot(1:10, (1:10)^2, type = "l", col = "blue", lwd = 2)
 })
 
 observeEvent(input$download_plot, {
  pdf("plot.pdf")
  print(output$plot)
  dev.off()
 })
}

shinyApp(ui, server)
