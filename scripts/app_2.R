library(shiny)

ui <- fluidPage(
 # Add custom CSS for the scrolling text area
 tags$head(
  tags$style(HTML("
      #scrollable_Text {
        height: 400px;  /* Set height to limit the vertical size */
        overflow-y: auto;  /* Enable vertical scrolling */
        white-space: pre-wrap;  /* Preserve whitespace */
        background-color: #f5f5f5;
        border: 1px solid #ccc;
        padding: 10px;
      }
    "))
 ),
 
 # Verbatim Text Output wrapped in a div with the custom CSS class
 div(verbatimTextOutput("scrollable_Text"), id = "scrollable_Text")
)

server <- function(input, output, session) {
 # Example long text to demonstrate scrolling
 long_text <- paste(rep("This is a very long line of text that will wrap and scroll.", 100), collapse = "\n")
 
 output$scrollable_Text <- renderText({
  long_text
 })
}

shinyApp(ui, server)
