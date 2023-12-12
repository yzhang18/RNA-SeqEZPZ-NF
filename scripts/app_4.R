library(shiny)

## Only run this example in interactive R sessions
 # Define UI
 ui <- fluidPage(
  actionButton("rmv", "Remove UI"),
  textInput("setup_txt2", "This is no longer useful")
 )
 
 # Server logic
 server <- function(input, output, session) {
  observeEvent(input$rmv, {
   removeUI(
    selector = "div:has(> #setup_txt2)"
   )
  })
 }
 
 # Complete app with UI and server components
 shinyApp(ui, server)
