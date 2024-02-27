library(shiny)
library(shinyBS)

ui <- fluidPage(
 sidebarLayout(
  sidebarPanel(
   textInput("my_text", "Enter some text"),
   bsButton("update_button", "Update Text")
  ),
  mainPanel(
   verbatimTextOutput("updated_text")
  )
 )
)

server <- function(input, output, session) {
 observeEvent(input$update_button, {
  updated_text <- paste("Updated Text:", input$my_text)
  updateTextInput(session, "my_text", value = updated_text)
 })
 
 output$updated_text <- renderPrint({
  input$my_text
 })
}

shinyApp(ui, server)
