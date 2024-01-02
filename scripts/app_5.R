library(shiny)

ui <- fluidPage(
 
 titlePanel("Text Area update"),
 
 sidebarLayout(
  sidebarPanel(
   actionButton('btnUpdate', 'Update')
  ),
  
  mainPanel(
   textAreaInput(inputId = 'testInput', 'Test', value = '', height = '400px')
  )
 )
)

server <- function(input, output, session) {
 observeEvent(input$btnUpdate, {
  updateTextAreaInput(session = session,
                      inputId = 'testInput',
                      value = 'abcd\ncde')
 })
}
shinyApp(ui = ui, server = server)
