library(shiny)

ui <- fluidPage(
 
 
 textInput(inputId = "id",
           label = 'Please enter your id'
 ),
 
 
 checkboxInput("agree", label = "I agree", value = FALSE),
 conditionalPanel(condition = "(input.submit_info >= 1) & ((input.id == '') || (!input.agree))",
                  
                  textOutput('error_msg')
 ),
 
 actionButton("submit_info", "Submit"),
 textOutput('success_msg')
 
 
)

server <- function(input, output) {
 
 output$error_msg <- renderText({
  shiny::validate(
   shiny::need(input$id != '', 'You must enter your id above to continue.'
   ),
   shiny::need(input$agree, "You must agree to continue")
  )
  
 })
 
 observeEvent(input$submit_info, {
  
  shiny::req(input$id)
  shiny::req(input$agree)
  output$success_msg <- renderText({"Success"})
  
 })
}

shinyApp(ui = ui, server = server)