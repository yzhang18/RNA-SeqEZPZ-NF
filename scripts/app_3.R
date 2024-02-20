library(shiny)

ui <- fluidPage(
 titlePanel("Dynamic ReactiveValues Example"),
 
 sidebarLayout(
  sidebarPanel(
   # Input for the value to be stored
   numericInput("value", "Enter a value", 0),
   textInput("name", "Enter a name for the value")
  ),
  
  mainPanel(
   # Output: value stored in reactiveValues
   textOutput("output")
  )
 )
)

server <- function(input, output, session) {
 # Initialize an empty reactiveValues
 values <- reactiveValues()
 
 # Update reactiveValues with dynamic name when input$value and input$name change
 observeEvent(input$name, {
  name <- input$name
  if (name != "") {
   values[[name]] <- input$value
  }
 })
 
 # Render the value stored in reactiveValues
 output$output <- renderText({
  name <- input$name
  if (name != "") {
   paste("Value stored in reactiveValues under name", name, ":", values[[name]])
  } else {
   "Enter a name to store the value"
  }
 })
}

shinyApp(ui, server)

