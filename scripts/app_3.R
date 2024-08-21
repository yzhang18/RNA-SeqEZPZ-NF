library(shiny)
library(shinydashboard)

ui <- fluidPage(
 titlePanel("Shiny App with Separated Rows Using Box"),
 
 # First row with content
 fluidRow(
  column(12,
         h3("This is the first row"),
         p("Some content for the first row goes here.")
  )
 ),
 
 # Box to separate rows
 fluidRow(
  div(style = "border-style: solid; border-color: black;",
          p("This box separates different sections of the app."))
  
 ),
 
 # Second row with content
 fluidRow(
  column(12,
         h3("This is the second row"),
         p("Some content for the second row goes here.")
  )
 )
)

server <- function(input, output, session) {
 # Server logic, if any
}

shinyApp(ui, server)
