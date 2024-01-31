library(shiny)

ui <- fluidPage(
 tags$head(
  tags$style(
   HTML("
        /* Set a consistent height for selectInput, textInput, and numericInput */
        .same-height {
          height: 30px; /* Set your desired height here */
        }
      ")
  )
 ),
 fluidRow(
  column(4, selectInput("select", "Select", choices = c("Option 1", "Option 2"), class = "same-height")),
  column(4, textInput("text", "Text", class = "same-height")),
  column(4, numericInput("numeric", "Numeric", value = 1, class = "same-height"))
 )
)

server <- function(input, output) {
 # Your server logic goes here
}

shinyApp(ui, server)
