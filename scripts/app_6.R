library(shiny)
library(dplyr)

# Sample list of data frames
data_list <- list(
 df1 = data.frame(
  Name = c("John", "Alice", "Bob", "David", "Eva"),
  Age = c(25, 30, 22, 35, 28),
  City = c("New York", "Los Angeles", "Chicago", "Houston", "San Francisco")
 ),
 df2 = data.frame(
  Name = c("Tom", "Mary", "Steve", "Anna", "Chris"),
  Age = c(22, 28, 35, 30, 40),
  City = c("Seattle", "Denver", "Atlanta", "Boston", "Dallas")
 )
)

ui <- fluidPage(
 numericInput("currentPage", "Current Page:", value = 1, min = 1, max = 1),
 numericInput("rowsPerPage", "Rows per page:", value = 2, min = 1),
 tableOutput("dataTable")
)

server <- function(input, output) {
 # Function to get paginated data
 get_paginated_data <- reactive({
  current_page <- input$currentPage
  rows_per_page <- input$rowsPerPage
  
  start_row <- (current_page - 1) * rows_per_page + 1
  end_row <- min(current_page * rows_per_page, sum(sapply(data_list, nrow)))
  
  data <- do.call(rbind, lapply(data_list, function(df) df))
  sliced_data <- slice(data, start_row:end_row)
  return(sliced_data)
 })
 
 # Render the paginated data table
 output$dataTable <- renderTable({
  as.data.frame(get_paginated_data())
 })
}

shinyApp(ui, server)
