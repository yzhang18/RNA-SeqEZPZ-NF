library(shiny)

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
 textInput("searchInput", "Search by Name:"),
 uiOutput("tables")
)

server <- function(input, output) {
 # Filter data frames based on the search input
 filtered_data <- reactive({
  search_term <- input$searchInput
  if (is.null(search_term) || search_term == "") {
   return(data_list)
  } else {
   filtered_list <- lapply(data_list, function(df) {
    df[grep(search_term, df$Name, ignore.case = TRUE), , drop = FALSE]
   })
   return(filtered_list)
  }
 })
 
 # Dynamically generate tables with titles
 output$tables <- renderUI({
  lapply(names(data_list), function(tab_name) {
   table_id <- paste0("table_", tab_name)
   tagList(
    h3(paste("Title for", tab_name)),
    tableOutput(table_id)
   )
  })
 })
 
 # Render the filtered data tables
 observe({
  for (tab_name in names(data_list)) {
   table_id <- paste0("table_", tab_name)
   output[[table_id]] <- renderTable({
    as.data.frame(filtered_data()[[tab_name]])
   })
  }
 })
}

shinyApp(ui, server)
