library(shiny)
library(shinyFiles)

allinhome <- list.files("~")

ui <- fluidPage(
 selectInput("2LevelFolder", 
             "Select a folder in home directory.",
             choices = allinhome),
 
 shinyFiles::shinyFilesButton("chooseFile", 
                              "Explore and choose file", 
                              "This is the title",
                              multiple = TRUE),
 
 verbatimTextOutput("path")
 
)

server <- function(input, output, session) {
 
 updateFileChoose <- function(folder) {
  print(folder)
  shinyFiles::shinyFileChoose(input = input, "chooseFile", 
                              roots = c(chosenFolder = file.path("~", folder)), 
                              session = session)
 }
 
 observe({
  updateFileChoose(input$`2LevelFolder`)
 })
 
 output$path <- renderPrint({
  input$chooseFile
 })
}

shinyApp(ui, server)