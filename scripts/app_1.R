library(shiny)
library(shinyFiles)

ui <- fluidPage(
 shinyDirButton("folder", "Select Folder", "Please select a folder to delete"),
 textOutput("selected_folder"),  # Display the selected folder path
 actionButton("delete", "Delete Folder"),
 textOutput("status")
)

server <- function(input, output, session) {
 
 # Set up the root directory for folder selection
 roots <- c(home = normalizePath("~"))
 shinyDirChoose(input, "folder", roots = roots, filetypes = c('', 'txt', 'csv'))
 
 # Observe the folder selection and display the selected path
 observe({
  req(input$folder)
  
  # Get the selected folder path
  folder_path <- parseDirPath(roots, input$folder)
  
  # Display the folder path
  output$selected_folder <- renderText({
   paste("Selected folder:", folder_path)
  })
 })
 
 # Observe the delete button click and delete the selected folder
 observeEvent(input$delete, {
  req(input$folder)
  
  # Get the selected folder path
  folder_path <- parseDirPath(roots, input$folder)
  
  if (dir.exists(folder_path)) {
   # Attempt to delete the folder
   unlink(folder_path, recursive = TRUE)
   
   if (!dir.exists(folder_path)) {
    output$status <- renderText("Folder deleted successfully.")
   } else {
    output$status <- renderText("Failed to delete the folder.")
   }
  } else {
   output$status <- renderText("Folder does not exist.")
  }
 })
}

shinyApp(ui, server)
