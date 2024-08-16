# app to select a folder, display folder size and delete with confirmation

library(shiny)
library(shinyFiles)

# Function to calculate directory size and return in human-readable format
get_dir_size <- function(path) {
 if (dir.exists(path)) {
  size <- sum(file.info(list.files(path, full.names = TRUE, recursive = TRUE))$size)
  # Convert size to human-readable format
  if (size < 1024) {
   size_str <- paste(size, "bytes")
  } else if (size < 1024^2) {
   size_str <- paste(round(size / 1024, 2), "KB")
  } else if (size < 1024^3) {
   size_str <- paste(round(size / 1024^2, 2), "MB")
  } else {
   size_str <- paste(round(size / 1024^3, 2), "GB")
  }
  return(size_str)
 } else {
  return("Directory does not exist.")
 }
}

ui <- fluidPage(
  shinyDirButton("folder", "Select Folder", "Please select a folder to delete"),
  textOutput("selected_folder"),  # Display the selected folder path
  textOutput("dir_size"),  # Display the selected folder path
  actionButton("delete", "Delete Folder"),
  textOutput("status")
)

server <- function(input, output, session) {
  
  # Set up the root directory for folder selection
  roots <- c(root = ("/filepath"))
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
    output$dir_size <- renderText({
     paste("Directory Size:", get_dir_size(folder_path), "bytes")
    })
  })
  
  # Observe the delete button click
  observeEvent(input$delete, {
   # Show a modal dialog asking for confirmation
   showModal(modalDialog(
    title = "Confirm Deletion",
    "Are you sure you want to delete this folder?",
    easyClose = FALSE,
    footer = tagList(
     modalButton("Cancel"),
     actionButton("confirm_delete", "Yes, Delete")
    )
   ))
  })
  
  # Handle the actual deletion when the user confirms
  observeEvent(input$confirm_delete, {
   # Close the modal dialog
   removeModal()
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
