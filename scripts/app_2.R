library(shinyFiles)

ui <- fluidPage(
 
 titlePanel("File Browser"),
 
 sidebarLayout(
  sidebarPanel(
   
   shinyFilesButton('files', label = 'Select', title = 'Please select a 
                       file', multiple = TRUE),
   verbatimTextOutput("filechosen")
  ),
  
  mainPanel(
  )
 )
)


server <- function(input, output) {
 
 shinyFileChoose(input, 'files', root = c(root = '/root'),
                 filetypes = c('',"gz"))
 
 file.lst <- reactiveValues(datapath="")
 react.file <- reactive(input$files)
 
 observeEvent(input$files,{
  file = react.file()
  new.path=as.character(parseFilePaths(c(root = "/root"),file)$datapath)
  print("new.path")
  print(new.path)
  print("file.lst$datapath")
  print(file.lst$datapath)
  file.lst$datapath <- c(file.lst$datapath,new.path)
  print(file.lst)
  })
 
 output$filechosen <- renderText({
  paste0(file.lst$datapath,"\n")
 })
 
}   
shinyApp(ui = ui, server = server)