library(shinyFiles)

ui <- fluidPage(
 
 titlePanel("File Browser"),
 
 sidebarLayout(
  sidebarPanel(
   
   shinyFilesButton('setup_grp1_r2_files', label = 'Group 1, R2 fastq files', 
                    title = 'Please select Group 1, R2 fastq file', multiple = TRUE),
   verbatimTextOutput("setup.grp1.r2.filepaths")
  ),
  
  mainPanel(
  )
 )
)


server <- function(input, output) {
 # browse from /filepath if it's specified otherwise /root
 if(file.exists("/filepath")){
  volumes=c(root="/filepath")
 }else{ 
  volumes=c(root="/root")
 }
 len.grp=reactiveVal(1)
 
 shinyFileChoose(input, 'setup_grp1_r2_files', root = volumes,
                 filetypes = c('',"gz"))
 
 file.lst <- reactiveValues(datapath="")
 react.file <- reactive(input$setup_grp1_r2_files)
 
 observeEvent(input$setup_grp1_r2_files,{
  file = react.file()
  new.path=as.character(parseFilePaths(root = volumes,file)$datapath)
  print("new.path")
  print(new.path)
  print("file.lst$datapath")
  print(file.lst$datapath)
  file.lst$datapath <- c(file.lst$datapath,new.path)
  print(file.lst)
 })
 
 output$setup.grp1.r2.filepaths <- renderText({
  paste0(file.lst$datapath,"\n")
 })
 
}   
shinyApp(ui = ui, server = server)