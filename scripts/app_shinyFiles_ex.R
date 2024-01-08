# example app to dynamically add/remove shinyFiles objects

# remove user library path to avoid confusion
if(length(.libPaths())>1) .libPaths(.libPaths()[-1])

library(shiny)
library(shinyjs)
library(shinyFiles)

css <- "
.nav li a.disabled {
background-color: #d8d8d8 !important;
color: #8a8a8a !important;
cursor: not-allowed !important;
border-color: #aaa !important;
}"

volumes=c(root="/")

ui <- fluidPage(
 shinyjs::useShinyjs(),
 shinyjs::inlineCSS(css),
 # horizontal vertical bar for r1 and r2 filepaths input text
 # tags to have horizontal scrolling bar
 # white-space: pre will make newline character works. without this, newline chr in text will be ignored
 tags$style(HTML("#setup_grp1_r1_filepaths { width: 300px; overflow-x: auto; white-space: pre;}")),
 tags$style(HTML("#setup_grp1_r2_filepaths { width: 300px; overflow-x: auto; white-space: pre;}")),
 tags$style(HTML("#setup_grp2_r1_filepaths { width: 300px; overflow-x: auto; white-space: pre;}")),
 tags$style(HTML("#setup_grp2_r2_filepaths { width: 300px; overflow-x: auto; white-space: pre;}")),
 tabsetPanel(
  tabPanel("Run Analysis", id="run_analysis", fluid = TRUE,
           mainPanel(
            # adding horizontal scrolling
            style = "overflow-x: auto;", # Apply CSS styles
            width = 9,
            fluidRow(
             column(3,tags$label("Path to R1 fastq files")),
             column(3,tags$label("Path to R2 fastq files"))
            ),
            hr(style="border-color:black;margin-top:0px;margin-bottom:0px"),
            fluidRow(style = "background-color:#f9f9f9",
                     column(3,shinyFilesButton("setup_grp1_r1_files", 
                                               label = "Select R1 fastq files",
                                               title="Please select Group 1 R1 fastq files",multiple=TRUE ),
                            textAreaInput('setup_grp1_r1_filepaths',label="",value="",width='150px',height='100px')),
                     column(3,shinyFilesButton("setup_grp1_r2_files", 
                                               label = "Select R2 fastq files",
                                               title="Please select Group 1 R2 fastq files",multiple=TRUE ),
                            textAreaInput('setup_grp1_r2_filepaths',label="",width='150px',height='100px')),
            ), #fluidRow
            fluidRow(
             column(3,shinyFilesButton("setup_grp2_r1_files", 
                                       label = "Select R1 fastq files",
                                       title="Please select Group 2 R1 fastq files",multiple=TRUE ),
                    textAreaInput('setup_grp2_r1_filepaths',label="",value="",width='150px',height='100px')),
             column(3,shinyFilesButton("setup_grp2_r2_files", 
                                       label = "Select R2 fastq files",
                                       title="Please select Group 2 R2 fastq files",multiple=TRUE ),
                    textAreaInput('setup_grp2_r2_filepaths',label="",value="",width='150px',height='100px')),
            ), #fluidRow
            # tag$div id doesn't work if there is a dot!
            tags$div(id = "setup-placeholder"),
            hr(style="border-color:black;margin-top:0px;margin-bottom:0px"),
            br(),
            column(6,actionButton("setup.remove.set", "Remove row", width = "100%"),
                   style="padding-bottom:20px"),
            column(6,actionButton("setup.insert.set", "Add row", width = "100%"),
                   style="padding-bottom:20px")
           )#mainPanel
            )# tabPanel run analysis
  )#tabSetPanel
 )#fluidpage

server <- function(input, output,session) {
 # reactiveVal to keep track of number of added row
 # start with 2 as the minimum # of row
 setup.value <- reactiveVal(2)
 # initialize inserted div
 setup.inserted.div <- c()
# inserting row
observeEvent(input$setup.insert.set, {
 setup.btn <- setup.value() +1
 # update the reactiveVal
 setup.value(setup.btn)
 # id for new div
 setup.div.id <- paste0("setup_div_",setup.btn)
 setup.id.r1.files <- paste0("setup_grp",setup.btn,"_r1_files")
 setup.id.r2.files <- paste0("setup_grp",setup.btn,"_r2_files")
 setup.id.r1.filepath <- paste0("setup_grp",setup.btn,"_r1_filepaths")
 setup.id.r2.filepath <- paste0("setup_grp",setup.btn,"_r2_filepaths")
 if(setup.btn %% 2 == 0) row.color="#FFFFFF"
 if(setup.btn %% 2 != 0) row.color="#f9f9f9"
 print(paste0("background-color:",row.color))
 insertUI(
  selector = "#setup-placeholder",
  ui = tags$div(
   tags$style(HTML(paste0("#setup_grp",setup.btn,"_r1_filepaths { width: 300px; overflow-x: auto; white-space: pre;}"))),
   tags$style(HTML(paste0("#setup_grp",setup.btn,"_r2_filepaths { width: 300px; overflow-x: auto; white-space: pre;}"))),
   fluidRow(style = paste0("background-color:",row.color),
            column(3,shinyFilesButton(setup.id.r1.files,
                                      label = "Select R1 fastq files",
                                      title=paste("Please select Group",setup.btn,"R1 fastq files"),
                                      multiple=TRUE ),
                   textAreaInput(setup.id.r1.filepath,label="",value="",width='150px',height='100px')),
            column(3,shinyFilesButton(setup.id.r2.files,
                                      label = "Select R2 fastq files",
                                      title=paste("Please select Group",setup.btn,"R1 fastq files"),
                                      multiple=TRUE ),
                   textAreaInput(setup.id.r2.filepath,label="",value="",width='150px',height='100px')),
            
   ),#fluidrow
   id = setup.div.id
  )
 )#insertUI
 setup.inserted.div <<- c(setup.inserted.div, setup.div.id)
}) #observeEvent input$setup.insert.set


observeEvent(input$setup.remove.set,{
 if(setup.value()>2){
  setup.btn <- setup.value()-1
  # update the reactiveVal
  setup.value(setup.btn)
  print(paste0("#",setup.inserted.div[length(setup.inserted.div)]))
  removeUI(
   selector = paste0("#",setup.inserted.div[length(setup.inserted.div)])
  )
  # remove the last values
  updateTextInput(
   session,
   paste0("setup_grp", length(setup.inserted.div)+1,"_r1_files"),
   NULL,
   ""
  )
  updateTextInput(
   session,
   paste0("setup_grp", length(setup.inserted.div)+1,"_r2_files"),
   NULL,
   ""
  )
  updateTextInput(
   session,
   paste0("setup.grp", length(setup.inserted.div)+1,".r1.filepath"),
   NULL,
   ""
  )
  updateTextInput(
   session,
   paste0("setup.grp", length(setup.inserted.div)+1,".r2.filepath"),
   NULL,
   ""
  )
  setup.inserted.div <<- setup.inserted.div[-length(setup.inserted.div)]
 }else{ # Do not remove row if there's only 2 rows
  setup.btn <- 2
  setup.value(setup.btn)
 }
})

updateFileChoose <- function(nsamples) {
 print(nsamples)
 for (i in 1:nsamples) {
  setup.grp.r1.name=paste0('setup_grp',i,'_r1_files')
  shinyFiles::shinyFileChoose(input, setup.grp.r1.name, root=volumes,
                              filetypes=c('', 'gz'),session=session)
  setup.grp.r2.name=paste0('setup_grp',i,'_r2_files')
  shinyFiles::shinyFileChoose(input, setup.grp.r2.name, root=volumes,
                              filetypes=c('', 'gz'),session=session)
 }
}

observe({
 nsamples <- setup.value()
 updateFileChoose(nsamples)
})


updateFilelst <- function(nsamples){ 
 # Initialize r1/r2 file list as reactiveValues to save previous choices
 x <- as.list(rep("",nsamples))
 names(x)=lapply(1:nsamples,function(i)paste0('grp',i))
 setup.grp.r1.file.lst <- do.call("reactiveValues",x)
 setup.grp.r2.file.lst <- do.call("reactiveValues",x)
 
 # getting all the r1 and r2 fastq files for all groups
 lapply(
  1:nsamples,
  function(i){
   observeEvent(input[[paste0('setup_grp',i,'_r1_files')]],{
    file = input[[paste0('setup_grp',i,'_r1_files')]]
    new.path=as.character(parseFilePaths(root = volumes,file)$datapath)
    # getting the existing files in the textbox
    setup.grp.r1.file.lst[[paste0('grp',i)]] <- input[[paste0('setup_grp',i,'_r1_filepaths')]]
    print("setup.grp.r1.file.lst[[paste0('grp',i)]]")
    print(setup.grp.r1.file.lst[[paste0('grp',i)]])
    setup.grp.r1.file.lst[[paste0('grp',i)]] <- c(setup.grp.r1.file.lst[[paste0('grp',i)]],new.path)
    print("setup.grp.r1.file.lst[[paste0('grp',i)]]")
    print(setup.grp.r1.file.lst[[paste0('grp',i)]])
    
    path.display=setup.grp.r1.file.lst[[paste0('grp',i)]]
    # remove the first entry if empty
    if(path.display[1]=="") path.display=path.display[-1]
    
    path.display=gsub('^/filepath/','',path.display)
    path.display=gsub('^/root/','',path.display)
    path.display <- paste0(path.display, collapse="\n")
    print("path.display")
    print(path.display)
    updateTextAreaInput(session, paste0('setup_grp',i,'_r1_filepaths'),value=path.display)
   })# observeEvent r1_files
   
   observeEvent(input[[paste0('setup_grp',i,'_r2_files')]],{
    file = input[[paste0('setup_grp',i,'_r2_files')]]
    new.path=as.character(parseFilePaths(root = volumes,file)$datapath)
    setup.grp.r2.file.lst[[paste0('grp',i)]] <- c(setup.grp.r2.file.lst[[paste0('grp',i)]],new.path)
    path.display.r2=setup.grp.r2.file.lst[[paste0('grp',i)]]
    # remove the first entry if empty
    if(path.display.r2[1]=="") path.display.r2=path.display.r2[-1]
    path.display.r2=gsub('^/filepath/','',path.display.r2)
    path.display.r2=gsub('^/root/','',path.display.r2)
    path.display.r2 <- paste0(path.display.r2, collapse="\n")
    print("path.display.r2")
    print(path.display.r2)
    updateTextAreaInput(session, paste0('setup_grp',i,'_r2_filepaths'),value=path.display.r2)
   })# observeEvent
  } #function
 )#lapply
}#updateFilelst


observe({
 # getting nsamples from add/remove button
 nsamples <- setup.value()
 updateFilelst(nsamples)
})

} #server

# app
shinyApp(ui = ui, server = server)