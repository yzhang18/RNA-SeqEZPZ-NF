# this app is based on compare_gene_list.R

# remove user library path to avoid confusion
if(length(.libPaths())>1) .libPaths(.libPaths()[-1])

library(shiny)
library(GeneOverlap)
library(gridExtra)
library(ggplot2)
library(gridGraphics)
library(eulerr)
library(reshape)
library(ggthemes)
library(ggpointdensity)
library(shiny)
library(RColorBrewer)
library(grid)
library(gridBase)
library(venn) # need to be 1.11
library(UpSetR)
library(DESeq2)
library(clusterProfiler) # need to be 3.8.1
library(msigdbr)
library(stringr)
library(ggplotify) #as.ggplot for upsetr
library(ggpolypath) #ggplot=true for venn
library(shinyFiles)
library(shinyjs)
library(ggrepel)
library(dplyr) #bind_rows, case_when
#library(plotly)
#library(DT)
library(purrr)
#library(tidyverse)

css <- "
.nav li a.disabled {
background-color: #d8d8d8 !important;
color: #8a8a8a !important;
cursor: not-allowed !important;
border-color: #aaa !important;
}"

#default hostfilepath to "/" if not exist
if(!exists("hostfilepath")) hostfilepath="/"
if(exists("hostfilepath") && substr(hostfilepath,nchar(hostfilepath),nchar(hostfilepath))!="/") 
 hostfilepath=paste0(hostfilepath,"/")

# default image directory used to run analysis out of container
#if(!exists("img.dir")) img.dir="~/Steve/virtual_server/rnaseq-singularity"

# ad-hoc variable to prevent error in initialization
log.lst=""
log.fail.lst=""


 # default checkboxes for venn.opts	
 venn.opts=data.frame(
  lbl= c("Numbers","Percentages","Labels","Legend")
 )
 venn.opts.lst=as.list(venn.opts$lbl)
 names(venn.opts.lst)=venn.opts$lbl
 
 # list of choices for msigdb species
 msigdb.species.lst <- as.list(msigdbr_species())$species_name
 
 ngroup=3
 tab0.ngroup=2
 
 colors=c("#F08080","#ADD8E6","#FFFACD")
 tab0.colors=c("#F08080","#ADD8E6","#FFFACD")
 tab3.colors=brewer.pal(12,"Paired")
 # color for pvalues (<2.2e-16,<=0.05,<0.5,>0.5)
 #col.pval=c("#fc8d59","#ffffbf","#91bfdb")
 col.pval = brewer.pal(11,"RdBu")[c(3,4,8,9)]
 # names for factorized pvalues
 pval.names=c("< 2.2e-16","< 0.05","< 0.5","> 0.5")
 # ad-hoc maximum plot. 
 max_plots=5

# arbitrary variables to prevent errors
 groups.lst=c("a","b","c")
 
# List of choices for genome
setup.genome.lst <- as.list(c("hg19","hg38","other"))
# ad-hoc max.nsamples if not exist set
if(!exists("max.nsamples")) max.nsamples=50

heatmap_sigf_overlap <- function(data,title=""){
 hm <- ggplot(data,aes(X1,X2,fill=pval.cat)) +
  # use thin border size 0.1 to separate tiles
  geom_tile(color="white",size=0.1)+
  # setting high and low color and color key name
  scale_fill_manual(values=col.pval,na.value = "transparent",drop=FALSE)+
  # draw nice square
  coord_equal() +
  # remove x and y labels
  labs(x=NULL, y=NULL,fill="p-values") +
  # removes a lot of chart junks including space around heatmap
  # from ggthemes package
  theme_tufte(base_family="Helvetica") +
  # remove tick marks on axes
  theme(axis.ticks=element_blank()) +
  # make text bigger
  theme(axis.text = element_text(size=14)) +
  # rotate x axis text
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  # make legend on top left
  theme(legend.position="top") +
  #theme(legend.justification = c(-1,0)) +
  # add pvalue as text
  geom_text(aes(label=round(jaccard,digits=2))) +
  # color key size
  theme(legend.key.width=unit(0.7,"cm")) +
  ggtitle(title)
 # plot venn and overlap heatmap together
 grid.arrange(hm,
              top=textGrob(""),
              bottom=textGrob(paste0("Color represents p-values of overlap.\n",
                                     "Number inside a cell is Jaccard similarity index.\n",
                                     "0 means no similarity, 1 means identical."),just=c("left"),
                              hjust=0,x=0.3,gp=gpar(fontsize=12)))	 
 # grob <- grobTree(textGrob(paste("Color represents p-values of overlap.\n",
 # "Number inside a cell is Jaccard similarity index.\n",
 # "0 means no similarity, 1 means identical."), x=0,  y=0, hjust=0,
 # gp=gpar(fontsize=13, fontface="italic")))
 # hm + annotation_custom(grob)
}

plot_euler <- function(s4,colors,cex,venn.opts,title){
 quantities=list(cex=cex);
 legend=FALSE
 if('Numbers' %in% venn.opts) quantities$type=c(quantities$type,"counts")
 if('Percentages' %in% venn.opts) quantities$type=c(quantities$type,"percent")
 if('Legend' %in% venn.opts) legend=list(side="bottom")
 
 if('Labels' %in% venn.opts ){ 
  labels=list(names(s4),cex=cex)
 }else {
  labels=""
 }
 if(is.null(venn.opts) | 
    sum(c('Numbers','Percentages') %in% venn.opts)==0) quantities=FALSE
 
 set.seed(1)
 plot(euler(s4,shape="ellipse"),
      quantities=quantities, labels=labels,
      fills=colors,edges=FALSE,main=title,
      adjust_labels = TRUE,legend=legend)
}

ui <- fluidPage(
 title = "Run full RNA-Seq analysis",
 shinyjs::useShinyjs(),
 shinyjs::inlineCSS(css),
 # horizontal vertical bar for r1 and r2 filepaths input text
 # tags to have horizontal scrolling bar
 # white-space: pre will make newline character works. without this, newline chr in text will be ignored
 tags$style(HTML("#setup_grp1_r1_filepaths { width: 300px; overflow-x: auto; white-space: pre;}")),
 tags$style(HTML("#setup_grp1_r2_filepaths { width: 300px; overflow-x: auto; white-space: pre;}")),
 tags$style(HTML("#setup_grp2_r1_filepaths { width: 300px; overflow-x: auto; white-space: pre;}")),
 tags$style(HTML("#setup_grp2_r2_filepaths { width: 300px; overflow-x: auto; white-space: pre;}")),
 tabsetPanel(id="inTabset",
  tabPanel("Run Analysis", id="run_analysis", fluid = TRUE,
           sidebarPanel(width=3,
                        shinyDirButton("setup_proj_dir", 
                                         label = "Select project folder",
                                         title="Please select a folder to run pipeline",multiple=FALSE),
                        verbatimTextOutput('setup.proj.dir.filepaths'),
                        br(),
                        #verbatimTextOutput("fileExists"),
                        conditionalPanel(condition="output.fileExists && !input.setup_proj_dir == '' ",
                               actionButton("setup.load.samples","Click to load existing samples.txt"),
                        br(),
                        br(),
                        ),
                        
                        selectInput(inputId = "setup.genome",label="Select genome",
                                    choices=setup.genome.lst,selected=1),
                        conditionalPanel(
                         condition= "input['setup.genome'] == 'other'",
                         textInput(inputId="setup.genome.name", 
                                   label = "Genome name",value=""),
                         shinyFilesButton("setup_genome_fa", 
                                          label = "Select genome fasta file",
                                          title="Please select genome fasta file",multiple=FALSE),
                         verbatimTextOutput('setup.genome.fa.filepaths'),
                         br(),
                         shinyFilesButton("setup_genome_gtf", 
                                          label = "Select genome GTF file",
                                          title="Please select genome GTF file",multiple=FALSE ),
                         verbatimTextOutput('setup.genome.gtf.filepaths'),
                         br()
                        ),#conditionPanel
                        textInput(inputId = "setup.time",label="Time limit",
                                  value="7-00:00:00",),
                        
                        checkboxInput(inputId = "setup.batch.adjust",label="Batch adjustment",
                                      value=TRUE),
                        
                        numericInput(inputId = "setup.ncpus.trim",label="# of CPUs for trimming",
                                     value=4,min=1),
                        
                        numericInput(inputId = "setup.ncpus.star",label="# of CPUs for STAR alignment",
                                     value=20,min=1),
                        checkboxInput(inputId = "setup.debug",label="Run debug",
                                      value=FALSE),
                        conditionalPanel(condition = "input.setup.proj.dir == ''",
                        textOutput('err_msg')),
                        #conditionalPanel(condition = "input.setup.genome == 'other' && (input.setup.genome.name == '' || input.setup_genome_fa == '' || input.setup_genome_gtf == '')",
                                         textOutput('err_msg_genome'),
                        actionButton(inputId="setup.run.analysis",label="Run full analysis",
                        )
           ),#sidebarPanel
           mainPanel(
            # adding horizontal scrolling
            style = "overflow-x: auto;", # Apply CSS styles
            width = 9,
                     fluidRow(
                         column(4,textInput(inputId="setup.email", 
                                         label = "Email address",value="" ))),
                     fluidRow(
                      column(2,tags$label("Group name"),style="padding-left:30px;padding-top:0px"),
                      column(2,tags$label("Control name")),
                      column(2,tags$label("Replicate name")),
                      column(3,tags$label("Path to R1 fastq files")),
                      column(3,tags$label("Path to R2 fastq files"))
                     ),
                     hr(style="border-color:black;margin-top:0px;margin-bottom:0px"),
                     fluidRow(style = "background-color:#f9f9f9",
                              column(2,
                                     column(2, 
                                            tags$label("1.", style = "padding-top: 30px;"),
                                            style = "padding-right: 0px;padding-left:0px"),
                                     column(10,
                                 
                                           textInput(inputId="setup.grp1.name", 
                                                         label = "",value=""),
                                           style ="padding-right:0px;padding-left:0px")
                              ),
                              column(2,textInput(inputId="setup.grp1.ctrl.name", 
                                                 label = "",value="" )),
                              column(2,textInput(inputId="setup.grp1.rep.name", 
                                                 label = "",value="" )),
                              br(),
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
                      column(2,column(2, 
                                      tags$label("2.", style = "padding-top: 30px;"),
                                      style = "padding-right: 0px;padding-left:0px"),
                             column(10,textInput(inputId="setup.grp2.name", 
                                                 label = "",value="" ),
                                    style ="padding-right:0px;padding-left:0px")
                      ),
                      column(2,textInput(inputId="setup.grp2.ctrl.name", 
                                         label = "",value="" )),
                      column(2,textInput(inputId="setup.grp2.rep.name", 
                                         label = "",value="" )),
                      br(),
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
                            style="padding-bottom:20px"),
           ),
           
  ),# tabPanel run analysis
  tabPanel("Log", id="Log", fluid = TRUE,
           sidebarPanel(width=3,
                        # Input: Choose a log file to view
                        actionButton("logtab.refresh.log.path","Refresh list"),
                        br(),
                        br(),
                        selectInput("logtab.log.path", "Choose a log file to view:",
                                    choices = log.lst,selected="run_rnaseq_full.out"),
                        br(),
                        # Input: Choose a failed log file to view
                        selectInput("logtab.fail.log.path", "Choose a failed log file to view:",
                                    choices = log.fail.lst,selected="")
                        ),
           mainPanel(width=9,
                     tags$div(
                      h5("Log file content:"),
                      tags$style(type="text/css", ".shiny-text-output {word-wrap: break-word;}"),
                     verbatimTextOutput("logtab.log.content"),
                      h5("Failed log content:"),
                     verbatimTextOutput("logtab.fail.log.content")))
                     
  ),# tabPanel Log
  tabPanel("QCs",id="qctab",fluid=TRUE,
            #uiOutput("multiqc")
           htmlOutput("multiqc")
  ), #tabPanel Outputs
  tabPanel("Outputs",id="outputtab",fluid=TRUE,
            htmlOutput("diff_report")
           ), #tabPanel Outputs
  tabPanel("Plots", value="tab3",fluid = TRUE,
           sidebarLayout(
            sidebarPanel(
             tags$style(type='text/css', 
                        ".selectize-input { 
					 font-size: 10pt; 
					 line-height: 10pt;} 
					 .selectize-dropdown {
					 font-size: 10pt; 
					 line-height: 10pt;}"),
             width=2,
             selectInput(
              inputId="tab3.grp1.name",
              label = ("Group 1"),
              selected=groups.lst[1],
              choices = groups.lst),
             textInput(inputId="tab3.grp1.plot.title", 
                       label = ("Group 1 label"), 
                       value = groups.lst[1]),
             
             tags$div(id = "placeholder"),
             
             splitLayout(
              actionButton("insert_set", "Insert", width = "100%"),
              actionButton("remove_set", "Remove", width = "100%")
             ),
             br(),
             br(),
             splitLayout(
              numericInput("pdf.width", label = HTML("PDF <br/> width (in)"), 
                           value = 11),
              numericInput("pdf.height", label = HTML("PDF <br/> height (in)"), 
                           value = 8)),
             actionButton(inputId='export', label="Save plots.pdf"),
             br(),
             br(),
             actionButton(inputId='exportGene', label="Export filtered gene lists")
            ),
            mainPanel(
             fluidRow(style = "background-color:#F5F5F5",
                      column(2,
                             checkboxGroupInput("tab3.venn.opts", label = ("Display on Venns"), 
                                                choices = venn.opts.lst,
                                                selected = venn.opts$lbl[1:3]),
                             
                             # slider bar font size
                             sliderInput("tab3.venn.cex", label = ("Venns' font size"), min = 0, 
                                         max = 6, value = 1,step = 0.1),
                             # slider bar pad
                             sliderInput("tab3.venn.pad", label = ("Venns' circle size"), min = 0, 
                                         max = 8, value = 1,step = 0.5),
                             # upset plot number of intersection
                             numericInput("tab3.nintersects.upset", label = ("Upset plot max intersections"), 
                                          value = 40,min=1,step=1),
                      ),
                      column(2,
                             textInput("tab3.color.grp1", label = HTML("Color for <br/> group 1"), value = tab3.colors[1]),
                             
                             tags$div(id = "placeholder-col")),
                      column(2,
                             numericInput("tab3.fdr.grp1", label = ("FDR cut-off for group 1"), 
                                          value = 0.05,step=0.01,min=0,max=1),
                             tags$div(id = "placeholder-fdr")),
                      column(2,
                             numericInput("tab3.fc.grp1", label = ("Fold-change for group 1"), 
                                          value = 1,min=1,step=0.5),
                             tags$div(id = "placeholder-fc")),
                      column(2,
                             numericInput("tab3.meanDiff.grp1", label = ("Difference for group 1"), 
                                          value = 0,min=0,step=1),
                             tags$div(id = "placeholder-meanDiff"))),
             br(),
             br(),
             # Table of diff genes
             fluidRow(style = "background-color:#F5F5F5",
              textInput("tab3.search.term", "Search by gene names separated by commas:"),
              numericInput("currentPage","Current Page",value=1,min=1),
              uiOutput("tables")
             ),
             # This is the dynamic UI for the plots
             br(),
             br(),
             fluidRow(
              # genes to highlight in volcano plot
              column(5,textInput("tab3.hilite.genes",
                                 label = "Enter official gene names separated by comma:",value="")),
              column(12,uiOutput("tab3.plots"))),
             br(),
             br(),
             fluidRow(
              conditionalPanel(
               condition="output.tab3.plot1 != null",
              plotOutput(outputId = "tab3.plot1")
             )
             ),
             br(),
             br(),
             fluidRow(
              conditionalPanel(
               condition="output.tab3.plot2 != null",
              column(12,plotOutput(outputId = "tab3.plot2"))
              )
             ),
             br(),
             br(),
             fluidRow(
             # species for msigdb
             selectInput(
              inputId="tab3.msigdb.species",
              label = ("Select your species"),
              selected="Homo sapiens",
              choices = msigdb.species.lst),
             # gene set enrichment pvalue numeric input
             numericInput("tab3.enrich.pval.co", label = ("FDR for enrichment"), 
                          value = 0.05,min=0,step=0.01),
             actionButton(
              inputId = "gen.go",
              label = "Generate Enrichment plots"
             ),
             br(),
             br(),
             conditionalPanel(
              condition= "input['gen.go'] >= 1",
              column(12,textInput("tab3.cap.enrich.terms", label = ("Enter terms to capitalize in enrichment plots"),
                        value="")),
               column(12,plotOutput(outputId = "tab3.plot3.1",height="auto")),
         
              br(),
               column(12,plotOutput(outputId = "tab3.plot3.2",height="auto")),
              br(),
               column(12,plotOutput(outputId = "tab3.plot3.3",height="auto")),
              br(),
               column(12,plotOutput(outputId = "tab3.plot3.4",height="auto")),
             )#fluidrow
              ), # conditionalpanel
             textOutput("tab3.text")
            )
           )
  )		
  ,selected="Run Analysis")#tabsetpanel
)#fluidpage



# server
server <- function(input, output,session) {

 shinyjs::disable(selector = 'a[data-value="Two Groups"')
 shinyjs::disable(selector = 'a[data-value="Three Groups"')
 #shinyjs::disable(selector = 'a[data-value="Plots"')

 #write.table(isolate(session$clientData$url_port),file="/mnt/outputs/port.txt",
 #	quote=FALSE,row.names=FALSE,col.names=FALSE)
 
 #### run analysis tab ######
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
  # ids for entries in each row
  setup.id.grp.name <- paste0("setup.grp",setup.btn,".name")
  setup.id.ctrl.name <- paste0("setup.grp",setup.btn,".ctrl.name")
  setup.id.rep.name <- paste0("setup.grp",setup.btn,".rep.name")
  setup.id.r1.files <- paste0("setup_grp",setup.btn,"_r1_files")
  setup.id.r2.files <- paste0("setup_grp",setup.btn,"_r2_files")
  setup.id.r1.filepath <- paste0("setup_grp",setup.btn,"_r1_filepaths")
  setup.id.r2.filepath <- paste0("setup_grp",setup.btn,"_r2_filepaths")
  if(setup.btn %% 2 == 0) row.color="#FFFFFF"
  if(setup.btn %% 2 != 0) row.color="#f9f9f9"
  insertUI(
   selector = "#setup-placeholder",
   ui = tags$div(
    tags$style(HTML(paste0("#setup_grp",setup.btn,"_r1_filepaths { width: 300px; overflow-x: auto; white-space: pre;}"))),
    tags$style(HTML(paste0("#setup_grp",setup.btn,"_r2_filepaths { width: 300px; overflow-x: auto; white-space: pre;}"))),
    fluidRow(style = paste0("background-color:",row.color),
             column(2,
                    column(2,
                           tags$label(paste0(setup.btn,"."), style = "padding-top: 30px;"),
                           style = "padding-right: 0px;padding-left:0px"),
                    column(10,textInput(inputId=setup.id.grp.name,
                                        label = "",value=""),
                                        style = "padding-right: 0px;padding-left:0px")
             ),
             column(2,textInput(inputId=setup.id.ctrl.name,
                                label = "",value="" )),
             column(2,textInput(inputId=setup.id.rep.name,
                                label = "",value="" )),
             br(),
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
   removeUI(
    selector = paste0("#",setup.inserted.div[length(setup.inserted.div)])
   )
   # remove the last values
   updateTextInput(
    session,
    paste0("setup.grp", length(setup.inserted.div)+2,".name"),
    NULL,
    ""
   )
   updateTextInput(
    session,
    paste0("setup.grp", length(setup.inserted.div)+2,"ctrl.name"),
    NULL,
    ""
   )
   updateTextInput(
    session,
    paste0("setup.grp", length(setup.inserted.div)+2,".rep.name"),
    NULL,
    ""
   )
   updateTextInput(
    session,
    paste0("setup_grp", length(setup.inserted.div)+2,"_r1_files"),
    NULL,
    ""
   )
   updateTextInput(
    session,
    paste0("setup_grp", length(setup.inserted.div)+2,"_r2_files"),
    NULL,
    ""
   )
   updateTextInput(
    session,
    paste0("setup.grp", length(setup.inserted.div)+2,".r1.filepath"),
    NULL,
    ""
   )
   updateTextInput(
    session,
    paste0("setup.grp", length(setup.inserted.div)+2,".r2.filepath"),
    NULL,
    ""
   )
   setup.inserted.div <<- setup.inserted.div[-length(setup.inserted.div)]
  }else{ # Do not remove row if there's only 1 row
   setup.btn <- 2
   setup.value(setup.btn)
  }
 })
 
  react.setup.grp.name <- reactive({
  grp.name <- sapply(grep("setup\\.grp\\d+\\.name", x = names(input), value = TRUE),
                     function(x) input[[x]])
  grp.name=as.vector(grp.name[order(names(grp.name))])
  grp.name[grp.name!=""]
 })
 
 react.setup.grp.ctrl.name <- reactive({
  grp.ctrl.name <- sapply(grep("setup\\.grp\\d+\\.ctrl\\.name", x = names(input), value = TRUE),
                          function(x) input[[x]])
  grp.ctrl.name=as.vector(grp.ctrl.name[order(names(grp.ctrl.name))])
  grp.ctrl.name[grp.ctrl.name!=""]
 })
 
 react.setup.grp.rep.name <- reactive({
  grp.rep.name <- sapply(grep("setup\\.grp\\d+\\.rep\\.name", x = names(input), value = TRUE),
                         function(x) input[[x]])
  # need to do order otherwise the new one will be the first in the vector
  grp.rep.name=as.vector(grp.rep.name[order(names(grp.rep.name))])
  grp.rep.name[grp.rep.name!=""]
 })
 
 react.setup.email <- reactive({
  input$setup.email
 })
 
 react.setup.genome.name <- reactive({
  input$setup.genome.name
 })
 
 react.setup.genome <- reactive({
  input$setup.genome
 })
 
 react.setup.genome.fa <- reactive({
  fa.file = input$setup_genome_fa
  fa.path=as.character(parseFilePaths(root = volumes,fa.file)$datapath)
  fa.path
 })
 
 react.setup.genome.gtf <- reactive({
  gtf.file = input$setup_genome_gtf
  gtf.path=as.character(parseFilePaths(root = volumes,gtf.file)$datapath)
  gtf.path
 })
 
 
 # browse from /filepath if it's specified otherwise /root
 if(file.exists("/filepath")){
  volumes=c(root="/filepath/")
 }else{ 
  volumes=c(root="/root/")
 }
 
 # shinyFileChoose doesn't work inside a reactive
 # I had to arbitrarily set a max number of samples a
 # samples.txt can have
 for (i in 1:max.nsamples) {
  setup.grp.r1.name=paste0('setup_grp',i,'_r1_files')
  shinyFileChoose(input, setup.grp.r1.name, root=volumes,
                  filetypes=c('', 'gz'),session=session)
  setup.grp.r2.name=paste0('setup_grp',i,'_r2_files')
  shinyFileChoose(input, setup.grp.r2.name, root=volumes,
                  filetypes=c('', 'gz'),session=session)
 }
 
 
 # Initialize r1/r2 file list as reactiveValues to save previous choices
 # need to use ad-hoc max.nsamples for shinyFile Choose
 x <- as.list(rep("",max.nsamples))
 names(x)=lapply(1:max.nsamples,function(i)paste0('grp',i))
 setup.grp.r1.file.lst <- do.call("reactiveValues",x)
 setup.grp.r2.file.lst <- do.call("reactiveValues",x)

  # getting all the r1 and r2 fastq files for all groups
  lapply(
  1:max.nsamples,
  function(i){
   observeEvent(input[[paste0('setup_grp',i,'_r1_files')]],{
    # getting new selected path and clean-up paths
    file = input[[paste0('setup_grp',i,'_r1_files')]]
    new.path=as.character(parseFilePaths(root = volumes,file)$datapath)
    new.path=gsub('^/filepath/','',new.path)
    new.path=gsub('^/root/','',new.path)
    # getting path in textAreaInput
    text.path=input[[paste0('setup_grp',i,'_r1_filepaths')]]
    text.path=unlist(strsplit(text.path,"\n"))
    # combine new.path and paths in textarea
    setup.grp.r1.file.lst[[paste0('grp',i)]] <- sort(unique(c(new.path,text.path)))
    path.display=gsub('^/filepath/','',setup.grp.r1.file.lst[[paste0('grp',i)]])
    path.display=gsub('^/root/','',path.display)
    # remove empty path
    if(sum(path.display=="")>0){
         rem.id=which(path.display=="")
         path.display=path.display[-rem.id]
    }
     
    path.display <- paste0(path.display, collapse="\n")
    updateTextAreaInput(session, paste0('setup_grp',i,'_r1_filepaths'),value=path.display)
    # output[[paste0('setup.grp',i,'.r1.filepaths')]]<- renderText({
    #  path.display=gsub('^/filepath/','',setup.grp.r1.file.lst[[paste0('grp',i)]])
    #  path.display=gsub('^/root/','',path.display)
    #  path.display <- paste0(path.display[-1], "\n")
    #  path.display
    # }) # renderText
   })# observeEvent r1_files

   observeEvent(input[[paste0('setup_grp',i,'_r2_files')]],{
    file = input[[paste0('setup_grp',i,'_r2_files')]]
    new.path=as.character(parseFilePaths(root = volumes,file)$datapath)
    new.path=gsub('^/filepath/','',new.path)
    new.path=gsub('^/root/','',new.path)
    text.path=input[[paste0('setup_grp',i,'_r2_filepaths')]]
    text.path=unlist(strsplit(text.path,"\n"))
    # combine new.path and paths in textarea
    setup.grp.r2.file.lst[[paste0('grp',i)]] <- sort(unique(c(new.path,text.path)))
    path.display=gsub('^/filepath/','',setup.grp.r2.file.lst[[paste0('grp',i)]])
    path.display=gsub('^/root/','',path.display)
    if(sum(path.display=="")>0){
     rem.id=which(path.display=="")
     path.display=path.display[-rem.id]
    }
    
    path.display <- paste0(path.display, collapse="\n")
    updateTextAreaInput(session, paste0('setup_grp',i,'_r2_filepaths'),value=path.display)
    # output[[paste0('setup.grp',i,'.r1.filepaths')]]<- renderText({
    #  path.display=gsub('^/filepath/','',setup.grp.r1.file.lst[[paste0('grp',i)]])
    #  path.display=gsub('^/root/','',path.display)
    #  path.display <- paste0(path.display[-1], "\n")
    #  path.display
    # }) # renderText
   })# observeEvent r1_files
    } #function
 )#lapply


    # # # get current content of textareainput
    # setup.grp.r1.file.lst[[paste0('grp',i)]] <- input[[paste0('setup_grp',i,'_r1_filepaths')]]
    # output[[paste0('setup.grp',i,'.r1.filepaths')]]<- renderText({
    #  path.display=gsub('^/filepath/','',setup.grp.r1.file.lst[[paste0('grp',i)]])
    #  path.display=gsub('^/root/','',path.display)
    #  path.display <- paste0(path.display[-1], "\n")
    #  path.display
    # }) # renderText
   
 # select project directory
 shinyDirChoose(input, "setup_proj_dir", root=volumes)
 
react.setup.proj.dir <- reactive({
  proj.dir = input$setup_proj_dir
  proj.path=parseDirPath(roots = volumes,proj.dir)
  proj.path=paste0(proj.path,"/")
 })

output$fileExists <- reactive({
 proj.dir <- react.setup.proj.dir ()
 file_exists=file.exists(file.path(proj.dir,"samples.txt"))
 file_exists
})
outputOptions(output, 'fileExists', suspendWhenHidden=FALSE)

 output$setup.proj.dir.filepaths <- renderText({
  projdir <- react.setup.proj.dir()
  path.display=gsub('^/filepath/','',projdir)
  path.display=gsub('^/root/','',path.display)
  path.display
 })
 
 observeEvent(input$setup.load.samples,{
  projdir <- react.setup.proj.dir()
  # load existing samples treating "NA" as regulating string
  df = read.table(file.path(projdir,"samples.txt"),na.strings=character(0),stringsAsFactors = FALSE)
  updateTextInput(session,inputId="setup.email",value=df[1,5])
  setup.btn <- setup.value()
  # number of rows to add
  add.row <- dim(df)[1]-setup.btn
  for (i in (1:dim(df)[1])){
   # ids for entries in each row
   setup.id.grp.name <- paste0("setup.grp",i,".name")
   setup.id.ctrl.name <- paste0("setup.grp",i,".ctrl.name")
   setup.id.rep.name <- paste0("setup.grp",i,".rep.name")
   setup.id.r1.files <- paste0("setup_grp",i,"_r1_files")
   setup.id.r2.files <- paste0("setup_grp",i,"_r2_files")
   setup.id.r1.filepath <- paste0("setup_grp",i,"_r1_filepaths")
   setup.id.r2.filepath <- paste0("setup_grp",i,"_r2_filepaths")
   # formatting path read from table
   # remove hostfilepath 
   setup.df.path.display.r1=gsub(hostfilepath,'',df[i,6])
   setup.df.path.display.r1 <- gsub(",","\n",setup.df.path.display.r1)
   setup.df.path.display.r2=gsub(hostfilepath,'',df[i,7])
   setup.df.path.display.r2 <- gsub(",","\n",setup.df.path.display.r2)
   # add new rows if needed
  if(i > (dim(df)[1]-add.row)){
   # update the reactiveVal
   setup.value(i)
  # id for new div
  setup.div.id <- paste0("setup_div_",i)

  if(i %% 2 == 0) row.color="#FFFFFF"
  if(i %% 2 != 0) row.color="#f9f9f9"
  insertUI(
   selector = "#setup-placeholder",
   ui = tags$div(
    tags$style(HTML(paste0("#setup_grp",i,"_r1_filepaths { width: 300px; overflow-x: auto; white-space: pre;}"))),
    tags$style(HTML(paste0("#setup_grp",i,"_r2_filepaths { width: 300px; overflow-x: auto; white-space: pre;}"))),
    fluidRow(style = paste0("background-color:",row.color),
             column(2,
                    column(2,
                           tags$label(paste0(i,"."), style = "padding-top: 30px;"),
                           style = "padding-right: 0px;padding-left:0px"),
                    column(10,textInput(inputId=setup.id.grp.name,
                                        label = "",value=df[i,1]),
                           style = "padding-right: 0px;padding-left:0px")
             ),
             column(2,textInput(inputId=setup.id.ctrl.name,
                                label = "",value=df[i,2] )),
             column(2,textInput(inputId=setup.id.rep.name,
                                label = "",value=df[i,3] )),
             br(),
             column(3,shinyFilesButton(setup.id.r1.files,
                                       label = "Select R1 fastq files",
                                       title=paste("Please select Group",i,"R1 fastq files"),
                                       multiple=TRUE ),
                    textAreaInput(setup.id.r1.filepath,label="",width='150px',height='100px',
                                  value=setup.df.path.display.r1)),
             column(3,shinyFilesButton(setup.id.r2.files,
                                       label = "Select R2 fastq files",
                                       title=paste("Please select Group",i,"R1 fastq files"),
                                       multiple=TRUE),
                    textAreaInput(setup.id.r2.filepath,label="",width='150px',height='100px',
                                  value=setup.df.path.display.r2)),
             
    ),#fluidrow
    id = setup.div.id
   )
  )#insertUI
  setup.inserted.div <<- c(setup.inserted.div, setup.div.id)
  } # if(i > (dim(df)[1]-add.row))
   # fill out rows
   updateTextInput(session,setup.id.grp.name,value = df[i,1])
   updateTextInput(session,setup.id.ctrl.name,value = df[i,2])
   updateTextInput(session,setup.id.rep.name,value = df[i,3])
   updateTextAreaInput(session,setup.id.r1.filepath,value = setup.df.path.display.r1)
   updateTextAreaInput(session,setup.id.r2.filepath,value = setup.df.path.display.r2)
  }# for (i in (1:dim(df)[1]))
 })
 
 # select other genome files
 shinyFileChoose(input, "setup_genome_fa", root=volumes,
                 filetypes=c('', 'gz','fasta','fa'))
 shinyFileChoose(input, "setup_genome_gtf", root=volumes,
                 filetypes=c('', 'gz','GTF','gtf'))
 
 observeEvent(input$setup_genome_fa,{
  fa.file = input$setup_genome_fa
  fa.path=as.character(parseFilePaths(root = volumes,fa.file)$datapath)
  output$setup.genome.fa.filepaths <- renderText({
   path.display=gsub('^/filepath/','',fa.path)
   path.display=gsub('^/root/','',path.display)
   path.display
  })
 })
 
 observeEvent(input$setup_genome_gtf,{
  file = input$setup_genome_gtf
  path=as.character(parseFilePaths(root = volumes,file)$datapath)
  output$setup.genome.gtf.filepaths <- renderText({
   path.display=gsub('^/filepath/','',path)
   path.display=gsub('^/root/','',path.display)
   path.display
  })
 })
 
 output$err_msg_genome <- renderText({
  genome.fa = react.setup.genome.fa()
  genome.name=react.setup.genome.name()
  genome=react.setup.genome()
  genome.gtf=react.setup.genome.gtf()
  if(genome=="other"){
  shiny::validate(
   shiny::need(genome.name != '', "You must enter a genome name"),
   shiny::need(genome.fa != '', 'You must select genome fasta file'),
   shiny::need(genome.gtf != '', "You must select a GTF file")
  )
  }
 })
 
 output$err_msg <- renderText({
  projdir=react.setup.proj.dir()
  grpname=react.setup.grp.name()
  ctrlname=react.setup.grp.ctrl.name()
  repname=react.setup.grp.rep.name()
  nsamples=setup.value()
 shiny::validate(
  shiny::need(projdir != '', 'You must select a project folder'),
  shiny::need(length(grpname) == nsamples, 'You must enter group name'),
  shiny::need(length(ctrlname) == nsamples, 'You must enter control name'),
  shiny::need(length(repname)==nsamples,'You must enter replicate name')
 )
 })
 
 observeEvent(input$setup.run.analysis,{
  # create samples.txt
   projdir=react.setup.proj.dir()
   grpname=react.setup.grp.name()
   ctrlname=react.setup.grp.ctrl.name()
   repname=react.setup.grp.rep.name()
   spikename=rep("NA",length(grpname))
   email=react.setup.email()
   nsamples=setup.value()
   genome=react.setup.genome()
   genome.name = react.setup.genome.name()
   genome.fa=react.setup.genome.fa()
   genome.gtf=react.setup.genome.gtf()
   if(email==""){
    email=rep("NA",length(grpname))
   }else{
    email=rep(email,length(grpname))
   }
   shiny::req(projdir!='')
   shiny::req(length(grpname)==nsamples)
   shiny::req(length(ctrlname)==nsamples)
   shiny::req(length(repname)==nsamples)
   
   # getting the r1 and r2 fastq for all groups
   
   # function to remove empty paths
   rem_empty_path <- function(input.name) {
    filt.path = input[[input.name]]
    filt.path = unlist(strsplit(filt.path,"\n"))
    print("filt.path")
    print(filt.path)
    empty.id = which(filt.path=="")
    print("empty.id")
    print(empty.id)
    if(length(empty.id)>0) filt.path = filt.path[-empty.id]
    print("filt.path after rem empty")
    print(filt.path)
   }
   print("input[[setup_grp1_r1_filepaths)]]")
   print(input[[paste0("setup_grp",1,"_r1_filepaths")]])
   r1_fastq_lst=lapply(1:length(grpname),function(x) rem_empty_path(paste0("setup_grp",x,"_r1_filepaths")))
   r2_fastq_lst=lapply(1:length(grpname),function(x) rem_empty_path(paste0("setup_grp",x,"_r2_filepaths")))
   print("r1_fastq_lst")
   print(r1_fastq_lst)
   # sort and unique to make sure the correct pairing of r1 and r2
   r1_fastq_lst=lapply(1:length(grpname),
          function(x) sort(unique(r1_fastq_lst[[x]])))
   r2_fastq_lst=lapply(1:length(grpname),
           function(x) sort(unique(r2_fastq_lst[[x]])))
   # changing to hostpath
   r1_fastq=lapply(1:length(grpname),
                   function(x) paste(paste0(hostfilepath,r1_fastq_lst[[x]]),collapse=","))
   r2_fastq=lapply(1:length(grpname),
                   function(x) paste(paste0(hostfilepath,r2_fastq_lst[[x]]),collapse=","))
   print("r1_fastq")
   print(r1_fastq)
   df <- data.frame(
    grpname=grpname,
    ctrlname=ctrlname,
    repname=repname,
    spikename=spikename,
    email=email,
    path_to_r1_fastq=unlist(r1_fastq),
    path_to_r2_fastq=unlist(r2_fastq))
   write.table(paste0("#Groupname\tControlname\tReplicatename\tspikename\temail",
                      "\tpath_to_r1_fastq\tpath_to_r2_fastq"),
               file=file.path(projdir,"samples.txt"),quote=FALSE,row.names=FALSE,col.names=FALSE)
   write.table(df,file=file.path(projdir,"samples.txt"),quote=FALSE,row.names=FALSE,col.names=FALSE,
               append = TRUE)
  #system("echo 'sbatch --help' > /hostpipe")
  # getting genome options
  if(genome=="hg38") options="genome=hg38"
  if(genome=="hg19") options="genome=hg19"
  if(genome=="other"){
   shiny::req(genome.fa != "")
   shiny::req(genome.gtf!="")
   shiny::req(genome!="")
  }
   fa.file=react.setup.genome.fa() 
   gtf.file=react.setup.genome.gtf ()
   # convert path to hostpath
   if(file.exists("/filepath")){
    hostfa=gsub('/filepath/',hostfilepath,fa.file)
    hostgtf=gsub('/filepath/',hostfilepath,gtf.file)
   }else{ 
    hostfa=gsub('/root/','/',fa.file)
    hostgtf=gsub('/root/','/',gtf.file)
   }
   if(input$setup.genome.name=="other"){
   options=paste0("genome=",input$setup.genome.name," ref_fa=",hostfa," ref_gtf=",hostgtf)
   }else{options=paste0("genome=",input$setup.genome)}
    
  if(input$setup.batch.adjust==TRUE){ batch_adjust="yes"
  }else{ batch_adjust="no"}
   
   if(input$setup.debug==TRUE){ run="debug"
   }else{ run=""}
  options=paste0(options," time=",input$setup.time," batch_adjust=",batch_adjust,
                 " ncpus_trim=",input$setup.ncpus.trim," ncpus_star=",input$setup.ncpus.star,
                 " run=",run)
  # converting project dir to host path
  if(file.exists("/filepath")){
   hostprojdir=gsub('/filepath/',hostfilepath,projdir)
  }else{ 
   hostprojdir=gsub('/root/','/',projdir)
  }
   print(paste0("echo 'cd ",projdir,"&& bash ",img.dir,"/scripts/run_rnaseq_full.sh ",options,
                " &> run_rnaseq_full.out' > /hostpipe"))
   #system(paste0("echo 'cd ",hostprojdir,"&& bash ",img.dir,"/scripts/run_rnaseq_full.sh ",options,
   #             " &> run_rnaseq_full.out' > /hostpipe"))
})
 
 #### log tab #####
 # update log file list
 observeEvent(input$logtab.refresh.log.path,{
  # list of files in log directory from most recent
  projdir <- react.setup.proj.dir()
  files=list.files(paste0(projdir,"outputs/logs"),pattern = "\\.txt$|\\.out$",full.names = TRUE)
  # Sort files, placing filenames starting with "run_" at the top
  sorted_files <- c(sort(files[startsWith(files, "run_")]), sort(files[!startsWith(files, "run_")]))
  log.lst = basename(sorted_files)
  if(file.exists(paste0(projdir,"outputs/logs/run_rnaseq_full.out"))){
   updateSelectInput(session,"logtab.log.path",choices=log.lst,selected = "run_rnaseq_full.out")
  }else{
   updateSelectInput(session,"logtab.log.path",choices=log.lst,selected="")
  }
     # display file content in UI
   output$logtab.log.content <- renderPrint({ 
    shiny::req(input$logtab.log.path!="")
    cat(react.logtab.log.content(),sep="\n")
   })
  # list only failed *.out file
  files=list.files(paste0(projdir,"outputs/logs"),pattern = "\\.out$",full.names = TRUE)
  # Filter files containing the word "FAILED"
  failed.files <- lapply(files, function(file) {
   if (any(grepl("FAILED", tail(readLines(file),2), ignore.case = FALSE))) {
    return(file)
   }
  })
  # Filter out NULL values (files without the word "FAILED")
  failed.files <- unlist(Filter(Negate(is.null), failed.files))
  if(is.null(failed.files)){
   log.fail.lst=""
  }else{
   # get the first FAILED
   file_info <- file.info(failed.files)
   sorted_files <- failed.files[order(as.numeric(file_info$mtime),decreasing=FALSE)]
   log.fail.lst = basename(sorted_files)
  }
  updateSelectInput(session,"logtab.fail.log.path",choices=log.fail.lst)
 
  })#observEvent refresh
 
 react.logtab.log.content <- reactive({
    projdir <- react.setup.proj.dir()
    path <- paste0(projdir,"outputs/logs/",input$logtab.log.path)
  content <- readLines(path)
  return(content)
 })
 
 react.logtab.fail.log.content <- reactive({
  projdir <- react.setup.proj.dir()
  if(input$logtab.fail.log.path=="NA") return("Nothing failed!")
  path <- paste0(projdir,"outputs/logs",input$logtab.fail.log.path)
  content <- readLines(path)
  return(content)
 })
 
 # display file content in UI
 output$logtab.fail.log.content <- renderPrint({
  shiny::req(input$logtab.fail.log.path!="")
  react.logtab.fail.log.content()
 })
 
 #### QC tab #####
 output$multiqc <- renderUI({
  projdir <- react.setup.proj.dir()
  # path for fastqc report html to be referred to in iframe
  addResourcePath("fastqc_rslt", paste0(projdir,"outputs/fastqc_rslt"))
  tags$iframe(seamless="seamless",
              src="fastqc_rslt/multiqc_report.html",
              width="100%",
              height="1000")
 })
 
 ## Differential genes analysis report ###
 output$diff_report <- renderUI({
  projdir <- react.setup.proj.dir()
  # path for differential report html to be referred to in iframe
  addResourcePath("diff_report_path",paste0(projdir,"outputs/diff_analysis_rslt/"))
  tags$iframe(seamless="seamless", 
              src= "diff_report_path/RNA-seq_differential_analysis_report.html",
              width="100%", 
              height=800)
 })
 
 #### processing plots tab #####
 
 react.tab3.data.loc <- reactive({
  projdir <- react.setup.proj.dir()
  data.loc = paste0(projdir,"outputs/")
 })
 
 react.tab3.out.loc <- reactive({
  shiny::req(dir.exists(data.loc))
  # create directory for shiny outputs
  dir.create(paste0(data.loc, "shiny_outputs"))
  out.loc=paste0(data.loc,"shiny_outputs/")
 })
 
 # initialize values for tab
react.tab3.rdata <- reactive({
 selected_tab <- input$inTabset
 projdir <- react.setup.proj.dir()
 data.loc <- react.tab3.data.loc()
 rdata.loc1=paste0(data.loc,"diff_analysis_rslt/RNA-seq_differential_analysis.RData")
 rdata.loc2=paste0(data.loc,"diff_analysis_rslt/RNA-seq differential analysis.RData")
 shiny::req(file.exists(rdata.loc1)||file.exists(rdata.loc2))
 # filename was changed. Make sure can read both
 if(file.exists(rdata.loc1)){
  load(rdata.loc1)
   } else if (file.exists(rdata.loc2)){
  load(rdata.loc2)
   }
 return(out.DESeq2)
})

 react.tab3.groups.lst <- reactive({
  out.DESeq2 <- react.tab3.rdata()
  groups= data.frame(
  var = names(out.DESeq2$results),
  names = names(out.DESeq2$results)
 )
 # List of choices for groupnames
 groups.lst <- as.list(groups$names)
 # Name it
 names(groups.lst) <- groups$var
 groups.lst
})

 observe({
  groups.lst <- react.tab3.groups.lst()
  updateSelectInput(session, "tab3.grp1.name",choices=groups.lst)
  updateTextInput(session, "tab3.grp1.plot.title",value=groups.lst[[1]])
 })

 # vals will contain all plots and table grobs
 vals.plot <- reactiveValues(venn.up1=NULL,venn.up2=NULL,venn.dwn1=NULL,
                             venn.dwn2=NULL,hm.up=NULL,hm.dwn=NULL,
                             upset.up=NULL,upset.dwn=NULL,msig.mf=NULL,
                             msig.bp=NULL,msig.cc=NULL,msig.curate=NULL)
 
 # adding groups
 value <- reactiveVal(1)
 inserted <- c()		
 inserted.col <- c()
 inserted.fdr <- c()
 inserted.fc <- c()
 inserted.meanDiff <- c()
 # reactive value to plot gene enrichment
 react.val.gen.go <- reactiveVal(0)
 # adjust reactive value as gen.go button is pressed
 observeEvent(input$gen.go, {
  gen.go <- react.val.gen.go()
  gen.go <- react.val.gen.go(1)
 })
 
 observeEvent(input$insert_set, {
  groups.lst <- react.tab3.groups.lst()
  btn <- value() +1
  value(btn)
  id <- paste0("txt1_", btn)
  insertUI(
   selector = "#placeholder",
   ui = tags$div(
    selectInput(inputId=paste0("tab3.grp",btn,".name"),
                label=(paste0("Group ",btn)), 
                selected=groups.lst[btn],
                choices=groups.lst),
    textInput(inputId=paste0("tab3.grp",btn,".plot.title"),
              label=(paste0("Group ",btn," label")), 
              value = groups.lst[btn]),
    id = id
   )
  )
  inserted <<- c(inserted, id)
  
  id.col <- paste0("txt2_",btn)
  insertUI(
   selector = "#placeholder-col",
   ui = tags$div(
    textInput(inputId=paste0("tab3.color.grp",btn), 
              label = (paste0("Color for group ",btn)), 
              value = tab3.colors[btn]),
    id = id.col
   )
  )
  inserted.col <<- c(inserted.col, id.col)
  
  id.fdr <- paste0("txt3_",btn)
  insertUI(
   selector = "#placeholder-fdr",
   ui = tags$div(
    numericInput(inputId=paste0("tab3.fdr.grp",btn), 
                 label = (paste0("FDR cut-off for group ",btn)), 
                 value = 0.05,step=0.01,min=0,max=1),
    id = id.fdr
   )
  )
  inserted.fdr <<- c(inserted.fdr, id.fdr)
  
  id.fc <- paste0("txt4_",btn)
  insertUI(
   selector = "#placeholder-fc",
   ui = tags$div(
    numericInput(inputId=paste0("tab3.fc.grp",btn), 
                 label = (paste0("Fold-change for group ",btn)), 
                 value = 1,min=1,step=0.5),
    id = id.fc
   )
  )
  inserted.fc <<- c(inserted.fc, id.fc)
  
  id.meanDiff <- paste0("txt5_",btn)
  insertUI(
   selector = "#placeholder-meanDiff",
   ui = tags$div(
    numericInput(inputId=paste0("tab3.meanDiff.grp",btn), 
                 label = (paste0("Difference for group ",btn)), 
                 value = 0,min=0,step=1),
    id = id.meanDiff
   )
  )
  inserted.meanDiff <<- c(inserted.meanDiff, id.meanDiff)
 })
 
 observeEvent(input$remove_set, {
  if(value()>1){
   btn <- value() -1
   value(btn)
   removeUI(
    selector = paste0("#",inserted[length(inserted)])
   )
   removeUI(
    selector = paste0("#",inserted.col[length(inserted.col)])
   )
   removeUI(
    selector = paste0("#",inserted.fdr[length(inserted.fdr)])
   )
   removeUI(
    selector = paste0("#",inserted.fc[length(inserted.fc)])
   )
   removeUI(
    selector = paste0("#",inserted.meanDiff[length(inserted.meanDiff)])
   )
   updateSelectInput(
    session,
    paste0("tab3.grp",length(inserted)+1 ,".name"),
    NULL,
    NA
   )
   updateTextInput(
    session,
    paste0("tab3.grp", length(inserted)+1,".plot.title"),
    NULL,
    NA
   )
   updateTextInput(
    session,
    paste0("tab3.color.grp",length(inserted.col)+1),
    NULL,
    NA
   )
   updateNumericInput(
    session,
    paste0("tab3.fdr.grp",length(inserted.fdr)+1),
    NULL,
    NA
   )
   updateNumericInput(
    session,
    paste0("tab3.fc.grp",length(inserted.fc)+1),
    NULL,
    NA
   )
   updateNumericInput(
    session,
    paste0("tab3.meanDiff.grp",length(inserted.meanDiff)+1),
    NULL,
    NA
   )
   
   inserted <<- inserted[-length(inserted)]
   inserted.col <<- inserted.col[-length(inserted.col)]
   inserted.fdr <<- inserted.fdr[-length(inserted.fdr)]
   inserted.fc <<- inserted.fc[-length(inserted.fc)]
   inserted.meanDiff <<- inserted.meanDiff[-length(inserted.meanDiff)]
  }else{ btn <- 1
  value(btn)}
 })
 
 #### processing three groups #######
 react.grp.name <- reactive({
  c(input$grp1.name,input$grp2.name,input$grp3.name)
 })
 react.grp.plot.title <- reactive({
  c(input$grp1.plot.title,input$grp2.plot.title,input$grp3.plot.title)
 })
 
 react.fc.cutoff <- reactive({
  c(input$fc.grp1,input$fc.grp2,input$fc.grp3)
 })
 
 react.venn.cex <- reactive({
  input$venn.cex
 })	
 
 react.venn.pad <- reactive({
  input$venn.pad
 })		
 
 react.fdr.cutoff <- reactive({
  c(input$fdr.grp1,input$fdr.grp2,input$fdr.grp3)
 })
 
 react.venn.opts <- reactive({input$venn.opts})
 
 react.result <- reactive({
  grp.name=react.grp.name()
  grp.plot.title=react.grp.plot.title()
  idx=match(grp.name,names(out.DESeq2$results))
  results=out.DESeq2$results[idx]
 })
 
 gs.RNASeq <- reactive({
  grp.name=react.grp.name()
  grp.plot.title=react.grp.plot.title()
  results=react.result()
  gs.RNASeq <- max(sapply(results,function(x) sum(!is.na(x$padj))))
 })
 
 grp.up <- reactive ({
  grp.name=react.grp.name()
  grp.plot.title=react.grp.plot.title()
  results=react.result()
  fdr.co <- react.fdr.cutoff()
  fc.cutoff <- react.fc.cutoff()
  # up-regulated genes
  grp.up=list()
  for(i in 1:ngroup){
   grp.up[[i]]=results[[grp.name[i]]][results[[grp.name[i]]]$padj<=fdr.co[i] &
                                       results[[grp.name[i]]]$log2FoldChange>(log2(fc.cutoff[i])) & 
                                       !is.na(results[[grp.name[i]]]$padj),]
  }
  grp.up
 })
 
 genes.lists.up <- reactive ({
  grp.up <- grp.up()
  grp.plot.title=react.grp.plot.title()
  genes.lists.up=lapply(grp.up,function(x) rownames(x))
  names(genes.lists.up)<-grp.plot.title
  genes.lists.up
 })
 
 gom.obj.up <- reactive({
  genes.lists.up <- genes.lists.up()
  gs.RNASeq <- gs.RNASeq()
  newGOM(genes.lists.up,genes.lists.up,genome.size=gs.RNASeq)
 })
 
 s4.up <- reactive({
  # plot euler venn
  s4<-genes.lists.up()
  names(s4)<-react.grp.plot.title()
  s4
 })
 
 grp.dwn <- reactive({
  grp.name=react.grp.name()
  grp.plot.title=react.grp.plot.title()
  results=react.result()
  fdr.co <- react.fdr.cutoff()
  fc.cutoff <- react.fc.cutoff()
  
  # down-regulated genes
  grp.dwn=list()
  for(i in 1:ngroup){
   grp.dwn[[i]]=results[[grp.name[i]]][results[[grp.name[i]]]$padj<=fdr.co[i] &
                                        results[[grp.name[i]]]$log2FoldChange<(-log2(fc.cutoff[i])) & 
                                        !is.na(results[[grp.name[i]]]$padj),]
  }
  grp.dwn
 })
 
 genes.lists.dwn <- reactive({
  grp.dwn <- grp.dwn()
  grp.plot.title=react.grp.plot.title()
  genes.lists.dwn <-
   lapply(grp.dwn,function(x) rownames(x))
  names(genes.lists.dwn)<-grp.plot.title
  genes.lists.dwn
 })
 
 gom.obj.dwn <- reactive({
  gs.RNASeq<-gs.RNASeq()
  genes.lists.dwn<-genes.lists.dwn()
  names(genes.lists.dwn)<-react.grp.plot.title()
  gom.obj.dwn <- newGOM(genes.lists.dwn,genes.lists.dwn,genome.size=gs.RNASeq)
  gom.obj.dwn
 })
 
 genes.lists <- reactive ({
  grp.name=react.grp.name()
  grp.plot.title=react.grp.plot.title()
  results=react.result()
  fdr.co <- react.fdr.cutoff()
  fc.cutoff <- react.fc.cutoff()
  grp.dwn <- grp.dwn ()
  grp.up <- grp.up ()
  # all regulated-genes
  grp=list()
  for(i in 1:ngroup){
   grp[[i]]=rbind(grp.up[[i]],grp.dwn[[i]])
  }
  
  genes.lists=lapply(grp,function(x) rownames(x))
  names(genes.lists)<-grp.plot.title
  genes.lists
 })
 
 gom.obj <- reactive({
  genes.lists <- genes.lists()
  gs.RNASeq <- gs.RNASeq()
  newGOM(genes.lists,genes.lists,genome.size=gs.RNASeq)
 })
 
 s4 <- reactive({
  # plot euler venn
  s4<-genes.lists()
  names(s4)<-react.grp.plot.title()
  s4
 })
 
 react.colors <- reactive({
  colors[1]=input$color.grp1
  colors[2]=input$color.grp2
  colors[3]=input$color.grp3
  colors
 })
 s4.dwn <- reactive({
  # plot euler venn
  s4<-genes.lists.dwn()
  names(s4)<-react.grp.plot.title()
  s4
 })
 
 data.up <- reactive({
  gom.obj.up <- gom.obj.up()
  grp.plot.title <- react.grp.plot.title()
  # calculate overlap
  data <- getMatrix(gom.obj.up,"pval")
  data <- melt(data,id=grp.plot.title)
  data <- data[c(2:3,6),]
  names(data)[names(data)=="value"]="pval"
  tmp <- getMatrix(gom.obj.up,"Jaccard")
  tmp <- melt(tmp,id=grp.plot.title)
  tmp <- tmp[c(2,3,6),]
  data$jaccard=tmp$value
  data$pval.cat=ifelse(data$pval==0,pval.names[1],ifelse(data$pval<=0.05,pval.names[2],
                                                         ifelse(data$pval<0.5,pval.names[3],pval.names[4])))
  data$pval.cat=factor(data$pval.cat,levels=pval.names)
  data
 })
 
 data.dwn <- reactive({
  gom.obj.dwn <- gom.obj.dwn()
  grp.plot.title <- react.grp.plot.title()
  ## heatmap down regulation
  data <- getMatrix(gom.obj.dwn,"pval")
  data <- melt(data,id=grp.plot.title)
  data <- data[c(2,3,6),]
  names(data)[names(data)=="value"]="pval"
  tmp <- getMatrix(gom.obj.dwn,"Jaccard")
  tmp <- melt(tmp,id=grp.plot.title)
  tmp <- tmp[c(2,3,6),]
  data$jaccard=tmp$value
  data$pval.cat=ifelse(data$pval==0,pval.names[1],ifelse(data$pval<=0.05,pval.names[2],
                                                         ifelse(data$pval<0.5,pval.names[3],pval.names[4])))
  data$pval.cat=factor(data$pval.cat,levels=pval.names)
  data
 })
 
 data <- reactive({
  gom.obj <- gom.obj()
  grp.plot.title <- react.grp.plot.title()
  # calculate overlap
  data <- getMatrix(gom.obj,"pval")
  data <- melt(data,id=grp.plot.title)
  data <- data[c(2:3,6),]
  names(data)[names(data)=="value"]="pval"
  tmp <- getMatrix(gom.obj,"Jaccard")
  tmp <- melt(tmp,id=grp.plot.title)
  tmp <- tmp[c(2,3,6),]
  data$jaccard=tmp$value
  data$pval.cat=ifelse(data$pval==0,pval.names[1],ifelse(data$pval<=0.05,pval.names[2],
                                                         ifelse(data$pval<0.5,pval.names[3],pval.names[4])))
  data$pval.cat=factor(data$pval.cat,levels=pval.names)
  data
 })
 
 output$plot1 <- renderPlot({
  venn.opts=react.venn.opts()
  colors=react.colors()
  cex = react.venn.cex()
  pad = react.venn.pad()
  set.seed(1)
  s4 <- s4.up()
  venn.up <- plot_euler(s4,colors,cex,venn.opts,title="Up-regulation")
  s4 <- s4.dwn()
  venn.dwn <- plot_euler(s4,colors,cex,venn.opts,title="Down-regulation")
  grid.arrange(venn.up,venn.dwn,ncol=2,padding=unit(pad,"line"),
               bottom="",right="",left="",top="")
 })
 
 output$plot2 <- renderPlot({
   pad = react.venn.pad()
  data=data.up()
  hm.up <- heatmap_sigf_overlap(data)
  data=data.dwn()
  hm.dwn <- heatmap_sigf_overlap(data)
  grid.arrange(hm.up,hm.dwn,ncol=2,padding=unit(pad,"line"),
               bottom="",right="",left="",top="")
 })
 
 output$plot3 <- renderPlot({	
  venn.opts=react.venn.opts()
  colors=react.colors()
  cex = react.venn.cex()
  set.seed(1)
  s4 <- s4()
  plot_euler(s4,colors,cex,venn.opts,title="All regulated genes")
 })
 
 output$plot4 <- renderPlot({
  data=data()
  heatmap_sigf_overlap(data)
 })
 
 ######## processing two groups ############
 
 react.tab0.grp.name <- reactive({
  c(input$tab0.grp1.name,input$tab0.grp2.name)
 })
 react.tab0.grp.plot.title <- reactive({
  c(input$tab0.grp1.plot.title,input$tab0.grp2.plot.title)
 })
 
 react.tab0.venn.cex	<- reactive({
  input$tab0.venn.cex
 })
 
 react.tab0.venn.pad	<- reactive({
  input$tab0.venn.pad
 })
 
 react.tab0.fc.cutoff <- reactive({
  c(input$tab0.fc.grp1,input$tab0.fc.grp2)
 })
 
 react.tab0.fdr.cutoff <- reactive({
  c(input$tab0.fdr.grp1,input$tab0.fdr.grp2)
 })
 
 react.tab0.venn.opts <- reactive({input$tab0.venn.opts})
 
 react.tab0.result <- reactive({
  tab0.grp.name=react.tab0.grp.name()
  tab0.grp.plot.title=react.tab0.grp.plot.title()
  idx=match(tab0.grp.name,names(out.DESeq2$results))
  tab0.results=out.DESeq2$results[idx]
 })
 
 tab0.gs.RNASeq <- reactive({
  tab0.results=react.tab0.result()
  gs.RNASeq <- max(sapply(tab0.results,function(x) sum(!is.na(x$padj))))
 })
 
 tab0.grp.up <- reactive ({
  grp.name=react.tab0.grp.name()
  grp.plot.title=react.tab0.grp.plot.title()
  results=react.tab0.result()
  fdr.co <- react.tab0.fdr.cutoff()
  fc.cutoff <- react.tab0.fc.cutoff()
  # up-regulated genes
  grp.up=list()
  for(i in 1:tab0.ngroup){
   grp.up[[i]]=results[[grp.name[i]]][results[[grp.name[i]]]$padj<=fdr.co[i] &
                                       results[[grp.name[i]]]$log2FoldChange>(log2(fc.cutoff[i])) & 
                                       !is.na(results[[grp.name[i]]]$padj),]
  }
  grp.up
 })
 
 tab0.genes.lists.up <- reactive ({
  grp.up <- tab0.grp.up()
  grp.plot.title=react.tab0.grp.plot.title()
  genes.lists.up=lapply(grp.up,function(x) rownames(x))
  names(genes.lists.up)<-grp.plot.title
  genes.lists.up
 })
 
 tab0.gom.obj.up <- reactive({
  genes.lists.up <- tab0.genes.lists.up()
  gs.RNASeq <- tab0.gs.RNASeq()
  newGOM(genes.lists.up,genes.lists.up,genome.size=gs.RNASeq)
 })
 
 tab0.s4.up <- reactive({
  # plot euler venn
  s4<-tab0.genes.lists.up()
  names(s4)<-react.tab0.grp.plot.title()
  s4
 })
 
 tab0.grp.dwn <- reactive({
  grp.name=react.tab0.grp.name()
  grp.plot.title=react.tab0.grp.plot.title()
  results=react.tab0.result()
  fdr.co <- react.tab0.fdr.cutoff()
  fc.cutoff <- react.tab0.fc.cutoff()
  
  # down-regulated genes
  grp.dwn=list()
  for(i in 1:tab0.ngroup){
   grp.dwn[[i]]=results[[grp.name[i]]][results[[grp.name[i]]]$padj<=fdr.co[i] &
                                        results[[grp.name[i]]]$log2FoldChange<(-log2(fc.cutoff[i])) & 
                                        !is.na(results[[grp.name[i]]]$padj),]
  }
  grp.dwn
 })
 
 tab0.genes.lists.dwn <- reactive({
  grp.dwn <- tab0.grp.dwn()
  grp.plot.title=react.tab0.grp.plot.title()
  genes.lists.dwn <-
   lapply(grp.dwn,function(x) rownames(x))
  names(genes.lists.dwn)<-grp.plot.title
  genes.lists.dwn
 })
 
 tab0.gom.obj.dwn <- reactive({
  gs.RNASeq<-tab0.gs.RNASeq()
  genes.lists.dwn<-tab0.genes.lists.dwn()
  names(genes.lists.dwn)<-react.tab0.grp.plot.title()
  gom.obj.dwn <- newGOM(genes.lists.dwn,genes.lists.dwn,
                        genome.size=gs.RNASeq)
  gom.obj.dwn
 })
 
 tab0.genes.lists <- reactive ({
  tab0.grp.name=react.tab0.grp.name()
  tab0.grp.plot.title=react.tab0.grp.plot.title()
  tab0.results=react.tab0.result()
  tab0.fdr.co <- react.tab0.fdr.cutoff()
  tab0.fc.cutoff <- react.tab0.fc.cutoff()
  tab0.grp.dwn <- tab0.grp.dwn ()
  tab0.grp.up <- tab0.grp.up ()
  # all regulated-genes
  tab0.grp=list()
  for(i in 1:tab0.ngroup){
   tab0.grp[[i]]=rbind(tab0.grp.up[[i]],tab0.grp.dwn[[i]])
  }
  
  tab0.genes.lists=lapply(tab0.grp,function(x) rownames(x))
  names(tab0.genes.lists)<-tab0.grp.plot.title
  tab0.genes.lists
 })
 
 tab0.gom.obj <- reactive({
  tab0.genes.lists <- tab0.genes.lists()
  tab0.gs.RNASeq <- tab0.gs.RNASeq()
  newGOM(tab0.genes.lists,tab0.genes.lists,genome.size=tab0.gs.RNASeq)
 })
 
 tab0.s4 <- reactive({
  # plot euler venn
  tab0.s4<-tab0.genes.lists()
  names(tab0.s4)<-react.tab0.grp.plot.title()
  tab0.s4
 })
 
 react.tab0.colors <- reactive({
  tab0.colors[1]=input$tab0.color.grp1
  tab0.colors[2]=input$tab0.color.grp2
  tab0.colors
 })
 tab0.s4.dwn <- reactive({
  # plot euler venn
  tab0.s4<-tab0.genes.lists.dwn()
  names(tab0.s4)<-react.tab0.grp.plot.title()
  tab0.s4
 })
 
 tab0.data.up <- reactive({
  gom.obj.up <- tab0.gom.obj.up()
  grp.plot.title <- react.tab0.grp.plot.title()
  # calculate overlap
  data <- getMatrix(gom.obj.up,"pval")
  data <- melt(data,id=grp.plot.title)
  data <- data[2,]
  names(data)[names(data)=="value"]="pval"
  tmp <- getMatrix(gom.obj.up,"Jaccard")
  tmp <- melt(tmp,id=grp.plot.title)
  tmp <- tmp[2,]
  data$jaccard=tmp$value
  data$pval.cat=ifelse(data$pval==0,pval.names[1],ifelse(data$pval<=0.05,pval.names[2],
                                                         ifelse(data$pval<0.5,pval.names[3],pval.names[4])))
  data$pval.cat=factor(data$pval.cat,levels=pval.names)
  data
 })
 
 tab0.data.dwn <- reactive({
  gom.obj.dwn <- tab0.gom.obj.dwn()
  grp.plot.title <- react.tab0.grp.plot.title()
  ## heatmap down regulation
  data <- getMatrix(gom.obj.dwn,"pval")
  data <- melt(data,id=grp.plot.title)
  data <- data[2,]
  names(data)[names(data)=="value"]="pval"
  tmp <- getMatrix(gom.obj.dwn,"Jaccard")
  tmp <- melt(tmp,id=grp.plot.title)
  tmp <- tmp[2,]
  data$jaccard=tmp$value
  data$pval.cat=ifelse(data$pval==0,pval.names[1],ifelse(data$pval<=0.05,pval.names[2],
                                                         ifelse(data$pval<0.5,pval.names[3],pval.names[4])))
  data$pval.cat=factor(data$pval.cat,levels=pval.names)
  data
 })
 
 tab0.data <- reactive({
  gom.obj <- tab0.gom.obj()
  grp.plot.title <- react.tab0.grp.plot.title()
  # calculate overlap
  data <- getMatrix(gom.obj,"pval")
  data <- melt(data,id=grp.plot.title)
  data <- data[2,]
  names(data)[names(data)=="value"]="pval"
  tmp <- getMatrix(gom.obj,"Jaccard")
  tmp <- melt(tmp,id=grp.plot.title)
  tmp <- tmp[2,]
  data$jaccard=tmp$value
  data$pval.cat=ifelse(data$pval==0,pval.names[1],ifelse(data$pval<=0.05,pval.names[2],
                                                         ifelse(data$pval<0.5,pval.names[3],pval.names[4])))
  data$pval.cat=factor(data$pval.cat,levels=pval.names)
  data
 })
 
 output$tab0.plot1 <- renderPlot({
  venn.opts=react.tab0.venn.opts()
  colors=react.tab0.colors()
  cex=react.tab0.venn.cex()
  pad=react.tab0.venn.pad()
  set.seed(1)
  s4 <- tab0.s4.up()
  venn.up <- plot_euler(s4=s4,colors=colors,cex=cex,venn.opts=venn.opts,title="Up-regulation")
  s4 <- tab0.s4.dwn()
  venn.dwn <- plot_euler(s4=s4,colors=colors,cex=cex,venn.opts=venn.opts,title="Down-regulation")
  grid.arrange(venn.up,venn.dwn,padding=unit(pad,"line"),ncol=2,
               top="",left="",right="",bottom="")
 })
 
 output$tab0.plot2 <- renderPlot({
  data=tab0.data.up()
  hm.up <- heatmap_sigf_overlap(data)
  data=tab0.data.dwn()
  hm.dwn <- heatmap_sigf_overlap(data)
  pad=react.tab0.venn.pad()
  grid.arrange(hm.up,hm.dwn,padding=unit(pad,"line"),ncol=2,
               top="",left="",right="",bottom="")
 })
 
 output$tab0.plot3 <- renderPlot({	
  venn.opts=react.tab0.venn.opts()
  colors=react.tab0.colors()
  cex=react.tab0.venn.cex()
  pad=react.tab0.venn.pad()
  set.seed(1)
  s4 <- tab0.s4()
  plot_euler(s4=s4,colors=colors,cex=cex,venn.opts=venn.opts,title="All regulated genes")
 })
 
 output$tab0.plot4 <- renderPlot({
  data=tab0.data()
  heatmap_sigf_overlap(data)
 })	
 
 #### processing > 3 groups	
 
 react.tab3.grp.name <- reactive({
  grp.name <- sapply(grep("tab3\\.grp.+\\.name", x = names(input), value = TRUE),
                     function(x) input[[x]])
  grp.name=as.vector(grp.name[order(names(grp.name))])
  grp.name[grp.name!="NA"]
 })
 react.tab3.grp.plot.title <- reactive({
  grp.plot.title <- sapply(grep("tab3\\.grp.+\\.plot\\.title", x = names(input), value = TRUE),
                           function(x) input[[x]])
  grp.plot.title=as.vector(grp.plot.title[order(names(grp.plot.title))])
  grp.plot.title[grp.plot.title!=""]
 })
 react.tab3.color <- reactive({
  colors <- sapply(grep("tab3\\.color.+", x = names(input), value = TRUE),
                   function(x) input[[x]])
  colors=as.vector(colors[order(names(colors))])
  colors[colors!=""]
 })
 react.tab3.fdr.cutoff <- reactive({
  fdr.cutoff <- sapply(grep("tab3\\.fdr.+", x = names(input), value = TRUE),
                       function(x) input[[x]])
  fdr.cutoff=as.vector(fdr.cutoff[order(names(fdr.cutoff))])
  fdr.cutoff[fdr.cutoff!="NA"]
 })  
 react.tab3.fc.cutoff <- reactive({
  fc.cutoff <- sapply(grep("tab3\\.fc.+", x = names(input), value = TRUE),
                      function(x) input[[x]])
  fc.cutoff=as.vector(fc.cutoff[order(names(fc.cutoff))])
  fc.cutoff[fc.cutoff!="NA"]
 })
 react.tab3.meanDiff.cutoff <- reactive({
  meanDiff.cutoff <- sapply(grep("tab3\\.meanDiff.+", x = names(input), value = TRUE),
                            function(x) input[[x]])
  meanDiff.cutoff=as.vector(meanDiff.cutoff[order(names(meanDiff.cutoff))])
  meanDiff.cutoff[meanDiff.cutoff!="NA"]
 }) 
 react.tab3.venn.opts <- reactive({input$tab3.venn.opts})
 
 react.tab3.venn.cex <- reactive({input$tab3.venn.cex})
 
 react.tab3.venn.pad <- reactive({input$tab3.venn.pad})
 
 react.tab3.msigdb.species <- reactive({
  input$tab3.msigdb.species
 })
 
 react.tab3.nintersects.upset <- reactive({input$tab3.nintersects.upset})
 
 react.tab3.enrich.pval.co <- reactive({input$tab3.enrich.pval.co})
 
 react.tab3.cap.enrich.terms<- reactive({
  terms=toupper(str_trim(unlist(strsplit(input$tab3.cap.enrich.terms,","))))
  names(terms)=str_to_title(terms)
  dict = c("Mrna" = "mRNA", "Dna" = "DNA", "Rna" = "RNA", 
           "Trna" = "tRNA", "Mirna" = "miRNA", "Rrna" = "rRNA",
           "Atp" = "ATP","Adp" = "ADP","Snorna" = "snoRNA",
           "Lncrna" = "lncRNA","Snrna" = "snRNA",
           "Ncrna" = "ncRNA")
  
  if(length(terms)>0){
   dict = c(dict,terms)
  }
  dict
  })
 
 react.tab3.result <- reactive({
  grp.name=react.tab3.grp.name()
  grp.plot.title=react.tab3.grp.plot.title()
  out.DESeq2 <- react.tab3.rdata()
  idx=match(grp.name,names(out.DESeq2$results))
  results=out.DESeq2$results[idx]
 })
 
 react.tab3.normCts <- reactive({
  out.DESeq2 <- react.tab3.rdata()
  normCts=counts(out.DESeq2$dds,normalized=TRUE)
 })
 
 tab3.gs.RNASeq <- reactive({
  grp.name=react.tab3.grp.name()
  grp.plot.title=react.tab3.grp.plot.title()
  results=react.tab3.result()
  gs.RNASeq <- max(sapply(results,function(x) sum(!is.na(x$padj))))
 })
 
 tab3.grp.up <- reactive ({
  grp.name=react.tab3.grp.name()
  grp.plot.title=react.tab3.grp.plot.title()
  results=react.tab3.result()
  normCts=react.tab3.normCts()
  fdr.co <- react.tab3.fdr.cutoff()
  fc.cutoff <- react.tab3.fc.cutoff()
  meanDiff.cutoff <- react.tab3.meanDiff.cutoff()
  
  # up-regulated genes
  grp.up=list()
  for(i in 1:length(grp.name)){
   
   # filter out by mean difference
   grps=strsplit(grp.name[i],"_vs_")
   treat.grp=grps[[1]][1]
   ref.grp=grps[[1]][2]
   treat.mean=rowMeans(normCts[,grep(treat.grp,colnames(normCts))])
   ref.mean=rowMeans(normCts[,grep(ref.grp,colnames(normCts))])
   diff.mean=treat.mean-ref.mean
   # grp.up[[i]]=results[[grp.name[i]]][results[[grp.name[i]]]$padj<=fdr.co[i] &
   # results[[grp.name[i]]]$log2FoldChange>(log2(fc.cutoff[i])) &
   # !is.na(results[[grp.name[i]]]$padj)& diff.mean >= meanDiff.cutoff[i],]
   row.pass=which(results[[grp.name[i]]]$padj<=fdr.co[i] &
                   results[[grp.name[i]]]$log2FoldChange>(log2(fc.cutoff[i])) &
                   !is.na(results[[grp.name[i]]]$padj)& diff.mean >= meanDiff.cutoff[i])
   grp.up[[i]]=results[[grp.name[i]]][row.pass,]
   grp.up[[i]]$row.pass = NULL
   grp.up[[i]]$norm.diff.mean = NULL
   if(length(row.pass)>0){
    grp.up[[i]]$row.pass = row.pass
    grp.up[[i]]$norm.diff.mean = diff.mean[row.pass]
   }
  }
  grp.up
 })
 
 tab3.grp.dwn <- reactive({
  grp.name=react.tab3.grp.name()
  grp.plot.title=react.tab3.grp.plot.title()
  results=react.tab3.result()
  normCts = react.tab3.normCts()
  fdr.co <- react.tab3.fdr.cutoff()
  fc.cutoff <- react.tab3.fc.cutoff()
  meanDiff.cutoff <- react.tab3.meanDiff.cutoff()
  
  # down-regulated genes
  grp.dwn=list()
  for(i in 1:length(grp.name)){
   # filter out by mean difference
   grps=strsplit(grp.name[i],"_vs_")
   treat.grp=grps[[1]][1]
   ref.grp=grps[[1]][2]
   treat.mean=rowMeans(normCts[,grep(treat.grp,colnames(normCts))])
   ref.mean=rowMeans(normCts[,grep(ref.grp,colnames(normCts))])
   diff.mean=ref.mean-treat.mean
   grp.dwn[[i]]=results[[grp.name[i]]][results[[grp.name[i]]]$padj<=fdr.co[i] &
                                        results[[grp.name[i]]]$log2FoldChange<(-log2(fc.cutoff[i])) &
                                        !is.na(results[[grp.name[i]]]$padj) & diff.mean >= meanDiff.cutoff[i],]
   row.pass=which(results[[grp.name[i]]]$padj<=fdr.co[i] &
                   results[[grp.name[i]]]$log2FoldChange<(-log2(fc.cutoff[i])) &
                   !is.na(results[[grp.name[i]]]$padj) & diff.mean >= meanDiff.cutoff[i])
   grp.dwn[[i]]=results[[grp.name[i]]][row.pass,]
   grp.dwn[[i]]$row.pass = NULL
   grp.dwn[[i]]$norm.diff.mean = NULL
   if(length(row.pass)>0){
    grp.dwn[[i]]$row.pass = row.pass
    grp.dwn[[i]]$norm.diff.mean = diff.mean[row.pass]
   }
  }#for(i in 1:length(grp.name))
  grp.dwn
 })
 
 tab3.genes.lists <- reactive ({
  grp.name=react.tab3.grp.name()
  grp.plot.title=react.tab3.grp.plot.title()
  results=react.tab3.result()
  grp.dwn <- tab3.grp.dwn ()
  grp.up <- tab3.grp.up ()
  
  # all regulated-genes
  grp=list()
  for(i in 1:length(grp.name)){
   grp[[i]]=rbind(grp.up[[i]],grp.dwn[[i]])
  }
  
  genes.lists=lapply(grp,function(x) rownames(x))
  names(genes.lists)<-grp.plot.title
  genes.lists
 })
 
 tab3.genes.lists.up <- reactive ({
  grp.name=react.tab3.grp.name()
  grp.plot.title=react.tab3.grp.plot.title()
  results=react.tab3.result()
  grp.up <- tab3.grp.up ()
  genes.lists.up=lapply(grp.up,function(x) rownames(x))
  names(genes.lists.up) <- grp.plot.title
  genes.lists.up
 })
 
 tab3.genes.lists.dwn <- reactive ({
  grp.name=react.tab3.grp.name()
  grp.plot.title=react.tab3.grp.plot.title()
  results=react.tab3.result()
  grp.dwn <- tab3.grp.dwn ()
  genes.lists.dwn=lapply(grp.dwn,function(x) rownames(x))
  names(genes.lists.dwn) <- grp.plot.title
  genes.lists.dwn
 })
 
 tab3.gom.obj <- reactive({
  genes.lists <- tab3.genes.lists()
  gs.RNASeq <- tab3.gs.RNASeq()
  
  newGOM(genes.lists,genes.lists,genome.size=gs.RNASeq)
 })
 
 tab3.gom.obj.up <- reactive({
  genes.lists <- tab3.genes.lists.up()
  gs.RNASeq <- tab3.gs.RNASeq()
  newGOM(genes.lists,genes.lists,genome.size=gs.RNASeq)
 })
 
 tab3.gom.obj.dwn <- reactive({
  genes.lists <- tab3.genes.lists.dwn()
  gs.RNASeq <- tab3.gs.RNASeq()
  newGOM(genes.lists,genes.lists,genome.size=gs.RNASeq)
 })
 
 tab3.data <- reactive({
  gom.obj <- tab3.gom.obj()
  grp.plot.title <- react.tab3.grp.plot.title()
  
  # calculate overlap
  data <- getMatrix(gom.obj,"pval")
  data <- melt(data,id=grp.plot.title)
  data <- data
  names(data)[names(data)=="value"]="pval"
  tmp <- getMatrix(gom.obj,"Jaccard")
  tmp <- melt(tmp,id=grp.plot.title)
  tmp <- tmp
  data$jaccard=tmp$value
  # categorizing pvalues
  data$pval.cat=ifelse(data$pval==0,pval.names[1],ifelse(data$pval<=0.05,pval.names[2],
                                                         ifelse(data$pval<0.5,pval.names[3],pval.names[4])))
  data$pval.cat=factor(data$pval.cat,levels=pval.names)
  data
 })
 
 tab3.data.up <- reactive({
  gom.obj <- tab3.gom.obj.up()
  grp.plot.title <- react.tab3.grp.plot.title()
  # calculate overlap
  data <- getMatrix(gom.obj,"pval")
  data <- melt(data,id=grp.plot.title)
  data <- data
  names(data)[names(data)=="value"]="pval"
  tmp <- getMatrix(gom.obj,"Jaccard")
  tmp <- melt(tmp,id=grp.plot.title)
  tmp <- tmp
  data$jaccard=tmp$value
  # categorizing pvalues
  data$pval.cat=ifelse(data$pval==0,pval.names[1],ifelse(data$pval<=0.05,pval.names[2],
                                                         ifelse(data$pval<0.5,pval.names[3],pval.names[4])))
  data$pval.cat=factor(data$pval.cat,levels=pval.names)
  data
 })
 
 tab3.data.dwn <- reactive({
  gom.obj <- tab3.gom.obj.dwn()
  grp.plot.title <- react.tab3.grp.plot.title()
  # calculate overlap
  data <- getMatrix(gom.obj,"pval")
  data <- melt(data,id=grp.plot.title)
  data <- data
  names(data)[names(data)=="value"]="pval"
  tmp <- getMatrix(gom.obj,"Jaccard")
  tmp <- melt(tmp,id=grp.plot.title)
  tmp <- tmp
  data$jaccard=tmp$value
  # categorizing pvalues
  data$pval.cat=ifelse(data$pval==0,pval.names[1],ifelse(data$pval<=0.05,pval.names[2],
                                                         ifelse(data$pval<0.5,pval.names[3],pval.names[4])))
  data$pval.cat=factor(data$pval.cat,levels=pval.names)
  data
 })
 
 tab3.s4 <- reactive({
  # input for euler venn
  s4 <-tab3.genes.lists()
  grp.plot.title <-react.tab3.grp.plot.title()
  
  names(s4) <- grp.plot.title
  s4
 })
 
 tab3.s4.up <- reactive({
  # input for euler venn
  s4 <-tab3.genes.lists.up()
  grp.plot.title <-react.tab3.grp.plot.title()
  
  names(s4) <- grp.plot.title
  s4
 })
 
 tab3.s4.dwn <- reactive({
  # input for euler venn
  s4 <-tab3.genes.lists.dwn()
  grp.plot.title <-react.tab3.grp.plot.title()
  names(s4) <- grp.plot.title
  s4
 })
 
 tab3.compare.df <- reactive({
  # input for clusterProfiler
  grp = tab3.grp.up()
  grp.name = react.tab3.grp.name()
  grp.plot.title =react.tab3.grp.plot.title()
  compare.df = NULL
  for(i in 1:length(grp.plot.title)){
   compare.df=rbind(compare.df,data.frame(SYMBOL=rownames(grp[[i]]),
                                          group1=rep(grp.plot.title[i],dim(grp[[i]])[1]),
                                          group2=rep("Up-regulated",dim(grp[[i]])[1])))
  }
  grp = tab3.grp.dwn()
  for(i in 1:length(grp.plot.title)){
   compare.df=rbind(compare.df,data.frame(SYMBOL=rownames(grp[[i]]),
                                          group1=rep(grp.plot.title[i],dim(grp[[i]])[1]),
                                          group2=rep("Down-regulated",dim(grp[[i]])[1])))
  }
  compare.df
 })
 
  # preprocess multiple search term
 react.search.term <- reactive({
  if(!is.null(input$tab3.search.term) && !input$tab3.search.term =="")
   paste(str_trim(unlist(strsplit(input$tab3.search.term,","))),collapse="|")
 })
 
 # Filter data based on the search input
 react.tab3.filtered.data <- reactive({
  search_term <- react.search.term()
  data_list <- react.tab3.expr.tbl()
  current_page <- input$currentPage
  if (is.null(search_term) || search_term == "") {
   rows_per_page <- 10
   start_row <- (current_page-1)*rows_per_page+1
   end_row <- min(current_page*rows_per_page,sum(sapply(data_list, nrow)))
   sliced_data_lst <- lapply(data_list, function(df) {
     slice(df, start_row:end_row)
    })
   return(sliced_data_lst)
  } else {
   filtered_list <- lapply(data_list, function(df) {
    df[grep(search_term, df$Genes, ignore.case = TRUE), , drop = FALSE]
   })
   return(filtered_list)
  }
 })

 # Dynamically generate tables with titles
 output$tables <- renderUI({
  data_list <- react.tab3.expr.tbl()
  lapply(names(data_list), function(tab_name) {
   table_id <- paste0("table_", tab_name)
   tagList(
    h4(paste("Gene expression for", tab_name)),
    tableOutput(table_id)
   )
  })
 })
 
 # Render the filtered data tables
 observe({
  data_list <- react.tab3.expr.tbl()
  filtered.data <- react.tab3.filtered.data()
  for (tab_name in names(data_list)) {
   # Need local so that each item gets its own number. Without it, the value
   # of tab_name in the renderTable() will be the same across all instances, because
   # of when the expression is evaluated.
   # https://gist.github.com/wch/5436415/
   local({
    tab_name <- tab_name
   table_id <- paste0("table_", tab_name)
   output[[table_id]] <- renderTable({
     as.data.frame(filtered.data[[tab_name]])
   })
   })#local
  }
 })

 # Insert the right number of plot output objects into the web page
 output$tab3.plots <- renderUI({
  grp.name=react.tab3.grp.name()
  genes.lists.dwn = tab3.genes.lists.dwn()
  genes.lists.up = tab3.genes.lists.up()
  plot_output_list<-list()
  for (i in 1:length(grp.name)){
   plot_output_list[[i]] <- plotOutput(paste0("tab3.volcano.grp",i),height="600px")
  }
  plot.i=1
  if(length(grp.name)<8){
   if(sum(unlist(sapply(genes.lists.up,function(x)length(x))))>0){
    plot_output_list[[length(grp.name)+plot.i]] <- plotOutput("tab3.venn.up")
    plot.i = plot.i +1
   }
   if(sum(unlist(sapply(genes.lists.dwn,function(x)length(x))))>0){
    plot_output_list[[length(grp.name)+plot.i]] <- plotOutput("tab3.venn.dwn")
    plot.i=plot.i+1
   }
   if(length(grp.name)>1)
    plot_output_list[[length(grp.name)+plot.i]] <- plotOutput("tab3.hm")
  }else{
   plot_output_list[[length(grp.name)+plot.i]] <- plotOutput("tab3.hm")
  }
  
  # Convert the list to a tagList - this is necessary for the list of items
  # to display properly.
  do.call(tagList, plot_output_list)
 })
 
 # # Call renderPlot for each one. Plots are only actually generated when they
 # # are visible on the web page.
 
 react.tab3.hilite.genes <- reactive({
  unlist(sapply(strsplit(input$tab3.hilite.genes,","),function(x)str_trim(x)))
 })
 
 react.tab3.expr.tbl <- reactive({
  print("input$inTabSet")
  print(input$inTabSet)
  req(input$inTabset=="tab3")
  grp.name=react.tab3.grp.name()
  grp.plot.title=react.tab3.grp.plot.title()
  results=react.tab3.result()
  normCts=react.tab3.normCts()
  fdr.co <- react.tab3.fdr.cutoff()
  fc.cutoff <- react.tab3.fc.cutoff()
  meanDiff.cutoff <- react.tab3.meanDiff.cutoff()
  expr.tbl <- list()
  for(i in 1:length(grp.name)){
   # get indiv grpname
   grps=strsplit(grp.name[i],"_vs_")
   treat.grp=grps[[1]][1]
   ref.grp=grps[[1]][2]
   treat.mean=rowMeans(normCts[,grep(treat.grp,colnames(normCts))])
   ref.mean=rowMeans(normCts[,grep(ref.grp,colnames(normCts))])
   diff.mean=treat.mean-ref.mean
   data=data.frame(
    Genes=rownames(results[[grp.name[i]]]),
    logFC=results[[grp.name[i]]]$log2FoldChange,
    FDR=results[[grp.name[i]]]$padj,
    diff.mean=diff.mean)
   data <- data %>% 
    mutate(
     Expression = 
      case_when(logFC >= log2(fc.cutoff[i]) & 
                 FDR <= fdr.co[i] & 
                 diff.mean >= meanDiff.cutoff[i] ~ "Up-regulated",
                logFC < -log2(fc.cutoff[i]) &
                 FDR <= fdr.co[i] & 
                 diff.mean <= -meanDiff.cutoff[i] ~ "Down-regulated",
                TRUE ~ "NS")
    )
   expr.tbl[[i]] <- data
  }
  names(expr.tbl)=grp.plot.title
  expr.tbl
 })
 
 observe({
 grp.name=react.tab3.grp.name()
 grp.plot.title=react.tab3.grp.plot.title()
 results=react.tab3.result()
 normCts=react.tab3.normCts()
 fdr.co <- react.tab3.fdr.cutoff()
 fc.cutoff <- react.tab3.fc.cutoff()
 meanDiff.cutoff <- react.tab3.meanDiff.cutoff()
 hilite.genes <- react.tab3.hilite.genes()
 for(i in 1:length(grp.name)){
  # Need local so that each item gets its own number. Without it, the value
  # of i in the renderPlot() will be the same across all instances, because
  # of when the expression is evaluated.
  # https://gist.github.com/wch/5436415/
  local({
   my_i <- i
 output[[paste0("tab3.volcano.grp",my_i)]] <- renderPlot({
   # get indiv group name
   grps=strsplit(grp.name[my_i],"_vs_")
   treat.grp=grps[[1]][1]
   ref.grp=grps[[1]][2]
   treat.mean=rowMeans(normCts[,grep(treat.grp,colnames(normCts))])
   ref.mean=rowMeans(normCts[,grep(ref.grp,colnames(normCts))])
   diff.mean=treat.mean-ref.mean
   data=data.frame(
    Genes=rownames(results[[grp.name[my_i]]]),
    logFC=results[[grp.name[my_i]]]$log2FoldChange,
    FDR=results[[grp.name[my_i]]]$padj,
    diff.mean=diff.mean)
   data <- data %>% 
    mutate(
     Expression = 
      case_when(logFC >= log2(fc.cutoff[my_i]) & 
                 FDR <= fdr.co[my_i] & 
                 diff.mean >= meanDiff.cutoff[my_i] ~ "Up-regulated",
                logFC < -log2(fc.cutoff[my_i]) &
                 FDR <= fdr.co[my_i] & 
                 diff.mean <= -meanDiff.cutoff[my_i] ~ "Down-regulated",
                TRUE ~ "NS")
    )
   # change FDR=0 so it can be graphed correctly
   data$logFDR=-log(data$FDR+min(c(data$FDR[data$FDR>0],1e-32),na.rm=TRUE),10)
   genes_show = hilite.genes
   genes_show_data <- dplyr::bind_rows(
    data %>%
     filter(Genes %in% genes_show)
   )
   xlim=c(min(data$logFC)*1.1,max(data$logFC)*1.1)
   ylim=c(0,max(data$logFDR,na.rm=TRUE)*1.1)
   p <- ggplot(data, aes(logFC, logFDR)) +
    geom_hline(yintercept=-log10(fdr.co[my_i]), col="#5d5d5d",linetype="dashed")+
    geom_point(aes(color = Expression), size = 4,alpha=4/5) +
    xlab(expression("log"[2]*"FC")) + 
    ylab(expression("-log"[10]*"FDR")) +
    scale_color_manual(values = c(tab3.colors[2], "gray50", tab3.colors[6])) +
    guides(colour = guide_legend(override.aes = list(size=2))) + 
    geom_label_repel(data = genes_show_data,mapping=aes(logFC,logFDR,label=Genes),
                     size=5,min.segment.length = 0) +
    theme_classic(base_size = 18) +
    theme(panel.grid.major=element_line(color = "#EBEBEB"),
          panel.grid.minor=element_line(color = "#EBEBEB"),
          plot.title=element_text(hjust=0.5))+
    ggtitle(paste0("Volcano plot for ",grp.plot.title[my_i]))
   if(fc.cutoff[my_i]>1){
    p <- p + 
     geom_vline(xintercept=c(-log2(fc.cutoff[my_i]), log2(fc.cutoff[my_i])), col="#5d5d5d",linetype="dashed")
   }
   p
 })#renderPlot  
  })#local
  }#for(i in 1:length(grp.name))
 })#observe

 output[["tab3.venn.up"]] <- renderPlot({
  venn.opts=react.tab3.venn.opts()
  grp.name=react.tab3.grp.name()
  colors=react.tab3.color()
  cex=react.tab3.venn.cex()
  pad=react.tab3.venn.pad()
  set.seed(1)
  s4 <- tab3.s4.up()
  # if all lists contains 0 genes then do not proceeds
  nonzero_s4=sum(sapply(s4,function(x) length(x)!=0))
  shiny::req(nonzero_s4>0)
  vals.plot$venn.up2<-plot_euler(s4=s4,colors=colors,cex=cex,
                                 venn.opts=venn.opts,title="")
  # remove label if requested
  if(!'Labels' %in% venn.opts )
   names(s4)=rep(" ",length(s4))
  vals.plot$venn.up1<-venn(s4,zcolor=colors,opacity=.8,box=FALSE,
                           ilcs = 0.8, sncs = 1,ggplot=TRUE)
  grid.arrange(vals.plot$venn.up1,vals.plot$venn.up2,ncol=2,
               padding=unit(pad,"line"),
               top=textGrob("Up-regulated Genes",gp=gpar(fontsize=20)))
 })
 
 output[["tab3.venn.dwn"]] <- renderPlot({
  venn.opts=react.tab3.venn.opts()
  grp.name=react.tab3.grp.name()
  colors=react.tab3.color()
  cex=react.tab3.venn.cex()
  pad=react.tab3.venn.pad()
  set.seed(1)
  s4 <- tab3.s4.dwn()
  # if all lists contains 0 genes then do not proceeds
  nonzero_s4=sum(sapply(s4,function(x) length(x)!=0))
  shiny::req(nonzero_s4>0)
  vals.plot$venn.dwn2<-plot_euler(s4=s4,colors=colors,cex=cex,
                                  venn.opts=venn.opts,title="")
  # remove label if requested
  if(!'Labels' %in% venn.opts )
   names(s4)=rep(" ",length(s4))
  vals.plot$venn.dwn1 <- venn(s4,zcolor=colors,opacity=.8,box=FALSE,ilcs = 0.8,
                              sncs = 1,ggplot=TRUE)
  grid.arrange(vals.plot$venn.dwn1,vals.plot$venn.dwn2,ncol=2,
               padding=unit(pad,"line"),
               top=textGrob("Down-regulated Genes",gp=gpar(fontsize=20)))
 })
 
 output[["tab3.hm"]] <- renderPlot({
  pad = react.venn.pad()
  data = tab3.data.up()
  grp.name = react.tab3.grp.name()
  shiny::req(length(grp.name)>1)
  vals.plot$hm.up <- heatmap_sigf_overlap(data,"Up-regulated genes")
  data=tab3.data.dwn()
  vals.plot$hm.dwn<- heatmap_sigf_overlap(data,"Down-regulated genes")
  grid.arrange(vals.plot$hm.dwn,vals.plot$hm.up,ncol=2)
 })
 
 output$tab3.plot1 <- renderPlot({
  s4 <- tab3.s4.up()
  grp.name <- react.tab3.grp.name()
  # only draw upset up-regulated if there are more than 1 group with non-empty list of up genes
  shiny::req(length(grp.name)>1 && sum(sapply(s4,function(x)length(x)!=0))>1)
  nintersects <- react.tab3.nintersects.upset()
  vals.plot$upset.up <- as.ggplot(upset(fromList(s4), order.by = "freq",
                                        nsets=length(grp.name), nintersects=nintersects, 
                                        mainbar.y.label = "# Up-regulated Genes", sets.x.label = "Significant Genes", 
                                        point.size=2.8,line.size=1,
                                        #c(intersection size title, intersection size tick labels, set size title, set size tick labels, set names, numbers above bars)
                                        text.scale = c(1.8, 1.8, 1.4, 1.8, 1.8, 1.5)))




  vals.plot$upset.up
 })

 
 output$tab3.plot2 <- renderPlot({
  s4 <- tab3.s4.dwn()
  grp.name <- react.tab3.grp.name()
  # only draw upset if there are more than two groups with non-empty list of diff genes
  shiny::req(length(grp.name)>1 && sum(sapply(s4,function(x)length(x)!=0))>1)
  nintersects <- react.tab3.nintersects.upset()
  vals.plot$upset.dwn <- as.ggplot(upset(fromList(s4), order.by = "freq",nsets=length(grp.name), 
                                         nintersects=nintersects,
                                         mainbar.y.label = "# Down-regulated Genes", sets.x.label = "Significant Genes", 
                                         point.size=2.8,line.size=1,
                                         #c(intersection size title, intersection size tick labels, set size title, set size tick labels, set names, numbers above bars)
                                         text.scale = c(1.8, 1.8, 1.4, 1.8, 1.8, 1.5)))
  vals.plot$upset.dwn
 })
 
 # gene sets from msigdb
 # # oncogenic signature gene
 # m_t2g.c6 <- msigdbr(species = msigdb.species, category = "C6") %>%
 #  dplyr::select(gs_name, gene_symbol)
 # # hallmark gene sets
 # m_t2g.h <- msigdbr(species = msigdb.species, category = "H") %>%
 #  dplyr::select(gs_name, gene_symbol)
 # m_t2g.c2.kegg <- msigdbr(species = msigdb.species, category = "C2",subcategory = "CP:KEGG") %>%
 #  dplyr::select(gs_name, gene_symbol)
 # m_t2g.c2.reactome <- msigdbr(species = msigdb.species, category = "C2",subcategory = "CP:REACTOME") %>%
 #  dplyr::select(gs_name, gene_symbol)
 # m_t2g.c2.biocarta <- msigdbr(species = msigdb.species, category = "C2",subcategory = "CP:BIOCARTA") %>%
 #  dplyr::select(gs_name, gene_symbol)
 
 observeEvent(input$gen.go,{
  output$tab3.plot3.1 <- renderPlot({
   withProgress(message="Generating enrichment plots",{
    enrich.pval.co <- react.tab3.enrich.pval.co()
    compare.df <- tab3.compare.df()
    grp.plot.title <- react.tab3.grp.plot.title()
    cap.enrich.terms <- react.tab3.cap.enrich.terms()
    # Using clusterProfiler to perform hypergeometric test on msigdb signatures
    msigdb.species <- react.tab3.msigdb.species()
    msig.gene.set <- msigdbr(species = msigdb.species, category = "C5",subcategory = "MF") %>%
     dplyr::select(gs_name, gene_symbol)
    msig.name ="MSigDB GO Molecular Function"

    formula_res <- compareCluster(SYMBOL~group1+group2, data=compare.df, fun="enricher",
                                  TERM2GENE=msig.gene.set,pvalueCutoff=enrich.pval.co,
                                  pAdjustMethod="BH")
    print("formula_res")
    print(formula_res)
    # re-arrange datasets using factor
    # and do pathway analysis using up- and down-regulated genes separately
    vals.plot$msig.mf <- dotplot(formula_res,x=~factor(group1),font.size=14,title=msig.name) + 
     facet_grid(~group2) +
     scale_y_discrete(labels=function(x) 
      str_wrap(str_replace_all(str_to_title(tolower(gsub("_"," ",gsub("GOMF_","",x)))),cap.enrich.terms), width=40)) +
     scale_x_discrete(labels=function(x) 
      str_wrap(x,width=10)) +
     scale_color_distiller(palette = 'Blues')
    vals.plot$msig.mf
   })#withProgress
  },height=1000)
 })#observeEvent
 
 observeEvent(input$gen.go,{
  output$tab3.plot3.2 <- renderPlot({
   withProgress(message="Generating enrichment plots",{
    enrich.pval.co <- react.tab3.enrich.pval.co()
    compare.df <- tab3.compare.df()
    grp.plot.title <- react.tab3.grp.plot.title()
    # Using clusterProfiler to perform hypergeometric test on msigdb signatures
    msigdb.species <- react.tab3.msigdb.species()
    msig.gene.set <- msigdbr(species = msigdb.species, category = "C5",subcategory = "BP") %>%
     dplyr::select(gs_name, gene_symbol)
    msig.name ="MSigDB GO Biological Process"
    cap.enrich.terms <- react.tab3.cap.enrich.terms()
    formula_res <- compareCluster(SYMBOL~group1+group2, data=compare.df, fun="enricher",
                                  TERM2GENE=msig.gene.set,pvalueCutoff=enrich.pval.co,
                                  pAdjustMethod="BH")
    vals.plot$msig.bp <- dotplot(formula_res,x=~factor(group1),font.size=14,title=msig.name) + 
     facet_grid(~group2) +
     scale_y_discrete(labels=function(x) 
      str_wrap(str_replace_all(str_to_title(tolower(gsub("_"," ",gsub("GOBP_","",x)))),cap.enrich.terms), width=40)) +
     scale_x_discrete(labels=function(x) str_wrap(x,width=10)) +
     scale_color_distiller(palette = 'Blues')
    vals.plot$msig.bp
   })#withProgress
  },height=1000)
 })#observeEvent
 
 observeEvent(input$gen.go,{
  output$tab3.plot3.3 <- renderPlot({
   withProgress(message="Generating enrichment plots",{
    enrich.pval.co <- react.tab3.enrich.pval.co()
    compare.df <- tab3.compare.df()
    grp.plot.title <- react.tab3.grp.plot.title()
    cap.enrich.terms <- react.tab3.cap.enrich.terms()
    # Using clusterProfiler to perform hypergeometric test on msigdb signatures
    msigdb.species <- react.tab3.msigdb.species()
    msig.gene.set <- msigdbr(species = msigdb.species, category = "C5",subcategory = "CC") %>%
     dplyr::select(gs_name, gene_symbol)
    msig.name ="MSigDB GO Cellular Component"
    formula_res <- compareCluster(SYMBOL~group1+group2, data=compare.df, fun="enricher",
                                  TERM2GENE=msig.gene.set,pvalueCutoff=enrich.pval.co,
                                  pAdjustMethod="BH")
    vals.plot$msig.cc <- dotplot(formula_res,x=~factor(group1),font.size=14,title=msig.name) + 
     facet_grid(~group2) +
     scale_y_discrete(labels=function(x) 
      str_wrap(str_replace_all(str_to_title(tolower(gsub("_"," ",gsub("GOCC_","",x)))),cap.enrich.terms), width=40)) +
     scale_x_discrete(labels=function(x) str_wrap(x,width=10)) +
     scale_color_distiller(palette = 'Blues')
    vals.plot$msig.cc
   }) #withProgress
  },height=1000)
 }) #observeEvent
 
 
 observeEvent(input$gen.go,{
  output$tab3.plot3.4 <- renderPlot({
   withProgress(message="Generating enrichment plots",{
    enrich.pval.co <- react.tab3.enrich.pval.co()
    compare.df <- tab3.compare.df()
    grp.plot.title <- react.tab3.grp.plot.title()
    cap.enrich.terms <- react.tab3.cap.enrich.terms()
    # Using clusterProfiler to perform hypergeometric test on msigdb signatures
    msigdb.species <- react.tab3.msigdb.species()
    msig.gene.set <- msigdbr(species = msigdb.species, category = "C2") %>%
     dplyr::select(gs_name, gene_symbol)
    msig.name ="MSigDB Curated Gene Sets"
    formula_res <- compareCluster(SYMBOL~group1+group2, data=compare.df, fun="enricher",
                                  TERM2GENE=msig.gene.set,
                                  pvalueCutoff=enrich.pval.co,pAdjustMethod="BH")
    vals.plot$msig.curate <- dotplot(formula_res,x=~factor(group1),font.size=14,title=msig.name) + facet_grid(~group2) +
     scale_y_discrete(labels=function(x) 
            str_wrap(str_replace_all(str_to_title(tolower(gsub("_"," ",x))),cap.enrich.terms), width=40)) +
     scale_x_discrete(labels=function(x) str_wrap(x,width=10)) +
     scale_color_distiller(palette = 'Blues')
    vals.plot$msig.curate
   }) # withProgress
  },height=1000)
 })#observeEvent
 
 ## clicking on the export button will generate a pdf file 
 ## containing all plots
 observeEvent( input$export,{
  if(input$export==0) return()
  gen.go <- react.val.gen.go()
  out.loc <- react.tab3.out.loc()
  withProgress(message="Saving pdf",{
   pad = react.tab3.venn.pad()
   file=file.path(out.loc,"plots.pdf")
   if(file.exists(file)) file.remove(file)
   pdf(file, onefile = TRUE,width=input$pdf.width,height=input$pdf.height)
   gridExtra::grid.arrange(vals.plot$venn.up1,vals.plot$venn.up2,ncol=2,
                           padding=unit(pad,"line"),
                           top=textGrob("Up-regulated Genes",gp=gpar(fontsize=20)),
                           vp=viewport(width=0.8, height=0.8))
   gridExtra::grid.arrange(vals.plot$venn.dwn1,vals.plot$venn.dwn2,ncol=2,
                           padding=unit(pad,"line"),
                           top=textGrob("Down-regulated Genes",gp=gpar(fontsize=20)),
                           vp=viewport(width=0.8, height=0.8))
   gridExtra::grid.arrange(vals.plot$hm.up,vals.plot$hm.dwn,ncol=2,padding=unit(pad,"line"),
                           top=textGrob("Significance of overlaps",gp=gpar(fontsize=20)),
                           vp=viewport(width=0.8, height=0.8))
   gridExtra::grid.arrange(vals.plot$upset.up,vp=viewport(width=0.8, height=0.8))
   gridExtra::grid.arrange(vals.plot$upset.dwn,vp=viewport(width=0.8, height=0.8))
   if(gen.go!=0){
    gridExtra::grid.arrange(vals.plot$msig.mf,vp=viewport(width=0.8, height=0.8))
    gridExtra::grid.arrange(vals.plot$msig.bp,vp=viewport(width=0.8, height=0.8))
    gridExtra::grid.arrange(vals.plot$msig.cc,vp=viewport(width=0.8, height=0.8))
    gridExtra::grid.arrange(vals.plot$msig.curate,vp=viewport(width=0.8, height=0.8))
   } #if(gen.go!=0)
   dev.off()
  })#withProgress
 }) # observeEvent
 
 ## clicking on the export button will generate a txt file 
 ## containing filtered gene lists
 observeEvent( input$exportGene,{
  if(input$exportGene==0) return()
  withProgress(message="Exporting gene lists....",{
   grp.name = react.tab3.grp.name()
   grp.up = tab3.grp.up()
   grp.dwn = tab3.grp.dwn()
   fc.co = react.tab3.fc.cutoff()
   fdr.co = react.tab3.fdr.cutoff()
   meanDiff.co = react.tab3.meanDiff.cutoff()
   projdir = react.setup.proj.dir()
   out.loc = react.tab3.out.loc()
   for (i in (1:length(grp.name))){
    # filter the complete table file
    in.file=paste0(gsub("_","",grp.name[i]),".complete.txt")
    df=read.table(paste0(projdir,"outputs/diff_analysis_rslt/tables/",in.file),
                  sep="\t",header=TRUE)
    df.filt=df[grp.up[[i]]$row.pass,]
    df.filt$norm.diff.mean = df[grp.up[[i]]$row.pass,"norm.diff.mean"]
    # write filter up txt
    out.up.file=paste0(gsub("_","",grp.name[i]),
                       paste0("_FC_",fc.co[i],"_FDR_",fdr.co[i],
                              "_meanDiff_",meanDiff.co[i]),".up.txt")
    write.table(df.filt,file=paste0(out.loc,"/",out.up.file),sep="\t",row.names = FALSE,
                quote=FALSE)
    df.filt=df[grp.dwn[[i]]$row.pass,]
    df.filt$norm.diff.mean = df[grp.dwn[[i]]$row.pass,"norm.diff.mean"]
    out.dwn.file=paste0(gsub("_","",grp.name[i]),
                        paste0("_FC_",fc.co[i],"_FDR_",fdr.co[i],
                               "_meanDiff_",meanDiff.co[i]),".down.txt")
    write.table(df.filt,file=paste0(out.loc,"/",out.dwn.file),sep="\t",row.names = FALSE,
                quote=FALSE)
   }
  })#withProgress
 }) # observeEvent
 
 output$tab3.text <- renderPrint({ grp.up = tab3.grp.up(); length(grp.up)})
}
# app
shinyApp(ui = ui, server = server)
