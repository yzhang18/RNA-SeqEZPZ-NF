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

# data directory
data.loc="/mnt/outputs"
# create directory for shiny outputs
dir.create(file.path(data.loc, "shiny_outputs"))
out.loc=file.path(data.loc,"shiny_outputs")

# filename was changed. Make sure can read both
if(file.exists("/mnt/outputs/diff_analysis_rslt/RNA-seq_differential_analysis.RData")){
	load("/mnt/outputs/diff_analysis_rslt/RNA-seq_differential_analysis.RData")
} else {
	load("/mnt/outputs/diff_analysis_rslt/RNA-seq differential analysis.RData")
}

groups= data.frame(
  var = names(out.DESeq2$results),
  names = names(out.DESeq2$results)
)

# List of choices for groupnames
groups.lst <- as.list(groups$names)
# Name it
names(groups.lst) <- groups$var
	
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
    tabsetPanel(
      tabPanel("Two Groups", fluid = TRUE,
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
				inputId="tab0.grp1.name",
				label = ("Group 1"),
				selected=names(out.DESeq2$results)[1],
				choices = groups.lst),
			textInput(inputId="tab0.grp1.plot.title", 
				label = ("Group 1 label"), 
				value = names(out.DESeq2$results)[1]),
			selectInput(
				inputId="tab0.grp2.name",
				label = ("Group 2"),
				selected=names(out.DESeq2$results)[2],
				choices = groups.lst),
			textInput(
				inputId="tab0.grp2.plot.title", 
				label = ("Group 2 label"), 
				value = names(out.DESeq2$results)[2])
			),
			mainPanel(
				fluidRow(style = "background-color:#F5F5F5",
					column(3,checkboxGroupInput("tab0.venn.opts", label = ("Display on Venns"), 
						choices = venn.opts.lst,
							selected = venn.opts$lbl[1:3]),
					# slider bar font size
					sliderInput("tab0.venn.cex", label = ("Venns' font size"), min = 0, 
						max = 6, value = 1,step = 0.1),
					# slider bar pad
					sliderInput("tab0.venn.pad", label = ("Venns' circle size"), min = 0, 
						max = 8, value = 1,step = 0.5)	
						),
					column(3,
					textInput("tab0.color.grp1", label = ("Color for group 1"), 
						value = tab0.colors[1]),
					textInput("tab0.color.grp2", label = ("Color for group 2"), 
						value = tab0.colors[2])),
					column(3,
					numericInput("tab0.fdr.grp1", label = ("FDR cut-off for group 1"), 
						value = 0.05,step=0.01,min=0,max=1),
					numericInput("tab0.fdr.grp2", label = ("FDR cut-off for group 2"),
						value = 0.05,step=0.01,min=0,max=1)),
					column(3,
					numericInput("tab0.fc.grp1", label = ("Fold-change for group 1"), 
						value = 1,step=0.5,min=1),
					numericInput("tab0.fc.grp2", label = ("Fold-change for group 2"), 
						value = 1,step=0.5,min=1))
				),
				br(),
				fluidRow(
					plotOutput(outputId = "tab0.plot1")),
				br(),
				fluidRow(
					plotOutput(outputId = "tab0.plot2")),
				br(),
				br(),
				fluidRow(column(12,plotOutput(outputId = "tab0.plot3"))),
				fluidRow(column(12,plotOutput(outputId = "tab0.plot4"))),
				br()
			)
		)
      ),
      tabPanel("Three Groups", fluid = TRUE,
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
				inputId="grp1.name",
				label = ("Group 1"),
				selected=names(out.DESeq2$results)[1],
				choices = groups.lst),
			textInput(inputId="grp1.plot.title", 
				label = ("Group 1 label"), 
				value = names(out.DESeq2$results)[1]),
			selectInput(
				inputId="grp2.name",
				label = ("Group 2"),
				selected=names(out.DESeq2$results)[2],
				choices = groups.lst),
			textInput(
				inputId="grp2.plot.title", 
				label = ("Group 2 label"), 
				value = names(out.DESeq2$results)[2]),
			selectInput(
				inputId="grp3.name",
				label = ("Group 3"),
				selected=names(out.DESeq2$results)[3],
				choices = groups.lst),
			textInput(
				inputId="grp3.plot.title", 
				label = ("Group 3 label"), 
				value = names(out.DESeq2$results)[3])
			),
			mainPanel(
				fluidRow(style = "background-color:#F5F5F5",column(3,
					checkboxGroupInput("venn.opts", label = ("Display on Venns"), 
						choices = venn.opts.lst,
							selected = venn.opts$lbl[1:3]),
					    
					# slider bar font size
					sliderInput("venn.cex", label = ("Venns' font size"), min = 0, 
						max = 6, value = 1,step = 0.1),
					# slider bar pad
					sliderInput("venn.pad", label = ("Venns' circle size"), min = 0, 
						max = 8, value = 1,step = 0.5)	
						),
					column(3,
					textInput("color.grp1", label = ("Color for group 1"), value = colors[1]),
					textInput("color.grp2", label = ("Color for group 2"), value = colors[2]),
					textInput("color.grp3", label = ("Color for group 3"), value = colors[3])),
					column(3,
					numericInput("fdr.grp1", label = ("FDR cut-off for group 1"), 
						value = 0.05,step=0.01,min=0,max=1),
					numericInput("fdr.grp2", label = ("FDR cut-off for group 2"), 
						value = 0.05,step=0.01,min=0,max=1),
					numericInput("fdr.grp3", label = ("FDR cut-off for group 3"), 
						value = 0.05,step=0.01,min=0,max=1)),
					column(3,
					numericInput("fc.grp1", label = ("Fold-change for group 1"), 
						value = 1,min=1,step=0.5),
					numericInput("fc.grp2", label = ("Fold-change for group 2"), 
						value = 1,min=1,step=0.5),
					numericInput("fc.grp3", label = ("Fold-change for group 3"), 
						value = 1,min=1,step=0.5))			
				),
				hr(),
				br(),
				fluidRow(
					plotOutput(outputId = "plot1")),
				br(),
				fluidRow(
					plotOutput(outputId = "plot2")),
				br(),
				br(),
				fluidRow(
					plotOutput(outputId = "plot3")),
				br(),
				fluidRow(
					plotOutput(outputId = "plot4")),
			    br()
			)
      )
    ),
	tabPanel("> 3 Groups", fluid = TRUE,
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
					selected=names(out.DESeq2$results)[1],
					choices = groups.lst),
				textInput(inputId="tab3.grp1.plot.title", 
					label = ("Group 1 label"), 
					value = names(out.DESeq2$results)[1]),
				selectInput(
					inputId="tab3.grp2.name",
					label = ("Group 2"),
					selected=names(out.DESeq2$results)[2],
					choices = groups.lst),
				textInput(inputId="tab3.grp2.plot.title", 
					label = ("Group 2 label"), 
					value = names(out.DESeq2$results)[2]),
					
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
				actionButton(inputId='export', label="Save plots.pdf")
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
					br(),br(),
					),
					column(2,
					textInput("tab3.color.grp1", label = ("Color for group 1"), value = tab3.colors[1]),
					textInput("tab3.color.grp2", label = ("Color for group 2"), value = tab3.colors[2]),
					tags$div(id = "placeholder-col")),
					column(2,
					numericInput("tab3.fdr.grp1", label = ("FDR cut-off for group 1"), 
						value = 0.05,step=0.01,min=0,max=1),
					numericInput("tab3.fdr.grp2", label = ("FDR cut-off for group 2"), 
						value = 0.05,step=0.01,min=0,max=1),
					tags$div(id = "placeholder-fdr")),
					column(2,
					numericInput("tab3.fc.grp1", label = ("Fold-change for group 1"), 
						value = 1,min=1,step=0.5),
					numericInput("tab3.fc.grp2", label = ("Fold-change for group 2"), 
						value = 1,min=1,step=0.5),
						tags$div(id = "placeholder-fc")),
				column(2,
				       numericInput("tab3.meanDiff.grp1", label = ("Counts difference for group 1"), 
				                    value = 0,min=0,step=1),
				       numericInput("tab3.meanDiff.grp2", label = ("Counts difference for group 2"), 
				                    value = 0,min=0,step=1),
				       tags$div(id = "placeholder-meanDiff"))),
				br(),
				br(),
				# This is the dynamic UI for the plots
				fluidRow(column(12,
				uiOutput("tab3.plots"))),
				br(),
				br(),
				fluidRow(
				column(12,plotOutput(outputId = "tab3.plot1"))
				),
				br(),
				br(),
				fluidRow(
				 column(12,plotOutput(outputId = "tab3.plot2"))
				),
				br(),
				br(),
				conditionalPanel(
				 condition= "input['gen.go'] >= 1",
				 fluidRow(
				  column(12,plotOutput(outputId = "tab3.plot3.1",height="auto"))
				 ),
				 br(),
				fluidRow(
				 column(12,plotOutput(outputId = "tab3.plot3.2",height="auto"))
				),
				br(),
				fluidRow(
				 column(12,plotOutput(outputId = "tab3.plot3.3",height="auto"))
				),
				br(),
				fluidRow(
				 column(12,plotOutput(outputId = "tab3.plot3.4",height="auto"))
				)
				), # conditionalpanel
				textOutput("tab3.text")
			)
		)
	)		
  )
)  

  
  # server
server <- function(input, output,session) {

#write.table(isolate(session$clientData$url_port),file="/mnt/outputs/port.txt",
#	quote=FALSE,row.names=FALSE,col.names=FALSE)

#### processing > 3 groups

 # vals will contain all plots and table grobs
 vals.plot <- reactiveValues(venn.up1=NULL,venn.up2=NULL,venn.dwn1=NULL,
                             venn.dwn2=NULL,hm.up=NULL,hm.dwn=NULL,
                             upset.up=NULL,upset.dwn=NULL,msig.mf=NULL,
                             msig.bp=NULL,msig.cc=NULL,msig.curate=NULL)
 
value <- reactiveVal(2)
inserted <- c()		
inserted.col <- c()
inserted.fdr <- c()
inserted.fc <- c()
inserted.meanDiff <- c()
# reactive value to plot gene enrichment
react.val.gen.go <- reactiveVal(0)
# adjust reactive value as gen.go button is pressed
observeEvent(input$gen.go, {
 gen.go <- reactiveVal(1)
 react.val.gen.go(gen.go)
 print("gen.go=")
 print(gen.go)
})

observeEvent(input$insert_set, {
    btn <- value() +1
	value(btn)
    id <- paste0("txt1_", btn)
    insertUI(
      selector = "#placeholder",
      ui = tags$div(
          selectInput(inputId=paste0("tab3.grp",btn,".name"),
			label=(paste0("Group ",btn)), 
			selected=names(out.DESeq2$results)[btn],
			choices=groups.lst),
          textInput(inputId=paste0("tab3.grp",btn,".plot.title"),
			label=(paste0("Group ",btn," label")), 
			value = names(out.DESeq2$results)[btn]),
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
	               label = (paste0("Counts difference for group ",btn)), 
	               value = 0,min=0,step=1),
	  id = id.meanDiff
	 )
	)
	inserted.meanDiff <<- c(inserted.meanDiff, id.meanDiff)
  })
 
   observeEvent(input$remove_set, {
	if(value()>2){
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
		  paste0("tab3.grp",length(inserted)+2 ,".name"),
		  NULL,
		  NA
		)
		updateTextInput(
		  session,
		  paste0("tab3.grp", length(inserted)+2,".plot.title"),
		  NULL,
		  NA
		)
		updateTextInput(
		  session,
		  paste0("tab3.color.grp",length(inserted.col)+2),
		  NULL,
		  NA
		)
		updateNumericInput(
		  session,
		  paste0("tab3.fdr.grp",length(inserted.fdr)+2),
		  NULL,
		  NA
		)
		updateNumericInput(
		  session,
		  paste0("tab3.fc.grp",length(inserted.fc)+2),
		  NULL,
		  NA
		)
   updateNumericInput(
    session,
    paste0("tab3.meanDiff.grp",length(inserted.meanDiff)+2),
    NULL,
    NA
   )

		inserted <<- inserted[-length(inserted)]
		inserted.col <<- inserted.col[-length(inserted.col)]
		inserted.fdr <<- inserted.fdr[-length(inserted.fdr)]
		inserted.fc <<- inserted.fc[-length(inserted.fc)]
		inserted.meanDiff <<- inserted.meanDiff[-length(inserted.meanDiff)]
		}else{ btn <- 2
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
		print(head(data))
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
		print(head(data))
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
	
	react.tab3.result <- reactive({
		grp.name=react.tab3.grp.name()
		grp.plot.title=react.tab3.grp.plot.title()
		
		idx=match(grp.name,names(out.DESeq2$results))
		results=out.DESeq2$results[idx]
	})
	
	react.tab3.normCts <- reactive({
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
			# grp.up[[i]]=results[[grp.name[i]]][results[[grp.name[i]]]$padj<=fdr.co[i] &
			# results[[grp.name[i]]]$log2FoldChange>(log2(fc.cutoff[i])) &
			# !is.na(results[[grp.name[i]]]$padj),]
			# cat("groupname=",grp.name[i],"\n")
			# cat("fdr.co=",fdr.co,"\n")
			# cat("dim(grp.up[[i]])=",dim(grp.up[[i]]),"\n")
			# cat("!isna.padj=",sum(!is.na(results[[grp.name[i]]]$padj)),"\n")
			
		 # filter out by mean difference
		 grps=strsplit(grp.name[i],"_vs_")
		 treat.grp=grps[[1]][1]
		 ref.grp=grps[[1]][2]
		 treat.mean=rowMeans(normCts[,grep(treat.grp,colnames(normCts))])
		 ref.mean=rowMeans(normCts[,grep(ref.grp,colnames(normCts))])
		 diff.mean=treat.mean-ref.mean
		 grp.up[[i]]=results[[grp.name[i]]][results[[grp.name[i]]]$padj<=fdr.co[i] &
		 results[[grp.name[i]]]$log2FoldChange>(log2(fc.cutoff[i])) &
		 !is.na(results[[grp.name[i]]]$padj)& diff.mean >= meanDiff.cutoff[i],]
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
		}
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
	
  # Insert the right number of plot output objects into the web page
  output$tab3.plots <- renderUI({
	grp.name=react.tab3.grp.name()
	
	plot_output_list<-list()
	if(length(grp.name)<8){
    plot_output_list[[1]] <- plotOutput("tab3.venn.up")
    plot_output_list[[2]] <- plotOutput("tab3.venn.dwn")
    plot_output_list[[3]] <- plotOutput("tab3.hm")
	}else{
	    plot_output_list[[1]] <- plotOutput("tab3.hm")
	}
    
    # Convert the list to a tagList - this is necessary for the list of items
    # to display properly.
    do.call(tagList, plot_output_list)
  })
  
  # # Call renderPlot for each one. Plots are only actually generated when they
  # # are visible on the web page.

  output[["tab3.venn.up"]] <- renderPlot({
		venn.opts=react.tab3.venn.opts()
		grp.name=react.tab3.grp.name()
		colors=react.tab3.color()
		cex=react.tab3.venn.cex()
		pad=react.tab3.venn.pad()
		set.seed(1)
		s4 <- tab3.s4.up()
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
  data=tab3.data.up()
  vals.plot$hm.up <- heatmap_sigf_overlap(data,"Up-regulated genes")
  data=tab3.data.dwn()
  vals.plot$hm.dwn<- heatmap_sigf_overlap(data,"Down-regulated genes")
  grid.arrange(vals.plot$hm.up,vals.plot$hm.dwn,ncol=2,padding=unit(pad,"line"),
               bottom="",right="",left="",top="")
      })
	
	output$tab3.plot1 <- renderPlot({
		s4 <- tab3.s4.up()
		grp.name <- react.tab3.grp.name()
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

	output$tab3.plot3.1 <- renderPlot({
	 enrich.pval.co <- react.tab3.enrich.pval.co()
	 compare.df <- tab3.compare.df()
	 grp.plot.title <- react.tab3.grp.plot.title()
	 # Using clusterProfiler to perform hypergeometric test on msigdb signatures
	 msigdb.species <- react.tab3.msigdb.species()
	 msig.gene.set <- msigdbr(species = msigdb.species, category = "C5",subcategory = "MF") %>%
	  dplyr::select(gs_name, gene_symbol)
	 msig.name ="MSigDB GO Molecular Function"
	 msig.gene.set$gs_name = gsub("GOMF_","",msig.gene.set$gs_name)
	 msig.gene.set$gs_name = gsub("_"," ",msig.gene.set$gs_name)
	 msig.gene.set$gs_name = tolower(msig.gene.set$gs_name)
	 msig.gene.set$gs_name = str_to_title(msig.gene.set$gs_name)
	 msig.gene.set$gs_name = str_replace_all(msig.gene.set$gs_name, 
	                                         c("Mrna" = "mRNA", "Dna" = "DNA", "Rna" = "RNA", 
	                                           "Trna" = "tRNA", "Mirna" = "miRNA", "Rrna" = "rRNA",
	                                           "Atp" = "ATP","Adp" = "ADP","Snorna" = "snoRNA",
	                                           "Lncrna" = "lncRNA","Snrna" = "snRNA"))
	 formula_res <- compareCluster(SYMBOL~group1+group2, data=compare.df, fun="enricher",
	                               TERM2GENE=msig.gene.set,pvalueCutoff=enrich.pval.co,
	                               pAdjustMethod="BH")
	 
	 # re-arrange datasets using factor
	 # and do pathway analysis using up- and down-regulated genes separately
	 vals.plot$msig.mf <- dotplot(formula_res,x=~factor(group1),font.size=14,title=msig.name) + 
	  facet_grid(~group2) +
	  scale_y_discrete(labels=function(x) str_wrap(x, width=40)) +
	  scale_x_discrete(labels=function(x) str_wrap(x,width=10)) +
	  scale_color_distiller(palette = 'Blues')
	 vals.plot$msig.mf
	},height=700)
	
	output$tab3.plot3.2 <- renderPlot({
	 enrich.pval.co <- react.tab3.enrich.pval.co()
	 compare.df <- tab3.compare.df()
	 grp.plot.title <- react.tab3.grp.plot.title()
	 # Using clusterProfiler to perform hypergeometric test on msigdb signatures
	 msigdb.species <- react.tab3.msigdb.species()
	 msig.gene.set <- msigdbr(species = msigdb.species, category = "C5",subcategory = "BP") %>%
	  dplyr::select(gs_name, gene_symbol)
	 msig.name ="MSigDB GO Biological Process"
	 msig.gene.set$gs_name = gsub("GOBP_","",msig.gene.set$gs_name)
	 msig.gene.set$gs_name = gsub("_"," ",msig.gene.set$gs_name)
	 msig.gene.set$gs_name = tolower(msig.gene.set$gs_name)
	 msig.gene.set$gs_name = str_to_title(msig.gene.set$gs_name)
	 msig.gene.set$gs_name = str_replace_all(msig.gene.set$gs_name, 
	                                         c("Mrna" = "mRNA", "Dna" = "DNA", "Rna" = "RNA", 
	                                           "Trna" = "tRNA", "Mirna" = "miRNA", "Rrna" = "rRNA",
	                                           "Atp" = "ATP","Adp" = "ADP","Snorna" = "snoRNA",
	                                           "Lncrna" = "lncRNA","Snrna" = "snRNA"))
	 formula_res <- compareCluster(SYMBOL~group1+group2, data=compare.df, fun="enricher",
	                               TERM2GENE=msig.gene.set,pvalueCutoff=enrich.pval.co,
	                               pAdjustMethod="BH")
	 vals.plot$msig.bp <- dotplot(formula_res,x=~factor(group1),font.size=14,title=msig.name) + 
	  facet_grid(~group2) +
	  scale_y_discrete(labels=function(x) str_wrap(x, width=40)) +
	  scale_x_discrete(labels=function(x) str_wrap(x,width=10)) +
	  scale_color_distiller(palette = 'Blues')
	 vals.plot$msig.bp
	},height=700)
	 
	output$tab3.plot3.3 <- renderPlot({
	 enrich.pval.co <- react.tab3.enrich.pval.co()
	 compare.df <- tab3.compare.df()
	 grp.plot.title <- react.tab3.grp.plot.title()
	 # Using clusterProfiler to perform hypergeometric test on msigdb signatures
	 msigdb.species <- react.tab3.msigdb.species()
	 msig.gene.set <- msigdbr(species = msigdb.species, category = "C5",subcategory = "CC") %>%
	  dplyr::select(gs_name, gene_symbol)
	 msig.name ="MSigDB GO Cellular Component"
	 msig.gene.set$gs_name = gsub("GOCC_","",msig.gene.set$gs_name)
	 msig.gene.set$gs_name = gsub("_"," ",msig.gene.set$gs_name)
	 msig.gene.set$gs_name = tolower(msig.gene.set$gs_name)
	 msig.gene.set$gs_name = str_to_title(msig.gene.set$gs_name)
	 msig.gene.set$gs_name = str_replace_all(msig.gene.set$gs_name, 
	                                         c("Mrna" = "mRNA", "Dna" = "DNA", "Rna" = "RNA", 
	                                           "Trna" = "tRNA", "Mirna" = "miRNA", "Rrna" = "rRNA",
	                                           "Atp" = "ATP","Adp" = "ADP","Snorna" = "snoRNA",
	                                           "Lncrna" = "lncRNA","Snrna" = "snRNA"))
	 formula_res <- compareCluster(SYMBOL~group1+group2, data=compare.df, fun="enricher",
	                               TERM2GENE=msig.gene.set,pvalueCutoff=enrich.pval.co,
	                               pAdjustMethod="BH")
	 vals.plot$msig.cc <- dotplot(formula_res,x=~factor(group1),font.size=14,title=msig.name) + 
	  facet_grid(~group2) +
	  scale_y_discrete(labels=function(x) str_wrap(x, width=40)) +
	  scale_x_discrete(labels=function(x) str_wrap(x,width=10)) +
	  scale_color_distiller(palette = 'Blues')
	 vals.plot$msig.cc
	},height=700)
	 
	 output$tab3.plot3.4 <- renderPlot({
	  # reset reactive value to plot gene enrichment
	  react.val.gen.go <- reactiveVal(0)
	  enrich.pval.co <- react.tab3.enrich.pval.co()
	  compare.df <- tab3.compare.df()
	  grp.plot.title <- react.tab3.grp.plot.title()
	  # Using clusterProfiler to perform hypergeometric test on msigdb signatures
	  msigdb.species <- react.tab3.msigdb.species()
	  msig.gene.set <- msigdbr(species = msigdb.species, category = "C2") %>%
	   dplyr::select(gs_name, gene_symbol)
	  msig.name ="MSigDB Curated Gene Sets"
	  msig.gene.set$gs_name = gsub("_"," ",msig.gene.set$gs_name)
	  msig.gene.set$gs_name = tolower(msig.gene.set$gs_name)
	  msig.gene.set$gs_name = str_to_title(msig.gene.set$gs_name)
	  msig.gene.set$gs_name = str_replace_all(msig.gene.set$gs_name, 
	                                          c("Mrna" = "mRNA", "Dna" = "DNA", "Rna" = "RNA", 
	                                            "Trna" = "tRNA", "Mirna" = "miRNA", "Rrna" = "rRNA",
	                                            "Atp" = "ATP","Adp" = "ADP","Snorna" = "snoRNA",
	                                            "Lncrna" = "lncRNA","Snrna" = "snRNA"))
	  formula_res <- compareCluster(SYMBOL~group1+group2, data=compare.df, fun="enricher",
	                                TERM2GENE=msig.gene.set,
	                                pvalueCutoff=enrich.pval.co,pAdjustMethod="BH")
	  vals.plot$msig.curate <- dotplot(formula_res,x=~factor(group1),font.size=14,title=msig.name) + facet_grid(~group2) +
	   scale_y_discrete(labels=function(x) str_wrap(x, width=40)) +
	   scale_x_discrete(labels=function(x) str_wrap(x,width=10)) +
	   scale_color_distiller(palette = 'Blues')
	  vals.plot$msig.curate
	 },height=700)

	 ## clicking on the export button will generate a pdf file 
	 ## containing all grobs
	 observeEvent( input$export,{
	  if(input$export==0) return()
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
	   if(input$gen.go!=0){
	    gridExtra::grid.arrange(vals.plot$msig.mf,vp=viewport(width=0.8, height=0.8))
	    gridExtra::grid.arrange(vals.plot$msig.bp,vp=viewport(width=0.8, height=0.8))
	    gridExtra::grid.arrange(vals.plot$msig.cc,vp=viewport(width=0.8, height=0.8))
	    gridExtra::grid.arrange(vals.plot$msig.curate,vp=viewport(width=0.8, height=0.8))
	   } #if(input$gen.go!=0)
	   dev.off()
	  })#withProgress
	 }) # observeEvent
	
	output$tab3.text <- renderPrint({ grp.up = tab3.grp.up(); length(grp.up)})
}
# app
shinyApp(ui = ui, server = server)
