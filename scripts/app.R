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
library(venn)
library(UpSetR)

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

# List of choices for selectInput
groups.lst <- as.list(groups$names)
# Name it
names(groups.lst) <- groups$var
	
# default checkboxes for venn.opts	
venn.opts=data.frame(
	lbl= c("Numbers","Percentages","Labels","Legend")
)
venn.opts.lst=as.list(venn.opts$lbl)
names(venn.opts.lst)=venn.opts$lbl

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

heatmap_sigf_overlap <- function(data){
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
		theme(legend.key.width=unit(0.7,"cm"))
	
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
              )
			  ),
			mainPanel(
				fluidRow(style = "background-color:#F5F5F5",column(3,
					checkboxGroupInput("tab3.venn.opts", label = ("Display on Venns"), 
						choices = venn.opts.lst,
							selected = venn.opts$lbl[1:3]),
					    
					# slider bar font size
					sliderInput("tab3.venn.cex", label = ("Venns' font size"), min = 0, 
						max = 6, value = 1,step = 0.1),
					# slider bar pad
					sliderInput("tab3.venn.pad", label = ("Venns' circle size"), min = 0, 
						max = 8, value = 1,step = 0.5)	
					),
					column(3,
					textInput("tab3.color.grp1", label = ("Color for group 1"), value = tab3.colors[1]),
					textInput("tab3.color.grp2", label = ("Color for group 2"), value = tab3.colors[2]),
					tags$div(id = "placeholder-col")),
					column(3,
					numericInput("tab3.fdr.grp1", label = ("FDR cut-off for group 1"), 
						value = 0.05,step=0.01,min=0,max=1),
					numericInput("tab3.fdr.grp2", label = ("FDR cut-off for group 2"), 
						value = 0.05,step=0.01,min=0,max=1),
					tags$div(id = "placeholder-fdr")),
					column(3,
					numericInput("tab3.fc.grp1", label = ("Fold-change for group 1"), 
						value = 1,min=1,step=0.5),
					numericInput("tab3.fc.grp2", label = ("Fold-change for group 2"), 
						value = 1,min=1,step=0.5),
						tags$div(id = "placeholder-fc"))),
				br(),
				br(),
				# This is the dynamic UI for the plots
				fluidRow(column(12,
				uiOutput("tab3.plots"))),
				br(),
				br(),
				fluidRow(
				column(12,plotOutput(outputId = "tab3.plot1"))
				)
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
value <- reactiveVal(2)
inserted <- c()		
inserted.col <- c()
inserted.fdr <- c()
inserted.fc <- c()

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

		inserted <<- inserted[-length(inserted)]
		inserted.col <<- inserted.col[-length(inserted.col)]
		inserted.fdr <<- inserted.fdr[-length(inserted.fdr)]
		inserted.fc <<- inserted.fc[-length(inserted.fc)]
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
		gs.RNASeq <- min(sapply(results,function(x) sum(!is.na(x$padj))))
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
			grp.up[[i]]=results[[grp.name[i]]][results[[grp.name[i]]]$padj<fdr.co[i] &
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
			grp.dwn[[i]]=results[[grp.name[i]]][results[[grp.name[i]]]$padj<fdr.co[i] &
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
		gs.RNASeq <- min(sapply(tab0.results,function(x) sum(!is.na(x$padj))))
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
			grp.up[[i]]=results[[grp.name[i]]][results[[grp.name[i]]]$padj<fdr.co[i] &
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
			grp.dwn[[i]]=results[[grp.name[i]]][results[[grp.name[i]]]$padj<fdr.co[i] &
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
	react.tab3.venn.opts <- reactive({input$tab3.venn.opts})
	
	react.tab3.venn.cex <- reactive({input$tab3.venn.cex})
	
	react.tab3.venn.pad <- reactive({input$tab3.venn.pad})
	
	react.tab3.result <- reactive({
		grp.name=react.tab3.grp.name()
		grp.plot.title=react.tab3.grp.plot.title()
		
		idx=match(grp.name,names(out.DESeq2$results))
		results=out.DESeq2$results[idx]
	})
	tab3.gs.RNASeq <- reactive({
		grp.name=react.tab3.grp.name()
		grp.plot.title=react.tab3.grp.plot.title()
		results=react.tab3.result()
		
		gs.RNASeq <- min(sapply(results,function(x) sum(!is.na(x$padj))))
	})
	tab3.grp.up <- reactive ({
		grp.name=react.tab3.grp.name()
		grp.plot.title=react.tab3.grp.plot.title()
		results=react.tab3.result()
		fdr.co <- react.tab3.fdr.cutoff()
		fc.cutoff <- react.tab3.fc.cutoff()
		
		# up-regulated genes
		grp.up=list()
		for(i in 1:length(grp.name)){
			grp.up[[i]]=results[[grp.name[i]]][results[[grp.name[i]]]$padj<fdr.co[i] &
			results[[grp.name[i]]]$log2FoldChange>(log2(fc.cutoff[i])) & 
			!is.na(results[[grp.name[i]]]$padj),]
		}
		grp.up
	})
	tab3.grp.dwn <- reactive({
		grp.name=react.tab3.grp.name()
		grp.plot.title=react.tab3.grp.plot.title()
		results=react.tab3.result()
		fdr.co <- react.tab3.fdr.cutoff()
		fc.cutoff <- react.tab3.fc.cutoff()

		# down-regulated genes
		grp.dwn=list()
		for(i in 1:length(grp.name)){
			grp.dwn[[i]]=results[[grp.name[i]]][results[[grp.name[i]]]$padj<fdr.co[i] &
			results[[grp.name[i]]]$log2FoldChange<(-log2(fc.cutoff[i])) & 
			!is.na(results[[grp.name[i]]]$padj),]
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
	tab3.gom.obj <- reactive({
		genes.lists <- tab3.genes.lists()
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
	tab3.s4 <- reactive({
		# input for euler venn
		s4 <-tab3.genes.lists()
		grp.plot.title <-react.tab3.grp.plot.title()
		
		names(s4) <- grp.plot.title
		s4
	})
  # Insert the right number of plot output objects into the web page
  output$tab3.plots <- renderUI({
	grp.name=react.tab3.grp.name()
	
	plot_output_list<-list()
	if(length(grp.name)<8){
    plot_output_list[[1]] <- plotOutput("tab3.venn")
    plot_output_list[[2]] <- plotOutput("tab3.hm")
	}else{
	    plot_output_list[[1]] <- plotOutput("tab3.hm")
	}
    
    # Convert the list to a tagList - this is necessary for the list of items
    # to display properly.
    do.call(tagList, plot_output_list)
  })
  
  # # Call renderPlot for each one. Plots are only actually generated when they
  # # are visible on the web page.

    output[["tab3.venn"]] <- renderPlot({
		venn.opts=react.tab3.venn.opts()
		grp.name=react.tab3.grp.name()
		colors=react.tab3.color()
		cex=react.tab3.venn.cex()
		pad=react.tab3.venn.pad()
		s4 <- tab3.s4()
		
		venn1<-grid.arrange(plot_euler(s4=s4,colors=colors,cex=cex,
			venn.opts=venn.opts,title="All regulated genes"),
			padding=unit(pad,"line"),top="",bottom="",right="",left="")
		layout(matrix(1:2, 1, byrow = TRUE))
		venn(s4,zcolor=colors,opacity=.8,box=FALSE,ilcs = 0.8, sncs = 1)
		frame()
		vps <- baseViewports()
		pushViewport(vps$inner, vps$figure, vps$plot)
		grid.draw(venn1)
		popViewport(2)
      })
      
    output[["tab3.hm"]] <- renderPlot({
        data=tab3.data()
		
		heatmap_sigf_overlap(data)
      })
	
	output$tab3.plot1 <- renderPlot({
		s4 <- tab3.s4()
		grp.name <- react.tab3.grp.name()
		
		upset(fromList(s4), order.by = "freq",nsets=length(grp.name), 
			mainbar.y.label = "Genes Overlaps", sets.x.label = "Significant Genes", 
			 mb.ratio = c(0.5, 0.5),
			#c(intersection size title, intersection size tick labels, set size title, set size tick labels, set names, numbers above bars)
			text.scale = c(1.4, 1.4, 1.4, 1.4, 1.8, 1.5))
	})
	
	output$tab3.text <- renderPrint({ cat(react.tab3.grp.name(),"n",react.tab3.grp.plot.title(),"\n",
	react.tab3.color(),"\n",names(react.tab3.result()),"\n",react.tab3.fdr.cutoff()) })
}
# app
shinyApp(ui = ui, server = server)
