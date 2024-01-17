addResourcePath("diff_report", "/mnt/outputs/diff_analysis_rslt")

ui <- fluidPage(
 tabPanel("Outputs",id="outputtab",fluid=TRUE,
 mainPanel(
 htmlOutput("diff_report")
)))

server <- function(input,output){
 output$diff_report <- renderUI({
  tags$iframe(seamless="seamless", 
              src= "diff_report/RNA-seq_differential_analysis_report.html",
              width=800, 
              height=800)
 })
}

shinyApp(ui, server)