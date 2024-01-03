library(shiny)
library(shinyjs)

css <- "
.nav li a.disabled {
background-color: #aaa !important;
color: #333 !important;
cursor: not-allowed !important;
border-color: #aaa !important;
}"

test=as.list(c("a","b"))
ui <- shinyUI(fluidPage(
 shinyjs::useShinyjs(),
 shinyjs::inlineCSS(css),
 navbarPage("Test",id="navbarPage",
            tabPanel("FirstTab", id = "first_tab",
                     sidebarLayout(
                      sidebarPanel(),
                      mainPanel()
                     )
            ),
            tabPanel("Secondtab", id = "second_tab",
                     sidebarLayout(
                      sidebarPanel(),
                      mainPanel(
                       selectInput(
                        inputId="tab0.grp1.name",
                        label = ("Group 1"),
                        selected=test[1],
                        choices = test),
                      )
                     )
            ),
            tabPanel("Third tab", id = "third_tab",
                     sidebarLayout(
                      sidebarPanel(),
                      mainPanel()
                     )
            )
 )
))

server <- shinyServer(function(input, output, session) {
 # disable tabs Exposure, Covariate, and Construct on page load
 #shinyjs::disable(selector = '.navbar-nav a[data-value="Secondtab"')
 shinyjs::disable(selector = '.navbar-nav a[data-value="Third tab"')
})

# Run the application
shinyApp(ui = ui, server = server)