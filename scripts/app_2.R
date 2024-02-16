library(shiny)
library(shinyBS)

# Define UI for application
ui <- fluidPage(
 radioButtons("entree", "Select an entree from the choices below",
              choiceNames = list(
               list("Surf + Turf",
                    bsButton("surf-info", label = "", icon = icon("info", lib = "font-awesome"), style = "default", size = "extra-small")), 
               list("Sushi Assortment",
                    bsButton("sushi-info", label = "", icon = icon("info", lib = "font-awesome"), style = "default", size = "extra-small")), 
               list("Pho",
                    bsButton("pho-info", label = "", icon = icon("info", lib = "font-awesome"), style = "default", size = "extra-small")), 
               list("Breakfast for dinner",
                    bsButton("break-info", label = "", icon = icon("info", lib = "font-awesome"), style = "default", size = "extra-small"))),
              choiceValues = list("surf", "sushi", "pho", "break")),
 
 bsPopover(
  id = "surf-info",
  title = "More information",
  content = HTML(paste0(
   "Ribeye steak, grilled jumbo shrimp, butter roasted potato medley, grilled asparagus."
  )),
  placement = "right",
  trigger = "hover",
  options = list(container = "body")
 ),
 bsPopover(
  id = "pho-info",
  title = "More information",
  content = HTML(paste0(
   "Rice noodles, braised oxtail, and anise-coriander scented broth."
  )),
  placement = "right",
  trigger = "hover",
  options = list(container = "body")
 ),
 bsPopover(
  id = "sushi-info",
  title = "More information",
  content = HTML(paste0(
   "20 piece assortment of nigiri and sashimi with tuna, salmon, eel, and yellowtail."
  )),
  placement = "right",
  trigger = "hover",
  options = list(container = "body")
 ),
 bsPopover(
  id = "break-info",
  title = "More information",
  content = HTML(paste0(
   "Sarahs Saturday special of chocolate Kodiak pancakes, soft scrambled eggs, and an assortment of seasonal fruit."
  )),
  placement = "right",
  trigger = "hover",
  options = list(container = "body")
 )
 
)


server <- function(input, output) {
 
}

# Run the application 
shinyApp(ui = ui, server = server)