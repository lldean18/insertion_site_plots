# install libraries if required
list.of.packages <- c("tidyverse", "data.table", "patchwork", "ggpattern", "shiny", "bslib")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)>0) {install.packages(new.packages)}
rm(list.of.packages, new.packages)

# load libraries
library(shiny)
library(bslib)
library(tidyverse)
library(patchwork)
library(ggpattern)


# Define UI for application that draws a histogram
ui <- page_sidebar(
  title = "Insertion site plot",
  sidebar = sidebar("Plotting options",
                    radioButtons(
                      inputId = "BINARY",
                      label = "Insertion bar display:",
                      choices = c("Scaled", "Binary")
                      ),
                    numericInput(
                      "MAX_PLOT_HEIGHT",
                      label = "To restrict the height of the y axis, enter a maximum plot height:",
                      value = NULL
                      ),
                    radioButtons(
                      inputId = "SEP_FWD_REV",
                      label = "Do you wish to separate features by strand?",
                      choices = c("no","yes")
                      ),
                    textInput(
                      "HIGHLIGHT_GENE",
                      label = "Enter the name of the gene to highlight:",
                      value = ""
                      ),
                    numericInput(
                      "PLOT_HEIGHT",
                      label = "Set the output plot height (cm)",
                      value = 10
                    ),
                    numericInput(
                      "PLOT_WIDTH",
                      label = "Set the output plot width (cm)",
                      value = 30
                    ),
                    selectInput("select", "Select option",
                                choices = list("Choice 1" = 1, "Choice 2" = 2, "Choice 3" = 3),
                                selected = 1
                    )),
  card(
    card_header("Upload your genome assembly here:"),
    fileInput("file", label = "must be in fasta format (.fa .fasta .fna)")
    ),
  card(
    card_header("Upload your annotation file here:"),
    fileInput("file", label = "must be in gff format (.gff)")
    ),
  card(
    card_header("Upload your insertion site data here:"),
    fileInput("file", label = "the .insert_site_plot.gz file ouput from bacteria_tradis")
  )
)

# Define server logic required to draw a histogram
server <- function(input, output) {

    # output$distPlot <- renderPlot({
    #     # generate bins based on input$bins from ui.R
    #     x    <- faithful[, 2]
    #     bins <- seq(min(x), max(x), length.out = input$bins + 1)
    # 
    #     # draw the histogram with the specified number of bins
    #     hist(x, breaks = bins, col = 'darkgray', border = 'white',
    #          xlab = 'Waiting time to next eruption (in mins)',
    #          main = 'Histogram of waiting times')
    # })
}

# Run the application 
shinyApp(ui = ui, server = server)
