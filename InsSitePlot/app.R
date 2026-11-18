# install packages if required
list.of.packages <- c("tidyverse", "data.table", "patchwork", "ggpattern", "shiny", "bslib")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)>0) {install.packages(new.packages)}
rm(list.of.packages, new.packages)

# load packages
library(shiny)
library(bslib)
library(tidyverse)
library(patchwork)
library(ggpattern)

# load APP source code
source("APP_file_prep_for_insertion_site_plot.R")


# Set up the user interface for the app
ui <- page_sidebar(
  title = "Insertion site plot",
  sidebar = sidebar("Plotting options",
                    radioButtons(
                      inputId = "BINARY",
                      label = "Insertion bar display:",
                      choices = c("red", "blue")
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
    fileInput("fasta_file", label = "must be in fasta format (.fa .fasta .fna)"),
    height = "600px"
    ),
  card(
    card_header("Upload your annotation file here:"),
    fileInput("gff_file", label = "must be in gff format (.gff)"),
    height = "600px"
    ),
  card(
    card_header("Upload your insertion site data here:"),
    fileInput("insertions_file", label = "the .insert_site_plot.gz file ouput from bacteria_tradis"),
    height = "700px"
  ),
  verbatimTextOutput("plot_status",
                     placeholder = T),
#  actionButton("plot", "Draw plot"),
  plotOutput("insertion_plot")
)




# setup the server logic
server <- function(input, output) {
  
  # make the insertions data frame from user input files
  insertion_sites <- reactive({
    #list required inputs before running
    req(input$fasta_file, input$gff_file, input$insertions_file)
    #prep insertion files
    prep_insertion_files(read.fasta(input$fasta_file$datapath),
                         read.gff(input$gff_file$datapath),
                         read.table(input$insertions_file$datapath))
  })
  

  output$plot_status <- renderText({
    sites <- insertion_sites()
    paste("Rows:", nrow(sites),
          "| Forward insertions:", sum(sites$fwd_insertions, na.rm = TRUE))
  })
  

  output$insertion_plot <- renderPlot({
    
    sites  <- insertion_sites()
    
    ggplot(sites, aes(fwd_insertions)) +
      geom_histogram(colour = input$BINARY) +
      labs(x = "Forward insertions")
  })
  
  
  
  
  
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

# Run the app
shinyApp(ui = ui, server = server)
