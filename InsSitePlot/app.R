# package for shiny app deployment
# install.packages('rsconnect')
# to deploy the app
# library(rsconnect)
# rsconnect::deployApp('InsSitePlot')

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
                    # input desired chromosome
                    textInput(
                      "chr",
                      label = "Chromosome to plot:",
                      value = "NZ_CP009273.1"
                    ),
                    # input start position to plot
                    numericInput(
                      "st",
                      label = "Start position:",
                      value = 3862000
                    ),
                    # input end position to plot
                    numericInput(
                      "en",
                      label = "End position:",
                      value = 3868000
                    ),
                    # input whether forward and reverse strand features should be separated on the y-axis
                    checkboxInput(
                      inputId = "SEP_FWD_REV",
                      label = "Separate features by strand?",
                      value = F
                    ),
                    # input whether the features should be rainbow colours
                    checkboxInput(
                      inputId = "RAINBOW",
                      label = "Plot rainbow feature colours?",
                      value = F
                    ),
                    # input whether insertion site y-axis should be scaled or binary
                    checkboxInput(
                      inputId = "BINARY",
                      label = "Only show insertion site locations? (not counts)",
                      value = F
                    ),
                    # input the colour for plotting forward insertion bars
                    selectInput("FWD_BAR_COLOUR",
                                label = "Select colour for forward-strand insertion bars", 
                                choices = list("blue" = "blue3", "red" = "red","darkred" = "darkred", "orange" = "orange2",
                                               "green" = "green4", "purple" = "purple4", "yellow" = "yellow2",
                                               "pink" = "deeppink2", "lightpink" = "pink", "black" = "black", "grey" = "grey30"), 
                                selected = "blue"),
                    # input the colour for plotting the reverse insertion bars
                    selectInput("REV_BAR_COLOUR",
                                label = "Select colour for reverse-strand insertion bars", 
                                choices = list("blue" = "blue3", "red" = "red","darkred" = "darkred", "orange" = "orange2",
                                               "green" = "green4", "purple" = "purple4", "yellow" = "yellow2",
                                               "pink" = "deeppink2", "lightpink" = "pink", "black" = "black", "grey" = "grey30"), 
                                selected = "pink"),
                    # input the maximum plot height
                    numericInput(
                      "MAX_PLOT_HEIGHT",
                      label = "To restrict the height of the y-axis, enter a maximum plot height:",
                      value = NULL
                      ),
                    # input the gene to highlight if any
                    textInput(
                      "HIGHLIGHT_GENE",
                      label = "Name of gene to highlight:",
                      value = ""
                      ),
                    # input the output file name
                    textInput(
                      "OUT_FILE",
                      label = "Name of output file:",
                      value = "insertion_site_plot.pdf"
                    ),
                    # input the desired output plot height
                    numericInput(
                      "PLOT_HEIGHT",
                      label = "Set the output plot height (cm)",
                      value = 10
                    ),
                    # input the desired output plot width
                    numericInput(
                      "PLOT_WIDTH",
                      label = "Set the output plot width (cm)",
                      value = 30
                    ),
                    # add save the plot button
                    downloadButton(
                      "PLOT_SAVE",
                      label = "Save the plot"
                      )
                    ),
  layout_columns(fill = F,
    card(
      card_header("Genome assembly upload:"),
      fileInput("fasta_file", label = "must be in fasta format (.fa .fasta .fna)")
    ),
    card(
      card_header("Annotation file upload:"),
      fileInput("gff_file", label = "must be in gff format (.gff)")
    ),
    card(
      card_header("Insertion site data upload:"),
      fileInput("insertions_file", label = "the .insert_site_plot.gz file ouput from bacteria_tradis")
    )
  ),

  plotOutput("insertion_plot")
)



##########################
# setup the server logic #
##########################
server <- function(input, output) {
  
  options(shiny.maxRequestSize=30*1024^2)
  
  # make the insertions data frame from user input files #
  ########################################################
  
  # # read in and prep fasta and inserions file
  # fastaR <- reactive({
  #   req(input$fasta_file)
  #   prep_insertion_sites(read.fasta(input$fasta_file$datapath),
  #                        read.table(input$insertions_file$datapath),
  #                        input$chr, input$st, input$en
  #                        )
  #   
  # })

  insertion_sites_R <- reactive({
    #list required inputs before running
    req(input$fasta_file, input$gff_file, input$insertions_file, input$chr, input$st, input$en)
    #prep insertion files
    prep_insertion_files(prep_insertion_sites(read.fasta(input$fasta_file$datapath),
                                              read.table(input$insertions_file$datapath),
                                              input$chr, input$st, input$en
                                              ),
                         read.gff(input$gff_file$datapath, input$chr, input$st, input$en))
  })
  
  
  
  
  #   insertion_sites_R <- reactive({
  #   #list required inputs before running
  #   req(input$fasta_file, input$gff_file, input$insertions_file)
  #   #prep insertion files
  #   prep_insertion_files(read.fasta(input$fasta_file$datapath),
  #                        read.gff(input$gff_file$datapath, input$chr, input$st, input$en),
  #                        read.table(input$insertions_file$datapath))
  # })
  
  # setup that doesn't depend on user input #
  ###########################################
  # Define arrow dimensions
  arrow_height <- 14
  arrow_head_length <- 0.2  # Fraction of feature length
  arrow_head_height <- 1.3   # Relative to arrow half-height
  
  
  # THE PLOT #
  ############
  plot_to_save <- reactive({
    
    # get the insertion sites from the input files
    insertion_sites <- insertion_sites_R()
    
    # filter data to the region for plotting
    my_region <- dplyr::filter(insertion_sites, chromosome == input$chr & position >= input$st & position <= input$en)
    
    # Summarise the genomic coordinates of each feature ###
    features <- if (all(is.na(my_region$feature_ID))) {
      tibble::tibble(
        feature_name = character(),
        feature_ID = character(),
        feature_type = character(),
        feature_strand = character(),
        feature_phase = character(),
        start = numeric(),
        end = numeric(),
        length = numeric(),
        body_height = numeric(),
        head_height = numeric(),
        head_length = numeric()
      )
    } else {
      my_region %>%
      dplyr::group_by(feature_name, feature_ID, feature_type, feature_strand, feature_phase) %>%
      dplyr::summarise(start = min(position), end = max(position), .groups = "drop") %>%
      dplyr::mutate(length = end - start, body_height = arrow_height * 0.6,
             head_height = arrow_height, head_length = length * arrow_head_length)
    }
    
    # ensure that feature ID is unique
    n <- nrow(dplyr::filter(features, duplicated(features$feature_ID)))
    features$feature_ID <- ifelse(duplicated(features$feature_ID),
                                  paste0(features$feature_ID, "_", seq_along(n)),
                                  features$feature_ID)
    rm(n)
    
    # remove the "cds-" string from the beginning of feature names
    features$feature_name <- sub("cds-", "", features$feature_name)
    
    # Determine whether features occur on both strands
    both_strands <- any(features$feature_strand == "+", na.rm = TRUE) &&
      any(features$feature_strand == "-", na.rm = TRUE)
    
    # Place features just below zero; separate strands only if both are present
    features <- features %>%
      dplyr::mutate(y_offset = if (both_strands && input$SEP_FWD_REV) {
        ifelse(feature_strand == "+", -arrow_height * 0.7, -arrow_height * 2.1)
      } else {-arrow_height * 0.7})
    
    # add feature fill colour
    features$fill_colour <- ifelse(tolower(features$feature_name) == tolower(input$HIGHLIGHT_GENE),
                                   {if (input$RAINBOW) "grey70" else "lightpink"},
                                   {if (input$RAINBOW) as.character(features$feature_ID) else "grey70"})
    
    # Prepare forward and reverse insertion bars
    if (input$BINARY) {
      fwd_bars <- my_region %>%
        dplyr::filter(fwd_insertions > 0) %>%
        dplyr::mutate(ymin = 0, ymax = 10)
      rev_bars <- my_region %>%
        dplyr::filter(rev_insertions > 0) %>%
        dplyr::mutate(ymin = 0, ymax = 10)
    } else {
      fwd_bars <- my_region %>%
        dplyr::mutate(ymin = 0, ymax = 0 + fwd_insertions)
      rev_bars <- my_region %>%
        dplyr::mutate(ymin = 0, ymax = 0 + rev_insertions)
    }
    
    # Calculate arrow polygon coordinates for CDS features and pseudo genes
    arrow_polygons <- if (all(is.na(features$feature_ID))) {
      tibble::tibble(
        feature_name = character(),
        feature_ID = character(),
        feature_type = character(),
        feature_strand = character(),
        fill_colour = character(),
        x = numeric(),
        y = numeric()
      )
    } else {
      features %>%
      dplyr::filter(feature_type == "CDS" | feature_type == "pseudogene") %>%
      rowwise() %>%
      do(data.frame(feature_name = .$feature_name,
                    feature_ID = .$feature_ID,
                    feature_type = .$feature_type,
                    feature_strand = .$feature_strand,
                    fill_colour = .$fill_colour,
                    x = if (.$feature_strand == "+") {
                      c(.$start, .$end - .$head_length, .$end - .$head_length,
                        .$end, .$end - .$head_length, .$end - .$head_length, .$start)
                    } else {
                      c(.$end, .$start + .$head_length, .$start + .$head_length,
                        .$start, .$start + .$head_length, .$start + .$head_length, .$end)
                    },
                    y = c(.$y_offset - .$body_height, .$y_offset - .$body_height,
                          .$y_offset - .$head_height, .$y_offset,
                          .$y_offset + .$head_height, .$y_offset + .$body_height,
                          .$y_offset + .$body_height)
      )) %>%
      ungroup() %>%
      dplyr::mutate(feature_ID = factor(feature_ID, levels = unique(feature_ID))) ###
    }
    
    # add calculation for body lines for pseudogenes
    body_lines <- if (nrow(dplyr::filter(arrow_polygons, feature_type == "pseudogene")) > 0 ) {
      arrow_polygons %>%
        dplyr::filter(feature_type == "pseudogene") %>%
        dplyr::group_by(feature_ID) %>%  # if feature_ID groups are intended
        dplyr::summarise(
          x_start = if (first(feature_strand) == "+") {
            min(x)
          } else { min(x[x != min(x)]) },
          x_end = if (first(feature_strand) == "+") {
            max(x[x != max(x)])
          } else { max(x) },
          y_max = max(y[y != max(y)]),
          y_min = min(y[y != min(y)]),
          .groups = "drop")
    }
    
    # prepare RNA triangles
    rna_triangles <- features %>%
      dplyr::filter(feature_type %in% c("rRNA", "tRNA", "ncRNA")) %>%
      dplyr::mutate(triangle_width = body_height,
             tip = ifelse(feature_strand == "+", end, start),
             base = ifelse(feature_strand == "+", start, end)) %>%
      rowwise() %>%
      reframe(feature_ID = feature_ID, feature_type = feature_type, x = c(base, base, tip),
              y = c(y_offset + triangle_width, y_offset - triangle_width, y_offset + triangle_width)) %>%
      ungroup()
    
    # set the feature colour order for rainbow arrows
    feature_order <- arrow_polygons %>%
      dplyr::distinct(feature_ID, x) %>%
      dplyr::arrange(x) %>%
      pull(feature_ID) %>%
      unique() %>%
      as.character()
    
    # INSERTIONS PLOT
    p_insertions <- 
      ggplot() +
      # Forward insertion bars
      geom_rect(data = fwd_bars,
                aes(xmin = position - ((input$en - input$st) /1000), xmax = position + ((input$en - input$st) /1000),
                    ymin = ymin, ymax = ymax), fill = input$FWD_BAR_COLOUR, alpha = 0.5) +
      #  Reverse insertion bars
      geom_rect(data = rev_bars,
                aes(xmin = position - ((input$en - input$st) /1000), xmax = position + ((input$en - input$st) /1000),
                    ymin = ymin, ymax = ymax), fill = input$REV_BAR_COLOUR, alpha = 0.5) +
      # Axis labels and theme
      labs(x = "Genomic position (bp)", y = {if (input$BINARY) "Insertions" else "Insertion count"}) +
      {if (!is.null(input$MAX_PLOT_HEIGHT)) coord_cartesian(ylim = c(0, input$MAX_PLOT_HEIGHT))} +
      scale_y_continuous(labels = {if (input$BINARY) NULL else abs}, limits = c(0, NA), expand = expansion(mult = c(0, 0.02))) +
      theme_classic() +
      theme(axis.title = element_text(size = 12),
            plot.title = element_text(size = 14, face = "bold"),
            axis.text.y = {if (input$BINARY) element_blank() else element_text()},
            axis.ticks.y = {if (input$BINARY) element_blank() else element_line()},
            axis.line.y = {if (input$BINARY) element_blank() else element_line()},
            plot.margin = margin(0, 0, 0, 0))
    
    # FEATURES PLOT
    p_features <- ggplot() +
      # CDS feature arrows
      geom_polygon(data = dplyr::filter(arrow_polygons, feature_type == "CDS"),
                   aes(x = x, y = y, group = feature_ID, fill = fill_colour, colour = fill_colour)) +
      # highlighted gene feature arrow if RAINBOW
      {if (input$RAINBOW)
        ggpattern::geom_polygon_pattern(data = dplyr::filter(arrow_polygons, feature_type == "CDS" &
                                             tolower(feature_name) == tolower(input$HIGHLIGHT_GENE)),
                             aes(x=x, y=y, group = feature_ID, fill = "grey"),
                             pattern = "gradient", pattern_orientation = "radial",
                             pattern_fill = "grey", pattern_fill2 = "black") } +
      # Pseudogene feature arrows
      ggpattern::geom_polygon_pattern(data = dplyr::filter(arrow_polygons, feature_type == "pseudogene"),
                           aes(x = x, y = y, group = feature_ID, fill = fill_colour),
                           pattern = "stripe", pattern_angle = 45, pattern_spacing = 0.05,
                           pattern_density = 0.05, colour = "black", pattern_colour = "grey40") +
      # crosses across pseudogene feature arrows
      {if (!is.null(body_lines))
        geom_segment(data = body_lines, aes(x = x_start, xend = x_end, y = y_min, yend = y_max),
                     inherit.aes = FALSE, colour = "black") } +
      {if (!is.null(body_lines))
        geom_segment(data = body_lines, aes(x = x_start, xend = x_end, y = y_max, yend = y_min),
                     inherit.aes = FALSE, colour = "black")} +
      # scale for arrow fill
      {if (input$RAINBOW) {
        scale_fill_manual(values = c(setNames(grDevices::colorRampPalette(
          c("goldenrod1", "sandybrown", "indianred2", "mediumpurple2", "skyblue3")
        )(length(feature_order)),
        feature_order), "grey70" = "grey70"), guide = "none")
      } else {
        scale_fill_identity(guide = "none")
      }} +
      # scale for arrow outline
      {if (input$RAINBOW) {
        scale_colour_manual(values = c(setNames(grDevices::colorRampPalette(
          c("goldenrod1", "sandybrown", "indianred2", "mediumpurple2", "skyblue3")
        )(length(feature_order)),
        feature_order), "grey70" = "black"), guide = "none")
      } else {
        scale_colour_identity(guide = "none")
      }} +
      
      # rRNA, ncRNA and tRNA triangles
      {if (any(rna_triangles$feature_type == "rRNA")) geom_polygon(data = dplyr::filter(rna_triangles, feature_type == "rRNA"),
                                                                   aes(x = x, y = y, group = feature_ID), fill = "grey30", colour = "grey30")} +
      {if (any(rna_triangles$feature_type == "ncRNA")) geom_polygon(data = dplyr::filter(rna_triangles, feature_type == "ncRNA"),
                                                                    aes(x = x, y = y, group = feature_ID), fill = "grey30", colour = "grey30")} +
      {if (any(rna_triangles$feature_type == "tRNA")) geom_polygon(data = dplyr::filter(rna_triangles, feature_type == "tRNA"),
                                                                   aes(x = x, y = y, group = feature_ID), fill = "grey30", colour = "grey30")} +
      # CDS feature IDs centered inside the arrows
      geom_text(data = features %>% filter(feature_type == "CDS"),
                aes(x = (start + end) / 2, y = y_offset, label = feature_name), fontface = "italic", size = 3, colour = "black") +
      # highlighted gene label with box behind
      geom_label(data = features %>% dplyr::filter(feature_type == "CDS" & tolower(feature_name) == tolower(input$HIGHLIGHT_GENE)), 
                 aes(x = (start + end) / 2, y = y_offset, label = feature_name, fill = fill_colour), fontface = "italic",
                 size = 3, colour = "black", linewidth = 0.2, label.padding = unit(0.1, "lines"), border.colour = NA) +
      # pseudogene feature IDs with box behind
      geom_label(data = features %>% dplyr::filter(feature_type == "pseudogene"), 
                 aes(x = (start + end) / 2, y = y_offset, label = feature_name, fill = fill_colour), fontface = "italic",
                 size = 3, colour = "black", linewidth = 0.2, label.padding = unit(0.1, "lines"), border.colour = NA) +
      # Axis labels and theme
      labs(x = "Genomic position (bp)", y = NULL) +
      scale_x_continuous(limits = c(input$st, input$en), expand = expansion(mult = 0)) +
      theme_classic() +
      theme(axis.title = element_text(size = 12),
            plot.title = element_text(size = 14, face = "bold"),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.line.y = element_blank(),
            legend.position = "None",
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.title.x = element_blank(),
            axis.line.x = element_blank(),
            plot.margin = margin(0, 0, 0, 0))
    
    # plot the two together
    p_insertions / p_features +
      plot_layout(heights = c(3, 1))
    
#    # save the plot
#    ggsave(input$OUT_FILE, device = "pdf", units = "cm",
#           height = input$PLOT_HEIGHT, input$width = PLOT_WIDTH)
    
    
  })

  # show the plot in the app
  output$insertion_plot <- renderPlot({
    plot_to_save()
  })
  
  # save the plot when the save plot button is pressed
  output$PLOT_SAVE <- downloadHandler(
    filename = function() {
      input$OUT_FILE
    },
    content = function(file) {
      ggsave(
        filename = file,
        plot = plot_to_save(),
        device = "pdf",
        units = "cm",
        height = input$PLOT_HEIGHT,
        width = input$PLOT_WIDTH
      )
    }
  )

}

# Run the app
shinyApp(ui = ui, server = server)
