###################################
### PRETTY INSERTION SITE PLOTS ###
###################################

# install libraries if required
list.of.packages <- c("tidyverse", "patchwork")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)>0) {install.packages(new.packages)}
rm(list.of.packages, new.packages)

# load libraries
library(tidyverse)
library(patchwork)

#### load the pre-processed insertion site data ####

#load("insertion_sites.RData") # for photorhabdus
load("insertion_sites_ecoli_annotated.RData") # for ecoli
insertion_sites <- insertion_sites_annotated


#### set user input variables ####

# user input parameters
BINARY = T # if false plot the insertion site count, if true, plot insertion site location only
BINARY = F

RAINBOW = T
RAINBOW = F

SEP_FWD_REV = T
SEP_FWD_REV = F

HIGHLIGHT_GENE = "None"
HIGHLIGHT_GENE = "argG"
HIGHLIGHT_GENE = "nadA"

MAX_PLOT_HEIGHT = 10

# The region you want to plot
chr <- "NZ_CP009273.1"
st <- 3862000
en <- 3868000
# trna
st <- 3310431
en <- 3318517
# rRNA
st <- 4153587
en <- 4159587
# ncRNA
st <- 775139
en <- 779139
# pseudogene
st <- 19000
en <- 21000

# features in ecoli to look at
#dnaA 3875686	3877089
#st <- 3872000
#en <- 3878000
#tRNA 3315431	3315517

#### prep data frames for plotting ####

# filter data to the region for plotting
my_region <- filter(insertion_sites, chromosome == chr & position >= st & position <= en)

# Define arrow dimensions
arrow_height <- 14
arrow_head_length <- 0.2  # Fraction of feature length
arrow_head_height <- 1.3   # Relative to arrow half-height

# # Summarise the genomic coordinates of each feature
features <- my_region %>%
  group_by(feature_name, feature_ID, feature_type, feature_strand, feature_phase) %>%
  summarise(start = min(position), end = max(position), .groups = "drop") %>%
  mutate(length = end - start, body_height = arrow_height * 0.6,
         head_height = arrow_height, head_length = length * arrow_head_length)

# ensure that feature ID is unique
n <- nrow(filter(features, duplicated(features$feature_ID)))
features$feature_ID <- ifelse(duplicated(features$feature_ID),
                              paste0(features$feature_ID, "_", seq_along(n)),
                              features$feature_ID)
rm(n)

# Determine whether features occur on both strands
both_strands <- any(features$feature_strand == "+", na.rm = TRUE) &&
  any(features$feature_strand == "-", na.rm = TRUE)

# Place features just below zero; separate strands only if both are present
features <- features %>%
  mutate(y_offset = if (both_strands && SEP_FWD_REV) {
    ifelse(feature_strand == "+", -arrow_height * 0.7, -arrow_height * 2.1)
  } else {-arrow_height * 0.7})

# Prepare forward and reverse insertion bars
if (BINARY) {
  fwd_bars <- my_region %>%
    filter(fwd_insertions > 0) %>%
    mutate(ymin = 0, ymax = 10)
  rev_bars <- my_region %>%
    filter(rev_insertions > 0) %>%
    mutate(ymin = 0, ymax = 10)
} else {
  fwd_bars <- my_region %>%
    mutate(ymin = 0, ymax = 0 + fwd_insertions)
  rev_bars <- my_region %>%
    mutate(ymin = 0, ymax = 0 + rev_insertions)
}

# # Calculate arrow polygon coordinates for CDS features
# arrow_polygons <- features %>%
#   filter(feature_type == "CDS") %>%
#   rowwise() %>%
#   do(data.frame(
#     feature_name = .$feature_name, feature_ID = .$feature_ID,
#     x = if (.$feature_strand == "+") {
#       c(.$start, .$end - .$head_length, .$end - .$head_length,
#         .$end, .$end - .$head_length, .$end - .$head_length, .$start)
#     } else {
#       c(.$end, .$start + .$head_length, .$start + .$head_length,
#         .$start, .$start + .$head_length, .$start + .$head_length, .$end)
#     },
#     y = c(.$y_offset - .$body_height, .$y_offset - .$body_height,
#           .$y_offset - .$head_height, .$y_offset,
#           .$y_offset + .$head_height, .$y_offset + .$body_height,
#           .$y_offset + .$body_height)
#   )) %>%
#   ungroup() %>%
#   mutate(feature_ID = factor(feature_ID, levels = unique(feature_ID)))
# # add arrow fill colour
# arrow_polygons$fill_colour <- ifelse(tolower(arrow_polygons$feature_name) == tolower(HIGHLIGHT_GENE),
#                                      if (RAINBOW) "grey70" else "lightpink",
#                                      if (RAINBOW) as.character(arrow_polygons$feature_ID) else "grey70")

# new with pseudogenes
# Calculate arrow polygon coordinates for CDS features
arrow_polygons <- features %>%
  filter(feature_type == "CDS" | feature_type == "pseudogene") %>%
  rowwise() %>%
  do(data.frame(feature_name = .$feature_name,
                feature_ID = .$feature_ID,
                feature_type = .$feature_type,
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
  mutate(feature_ID = factor(feature_ID, levels = unique(feature_ID)))

# add arrow fill colour
arrow_polygons$fill_colour <- ifelse(tolower(arrow_polygons$feature_name) == tolower(HIGHLIGHT_GENE),
                                     if (RAINBOW) "grey70" else "lightpink",
                                     if (RAINBOW) as.character(arrow_polygons$feature_ID) else "grey70")


# # prepare RNA triangles
rna_triangles <- features %>%
  filter(feature_type %in% c("rRNA", "tRNA", "ncRNA")) %>%
  mutate(triangle_width = body_height,
         tip = ifelse(feature_strand == "+", end, start),
         base = ifelse(feature_strand == "+", start, end)) %>%
  rowwise() %>%
  reframe(feature_ID = feature_ID, feature_type = feature_type, x = c(base, base, tip),
          y = c(y_offset + triangle_width, y_offset - triangle_width, y_offset + triangle_width)) %>%
  ungroup()

#### draw the plot with patchwork ####
 
p_insertions <- ggplot() +
 # Forward insertion bars
 geom_rect(data = fwd_bars,
           aes(xmin = position - ((en - st) /1000), xmax = position + ((en - st) /1000),
               ymin = ymin, ymax = ymax), fill = "grey30", alpha = 0.5) +
#  Reverse insertion bars
 geom_rect(data = rev_bars,
           aes(xmin = position - ((en - st) /1000), xmax = position + ((en - st) /1000),
               ymin = ymin, ymax = ymax), fill = "deeppink", alpha = 0.5) +
 # Axis labels and theme
 labs(x = "Genomic position (bp)", y = if (BINARY) "Insertions" else "Insertion count") +
  {if (exists("MAX_PLOT_HEIGHT")) coord_cartesian(ylim = c(0, MAX_PLOT_HEIGHT))} +
 scale_y_continuous(labels = if (BINARY) NULL else abs, limits = c(0, NA), expand = expansion(mult = c(0, 0.02))) +
 theme_classic() +
 theme(axis.title = element_text(size = 12),
       plot.title = element_text(size = 14, face = "bold"),
       axis.text.y = if (BINARY) element_blank() else element_text(),
       axis.ticks.y = if (BINARY) element_blank() else element_line(),
       axis.line.y = if (BINARY) element_blank() else element_line())
p_insertions


p_features <- ggplot() +
  # CDS feature arrows
  geom_polygon(data = arrow_polygons, aes(x = x, y = y, group = feature_ID,
                                          fill = fill_colour, colour = fill_colour)) +
  
  {if (RAINBOW) {
    scale_fill_manual(values = c(setNames(grDevices::colorRampPalette(
      c("yellow", "orange1", "red1", "purple1", "lightblue")
    )(length(unique(arrow_polygons$feature_ID))),
    unique(as.character(arrow_polygons$feature_ID))), "grey70" = "grey70"), guide = "none")
  } else {
    scale_fill_identity(guide = "none")
  }} +
  
  {if (RAINBOW) {
    scale_colour_manual(values = c(setNames(grDevices::colorRampPalette(
      c("yellow", "orange1", "red1", "purple1", "lightblue")
    )(length(unique(arrow_polygons$feature_ID))),
    unique(as.character(arrow_polygons$feature_ID))), "grey70" = "grey70"), guide = "none")
  } else {
    scale_colour_identity(guide = "none")
  }} +
  
  # rRNA, ncRNA and tRNA triangles
  {if (any(rna_triangles$feature_type == "rRNA")) geom_polygon(data = filter(rna_triangles, feature_type == "rRNA"),
                                                               aes(x = x, y = y, group = feature_ID), fill = "purple4", colour = "purple4")} +
  {if (any(rna_triangles$feature_type == "ncRNA")) geom_polygon(data = filter(rna_triangles, feature_type == "ncRNA"),
                                                                aes(x = x, y = y, group = feature_ID), fill = "black", colour = "black")} +
  {if (any(rna_triangles$feature_type == "tRNA")) geom_polygon(data = filter(rna_triangles, feature_type == "tRNA"),
                                                               aes(x = x, y = y, group = feature_ID), fill = "green3", colour = "green3")} +
  #  # Feature IDs centered inside the arrows
  geom_text(data = features %>% filter(feature_type == "CDS"),
            aes(x = (start + end) / 2, y = y_offset, label = feature_name), fontface = "italic", size = 3, colour = "black") +
  # Axis labels and theme
  labs(x = "Genomic position (bp)", y = NULL) +
#  scale_y_continuous(labels = NULL, limits = c(-arrow_height * 2.8, 0), expand = expansion(mult = 0)) +
  scale_x_continuous(limits = c(st, en), expand = expansion(mult = 0)) +
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
        axis.line.x = element_blank())
p_features
  
# plot the two together
p_insertions / p_features +
  plot_layout(heights = c(3, 1))


