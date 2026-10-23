###################################
### PRETTY INSERTION SITE PLOTS ###
###################################

# install libraries if required
list.of.packages <- c("tidyverse")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)>0) {install.packages(new.packages)}
rm(list.of.packages, new.packages)

# load libraries
library(tidyverse)

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

# The region you want to plot
chr <- "NZ_CP009273.1"
st <- 3862000
en <- 3868000

# features in ecoli to look at
#dnaA 3875686	3877089
#st <- 3872000
#en <- 3878000
#tRNA 3315431	3315517

#### prep data frames for plotting ####

# filter data to the region for plotting
my_region <- filter(insertion_sites, chromosome == chr & position >= st & position <= en)

# # Summarise the genomic coordinates of each feature
# features <- my_region %>%
#   group_by(feature_name, feature_ID, feature_type, feature_strand) %>%
#   summarise(start = min(position), end = max(position), .groups = "drop")

features <- my_region %>%
  group_by(feature_name, feature_ID, feature_type, feature_strand) %>%
  summarise(start = min(position), end = max(position), .groups = "drop") %>%
  mutate(
    length = end - start,
    body_height = arrow_height * 0.6,
    head_height = arrow_height,
    head_length = length * arrow_head_length)

# Determine whether features occur on both strands
both_strands <- any(features$feature_strand == "+", na.rm = TRUE) &&
  any(features$feature_strand == "-", na.rm = TRUE)

# Place features just below zero; separate strands only if both are present
features <- features %>%
  mutate(y_offset = if (both_strands) {
    ifelse(feature_strand == "+", -arrow_height * 0.7, -arrow_height * 2.1)
  } else {
    -arrow_height * 0.7
  })

# Define arrow dimensions
arrow_height <- 14
arrow_head_length <- 0.2  # Fraction of feature length
arrow_head_height <- 1.3   # Relative to arrow half-height

# # Prepare forward and reverse insertion bars
# if (BINARY) {
#   fwd_bars <- my_region %>% filter(fwd_insertions > 0) %>% mutate(ymin = 0, ymax = (arrow_height * 0.6) * 3)
#   rev_bars <- my_region %>% filter(rev_insertions > 0) %>% mutate(ymin = -(arrow_height * 0.6) * 3, ymax = 0)
# } else {
#   fwd_bars <- my_region %>% mutate(ymin = (arrow_height * 0.6), ymax = (arrow_height * 0.6) + fwd_insertions)
#   rev_bars <- my_region %>% mutate(ymin = -(arrow_height * 0.6) - rev_insertions, ymax = -(arrow_height * 0.6))
# }

# Determine the vertical extent of the feature arrows
if (BINARY) {
  feature_top <- max(features$y_offset + features$head_height, na.rm = TRUE)
  feature_bottom <- min(features$y_offset - features$head_height, na.rm = TRUE)
} else {
feature_top <- 0
}

# Prepare forward and reverse insertion bars
if (BINARY) {
  fwd_bars <- my_region %>%
    filter(fwd_insertions > 0) %>%
    mutate(ymin = feature_top, ymax = feature_top + (arrow_height * 0.6) * 3)
  rev_bars <- my_region %>%
    filter(rev_insertions > 0) %>%
    mutate(ymin = feature_bottom - (arrow_height * 0.6) * 3, ymax = feature_bottom)
} else {
  fwd_bars <- my_region %>%
    mutate(ymin = feature_top, ymax = feature_top + fwd_insertions)
  rev_bars <- my_region %>%
#    mutate(ymin = feature_bottom - rev_insertions, ymax = feature_bottom)
    mutate(ymin = feature_top, ymax = feature_top + rev_insertions)
}

# # Calculate arrow polygon coordinates for CDS features only
arrow_polygons <- features %>%
  filter(feature_type == "CDS") %>%
  rowwise() %>%
  do(data.frame(
    feature_name = .$feature_name, feature_ID = .$feature_ID,
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

# arrow_polygons <- features %>%
#   filter(feature_type == "CDS") %>%
#   rowwise() %>%
#   do(data.frame(feature_name = .$feature_name, feature_ID = .$feature_ID,
#     x = if (.$feature_strand == "+") {
#       c(.$start, .$end - .$head_length, .$end - .$head_length,
#         .$end, .$end - .$head_length, .$end - .$head_length, .$start)
#     } else {
#       c(.$end, .$start + .$head_length, .$start + .$head_length,
#         .$start, .$start + .$head_length, .$start + .$head_length, .$end)
#     },
#     y = c(
#       .$y_offset - .$body_height,
#       .$y_offset - .$body_height,
#       .$y_offset - .$head_height,
#       .$y_offset,
#       .$y_offset + .$head_height,
#       .$y_offset + .$body_height,
#       .$y_offset + .$body_height
#     )
#   )) %>%
#   ungroup()

# # Prepare rRNA and tRNA triangles
# rna_triangles <- features %>%
#   filter(feature_type %in% c("rRNA", "tRNA", "ncRNA")) %>%
#   mutate(mid = (start + end) / 2,
#          triangle_width = ifelse(feature_type == "rRNA", arrow_height * 0.3, ifelse(feature_type == "tRNA", arrow_height * 0.2, arrow_height * 0.15)),
#          triangle_length = (end - start),
#          tip = ifelse(feature_strand == "+", end, start),
#          base = ifelse(feature_strand == "+", end - triangle_length, start + triangle_length),
#          ybase = ifelse(feature_type == "rRNA", triangle_width, triangle_width)) %>%
#   rowwise() %>%
#   do(data.frame(feature_ID = .$feature_ID, feature_type = .$feature_type,
#                 x = c(.$base, .$tip, .$base),
#                 y = c(-.$ybase, 0, .$ybase))) %>%
#   ungroup()

rna_triangles <- features %>%
  filter(feature_type %in% c("rRNA", "tRNA", "ncRNA")) %>%
  mutate(triangle_width = ifelse(feature_type == "rRNA", arrow_height * 0.3,
      ifelse(feature_type == "tRNA", arrow_height * 0.2, arrow_height * 0.15)),
    triangle_length = end - start, tip = ifelse(feature_strand == "+", end, start),
    base = ifelse(feature_strand == "+", end - triangle_length, start + triangle_length)) %>%
  rowwise() %>%
  do(data.frame(feature_ID = .$feature_ID, feature_type = .$feature_type, x = c(.$base, .$tip, .$base),
                y = c(.$y_offset - .$triangle_width, .$y_offset, .$y_offset + .$triangle_width))) %>%
  ungroup()

# Plot
ggplot() +
  # Forward insertion bars
  geom_rect(data = fwd_bars,
    aes(xmin = position - ((en - st) /1000), xmax = position + ((en - st) /1000),
        ymin = ymin, ymax = ymax), fill = "grey30", alpha = 0.5) +
  # Reverse insertion bars
  geom_rect(data = rev_bars,
    aes(xmin = position - ((en - st) /1000), xmax = position + ((en - st) /1000),
        ymin = ymin, ymax = ymax), fill = "deeppink", alpha = 0.5) +
  # CDS feature arrows
  geom_polygon(data = arrow_polygons, aes(x = x, y = y, group = feature_ID,
                                          fill = if (RAINBOW) feature_ID else "grey70"),
                                          colour = if (RAINBOW) NA else "grey70") +
  {if (RAINBOW) {
    scale_fill_manual(values = setNames(grDevices::colorRampPalette(
          c("yellow", "orange", "red", "purple", "blue")
        )(length(unique(arrow_polygons$feature_ID))),
        unique(arrow_polygons$feature_ID)), guide = "none")
  } else {
    scale_fill_identity(guide = "none")
  }} +
    # rRNA, ncRNA and tRNA triangles
  {if (any(rna_triangles$feature_type == "rRNA")) geom_polygon(data = filter(rna_triangles, feature_type == "rRNA"),
                                                               aes(x = x, y = y, group = feature_ID), fill = "purple4", colour = "purple4")} +
  {if (any(rna_triangles$feature_type == "ncRNA")) geom_polygon(data = filter(rna_triangles, feature_type == "ncRNA"),
                                                                aes(x = x, y = y, group = feature_ID), fill = "black", colour = "black")} +
  {if (any(rna_triangles$feature_type == "tRNA")) geom_polygon(data = filter(rna_triangles, feature_type == "tRNA"),
                                                               aes(x = x, y = y, group = feature_ID), fill = "#00FF00", colour = "#00FF00")} +
#  # Feature IDs centered inside the arrows
  geom_text(data = features %>% filter(feature_type == "CDS"),
            aes(x = (start + end) / 2, y = y_offset, label = feature_name), fontface = "italic", size = 3, colour = "black") +
      # Axis labels and theme
  labs(x = "Genomic position (bp)", y = if (BINARY) NULL else "Insertion count") +
#  scale_y_continuous(labels = if (BINARY) NULL else abs) +
  theme_classic() +
  theme(axis.title = element_text(size = 12),
        plot.title = element_text(size = 14, face = "bold"),
        axis.text.y = if (BINARY) element_blank() else element_text(),
        axis.ticks.y = if (BINARY) element_blank() else element_line(),
        axis.line.y = if (BINARY) element_blank() else element_line())







