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

load("insertion_sites.RData")
load("insertion_sites_ecoli.RData")


#### set user input variables ####

# user input parameters
BINARY = T # if false plot the insertion site count, if true, plot insertion site location only
BINARY = F

# The region you want to plot
chr <- "NZ_CP009273.1"
st <- 5000
en <- 10000

#### prep data frames for plotting ####

# filter data to the region for plotting
my_region <- filter(insertion_sites, chromosome == chr & position >= st & position <= en)

# Summarise the genomic coordinates of each feature
features <- my_region %>%
  filter(!is.na(feature_name)) %>%
  group_by(feature_name, feature_ID, feature_type, feature_strand) %>%
  summarise(start = min(position), end = max(position), .groups = "drop")

# Define arrow dimensions
arrow_height <- 14
arrow_head_length <- 0.2  # Fraction of feature length
arrow_head_height <- 1.3   # Relative to arrow half-height

# Prepare forward and reverse insertion bars
if (BINARY) {
  fwd_bars <- my_region %>% filter(fwd_insertions > 0) %>% mutate(ymin = 0, ymax = (arrow_height * 0.6) * 3)
  rev_bars <- my_region %>% filter(rev_insertions > 0) %>% mutate(ymin = -(arrow_height * 0.6) * 3, ymax = 0)
} else {
  fwd_bars <- my_region %>% mutate(ymin = (arrow_height * 0.6), ymax = (arrow_height * 0.6) + fwd_insertions)
  rev_bars <- my_region %>% mutate(ymin = -(arrow_height * 0.6) - rev_insertions, ymax = -(arrow_height * 0.6))
}

# Calculate arrow polygon coordinates for CDS features only
arrow_polygons <- features %>%
  filter(feature_type == "CDS") %>%
  mutate(length = end - start, head_length = length * arrow_head_length,
         body_height = arrow_height * 0.6, head_height = arrow_height) %>%
  rowwise() %>%
  do(data.frame(feature_name = .$feature_name, feature_ID = .$feature_ID,
                x = if (.$feature_strand == "+") c(.$start, .$end - .$head_length, .$end - .$head_length, .$end,
                                                   .$end - .$head_length, .$end - .$head_length, .$start) else c(.$end, .$start + .$head_length, .$start + .$head_length, .$start,
                                                                                                                 .$start + .$head_length, .$start + .$head_length, .$end),
                y = c(-.$body_height, -.$body_height, -.$head_height, 0,
                      .$head_height, .$body_height, .$body_height))) %>%
  ungroup()

# Prepare rRNA and tRNA triangles
rna_triangles <- features %>%
  filter(feature_type %in% c("rRNA", "tRNA", "ncRNA")) %>%
  mutate(mid = (start + end) / 2,
         triangle_width = ifelse(feature_type == "rRNA", arrow_height * 0.3, ifelse(feature_type == "tRNA", arrow_height * 0.2, arrow_height * 0.15)),
         triangle_length = (end - start),
         tip = ifelse(feature_strand == "+", end, start),
         base = ifelse(feature_strand == "+", end - triangle_length, start + triangle_length),
         ybase = ifelse(feature_type == "rRNA", triangle_width, triangle_width)) %>%
  rowwise() %>%
  do(data.frame(feature_ID = .$feature_ID, feature_type = .$feature_type,
                x = c(.$base, .$tip, .$base),
                y = c(-.$ybase, 0, .$ybase))) %>%
  ungroup()


# Plot
ggplot() +
  # Forward insertion bars
  geom_rect(data = fwd_bars,
    aes(xmin = position - ((en - st) /400), xmax = position + ((en - st) /400),
        ymin = ymin, ymax = ymax), fill = "grey30") +
  # Reverse insertion bars
  geom_rect(data = rev_bars,
    aes(xmin = position - ((en - st) /400), xmax = position + ((en - st) /400),
        ymin = ymin, ymax = ymax), fill = "deeppink") +
  # CDS feature arrows
  geom_polygon(data = arrow_polygons,
               aes(x = x, y = y, group = feature_ID), fill = "grey70", colour = "grey70") +
  # rRNA, ncRNA and tRNA triangles
  {if (any(rna_triangles$feature_type == "rRNA")) geom_polygon(data = filter(rna_triangles, feature_type == "rRNA"),
                                                               aes(x = x, y = y, group = feature_ID), fill = "purple4", colour = "purple4")} +
  {if (any(rna_triangles$feature_type == "ncRNA")) geom_polygon(data = filter(rna_triangles, feature_type == "ncRNA"),
                                                                aes(x = x, y = y, group = feature_ID), fill = "black", colour = "black")} +
  {if (any(rna_triangles$feature_type == "tRNA")) geom_polygon(data = filter(rna_triangles, feature_type == "tRNA"),
                                                               aes(x = x, y = y, group = feature_ID), fill = "#00FF00", colour = "#00FF00")} +
  # geom_polygon(data = filter(rna_triangles, feature_type == "rRNA"),
  #              aes(x = x, y = y, group = feature_ID), fill = "purple4", colour = "purple4") +
  # geom_polygon(data = filter(rna_triangles, feature_type == "ncRNA"),
  #              aes(x = x, y = y, group = feature_ID), fill = "black", colour = "black") +
  # geom_polygon(data = filter(rna_triangles, feature_type == "tRNA"),
  #              aes(x = x, y = y, group = feature_ID), fill = "#00FF00", colour = "#00FF00") +
  # Feature IDs centered inside the arrows
  geom_text(data = features %>% filter(feature_type == "CDS"),
            aes(x = (start + end) / 2, y = 0, label = feature_name), fontface = 'italic', size = 3, colour = "black") +
    # Axis labels and theme
  labs(x = "Genomic position (bp)", y = if (BINARY) NULL else "Insertion count") +
  scale_y_continuous(labels = if (BINARY) NULL else abs) +
  theme_classic() +
  theme(axis.title = element_text(size = 12),
        plot.title = element_text(size = 14, face = "bold"),
        axis.text.y = if (BINARY) element_blank() else element_text(),
        axis.ticks.y = if (BINARY) element_blank() else element_line(),
        axis.line.y = if (BINARY) element_blank() else element_line())







