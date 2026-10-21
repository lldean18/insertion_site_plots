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

#### read in the files ####

# assembly
assembly <- readLines("323630L_Photorhabduskhanii.fna")
# insertion sites
insertions <- read.table("combined.insert_site_plot.gz")
# annotation
gff_lines <- readLines("323630L_Photorhabduskhanii.gff") # read each line separately
fasta_start <- grep("^##FASTA", gff_lines) # Find the start of the FASTA section
if (length(fasta_start) > 0) {
  gff_lines <- gff_lines[seq_len(fasta_start[1] - 1)] # Remove the FASTA section if present
}
gff <- read.table(text = gff_lines, sep = "\t", header = FALSE, comment.char = "#",
                  stringsAsFactors = FALSE) # Read the annotation into a data frame
colnames(gff) <- c("chromosome","source","feature_type","start","end","score",
                   "strand","phase","attributes") # Assign column names
rm(fasta_start, gff_lines) # clean up


#### Process the fasta assembly ####
# Find the chromosome headers
header_lines <- grep("^>", assembly)
# Extract chromosome names (removing > and any description after the first space)
chromosome_names <- sub("^>([^[:space:]]+).*", "\\1", assembly[header_lines])

# Calculate sequence lengths, accounting for sequences split over multiple lines
chromosome_lengths <- integer(length(header_lines))
for (i in seq_along(header_lines)) {
  start <- header_lines[i] + 1
  end <- if (i < length(header_lines)) header_lines[i + 1] - 1 else length(assembly)
  
  chromosome_lengths[i] <- sum(nchar(assembly[start:end]))
}

# Create the data frame
genome_positions <- data.frame(
  chromosome = rep(chromosome_names, chromosome_lengths),
  position = unlist(lapply(chromosome_lengths, seq_len), use.names = FALSE))

# tidy up
rm(assembly, chromosome_lengths, chromosome_names, end, start, header_lines, i)

#### Merge the assembly and insertion site info ####
insertion_sites <- cbind(genome_positions, insertions)
colnames(insertion_sites) <- c("chromosome","position","fwd_insertions","rev_insertions")
# tidy up
rm(genome_positions, insertions)

#### Process the GFF annotation ####

# remove columns containing the whole contigs as regions
gff <- gff[gff$feature_type != "region", ]

# extract feature names and ID's
gff$feature_name <- sub(".*Name=([^;]+).*", "\\1", gff$attributes)
gff$feature_ID <- sub(".*ID=([^;]+).*", "\\1", gff$attributes)
# Set the name or ID to NA if no attribute exists
gff$feature_name[!grepl("Name=", gff$attributes)] <- NA
gff$feature_ID[!grepl("ID=", gff$attributes)] <- NA

# add the feature names and ID's to the insertion sites data frame
# Initialize new columns
insertion_sites$feature_name <- NA_character_
insertion_sites$feature_ID <- NA_character_
insertion_sites$feature_strand <- NA_character_
insertion_sites$feature_type <- NA_character_
insertion_sites$feature_phase <- NA_character_

# add the annotation information to the insertion sites data
for (line in seq_len(nrow(gff))) {
  annotation <- gff[line, ] # Get the current annotation line
  # add feature names
  insertion_sites$feature_name <- ifelse(insertion_sites$chromosome == annotation$chromosome &
           insertion_sites$position >= annotation$start &
           insertion_sites$position <= annotation$end, annotation$feature_name, insertion_sites$feature_name)
  # add feature IDs
  insertion_sites$feature_ID <- ifelse(insertion_sites$chromosome == annotation$chromosome &
                                           insertion_sites$position >= annotation$start &
                                           insertion_sites$position <= annotation$end,
                                       annotation$feature_ID, insertion_sites$feature_ID)
  # add feature strand
  insertion_sites$feature_strand <- ifelse(insertion_sites$chromosome == annotation$chromosome &
                                         insertion_sites$position >= annotation$start &
                                         insertion_sites$position <= annotation$end,
                                       annotation$strand, insertion_sites$feature_strand)
  # add feature type
  insertion_sites$feature_type <- ifelse(insertion_sites$chromosome == annotation$chromosome &
                                             insertion_sites$position >= annotation$start &
                                             insertion_sites$position <= annotation$end,
                                           annotation$feature_type, insertion_sites$feature_type)
  # add feature phase
  insertion_sites$feature_phase <- ifelse(insertion_sites$chromosome == annotation$chromosome &
                                           insertion_sites$position >= annotation$start &
                                           insertion_sites$position <= annotation$end,
                                         annotation$phase, insertion_sites$feature_phase)
  # Report progress every 100 annotations (untested)
  if (line %% 100 == 0) {
    cat("Processed", line, "of", nrow(gff), "annotations (", round(100 * line / nrow(gff), 1), "%)\n")
  }
}

rm(annotation, line, gff) # clean up

#### save the plotting data frame at this point as it took a while to make ####

save(insertion_sites, file = "insertion_sites.RData")
#load("insertion_sites.RData")

#### Draw the plot ####

# user input parameters
BINARY = T # if false plot the insertion site count, if true, plot insertion site location only
BINARY = F

# first set the region you want to plot
chr <- "contig_2"
st <- 1
en <- 10000
my_region <- filter(insertion_sites, chromosome == chr & position >= st & position <= en)


###################



# # Summarise the genomic coordinates of each feature
# features <- my_region %>%
#   filter(!is.na(feature_name)) %>%
#   group_by(feature_name, feature_ID) %>%
#   summarise(start = min(position), end = max(position), .groups = "drop")

# Summarise the genomic coordinates of each feature
features <- my_region %>%
  filter(!is.na(feature_name)) %>%
  group_by(feature_name, feature_ID, feature_type, feature_strand) %>%
  summarise(start = min(position), end = max(position), .groups = "drop")


# Set the scale factor for insertion bar heights
insertion_scale <- 1

# Define arrow dimensions
arrow_height <- 14
arrow_head_length <- 0.2  # Fraction of feature length
arrow_head_height <- 1.3   # Relative to arrow half-height

# Prepare forward and reverse insertion bars
if (BINARY) {
  fwd_bars <- my_region %>% filter(fwd_insertions > 0) %>% mutate(ymin = 0, ymax = (arrow_height * 0.6) * 3)
  rev_bars <- my_region %>% filter(rev_insertions > 0) %>% mutate(ymin = -(arrow_height * 0.6) * 3, ymax = 0)
} else {
  fwd_bars <- my_region %>% mutate(ymin = (arrow_height * 0.6), ymax = (arrow_height * 0.6) + fwd_insertions * insertion_scale)
  rev_bars <- my_region %>% mutate(ymin = -(arrow_height * 0.6) - rev_insertions * insertion_scale, ymax = -(arrow_height * 0.6))
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
  filter(feature_type %in% c("rRNA", "tRNA")) %>%
  mutate(mid = (start + end) / 2,
         triangle_width = ifelse(feature_type == "rRNA", arrow_height * 0.4, arrow_height * 0.2),
         triangle_length = ifelse(feature_type == "rRNA", (end - start), (end - start)),
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
        ymin = ymin, ymax = ymax), fill = "grey20") +
  # Reverse insertion bars
  geom_rect(data = rev_bars,
    aes(xmin = position - ((en - st) /400), xmax = position + ((en - st) /400),
        ymin = ymin, ymax = ymax), fill = "deeppink") +
  # CDS feature arrows
  geom_polygon(data = arrow_polygons,
               aes(x = x, y = y, group = feature_ID), fill = "grey70", colour = "grey70") +
  # rRNA and tRNA triangles
  geom_polygon(data = filter(rna_triangles, feature_type == "rRNA"),
               aes(x = x, y = y, group = feature_ID), fill = "purple", colour = "black") +
  geom_polygon(data = filter(rna_triangles, feature_type == "tRNA"),
               aes(x = x, y = y, group = feature_ID), fill = "#00FF00", colour = "#00FF00") +
  # Feature IDs centered inside the arrows
  geom_text(data = features %>% filter(feature_type == "CDS"),
            aes(x = (start + end) / 2, y = 0, label = feature_ID), size = 3, colour = "black") +
    # Axis labels and theme
  labs(x = "Genomic position (bp)", y = if (BINARY) NULL else "Insertion count") +
  scale_y_continuous(labels = if (BINARY) NULL else abs) +
  theme_classic() +
  theme(axis.title = element_text(size = 12),
        plot.title = element_text(size = 14, face = "bold"),
        axis.text.y = if (BINARY) element_blank() else element_text(),
        axis.ticks.y = if (BINARY) element_blank() else element_line(),
        axis.line.y = if (BINARY) element_blank() else element_line())







