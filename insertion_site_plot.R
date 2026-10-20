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

# first set the region you want to plot
chr <- "contig_1"
st <- 86
en <- 4000
my_region <- filter(insertion_sites, chromosome == chr & position >= st & position <= en)

# Summarise the genomic coordinates of each feature within the region
features <- my_region %>%
  filter(!is.na(feature_name)) %>%
  group_by(feature_name, feature_ID) %>%
  summarise(
    start = min(position),
    end = max(position),
    .groups = "drop"
  )

# Calculate the maximum insertion count to scale the plot
max_insertions <- max(
  my_region$fwd_insertions,
  my_region$rev_insertions,
  na.rm = TRUE
)

# the plot
ggplot() +
  
  # Forward insertions: bars extending upwards
  geom_col(
    data = my_region,
    aes(x = position, y = fwd_insertions),
    width = 1,
    fill = "steelblue"
  ) +
  
  # Reverse insertions: bars extending downwards
  geom_col(
    data = my_region,
    aes(x = position, y = -rev_insertions),
    width = 1,
    fill = "firebrick"
  ) +
  
  # Feature arrows
  geom_segment(
    data = features,
    aes(
      x = start,
      xend = end,
      y = 0,
      yend = 0
    ),
    arrow = arrow(
      length = unit(2, "mm"),
      type = "closed"
    ),
    linewidth = 0.5,
    colour = "black"
  ) +
  
  # Feature IDs centred on the arrows
  geom_text(
    data = features,
    aes(
      x = (start + end) / 2,
      y = 0,
      label = feature_ID
    ),
    size = 3,
    vjust = -0.8
  ) +
  
  # Ensure the arrows and insertion bars share the same zero baseline
  geom_hline(
    yintercept = 0,
    linewidth = 0.3,
    colour = "black"
  ) +
  
  # Labels and theme
  labs(
    x = "Genomic position (bp)",
    y = "Insertion count",
    title = "Insertion sites across the genome"
  ) +
  
  scale_y_continuous(
    labels = abs
  ) +
  
  theme_classic() +
  
  theme(
    axis.text.y = element_text(colour = "black"),
    axis.title = element_text(size = 12),
    plot.title = element_text(size = 14, face = "bold")
  )



###################



# Summarise the genomic coordinates of each feature
features <- my_region %>%
  filter(!is.na(feature_name)) %>%
  group_by(feature_name, feature_ID) %>%
  summarise(
    start = min(position),
    end = max(position),
    .groups = "drop"
  )

# Set the half-height of the feature arrows
arrow_height <- 6

# Set the scale factor for insertion bar heights
insertion_scale <- 1

# Prepare forward and reverse insertion bars
fwd_bars <- my_region %>%
  mutate(
    ymin = arrow_height,
    ymax = arrow_height + fwd_insertions * insertion_scale
  )

rev_bars <- my_region %>%
  mutate(
    ymin = -arrow_height - rev_insertions * insertion_scale,
    ymax = -arrow_height
  )

# Plot
ggplot() +
  # Forward insertion bars
  geom_rect(data = fwd_bars,
    aes(xmin = position - 1, xmax = position + 1, ymin = ymin, ymax = ymax),
    fill = "steelblue") +
  # Reverse insertion bars
  geom_rect(data = rev_bars,
    aes(xmin = position - 0.5, xmax = position + 0.5, ymin = ymin, ymax = ymax),
    fill = "firebrick") +
  # Thick feature arrows
  geom_segment(data = features, aes(x = start, xend = end, y = 0, yend = 0),
    arrow = arrow(length = unit(2, "mm"), type = "closed"),
    linewidth = arrow_height, colour = "grey70", lineend = "butt", linejoin = "mitre") +
  # Feature IDs centred inside the arrows
  geom_text(data = features,
    aes(x = (start + end) / 2, y = 0, label = feature_ID), size = 3, colour = "black") +
  # Axis labels and theme
  labs(x = "Genomic position (bp)", y = "Insertion count") +
  scale_y_continuous(labels = abs) +
  # theme
  theme_classic() +
  theme(axis.title = element_text(size = 12),
        plot.title = element_text(size = 14, face = "bold"))








