######################################
# prep files for insertion site plot #
######################################

# install libraries if required
list.of.packages <- c("tidyverse", "data.table")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)>0) {install.packages(new.packages)}
rm(list.of.packages, new.packages)

# load libraries
library(tidyverse)

#### read in the files ####

# assembly
assembly <- readLines("GCF_000750555.1_ASM75055v1_genomic.fna")
# insertion sites
insertions <- read.table("BWtacXpress1_EKDL260002324-1A_23GK55LT4_L5_1.fq.gz.NZ_CP009273.1.insert_site_plot.gz")
# annotation
gff_lines <- readLines("GCF_000750555.1_ASM75055v1_genomic.gff") # read each line separately
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
gff <- filter(gff, feature_type != "region")

# extract feature names and ID's
gff$feature_name <- "Unknown"
gff$feature_name <- ifelse(grepl("gene=", gff$attributes),
                           sub(".*gene=([^;]+).*", "\\1", gff$attributes),
                           ifelse(grepl("Name=", gff$attributes),
                                  sub(".*Name=([^;]+).*", "\\1", gff$attributes),
                           sub(".*ID=([^;]+).*", "\\1", gff$attributes)))
#any(grepl("Unknown", gff$feature_name)) # check if any unknowns remain

gff$feature_ID <- "Unknown"
gff$feature_ID <- ifelse(grepl("ID=", gff$attributes),
                         sub(".*ID=([^;]+).*", "\\1", gff$attributes),
                         gff$feature_ID)
#any(grepl("Unknown", gff$feature_ID)) # check if any unknowns remain

# load library
library(data.table)

# rename columns to match
names(gff)[names(gff) == "strand"] <- "feature_strand"
names(gff)[names(gff) == "phase"] <- "feature_phase"

# Convert to data.table
insertion_sites_dt <- as.data.table(insertion_sites)
gff_dt <- as.data.table(gff)

# Preserve original row numbers
insertion_sites_dt[, site_id := .I]
gff_dt[, annotation_id := .I]

# Represent insertion sites as single-base intervals
insertion_sites_dt[, `:=`(
  start = position,
  end = position
)]

# Keep the required GFF columns
gff_dt <- gff_dt[, .(
  chromosome,
  start,
  end,
  annotation_id,
  feature_name,
  feature_ID,
  feature_strand,
  feature_type,
  feature_phase
)]

# Index the GFF annotations
setkey(gff_dt, chromosome, start, end)

# Find overlaps between insertion sites and GFF annotations
overlaps <- foverlaps(
  insertion_sites_dt,
  gff_dt,
  by.x = c("chromosome", "start", "end"),
  by.y = c("chromosome", "start", "end"),
  type = "within",
  nomatch = NA
)

# Select the last matching annotation for each insertion site
# Keep all overlapping annotations, preserving NA values
annotations <- overlaps[
  !is.na(annotation_id),
  .(
    site_id,
    feature_name,
    feature_ID,
    feature_strand,
    feature_type,
    feature_phase
  )
]

# Join annotations to the original insertion site information
insertion_sites_annotated <- merge(
  as.data.table(insertion_sites)[, .(site_id = .I, .SD)],
  annotations,
  by = "site_id",
  all.x = TRUE,
  sort = FALSE,
  suffixes = c("", "_annotation")
)

# Convert back to a data.frame if desired
insertion_sites_annotated <- as.data.frame(insertion_sites_annotated)

# get rid of the additional string in column names
insertion_sites_annotated <- insertion_sites_annotated %>% rename_with(~str_remove(., '.SD.'))

# cleanup
rm(insertion_sites, insertion_sites_dt, gff_dt, gff, overlaps, annotations)


#### save the plotting data frame at this point as it took a while to make ####

save(insertion_sites_annotated, file = "insertion_sites_ecoli_annotated.RData")


