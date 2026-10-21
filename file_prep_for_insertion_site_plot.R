######################################
# prep files for insertion site plot #
######################################

# install libraries if required
list.of.packages <- c("tidyverse")
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
gff <- gff[gff$feature_type != "region", ]

# extract feature names and ID's
gff$feature_name <- ifelse(grepl("gene=", gff$attributes), sub(".*gene=([^;]+).*", "\\1", gff$attributes), sub(".*Name=([^;]+).*", "\\1", gff$attributes))
#gff$feature_name <- sub(".*Name=([^;]+).*", "\\1", gff$attributes)
gff$feature_ID <- sub(".*ID=([^;]+).*", "\\1", gff$attributes)
# Set the name or ID to NA if no attribute exists
gff$feature_name[!grepl("Name=", gff$attributes)] <- NA
gff$feature_ID[!grepl("ID=", gff$attributes)] <- NA

# for ecoli needed to remove rows where the source is protein homology because they are all duplicates of genes


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

save(insertion_sites, file = "insertion_sites_ecoli.RData")
#load("insertion_sites.RData")

