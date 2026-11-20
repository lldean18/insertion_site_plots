######################################
# prep files for insertion site plot #
######################################

#### read in the files ####

#fasta_file = "GCF_000750555.1_ASM75055v1_genomic.fna"
#fasta_file = "323630L_Photorhabduskhanii.fna"
read.fasta <- function(fasta_file) {
  
  # read in assembly file
  assembly <- readLines(fasta_file)
  
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
  #rm(assembly, chromosome_lengths, chromosome_names, end, start, header_lines, i)
  return(genome_positions)
  
}
#genome_positions <- read.fasta(fasta_file)

#gff_file = "GCF_000750555.1_ASM75055v1_genomic.gff"
#gff_file = "323630L_Photorhabduskhanii.gff"
read.gff <- function(gff_file, chr, st, en) {
  
  # read in the file
  gff_lines <- readLines(gff_file) # read each line separately
  fasta_start <- grep("^##FASTA", gff_lines) # Find the start of the FASTA section
  if (length(fasta_start) > 0) {
    gff_lines <- gff_lines[seq_len(fasta_start[1] - 1)] # Remove the FASTA section if present
  }
  gff <- read.table(text = gff_lines, sep = "\t", header = FALSE, comment.char = "#",
                    stringsAsFactors = FALSE) # Read the annotation into a data frame
  colnames(gff) <- c("chromosome","source","feature_type","start","end","score",
                     "strand","phase","attributes") # Assign column names
  #rm(fasta_start, gff_lines) # clean up
  
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
  
  # rename columns to match
  names(gff)[names(gff) == "strand"] <- "feature_strand"
  names(gff)[names(gff) == "phase"] <- "feature_phase"
  
  # filter only to the region of interest
  gff <- dplyr::filter(gff, chromosome == chr & start >= st & end <= en)
  
  return(gff)
  
}
#gff <- read.gff(gff_file)
#gff <- read.gff(gff_file, "contig_1", 4000, 10000)

#insertions <- read.table("BWtacXpress1_EKDL260002324-1A_23GK55LT4_L5_1.fq.gz.NZ_CP009273.1.insert_site_plot.gz")
#insertions <- read.table("combined.insert_site_plot.gz")

#####################################################
########## TRYING TO USE LESS MEMORY ################
### process insertion sites only then filter
prep_insertion_sites <- function(genome_positions, insertions, chr, st, en) {
  
  # Merge the assembly and insertion site info
  insertion_sites <- cbind(genome_positions, insertions)
  colnames(insertion_sites) <- c("chromosome","position","fwd_insertions","rev_insertions")
  
  # filter the insertion sites to the region of interest
  insertion_sites_flt <- dplyr::filter(insertion_sites, chromosome == chr & position >= st & position <= en)
  
  return(insertion_sites_flt)
  
}

# testing
#test <- prep_insertion_sites(genome_positions, insertions, "contig_1", 4000, 10000)

### process the files ###
prep_insertion_files <- function(insertion_sites, gff) {
  
  # Convert to data.table
  insertion_sites_dt <- data.table::as.data.table(insertion_sites)
  gff_dt <- data.table::as.data.table(gff)
  
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
  data.table::setkey(gff_dt, chromosome, start, end)
  
  # Find overlaps between insertion sites and GFF annotations
  overlaps <- data.table::foverlaps(
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
    data.table::as.data.table(insertion_sites)[, .(site_id = .I, .SD)],
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
  
  return(insertion_sites_annotated)
  
}

# testing
#test2 <- prep_insertion_files(test, gff)

#####################################################








# ### process the files ###
# prep_insertion_files <- function(genome_positions, gff, insertions) {
#   
#   # Merge the assembly and insertion site info
#   insertion_sites <- cbind(genome_positions, insertions)
#   colnames(insertion_sites) <- c("chromosome","position","fwd_insertions","rev_insertions")
#   
#   # Convert to data.table
#   insertion_sites_dt <- data.table::as.data.table(insertion_sites)
#   gff_dt <- data.table::as.data.table(gff)
#   
#   # Preserve original row numbers
#   insertion_sites_dt[, site_id := .I]
#   gff_dt[, annotation_id := .I]
#   
#   # Represent insertion sites as single-base intervals
#   insertion_sites_dt[, `:=`(
#     start = position,
#     end = position
#   )]
#   
#   # Keep the required GFF columns
#   gff_dt <- gff_dt[, .(
#     chromosome,
#     start,
#     end,
#     annotation_id,
#     feature_name,
#     feature_ID,
#     feature_strand,
#     feature_type,
#     feature_phase
#   )]
#   
#   # Index the GFF annotations
#   data.table::setkey(gff_dt, chromosome, start, end)
#   
#   # Find overlaps between insertion sites and GFF annotations
#   overlaps <- data.table::foverlaps(
#     insertion_sites_dt,
#     gff_dt,
#     by.x = c("chromosome", "start", "end"),
#     by.y = c("chromosome", "start", "end"),
#     type = "within",
#     nomatch = NA
#   )
#   
#   # Select the last matching annotation for each insertion site
#   # Keep all overlapping annotations, preserving NA values
#   annotations <- overlaps[
#     !is.na(annotation_id),
#     .(
#       site_id,
#       feature_name,
#       feature_ID,
#       feature_strand,
#       feature_type,
#       feature_phase
#     )
#   ]
#   
#   # Join annotations to the original insertion site information
#   insertion_sites_annotated <- merge(
#     data.table::as.data.table(insertion_sites)[, .(site_id = .I, .SD)],
#     annotations,
#     by = "site_id",
#     all.x = TRUE,
#     sort = FALSE,
#     suffixes = c("", "_annotation")
#   )
#   
#   # Convert back to a data.frame if desired
#   insertion_sites_annotated <- as.data.frame(insertion_sites_annotated)
#   
#   # get rid of the additional string in column names
#   insertion_sites_annotated <- insertion_sites_annotated %>% rename_with(~str_remove(., '.SD.'))
#   
#   return(insertion_sites_annotated)
#   
# }
#insertion_sites_annotated <- prep_insertion_files(genome_positions, gff, insertions)



