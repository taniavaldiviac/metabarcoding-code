#Code to construct the RDF of chordata from Fasta file, to fasta file with taxonomy to be used in dada2

# environment
library(insect)
library(tidyverse)
library(Biostrings)
library(rentrez)
library(purrr)
library(dplyr)
library(taxize)
library(dplyr)
library(stringr)


fasta_path <- "metadata/Amplicons_12S_virtualPCR_raw.fasta"
df <- fasta_to_df(fasta_path)

# Extract accession number from the Title column
df <- df %>%
  mutate(accession = str_extract(Title, "^\\S+"))

# Now merge with sp_names
df_merged <- df %>%
  left_join(sp_names, by = c("accession" = "accession_numbers"))

# # Now, let's process the df dataframe
 df_processed <- df %>%
   # Extract accession number from the Title column
   mutate(accession = str_extract(Title, "^\\S+"))

 # Create new titles combining clean species name and accession number
df_merged <- df_merged %>%
  mutate(NewTitle = ifelse(!is.na(sp_names), 
                           paste(sp_names, accession, sep = "_"), 
                           Title))
write.csv(df_merged, "./df_merged.csv", row.names = FALSE)
# Remove rows where clean_sp_names is NA (i.e., no match found or filtered out)
df_merged <- df_merged %>%
  filter(!is.na(sp_names) & sp_names != "Homo sapiens")

# Check the result
cat("Number of rows after cleaning:", nrow(df_merged), "\n")
cat("Number of unique species:", length(unique(df_merged$sp_names)), "\n")

#Number of rows after cleaning: 38362
#Number of unique species: 11566 

# Print first few rows of the merged dataframe
print(head(df_merged))

# Optional: write the cleaned data to a CSV file
write.csv(df_merged, "./metadata/cleaned_merged_data.csv", row.names = FALSE)


species.list <- data.frame(species = (df_merged$sp_names))

# Create a small subset of species.list
# set.seed(123)
# subset_species <- sample(species.list$species, 5)
# subset_species_df <- data.frame(species = subset_species)
# subset_species_df <- as.character(subset_species_df$species)
# tax.db <- classification(subset_species_df, db = "itis", rows = 1) 

Sys.setenv(ENTREZ_KEY = "20b168480190bbac6a65958eab3cdb9fd209")

#load("18sep.RData")

species.list.ch <- as.character(species.list$species)

# Get the total number of species
total_species <- length(species.list.ch)

# Define the batch size
batch_size <- 1600  # You can adjust this value as needed

# Create a directory to store batch results if it doesn't exist
dir.create("./metadata/batch_results", showWarnings = FALSE)

# Loop through the species list in batches
for (i in seq(1, total_species, by = batch_size)) {
  # Define the start and end of the current batch
  batch_end <- min(i + batch_size - 1, total_species)  # Ensure not to exceed total species
  
  # Extract the current batch of species
  species_batch <- species.list.ch[i:batch_end]
  
  # Try running the classification function on the batch and handle errors if needed
  tryCatch({
    batch_results <- classification(species_batch, db = "ncbi", rows = 1)
    
    # Save the batch results
    batch_filename <- paste0("batch_results/batch_", i, "_to_", batch_end, ".RData")
    save(batch_results, file = batch_filename)
    
    cat("Saved batch", i, "to", batch_end, "as", batch_filename, "\n")
    
  }, error = function(e) {
    message(paste("Error on batch", i, "to", batch_end, ":", e))
  })
}

#Get taxonomy information
library(dplyr)
library(purrr)
library(tidyr)

# Obtener la lista de archivos de batch en el directorio
batch_files <- list.files("./batch_results", pattern = "^batch_.*\\.RData$", full.names = TRUE)

# Ordenar los archivos por nombre para asegurarnos de que el primero sea realmente el primer batch
batch_files <- sort(batch_files)


extract_taxonomy <- function(species_result, species_name) {
  if (is.null(species_result) || length(species_result) == 0 || !is.data.frame(species_result)) {
    return(data.frame(species = species_name, 
                      kingdom = NA, phylum = NA, class = NA, order = NA, family = NA, genus = NA,
                      stringsAsFactors = FALSE))
  }
  
  # Ensure 'rank' and 'name' columns exist
  if (!all(c("rank", "name") %in% names(species_result))) {
    return(data.frame(species = species_name, 
                      kingdom = NA, phylum = NA, class = NA, order = NA, family = NA, genus = NA,
                      stringsAsFactors = FALSE))
  }
  
  df <- species_result %>%
    filter(rank %in% c("kingdom", "phylum", "class", "order", "family", "genus")) %>%
    distinct(rank, .keep_all = TRUE) %>%
    select(rank, name) %>%
    pivot_wider(names_from = rank, values_from = name)
  
  df$species <- species_name
  
  return(df)
}

# Loop through the batch files and process each one
for (batch_file in batch_files) {
  cat("Processing:", batch_file, "\n")
  
  # Load the batch data
  load(batch_file)
  
  # Convert the batch results to a data frame
  batch_df <- tryCatch({
    imap_dfr(batch_results, ~extract_taxonomy(.x, .y))
  }, error = function(e) {
    cat("Error processing batch:", batch_file, "\n")
    cat("Error message:", e$message, "\n")
    data.frame()
  })
  
  # Save the batch data frame as a CSV file
  save_file_name <- gsub("batch_results/", "", batch_file)
  save_file_name <- gsub("\\.RData$", ".csv", save_file_name)
  write.csv(batch_df, file.path("batch_results", save_file_name), row.names = FALSE)
  
  cat("Saved batch data frame as:", save_file_name, "\n")
}

# Combine all the batch data frames into a single data frame
all_batch_files <- list.files("batch_results", pattern = "\\.csv$", full.names = TRUE)
all_batch_df <- bind_rows(lapply(all_batch_files, read.csv))

# Show the first few rows of the combined data frame
print(head(all_batch_df))

# Save the combined data frame as a CSV file
write.csv(all_batch_df, "all_taxonomy_results.csv", row.names = FALSE)

# Print a summary
cat("Total number of species in the combined data frame:", nrow(all_batch_df), "\n")
cat("Combined data frame saved as 'all_taxonomy_results.csv'\n")


library(dplyr)
library(Biostrings)
library(data.table)

# Load the data
taxonomy_data <- fread("all_taxonomy_results.csv")
df_merged <- fread("df_merged.csv")  # Replace with actual path

# Convert taxonomy_data to a named list for faster lookups
taxonomy_list <- taxonomy_data[, .(
  taxonomy = paste(kingdom, phylum, class, order, family, genus, species, sep = ";")
), by = species]
taxonomy_dict <- setNames(taxonomy_list$taxonomy, taxonomy_list$species)

# Function to create header
create_header <- function(species, accession) {
  taxonomy <- taxonomy_dict[species]
  headers <- paste0(taxonomy)
  headers[is.na(taxonomy)] <- NA
  headers
}

# Process data
result <- df_merged[, .(
  header = create_header(sp_names, accession),
  sequence = Sequence
)]

# Remove rows with NA headers
result <- result[!is.na(header)]

# Function to write FASTA file with sequences on a single line
write_fasta_single_line <- function(headers, sequences, file_path) {
  stopifnot(length(headers) == length(sequences))
  
  con <- file(file_path, "w")
  on.exit(close(con))
  
  for (i in seq_along(headers)) {
    writeLines(paste0(">", headers[i]), con)
    writeLines(sequences[i], con)
  }
}

# Write to FASTA file with sequences on a single line
write_fasta_single_line(result$header, result$sequence, "taxonomy_sequences.fasta")

cat("Created FASTA file with taxonomy information: taxonomy_sequences.fasta\n")
cat("Total sequences written:", nrow(result), "\n")

# Optional: Verify the format
cat("First few lines of the created FASTA file:\n")
system("head -n 4 taxonomy_sequences.fasta")

