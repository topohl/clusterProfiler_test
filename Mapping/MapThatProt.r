#' UniProt ID Mapping for ClusterProfiler Analysis
#'
#' This script processes gene data with UniProtKB IDs (e.g., "GENE_MOUSE") and maps them
#' to UniProt Accessions using an official UniProt ID mapping file.
#'
#' @details
#' The workflow includes:
#' \enumerate{
#'   \item Setting up the environment and loading necessary libraries.
#'   \item Defining project directories.
#'   \item Checking for the existence of input data and the mapping file.
#'   \item Preprocessing gene identifiers by extracting the primary ID when multiple are provided.
#'   \item Reading and filtering the official UniProt mapping file (MOUSE_10090_idmapping.dat) to obtain an Entry Name to Accession map.
#'   \item Merging the mapping data with the gene information.
#'   \item Saving the mapped results for downstream analysis.
#' }
#'
#' @note Ensure that the mapping file is downloaded, unzipped, and placed in the specified directory before running this script.
#'
#' @author Tobias Pohl

# -----------------------------------------------------
# Load Libraries
# -----------------------------------------------------

if (!requireNamespace("pacman", quietly = TRUE))
    install.packages("pacman")

# Load required libraries with pacman
pacman::p_load(dplyr, stringr, tidyr, purrr, readr)

# -----------------------------------------------------
# Define Paths and Project Directories
# -----------------------------------------------------

working_dir <- "/Users/tobiaspohl/Documents/clusterProfiler"
results_dir <- file.path(working_dir, "Results")
setwd(working_dir)

# Specify cell types which determine dataset file name
cell_types <- c("TESTmcherryG1", "mcherryG2")
file_name <- paste0(paste(cell_types, collapse = "_"), ".csv")
data_path <- file.path(working_dir, "Datasets", file_name)

# IMPORTANT: Set path for the downloaded UniProt mapping file
uniprot_mapping_file_path <- file.path(working_dir, "Datasets", "MOUSE_10090_idmapping.dat")

# -----------------------------------------------------
# Verify Input Files
# -----------------------------------------------------

if (!file.exists(data_path)) {
    stop("Input data file does not exist: ", data_path)
}
if (!file.exists(uniprot_mapping_file_path)) {
    stop("UniProt mapping file not found at: ", uniprot_mapping_file_path,
         "\nDownload and unzip from: ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/MOUSE_10090_idmapping.dat.gz")
}

# -----------------------------------------------------
# Preprocess Gene Data
# -----------------------------------------------------

cat("Reading gene data from:", data_path, "\n")
df <- readr::read_csv(data_path, show_col_types = FALSE) %>%
    mutate(gene_symbol = as.character(gene_symbol)) %>%
    # In the case of semicolon-separated IDs, retain the primary (first) identifier
    mutate(gene_symbol = map_chr(str_split(gene_symbol, ";"), ~str_trim(.x[1]))) %>%
    filter(gene_symbol != "" & !is.na(gene_symbol)) %>%
    distinct(gene_symbol, .keep_all = TRUE)

cat("Total genes loaded:", nrow(df), "\n")
if (nrow(df) == 0) {
    stop("No valid gene symbols were identified after preprocessing.")
}
cat("Preview of gene data BEFORE mapping:\n")
print(head(df))

# -----------------------------------------------------
# Parse UniProt Mapping File
# -----------------------------------------------------

cat("Loading UniProt mapping file from:", uniprot_mapping_file_path, "\n")
uniprot_mapping <- readr::read_tsv(
    uniprot_mapping_file_path,
    col_names = c("UniProt_Accession", "Type", "Value"),
    col_types = "ccc",  # reading all columns as characters
    quote = ""         # avoid handling quotes in the identifiers
)

cat("Total entries in UniProt mapping file:", nrow(uniprot_mapping), "\n")

# Extract mappings for UniProtKB-ID to Accession
entry_name_to_accession <- uniprot_mapping %>%
    dplyr::filter(Type == "UniProtKB-ID") %>%
    dplyr::select(UniProt_Accession, UniProtKB_ID = Value) %>% 
    dplyr::distinct(UniProtKB_ID, .keep_all = TRUE)

cat("Total unique UniProtKB-ID mappings identified:", nrow(entry_name_to_accession), "\n")
if (nrow(entry_name_to_accession) == 0) {
    stop("Failed to obtain UniProtKB-ID mappings from the file.")
}

# -----------------------------------------------------
# Map Gene IDs to UniProt Accessions
# -----------------------------------------------------

cat("Mapping gene identifiers to UniProt Accessions...\n")
df_mapped <- df %>%
    left_join(entry_name_to_accession, by = c("gene_symbol" = "UniProtKB_ID")) %>%
    # Replace the original gene_symbol with the Accession if mapping exists
    mutate(
        gene_symbol = if_else(
            !is.na(UniProt_Accession) & UniProt_Accession != "",
            UniProt_Accession,
            gene_symbol
        )
    ) %>%
    # Retain and order the columns for downstream analysis (adjust if necessary)
    dplyr::select(gene_symbol, pval, padj, log2fc)

cat("Preview of gene data AFTER mapping:\n")
print(head(df_mapped))

cat("Final gene count:", nrow(df_mapped), "\n")

# -----------------------------------------------------
# Save Mapped Results
# -----------------------------------------------------

results_file <- file.path(results_dir, paste0(tools::file_path_sans_ext(file_name), "_mapped.csv"))
readr::write_csv(df_mapped, results_file)
cat("Mapped results saved to:", results_file, "\n")

cat("Pipeline completed successfully.\n")
