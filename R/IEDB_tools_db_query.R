#' @title IEDB tools DB Query
#'
#' @description
#' This script processes information from IEDB database and its tools
#' 
#'
#' @details
#' You can retrieve info from IEDB thorugh a query approach. The script is set up to needing exclusively species scientific 
#'names. Some functions were created for epitope prediction of MHC I and MHC II but never properly used (beta version)
#' 
#'
#' @author Pedro Granjo
#' @date 5-JAN-2025
#'
#'
NULL  

############################# Install Packages if necessary################################################
# List of required libraries
libraries <- c(
  "httr", "jsonlite", "tibble", "dplyr", "stringr"
)

# Packages to install (those not already installed)
to_install <- setdiff(libraries, rownames(installed.packages()))

# Install missing packages
if (length(to_install) > 0) {
  install.packages(to_install, repos = "http://cran.us.r-project.org")
}

# Load required libraries
lapply(libraries, library, character.only = TRUE)

##################################################################################################
######################################### IEDB Query #############################################
##################################################################################################

#' @title Verify Validity of Query Parameters
#' @description Checks if the provided query parameters are valid based on a metadata table and general terms.
#' @param parameters A list of query parameters to be checked.
#' @param table_search A string representing the table being queried.
#' @param metadata_table A list of a single or multiple tables metadata for example their columns names, troubleshooting the possibility of invalid parameters.
#' @return A list with a boolean indicating if the parameters are valid, and a message or invalid columns if not.
verify_valid_parameters <- function(parameters, table_search, metadata_table) {
  
  # General terms that are always considered valid
  general_terms <- c("order", "offset", "limit", "Range", "Range_Unit", "Prefer", "select")
  
  # Get the valid parameters from metadata and general terms
  comparatives <- c(general_terms, metadata_table[[table_search]][["valid_parameters"]])
  
  # Check if all parameters are in the valid list
  if (!"order" %in% names(parameters)) {
    # If "order" parameter is missing
    return(list(boolean = FALSE, message = "Missing required 'order' parameter."))
  }
  
  if (all(names(parameters) %in% comparatives)) {
    # All parameters are valid
    return(list(boolean = TRUE))
  } else {
    # Find the invalid parameters
    invalid_columns <- setdiff(names(parameters), comparatives)
    return(list(boolean = FALSE, invalid_columns = invalid_columns, message = "Invalid parameters found."))
  }
}

#' @title Query IEDB Database
#' @description Sends a query to the IEDB database API and retrieves results in a paginated format. The maximum page size is 10,000
#' @param endpoint A string specifying the table name that we aim to query (e.g., 'antigen_search or reference_search').
#' @param query_params A named list of query parameters.
#' @param base_uri The base URI for the API (default: 'https://query-api.iedb.org/').
#' @param page_size The number of records per page (default: 10000).
#' @param metadata_table A list of a single or multiple tables metadata for example their columns names, troubleshooting the possibility of invalid parameters.
#' @return A tibble containing the query results.
iq_query <- function(endpoint, query_params, base_uri = 'https://query-api.iedb.org/', page_size = 10000,metadata_table) {
  
  # Check if query_params is a named list
  if (!is.list(query_params) || is.null(names(query_params))) {
    stop("query_params should be a named list of parameters.")
  }
  
  assesser <- verify_valid_parameters(query_params,endpoint, metadata_table)
  
  if(!(unlist(assesser[["boolean"]]))){
    query_params <- tables_metadata$endpoint[["default_values"]]
    message(paste("Warning: you are using invalid parameters into your query parameter the following columns are invalid:", assesser[["invalid_columns"]]) )
    message("Default values that will return the whole table will be used")
  }
  
  if(endpoint == "reference_search"){
    page_size = 5000
  }
  
  get_text <- 'NA'  # To initiate the while loop
  
  # Create a final tibble format to compile the results
  final_tbl <- tibble()
  
  url <- paste0(base_uri, endpoint)
  
  offset <- 0
  query_params[['offset']] <- format(offset, scientific = F)
  
  
  while (get_text != '[]') {
    #message(paste0("Fetching offset: ", query_params[['offset']]))
    
    tryCatch({
      # Attempt API GET request
      get_1 <- GET(url, query = query_params)
      
      
      if (http_error(get_1)) {
        stop("HTTP error: ", status_code(get_1), " - ", http_status(get_1)$message)
      }
      
      # Parse the content as text
      get_text <- content(get_1, as = 'text', encoding = "UTF-8")
      
      # Convert JSON to tibble (if valid)
      resp_tbl <- tibble::as_tibble(fromJSON(get_text), .name_repair = "unique")
      
      # Add to final tibble
      final_tbl <- bind_rows(final_tbl, resp_tbl)
      rownames(final_tbl) <- NULL
      
      # Update offset and query_params for the next batch
      offset <- offset + page_size
      
      query_params[['offset']] <- format(offset, scientific = FALSE)
      
      Sys.sleep(1)  # sleep for 1 second between calls so as not to overload the server
      
    }, error = function(e) {
      warning("An error occurred: ", conditionMessage(e))
      break  # Exit on error
    }, warning = function(w) {
      message("Warning: ", conditionMessage(w))
    })
  }
  
  # Return the final tibble
  return(final_tbl)
}


#' @title Query Species Data
#' @description Generates query parameters for a specific species and retrieves data from the IEDB database.
#' @param species_name A string specifying the specy name for example Escherichia coli. The specie name has to be according the IEDB nomemclature
#' @param metadata_table A data frame or list containing metadata for query validation.
#' @return A tibble containing the results for the specified species.
species_query <- function(species_name, metadata_table) {
  # Create query parameters dynamically using the species_name
  params <- list(
    parent_source_antigen_iri = 'like.UNIPROT*',
    select = 'parent_source_antigen_id, parent_source_antigen_iri, parent_source_antigen_source_org_iri, parent_source_antigen_names, parent_source_antigen_source_org_name, bcell_ids,tcell_ids,elution_ids',
    offset = 0,
    parent_source_antigen_source_org_name = paste0('like.*', species_name, '*'),  # Insert species name dynamically
    host_organism_names = "cs.{Homo sapiens (human)}",  # Case-insensitive matching
    order = 'parent_source_antigen_id'
  )
  
  # Execute the query using the iq_query function
  result <- iq_query('antigen_search', params, metadata_table = metadata_table)
  
  # Return the result
  return(result)
}



#' @title Split Vector into Batches
#' @description Divides a vector into smaller parts of a specified size.
#' @param ids A vector composed of elements that deemed to be splitted.
#' @param batch_size An integer specifying the size of each batch.
#' @return A list of vectors representing the split batches.
split_into_batches <- function(ids, batch_size) {
  split(as.vector(ids), ceiling(seq_along(as.vector(ids)) / batch_size))
}


#' @title Query Cell Types in Batches
#' @description Queries a table for cell type data in batches and handles timeout constraints.
#' @param id_batches A list of ID batches to query.
#' @param endpoint The endpoint (table name within iedb api) to query (e.g., 'tcell_search').
#' @param v_filtration vertical filtering - It will filter the fields and retrieve them.
#' @param order_f The field by which results are ordered.
#' @param id_field The name of the ID field to be queried.
#' @param rowfiltration Additional filters for the query (default: empty vector).
#' @param timeout The time limit in seconds (default: 4800).
#' @return A list containing processed data from iedb api and any leftover IDs (unprocessed ids due to timeout).
query_celltypes <- function(id_batches, endpoint, v_filtration, order_f, id_field, rowfiltration = c()) {
  
  processed_data <- list()  # To store processed data
  
  # Custom lapply using a for-loop to handle timeout
  for (i in seq_along(id_batches)) {
    
    # Process the current batch
    id_batch <- id_batches[[i]]
    single_str <- paste(id_batch, collapse = ",")
    
    # Prepare query parameters
    params <- list(
      select = v_filtration,
      order = order_f
    )
    params[[id_field]] <- paste0("in.(", single_str, ")")
    
    if (length(rowfiltration) > 0) {
      params <- c(params, rowfiltration)  # Add any additional filters
    }
    
    # Query the endpoint
    result <- iq_query(endpoint, params, metadata_table = tables_metadata)
    
    # Append results to processed data if not NULL
    if (!is.null(result)) {
      processed_data <- append(processed_data, list(result))
    }
  }
  
  # Combine all processed data into a single tibble/data.frame
  combined_data <- bind_rows(processed_data)
  
  # Return both processed data and any leftover IDs
  return(list(processed_data = combined_data))
}

#' @title Clean ID Vector
#' @description Cleans an ID vector by removing NA values.
#' @param ids A vector of IDs.
#' @param id_type A string indicating the type of IDs (e.g., "T-Cell").
#' @param species_name A string specifying the species for context.
#' @return A cleaned vector without NA values.
clean_ids <- function(ids, id_type, species_name) {
  # Count and remove NA values
  removed_na_count <- sum(is.na(ids))
  cleaned_ids <- ids[!is.na(ids)]
  
  # Print the number of NA values removed, if any
  if (removed_na_count > 0) {
    message(paste("Removed", removed_na_count, "NA", id_type, "IDs for species:", species_name))
  }
  
  return(cleaned_ids)
}

#' @title Ensure Non-Empty Data Frame
#' @description Ensures that a data frame has at least one row. Returns a placeholder if empty.
#' @param df A data frame to check.
#' @param required_col A string specifying the name of the placeholder column (default: 'iedb_assay_id').
#' @return The original data frame or a tibble with one empty column.
ensure_non_empty <- function(df, required_col = "iedb_assay_id") {
  if (nrow(df) == 0) {
    return(tibble(!!required_col := character()))
  }
  return(df)
}

#' @title Retrieve Data for a Species
#' @description Retrieves data for a given species using an IEDB query.
#' @param species_name A string specifying the species name.
#' @param tables_metadata A data frame or list containing metadata for query validation.
#' @return A tibble with species data from antigen_search table or NULL if no data is found.
get_species_data <- function(species_name, tables_metadata) {
  species_data <- species_query(species_name, tables_metadata)
  if (nrow(species_data) == 0) {
    
    message(paste("Species", species_name, "returned an empty tibble."))
    return(NULL)
  }
  print(paste("Species", nrow(species_data)))
  return(species_data)
}


#' @title Extract IDs from Data
#' @description Extracts T-Cell, B-Cell, and MHC IDs from species data or leftover IDs.
#' @param species_data A data frame with species data.
#' @param species_name A string specifying the species name.
#' @return A list containing cleaned T-Cell, B-Cell, and MHC IDs.
extract_ids <- function(species_data, species_name) {
  tcell_ids <- clean_ids(unique(unlist(species_data$tcell_ids)), "T-Cell", species_name)
  bcell_ids <- clean_ids(unique(unlist(species_data$bcell_ids)), "B-Cell", species_name)
  mhc_ids <- clean_ids(unique(unlist(species_data$elution_ids)), "MHC", species_name)
  return(list(tcell_ids = tcell_ids, bcell_ids = bcell_ids, mhc_ids = mhc_ids))
}





#' @title Query All Cell Types
#' @description Queries cell type data for T-Cell, B-Cell, and MHC types, respecting a timeout.
#' @param ids A list of IDs (tcell, bcell and mhc) for querying.
#' @return A list containing processed cell data and any leftover IDs.
query_all_celltypes <- function(ids) {
  cell_filtr <- list()
  
  if(length(ids[["mhc_ids"]]) > 0) {
      select_table <- 'elution_id, reference_id, pubmed_id, reference_authors, reference_titles, reference_dates, journal_name, qualitative_measure, assay_names, linear_sequence, linear_sequence_length, mhc_allele_name'
      query_res <- query_celltypes(split_into_batches(ids[["mhc_ids"]], batch_size = 200), 'mhc_search', select_table,'elution_id', 'elution_id', list(host_organism_name = "like.*Homo sapiens (human)*"))
      cell_filtr$mhc <- query_res[["processed_data"]]
    }

  
  
  if(length(ids[["tcell_ids"]]) > 0) {
      select_table <- 'tcell_id, reference_id, pubmed_id, reference_authors, reference_titles, reference_dates, journal_name, qualitative_measure, assay_names, linear_sequence, linear_sequence_length'
      query_res <- query_celltypes(split_into_batches(ids[["tcell_ids"]], batch_size = 200), 'tcell_search',select_table ,'tcell_id', 'tcell_id', list(host_organism_name = "like.*Homo sapiens (human)*"))
      cell_filtr$tcell <- query_res[["processed_data"]]
    }

  
  if(length(ids[["bcell_ids"]]) > 0){
      select_table <- 'bcell_id, reference_id, pubmed_id, reference_authors, reference_titles, reference_dates, journal_name, qualitative_measure, assay_names, linear_sequence, linear_sequence_length'
      query_res <- query_celltypes(split_into_batches(ids[["bcell_ids"]], batch_size = 200), 'bcell_search',select_table, 'bcell_id', 'bcell_id', list(host_organism_name = "like.*Homo sapiens (human)*"))
      cell_filtr$bcell <- query_res[["processed_data"]]
    }

  
  return(list(processed_data = cell_filtr))
}


#' @title Combine Cell Type Data
#' @description Combines data from multiple cell type queries into a unified data frame.
#' @param cell_filtr A list of data frames for each cell type.
#' @return A combined data frame or NULL if no data is available.
combine_cell_data <- function(cell_filtr) {
  # Define the specific ID columns we are looking for
  id_column_map <- list(
    tcell_id = "iedb_assay_id",
    bcell_id = "iedb_assay_id",
    elution_id = "iedb_assay_id"
  )
  
  # Rename the ID columns to "iedb_assay_id"
  cell_filtr <- lapply(cell_filtr, function(df) {
    if (!is.null(df)) {
      # Find the column names that match the specific ID columns and rename them
      for (old_name in names(id_column_map)) {
        if (old_name %in% colnames(df)) {
          colnames(df)[colnames(df) == old_name] <- id_column_map[[old_name]]
        }
      }
    }
    return(df)
  })
  
  # Combine the processed data frames
  df_combined <- bind_rows(lapply(cell_filtr, function(x) if (!is.null(x)) x))
  # Return NULL if no rows were combined, else return the combined dataframe
  if (nrow(df_combined) == 0) return(NULL)
  return(df_combined)
}



#' @title Generate Reference Table
#' @description Creates a reference table by filtering out assay-specific columns.
#' @param df_combined A combined data frame of assay data.
#' @return A data frame with distinct reference data.
generate_reference_table <- function(df_combined) {
  reference_table <- df_combined[, !(colnames(df_combined) %in% c("iedb_assay_id", "qualitative_measure", "assay_names", "mhc_allele_name", "linear_sequence", "linear_sequence_length"))]
  return(distinct(reference_table))
}



combine_assay_info <- function(species_data,
                               df_combined,
                               start_index = NULL, start_time = NULL, timeout = NULL, all_ids = NULL, run_id = NULL
                               ) {
  ids_df <- species_data %>%
      select(parent_source_antigen_iri, bcell_ids, tcell_ids, elution_ids) %>%
      rowwise() %>%
      mutate(valid_ids = list(unlist(c(bcell_ids, tcell_ids, elution_ids))[unlist(c(bcell_ids, tcell_ids, elution_ids)) %in% df_combined$iedb_assay_id])) %>%
      ungroup() %>%
      filter(lengths(valid_ids) > 0)

  species_data_cleaned <- species_data[, !(colnames(species_data) %in% c("bcell_ids", "tcell_ids", "elution_ids"))]
  
  # Initialize a list to store the results from processing each row of ids_df
  results_list <- list()
  processed_count <- 0
  
  # Process each row of ids_df using apply (or rowwise/mutate)
  
  processed_results <- apply(ids_df, 1, function(row) {
    antigen_iri <- row[["parent_source_antigen_iri"]]
    row_ids <- row[["valid_ids"]]
    
    row_result <- process_assay_row(row_ids, df_combined)
    row_result$parent_source_antigen_iri <- antigen_iri
    
    return(row_result)
  })
  
  final_combined_info <- bind_rows(processed_results)
  
  combined_info <- merge(final_combined_info, species_data_cleaned, by = "parent_source_antigen_iri", all.x = TRUE)
  
  saveRDS(combined_info, file = "combined_info.rds")
  
  message("All assays processed and combined.")
  
  return(
    combined_info = combined_info)
}
  

#' @title Process Row IDs for Assay Data
#' @description Processes a vector of row IDs to extract assay-related information from a combined data frame.
#' @param row_ids A vector of IDs to process.
#' @param df_combined A data frame containing combined assay data.
#' @return A tibble with collapsed information about assays, including IDs, references, qualitative measures, and sequences per HCP antigen
process_assay_row <- function(row_ids, df_combined) {
  # Initialize empty tibble to store results
  df_info <- tibble(
    iedb_assay_id = character(), 
    reference_id = character(), 
    qualitative_measure = character(), 
    assay_names = character(),
    linear_sequence = character(),
    linear_sequence_length = character(),
    mhc_allele_name = character()
  )
  
  # Initialize variables to collect and collapse multiple values
  iedb_assay_id <- character()
  collected_reference_ids <- character()
  collected_qualitative_measure <- character()
  collected_assay_names <- character()
  collected_linear_sequence <- character()
  collected_linear_sequence_length <- character()
  collected_mhc_allele_names <- character()
  
  # Iterate over each row ID in the current antigen row
  sapply(seq_along(row_ids), function(i) {
    iedb_id <- row_ids[i]
    
    # Filter the current row in df_combined for the current iedb_id
    inst_row <- df_combined[df_combined$iedb_assay_id == iedb_id, ]
    
    # Extract and append relevant fields
    reference_id <- ifelse(
      inst_row$qualitative_measure == "Positive-Low", paste(inst_row$reference_id, "PL", sep ="_"),
      ifelse(inst_row$qualitative_measure == "Positive-High", paste(inst_row$reference_id, "PH", sep ="_"),
             ifelse(inst_row$qualitative_measure == "Positive", paste(inst_row$reference_id, "P", sep ="_"),
                    ifelse(inst_row$qualitative_measure == "Positive-Intermediate", paste(inst_row$reference_id, "PI", sep ="_"),
                           inst_row$reference_id)))
    )
    qualitative_measure <- inst_row$qualitative_measure
    assay_names <- inst_row$assay_names
    mhc_allele_name <- inst_row$mhc_allele_name
    linear_sequence <- inst_row$linear_sequence
    linear_sequence_length <- inst_row$linear_sequence_length
    
    iedb_assay_id <<- c(iedb_assay_id, iedb_id)
    collected_reference_ids <<- c(collected_reference_ids, reference_id)
    collected_qualitative_measure <<- c(collected_qualitative_measure, qualitative_measure)
    collected_assay_names <<- c(collected_assay_names, assay_names)
    collected_linear_sequence <<- c(collected_linear_sequence, linear_sequence)
    collected_linear_sequence_length <<- c(collected_linear_sequence_length, linear_sequence_length)
    
    if (!is.na(mhc_allele_name) && !is.null(mhc_allele_name)) {
      collected_mhc_allele_names <<- c(collected_mhc_allele_names, paste(mhc_allele_name, iedb_id, sep = "_"))
    }
  })
  
  # Collapse the collected information into a single row, concatenating multiple values with commas
  collapsed_iedb_assay_id <- paste(iedb_assay_id, collapse = ", ")
  collapsed_reference_ids <- paste(unique(collected_reference_ids), collapse = ", ")
  collapsed_qualitative_measure <- paste(collected_qualitative_measure, collapse = ", ")
  collapsed_assay_names <- paste(collected_assay_names, collapse = ", ")
  collapsed_mhc_allele_names <- paste(collected_mhc_allele_names, collapse = ", ")
  collapsed_linear_sequence <- paste(collected_linear_sequence, collapse = ", ")
  collapsed_linear_sequence_length <- paste(collected_linear_sequence_length, collapse = ", ")
  
  # Add the collapsed information to df_info
  df_info <- rbind(df_info, data.frame(
    iedb_assay_id = collapsed_iedb_assay_id,
    reference_id = collapsed_reference_ids,
    qualitative_measure = collapsed_qualitative_measure,
    assay_names = collapsed_assay_names,
    mhc_allele_name = collapsed_mhc_allele_names,
    linear_sequence = collapsed_linear_sequence,
    linear_sequence_length = collapsed_linear_sequence_length,
    stringsAsFactors = FALSE
  ))
  
  return(df_info)
}


#' @title Merge Partial Combined Information
#' @description Merges new partial combined information with previously saved results, if available.
#' @param new_partial_combined_info A data frame with new partial combined information.
#' @param partial_results A data frame with previously saved partial results (optional).
#' @return A merged data frame combining the new and existing partial information.
merge_partial_combined_info <- function(new_partial_combined_info, partial_results) {
  
  if (!is.null(partial_results)) {
    
    partial_combined_info <- partial_results
    
    message("Merging new partial_combined_info with previously saved partial results")
    
    # Combine the new partial results with the old partial_combined_info
    merged_partial_combined_info <- bind_rows(partial_combined_info, new_partial_combined_info)

    return(merged_partial_combined_info)
  } else {
    message("No previous partial_combined_info results found. Proceeding with new data only.")
    return(new_partial_combined_info)  # No partial results, return the new data
  }
}

#' @title Process Accession Column
#' @description Extracts and processes the accession part from a column in a combined data frame.
#' @param combined_info A data frame containing assay data with a 'parent_source_antigen_iri' column.
#' @return The data frame with a processed 'Accession' column.
process_accession_column <- function(combined_info) {
  Process <- lapply(str_split(combined_info$parent_source_antigen_iri, ":"), function(accession) {
    return(accession[2])
  })
  
  colnames(combined_info)[colnames(combined_info) == "parent_source_antigen_iri"] <- "Accession"
  combined_info$Accession <- do.call(rbind, Process)
  return(combined_info)
}


#' @title Process Information for a Species
#' @description Processes assay data for a specified species, handling IDs, queries.
#' @param species_name The name of the species to process.
#' @param tables_metadata Metadata for table validation (default: NULL).
#' @param df_combined A data frame of previously combined data (default: NULL) from a previous run of this function
#' @return A list containing the reference table, combined info, and leftover IDs.
process_species_info <- function(species_name, 
                                 tables_metadata =NULL) {

  
  # 1. Fetch species_data 
  species_data <- get_species_data(species_name, tables_metadata)
  
  if (is.null(species_data)) {
    return(list(reference_table = NULL, combined_info = NULL, empty_species = species_name))
  }
  
  message(paste(species_name, "Fetching of species tables was a success!"))
  
  # 2. Extract IDs (directly from species_data)
  ids <- extract_ids(species_data, species_name)
  
  # 3. Query the data for each cell type
  query_result <- query_all_celltypes(ids)
  cell_filtr <- query_result[["processed_data"]]
  
  # 4. Combine newly processed data with the existing df_combined (if provided)
  df_combined <- combine_cell_data(cell_filtr)
  
  # 5. Call `combine_assay_info` and pass in processed_ids
  combined_info_result <- combine_assay_info(
    species_data, df_combined= df_combined)
  
  # 6. Continue processing the combined data if all was successful
  combined_info <- combined_info_result
  combined_info <- process_accession_column(combined_info)
  
  #Generate the reference table from combined_info
  reference_table <- generate_reference_table(df_combined)
  
  # 7. Return the final combined data
  return(list(
    reference_table = reference_table, 
    combined_info = combined_info,
    df_combined = df_combined
  ))
}



#' @title Process Multiple Species
#' @description Processes information for multiple species and combines the results.
#' @param species_names A vector of species names to process.
#' @param tables_metadata Metadata for table validation.
#' @return A list containing the final combined reference table and combined assay data for all species.
process_multiple_species <- function(species_names, tables_metadata) {
  # Initialize empty lists to store results for all species
  all_reference_tables <- list()
  all_combined_info <- list()
  
  # Loop through each species name and apply the process_species_info function
  for (species_name in species_names) {
    species_result <- process_species_info(species_name, tables_metadata)
    
    # Append the reference table and combined info to the lists
    all_reference_tables[[species_name]] <- species_result$reference_table
    all_combined_info[[species_name]] <- species_result$combined_info
  }
  
  # Combine all reference tables and combined info into single tables
  final_reference_table <- bind_rows(all_reference_tables)
  final_combined_info <- bind_rows(all_combined_info)
  
  # Return the final combined tables
  return(list(final_reference_table = final_reference_table, final_combined_info = final_combined_info))
}


############################################################################################################
### Data Integration

#' @title Check Qualitative Measure
#' @description Checks for matches of positive values within the 'qualitative_measure' column of a data frame.
#' @param data A data frame containing a 'qualitative_measure' column.
#' @param positive_values A vector of strings to match against the 'qualitative_measure' column.
#' @return A logical vector indicating matches.
check_qualitative_measure <- function(data, positive_values) {
  boolean <- c()  # Initialize empty vector to store results
  
  # Apply the matching logic to each row of the specified column
  apply(as.data.frame(data[["qualitative_measure"]]), 1, function(x) {
    boolean <<- c(boolean, any(sapply(positive_values, grepl, x)))
  })
  
  return(boolean)
}


#' @title Clean Protein Name
#' @description Cleans a protein name by extracting the part before a specified pattern (e.g., '(Uni').
#' @param protein_string A string representing the protein name.
#' @return The cleaned protein name.
clean_protein_name <- function(protein_string){
  target_str <- protein_string
  return(str_split(target_str, fixed("(Uni"))[[1]][1])
}

#' @title Process Accession and Protein Name
#' @description Processes data to subset rows based on qualitative measures, clean protein names, and return a data frame with 'Accession' and 'Protein_name'.
#' @param data A data frame containing assay data with columns 'qualitative_measure' and 'parent_source_antigen_names'.
#' @param positive_values A vector of strings to match within the 'qualitative_measure' column.
#' @return A data frame with 'Accession' and cleaned 'Protein_name' columns.
process_accession_and_protein <- function(data, positive_values) {
  boolean_vector <- check_qualitative_measure(data,positive_values)
  # Subset the data where boolean_vector is TRUE
  subset_data <- data[boolean_vector, ]
  name_list <- lapply(subset_data$parent_source_antigen_names, function(x){
    if(length(x) > 1){
      cleared_name <- c()
      for(i in x){
        cleared_name <- c(cleared_name, clean_protein_name(i))
      }
      final_name <- paste(cleared_name, collapse = ", other name = ")
    } else{
      clean_protein_name(x)
    }
  })
  protein_names <- do.call(rbind,name_list)
  # Create a new data frame with Accession and Protein_name columns
  new_accession <- data.frame(Accession = unlist(subset_data$Accession), Protein_name =  protein_names, stringsAsFactors = FALSE)
  return(new_accession)
}


##############################################################################################################
######################################## MHC Elution binding tools ###########################################
##############################################################################################################
########################### MHC I


#' @title Prepare MHC I Sequences
#' @description Prepares and formats MHC I sequences for prediction, organizing them into FASTA-like format.
#' @param names_vector A vector of protein names.
#' @param sequence_vector A vector of corresponding sequences.
#' @return A list of formatted sequence batches.
prepare_mhci_sequences <- function(names_vector, sequence_vector){
  #Remove the rows that don't have an available sequence
  processed_names_vector <- names_vector[!(sequence_vector == "Not Available")]
  processed_sequence_vector <- sequence_vector[!(sequence_vector == "Not Available")]
  
  
  small_df <- data.frame(names = processed_names_vector, sequences = processed_sequence_vector, stringsAsFactors = FALSE)
  
  #I aggreaged the acession numbers into a single string for sequences that at Uniprot are the same (most likely Homologous)
  merged_df <- aggregate(names ~ sequences, data = small_df , FUN = function(x) paste(x, collapse = ", "))
  
  fasta_like_input <- c()
  
  for(i in seq_along(merged_df$sequences)){
    if(grepl("X",merged_df$sequences[i])){
      #The app doesn't like X inside the aminoacid sequence
      string_split <- strsplit(merged_df$sequences[i], split = 'X')[[1]]  # Get the first element as the split list
      cleaned_strings <- string_split[string_split != ""]  # Remove empty parts
      
      for(j in seq_along(cleaned_strings)){
        fasta_like_input <- c(fasta_like_input, paste0(">Uniprot:", merged_df$names[i], "\n", cleaned_strings[j]))
      }
      
      #fasta_like_input <- c(fasta_like_input, c(paste0(">","Uniprot:",  merged_df$names[i],"\n", cleaned_strings)))
      
    } else{
      fasta_like_input <- c(fasta_like_input, paste0(">","Uniprot:",  merged_df$names[i],"\n", merged_df$sequences[i]))
    }
    
  }
  
  sequences <- split_into_batches(fasta_like_input, batch_size = 150)
  results <- lapply(sequences , function(sequence)
    single_str <- paste(sequence, collapse ="\n"))
  
  return(results)
  
}


#' @title Split Alleles into Groups
#' @description Splits a string of alleles into groups of a specified size.
#' @param alleles_str A string containing alleles separated by commas.
#' @param group_size An integer specifying the group size.
#' @return A list of allele groups.
split_alleles <- function(alleles_str, group_size) {
  alleles <- unlist(strsplit(alleles_str, ","))
  split(alleles, ceiling(seq_along(alleles) / group_size))
}

#' @title MHC I Prediction
#' @description Submits a sequence for MHC I prediction using the IEDB API and returns the API response.
#' @param input_sequence_text A string of input sequences in FASTA format.
#' @param email User's email address for job notifications.
#' @param pipeline_title A title for the prediction job.
#' @param peptide_length_range A numerical list specifying the range of peptide lengths.
#' @param alleles A vector of alleles for the prediction.
#' @return The response from the API containing the result ID.
mhci_prediction <- function(input_sequence_text, email, pipeline_title, peptide_length_range, alleles) {
  
  api_url <- "https://api-nextgen-tools.iedb.org/api/v1/pipeline"
  
  # Define the pipeline request data
  request_data <- list( 
    pipeline_title = pipeline_title,  # Custom job title
    email = email,  # User's email
    run_stage_range = list(1, 1),  # Running stage range
    stages = list(
      list(
        stage_number = 1,
        tool_group = "mhci",  # Specify the tool group as mhci
        input_sequence_text = input_sequence_text,  # Input sequence
        input_parameters = list(
          alleles = alleles,  # Alleles for testing
          peptide_length_range = peptide_length_range,  # Peptide length range
          predictors = list(
            list(
              type = "binding",
              method = "netmhcpan_el"  # Predictor method
            ), list( type= "immunogenicity",
                     method = "immunogenicity"
            )
          )
        )
      )
    )
  )
  
  # Convert request body to JSON format
  request_body <- toJSON(request_data, auto_unbox = TRUE)
  
  # Send POST request to the API
  response <- POST(
    url = api_url,
    body = request_body,
    encode = "json",  # Send as JSON
    content_type_json()
  )
  
  if (http_status(response)$category == "Success") {
    print("Pipeline submission succeeded!")
    
    api_response <- content(response, as = "parsed", type = "application/json")
    return(api_response)
    
  } else {
    print("Pipeline submission failed.")
    print(response)
    return(NULL)
  }
}

#' @title Retrieve MHC-I Pipeline Results
#'
#' @description
#' Retrieves pipeline results using `result_id` from an API response.
#'
#' @param api_response A list containing the `result_id` for the requested pipeline result.
#' @return Parsed result data if successful; otherwise, NULL on failure.
mhci_pipeline_results <- function(api_response) {
  
  # Extract the result_id from the api_response
  result_id <- api_response$result_id
  
  # Construct the URL to get the results
  result_url <- paste0("https://api-nextgen-tools.iedb.org/api/v1/results/", result_id)
  
  # Send GET request to the API
  response <- GET(
    url = result_url,
    add_headers("accept" = "application/json")
  )
  
  # Check response status
  if (http_status(response)$category == "Success") {
    print("Pipeline result retrieval succeeded!")
    
    # Parse and return the result data
    result_data <- content(response, as = "parsed", encoding = "UTF-8")
    return(result_data)
    
  } else {
    # If request failed, print the error message
    print("Pipeline result retrieval failed.")
    print(content(response, as = "text"))
    
    # Return NULL to indicate failure
    return(NULL)
  }
}


#' @title Process MHC-I Pipeline Results
#'
#' @description
#' Processes and organizes semi-structured MHC-I pipeline results into a more structured format.
#'
#' @param result_data A list containing raw result data from the MHC-I pipeline API.
#' @return A list with cleaned data frames for each result type, including `peptide_table`, 
#'   `allele_distance`, and `input_sequence_table`.
process_mhci_results <- function(result_data) {
  global_results <- list()
  
  # "peptide_table" data processing
  if (result_data$data$results[[1]]$type == "peptide_table") {
    peptide_data <- do.call(rbind, result_data$data$results[[1]]$table_data)
    peptide_columns <- do.call(bind_rows, result_data$data$results[[1]]$table_columns)
    future_colnames <- unlist(distinct(peptide_columns[, 1]))
    future_colnames[13] <- "immunogenicity"
    colnames(peptide_data) <- future_colnames
    
    peptide_data <- as.data.frame(peptide_data)
    filtered_peptide_data <- peptide_data[peptide_data$score > 0 & peptide_data$percentile == 1 ,]
    
    # Store as "peptide" in global results
    global_results[["peptide_table"]] <- filtered_peptide_data
  }
  
  # 2. "netmhcpan_allele_distance" data processing
  if (result_data$data$results[[2]]$type == "netmhcpan_allele_distance") {
    allele_distance_data <- do.call(rbind, result_data$data$results[[2]]$table_data)
    allele_distance_columns <- do.call(bind_rows, result_data$data$results[[2]]$table_columns)
    colnames(allele_distance_data) <- as.vector(unlist(distinct(allele_distance_columns[, 1])))
    rownames(allele_distance_data) <- paste("allele", 1:nrow(allele_distance_data), sep = "_")
    
    # Store as "allele_distance" in global results
    global_results[["allele_distance"]] <- as.data.frame(allele_distance_data)
  }
  
  # 3. "input_sequence_table" data processing
  if (result_data$data$results[[3]]$type == "input_sequence_table") {
    sequence_data <- do.call(rbind, result_data$data$results[[3]]$table_data)
    sequence_columns <- do.call(bind_rows, result_data$data$results[[3]]$table_columns)
    colnames(sequence_data) <- as.vector(unlist(distinct(sequence_columns[, 1])))
    sequence_data <- as.data.frame(sequence_data)
    
    sequence_numbers <- unlist(filtered_peptide_data$sequence_number)
    filtered_input_sequence_table <- as.data.frame(sequence_data)[sequence_data$sequence_number %in% sequence_numbers, ]
    
    filtered_input_sequence_table$sequence_id <- paste(
      filtered_input_sequence_table$sequence_number, 
      filtered_input_sequence_table$sequence_name, 
      sep = "_"
    )
    # Store as "input_sequence" in global results
    global_results[["input_sequence_table"]] <- filtered_input_sequence_table
  }
  
  if (!is.null(filtered_peptide_data) && !is.null(filtered_input_sequence_table)) {
    # Ensure columns are atomic before merging
    filtered_input_sequence_table$sequence_number <- unlist(filtered_input_sequence_table$sequence_number)
    filtered_peptide_data$sequence_number <- unlist(filtered_peptide_data$sequence_number)
    
    # Merge sequence_id from filtered_input_sequence_table into filtered_peptide_data
    filtered_peptide_data <- merge(filtered_peptide_data, 
                                   filtered_input_sequence_table[, c("sequence_number", "sequence_id")], 
                                   by = "sequence_number", all.x = TRUE)
    
    # Update global results with the new filtered_peptide_data containing sequence_id
    global_results[["peptide_table"]] <- filtered_peptide_data
  }
  # Return the global results
  return(global_results)
}




#' @title Submit MHCI Pipeline Predictions
#' 
#' @description
#' Submits MHCI pipeline predictions for multiple sequences and allele groups, storing 
#' responses for each combination.
#' 
#' @param input_sequence A list of sequences to submit for MHCI prediction.
#' @param alleles_groups A list of allele groups, where each element is a vector of alleles.
#' @param email A character string containing the user’s email.
#' @param title A title for the pipeline submission.
#' @param peptide_range A character string specifying the peptide range.
#' @return A list of API responses, keyed by sequence and allele group indices.
#' @export
submit_mhci_predictions <- function(input_sequence, alleles_groups, email, title, peptide_range) {
  api_responses <- list()  # Initialize storage for responses
  
  for (i in seq_along(input_sequence)) {
    for (j in seq_along(alleles_groups)) {
      alleles_subset <- paste(alleles_groups[[j]], collapse = ",")
      
      # Submit MHCI prediction for the current sequence and allele group
      api_response <- mhci_prediction(input_sequence[[i]], email, title, peptide_range, alleles_subset)
      
      # Store response
      api_responses[[paste0("sequence_", i, "_alleles_", j)]] <- api_response
    }
  }
  
  return(api_responses)
}


#' @title Process Stored MHCI API Responses
#'
#' @description
#' Iterates through API responses, retrieves the results, checks their status, and processes
#' them if complete. Handles pending and failed requests.
#'
#' @param api_responses A list of stored API responses from `submit_mhci_predictions`.
#' @param input_sequence A list of sequences submitted for MHCI prediction.
#' @param alleles_groups A list of allele groups submitted for MHCI prediction.
#' @return A list of processed results, including status checks for each sequence and allele group.
#'
#' @examples
#' results_data <- process_mhci_responses(api_responses, input_sequence, alleles_groups)
#'
#' @export
process_mhci_responses <- function(api_responses, input_sequence, alleles_groups) {
  results_data <- list()  # Initialize storage for results
  
  for (i in seq_along(input_sequence)) {
    for (j in seq_along(alleles_groups)) {
      api_response <- api_responses[[paste0("sequence_", i, "_alleles_", j)]]
      
      # Retrieve results and check status
      result_data <- mhci_pipeline_results(api_response)
      
      if (result_data$status == "done") {
        print("Good to Go!")
        if (length(result_data$data$results[[1]]$table_data) != 0) {
          results_data[[paste0("sequence_", i, "_alleles_", j)]] <- process_mhci_results(result_data)
        } else {
          message(paste("input_sequence", i, "and allele group", j, "didn't work"))
        }
      } else if (result_data$status == "pending") {
        message("Still Pending - Wait a couple more minutes")
        results_data[[paste0("sequence_", i, "_alleles_", j)]] <- result_data
      }
    }
  }
  
  return(results_data)
}



#' @title Combine MHCI Results Data
#'
#' @description
#' Combines processed MHCI results, merging peptide, allele, and sequence tables.
#'
#' @param results_data A list of processed MHCI results from `process_mhci_responses`.
#' @return A list containing combined `peptide_table`, `allele_distance`, and `input_sequence_table`.
#'
#' @examples
#' combined_results <- combine_mhci_results(results_data)
#' print(combined_results$combined_peptide_table)
#'
#' @export
combine_mhci_results <- function(results_data) {
  combined_peptide_table <- data.frame()  # Initialize empty data frames
  combined_allele_distance <- data.frame()
  combined_input_sequence_table <- data.frame()
  
  count <- 1  # Counter for tracking runs
  
  lapply(results_data, function(run) {
    # Combine peptide_table
    if (!is.null(run[["peptide_table"]]) && nrow(run[["peptide_table"]]) > 0) {
      combined_peptide_table <<- rbind(combined_peptide_table, run[["peptide_table"]])
    } else {
      print(paste("Houston, we have a problem with peptide_table in run", count))
    }
    
    # Combine allele_distance
    if (!is.null(run[["allele_distance"]]) && nrow(run[["allele_distance"]]) > 0) {
      combined_allele_distance <<- rbind(combined_allele_distance, run[["allele_distance"]])
    } else {
      print(paste("Houston, we have a problem with allele_distance in run", count))
    }
    
    # Combine input_sequence_table
    if (!is.null(run[["input_sequence_table"]]) && nrow(run[["input_sequence_table"]]) > 0) {
      combined_input_sequence_table <<- rbind(combined_input_sequence_table, run[["input_sequence_table"]])
    } else {
      print(paste("Houston, we have a problem with input_sequence_table in run", count))
    }
    
    count <<- count + 1  # Increment count
  })
  
  # Return combined data frames in a list
  return(list(
    combined_peptide_table = combined_peptide_table,
    combined_allele_distance = combined_allele_distance,
    combined_input_sequence_table = combined_input_sequence_table
  ))
}

########################### MHC II

#' @title MHC-II Prediction Request
#'
#' @description
#' Sends a POST request to the MHC-II prediction API with specified parameters.
#'
#' @param sequence_text A character string representing the protein sequence for prediction.
#' @param alleles A character string of comma-separated alleles for prediction.
#' @param method Prediction method (check iedb mhcii prediction tools online) (default is "recommended").
#' @param method_version Optional character string for method version.
#' @param length Optional string for peptide length. You can add multiple peptide length separated by commas
#' @param email_address Optional email address to receive prediction results.
#' @return A character string containing the prediction results or an error message if failed.
mhcii_prediction <- function(
    sequence_text,       
    alleles,             
    method = "recommended",  
    method_version = NULL,
    length = NULL,        
    email_address = NULL 
) {
  
  # API URL is constant
  api_url <- "http://tools-cluster-interface.iedb.org/tools_api/mhcii/"
  
  # Construct the API request body
  body <- list(
    method = method,
    sequence_text = sequence_text,
    allele = alleles
  )
  
  # Add method version if provided
  if (!is.null(method_version)) {
    body$method <- paste0(method, "-", method_version)
  }
  
  # Add peptide length if provided
  if (!is.null(length)) {
    body$length <- length
  }
  
  # Add email address if provided
  if (!is.null(email_address)) {
    body$email_address <- email_address
  }
  
  # Perform the POST request to the API
  response <- POST(
    url = api_url,
    body = body,
    encode = "form"  # Form-encoded data
  )
  
  # Check if the request succeeded
  if (http_status(response)$category == "Success") {
    # Return the response content as text
    content <- content(response, as = "text")
    return(content)
  } else {
    # Return error if failed
    return(paste("Error:", http_status(response)$message))
  }
}


#' @title Process a Single Sequence with Accession for MHC-II
#'
#' @description
#' Processes a single sequence for MHC-II prediction, adding the result to a data frame with the accession.
#'
#' @param accession A character string for the protein accession number.
#' @param sequence The protein sequence as a character string.
#' @param alleles A character string of comma-separated alleles for prediction.
#' @param method Prediction method (check iedb mhcii prediction tools online) (default is "recommended").
#' @param method_version Optional character string for method version.
#' @param length Optional string for peptide length. You can add multiple peptide length separated by commas
#' @param email_address Optional email address to receive prediction results.
#' @return A list with the processed result data frame or NULL if no result was returned.
process_single_sequence_with_accession <- function(accession, sequence, alleles, method = "recommended", method_version = NULL, length = NULL, email_address = NULL) {
  
  # Call the mhcii_prediction function for the single sequence
  result <- mhcii_prediction(
    sequence = sequence,   # Provide the single sequence
    alleles = alleles,                  # Alleles list
    method = method,                    # Prediction method
    method_version = method_version,    # Method version (optional)
    length = length,                    # Peptide lengths (optional)
    email_address = email_address       # Email to receive results (optional)
  )
  
  # Check if the result is empty or NULL
  if (is.null(result) || result == "") {
    # Return an empty data frame with a warning or log
    warning(paste("No result returned for accession:", accession))
    return(list(empty_df = c(accession)))
  }
  
  # Split the result into lines
  result_lines <- strsplit(result, "\n")[[1]]
  
  # Check if there are enough lines to form a valid table
  if (length(result_lines) < 2) {
    # If the result doesn't contain enough rows, return an empty data frame
    return(list(empty_df = c(accession)))
  }
  
  # Process the first line as column names
  column_names <- strsplit(result_lines[1], "\t")[[1]]
  
  #print(column_names)
  
  # Process the remaining lines as data rows
  data_rows <- lapply(result_lines[-1], function(x) strsplit(x, "\t")[[1]])
  
  # Create a data frame with the result
  result_df <- as.data.frame(do.call(rbind, data_rows), stringsAsFactors = FALSE)
  
  # Assign the column names
  colnames(result_df) <- column_names
  
  # Add the accession number as a new column to each row
  result_df$Accession <- accession
  
  # Return the processed data frame
  return(list(single_accession_immuno_df = result_df, empty_df = NULL))
}



#' @title Process Multiple Sequences for MHC-II
#'
#' @description
#' Iterates over a data frame of sequences, processes each sequence for MHC-II prediction.
#'
#' @param df A data frame with columns `accession` and `sequence`.
#' @param alleles Character string of comma-separated alleles.
#' @param method The prediction method to use (default: "recommended").
#' @param method_version Optional method version.
#' @param length Optional peptide length. (string format)
#' @param email_address Optional email address to receive results.
#' @param past_result_df Optional data frame to combine with current results.
#' @param timeout Numeric value in seconds for the maximum time before stopping the process.
#' @return A list with processed results as a data frame and a vector of missing accessions.
process_sequences_mhcii_df <- function(df, alleles, method = "recommended", method_version = NULL, length = NULL, email_address = NULL, past_result_df = NULL, timeout = 700) {
  start_time <- Sys.time()
  
  # Initialize an empty list to store results and empty data frame records
  results_list <- list()
  empty_df <- c()
  
  # Loop over the rows of the dataframe
  for (i in 1:nrow(df)) {
    
    # Extract accession and sequence for the current row
    accession <- df[i, "accession"]
    sequence <- df[i, "sequence"]
    
    # Call the single sequence processing function
    result <- process_single_sequence_with_accession(
      accession = accession,
      sequence = sequence,
      alleles = alleles,
      method = method,
      method_version = method_version,
      length = length,
      email_address = email_address
    )
    
    if(!is.null(result$empty_df)){
      empty_df <- c(empty_df, result[["empty_df"]])
    } else{
      results_list[[i]] <- result[["single_accession_immuno_df"]]
    }
    
    
    if (i %% 50 == 0) {
      Sys.sleep(60)
    }
    
    
    current_time <- Sys.time()
    elapsed_time <- as.numeric(difftime(current_time, start_time, units = "secs"))
    
    if (elapsed_time > timeout) {
      # If the timeout is exceeded, return the current results so far
      if (is.null(past_result_df)) {
        results_list <- do.call(rbind, results_list)
        return(list(current_result_list = results_list, all_empty_df = empty_df))
      } else {
        
        results_df <- do.call(rbind,results_list)
        combined_results <-rbind(past_result_df,results_df)
        return(list(current_result_list = combined_results, all_empty_df = empty_df))
      }
    }
    
    
  }
  
  saveRDS(results_list, file = "notcombined_results.rds")
  message("Houston we saved the document!")
  # Combine all the results into a single data frame
  combined_results <- do.call(rbind, results_list)
  
  # Return the combined results and any empty data frames
  return(list(mhcii_results = combined_results, accession_research = empty_df))
}



#' @title Filter Immunogenic Peptides
#'
#' @description
#' Filters peptides for immunogenicity based on score and immunogenicity thresholds.
#'
#' @param maindatabase A data frame with columns including `Accession` and `Critical_for_short`.
#' @param final_predictive_peptide_table A data frame containing predicted peptides.
#' @param score_margin Numeric score margin for filtering (default: 0.10).
#' @param immune_margin Numeric immunogenicity margin for filtering (default: 0.10).
#' @return A data frame of peptides that meet the immunogenicity criteria.
filter_immunogenic_peptides <- function(maindatabase, final_predictive_peptide_table, score_margin = 0.10, immune_margin = 0.10) {
  
  # Identify immunogenic HCPs based on "Immunogenic" value in `Critical_for_short`
  Immuno_accession <- maindatabase[maindatabase$Critical_for_short == "Immunogenic", "Accession"]
  
  # Filter entries in `final_predictive_peptide_table` with matching accessions
  int_val <- final_predictive_peptide_table[final_predictive_peptide_table$Accession %in% Immuno_accession, ]
  
  # Calculate threshold values with 10% added margin
  median_threshold_score <- median(as.numeric(int_val$score)) * (1 + score_margin)
  
  # Check if `immunogenicity` column exists
  if ("immunogenicity" %in% colnames(final_predictive_peptide_table)) {
    median_threshold_immune_score <- median(as.numeric(int_val$immunogenicity)) * (1 + immune_margin)
    
    # Filter based on both score and immunogenicity
    bench_marked_immuno_table <- final_predictive_peptide_table[
      as.numeric(final_predictive_peptide_table$immunogenicity) > median_threshold_immune_score &
        as.numeric(final_predictive_peptide_table$score) > median_threshold_score, 
    ]
  } else {
    # Filter based only on score if `immunogenicity` column is not available
    bench_marked_immuno_table <- final_predictive_peptide_table[
      as.numeric(final_predictive_peptide_table$score) > median_threshold_score, 
    ]
  }
  
  # Return the filtered immunogenic peptide table
  return(bench_marked_immuno_table)
}



#############################################################################################################
######################################### Antibody Prediction ###############################################
#############################################################################################################

#' @title Antibody Epitope Prediction
#'
#' @description
#' Sends a POST request to the antibody epitope prediction API using the specified method.
#'
#' @param sequence_text A character string of the protein sequence for epitope prediction.
#' @param method A character string specifying the prediction method (default: "Bepipred-2.0").
#' @return A character string with prediction results or an error message.
antibody_epitope_prediction <- function(
    sequence_text,  
    method = "Bepipred-2.0" 
) {
  
  # API URL for Antibody Epitope Prediction
  api_url <- "http://tools-cluster-interface.iedb.org/tools_api/bcell/"
  
  # Construct the API request body
  body <- list(
    method = method,
    sequence_text = sequence_text
  )
  
  # Perform the POST request to the API
  response <- POST(
    url = api_url,
    body = body,
    encode = "form"  # Form-encoded data
  )
  
  # Check if the request succeeded
  if (http_status(response)$category == "Success") {
    # Return the response content as text
    content <- content(response, as = "text")
    return(content)
  } else {
    # Return error if failed
    return(paste("Error:", http_status(response)$message))
  }
}


#' @title Process a Single Antibody Sequence
#'
#' @description
#' Processes a single sequence for antibody epitope prediction and returns the result in a data frame.
#'
#' @param accession A character string for the protein accession number.
#' @param sequence A character string of the protein sequence.
#' @param method The prediction method to use.
#' @return A list with the processed result data frame or NULL if no result was returned.
process_antibody_sequence <- function(accession, sequence, method = "Emini") {
  
  # Call the antibody_epitope_prediction function for the single sequence
  result <- antibody_epitope_prediction(
    sequence_text = sequence,  # Provide the single sequence
    method = method            # Prediction method
  )
  
  # Check if the result is empty or NULL
  if (is.null(result) || result == "") {
    # Return an empty data frame with a warning or log
    warning(paste("No result returned for accession:", accession))
    return(list(empty_df = c(accession)))
  }
  
  # Split the result into lines
  result_lines <- strsplit(result, "\n")[[1]]
  
  # Check if there are enough lines to form a valid table
  if (length(result_lines) < 2) {
    # If the result doesn't contain enough rows, return an empty data frame
    return(list(empty_df = c(accession)))
  }
  
  # Process the first line as column names
  column_names <- strsplit(result_lines[1], "\t")[[1]]
  
  # Process the remaining lines as data rows
  data_rows <- lapply(result_lines[-1], function(x) strsplit(x, "\t")[[1]])
  
  # Create a data frame with the result
  result_df <- as.data.frame(do.call(rbind, data_rows), stringsAsFactors = FALSE)
  
  # Assign the column names
  colnames(result_df) <- column_names
  
  # Add the accession number as a new column to each row
  result_df$Accession <- accession
  
  # Return the processed data frame
  return(list(single_accession_immuno_df = result_df, empty_df = NULL))
}



#' @title Process Multiple Antibody Sequences
#'
#' @description
#' Processes multiple sequences for antibody prediction and organizes results within a specified timeout.
#'
#' @param df A data frame with columns `accession` and `sequence`.
#' @param method The prediction method to use (default: "Bepipred-2.0").
#' @param past_result_df Optional data frame to combine with current results.
#' @param timeout Numeric value for timeout duration in seconds.
#' @return A list with processed results and accessions with empty results.
process_sequences_antibodies_df <- function(df, method = "Bepipred-2.0", past_result_df = NULL, timeout = 7200) {
  start_time <- Sys.time()
  
  # Initialize an empty list to store results and empty data frame records
  results_list <- list()
  empty_df <- c()
  
  # Loop over the rows of the dataframe
  for (i in 1:nrow(df)) {
    
    # Extract accession and sequence for the current row
    accession <- df[i, "accession"]
    sequence <- df[i, "sequence"]
    
    # Call the single sequence processing function
    result <- process_antibody_sequence(
      accession = accession,
      sequence = sequence,
      method = method
    )
    
    if (!is.null(result$empty_df)) {
      empty_df <- c(empty_df, result[["empty_df"]])
    } else {
      results_list[[i]] <- result[["single_accession_immuno_df"]]
    }
    
    # Add a 60-second pause every 50 iterations
    if (i %% 50 == 0) {
      Sys.sleep(60)
    }
    
    current_time <- Sys.time()
    elapsed_time <- as.numeric(difftime(current_time, start_time, units = "secs"))
    
    if (elapsed_time > timeout) {
      # If the timeout is exceeded, return the current results so far
      results_df <- do.call(rbind, results_list)
      
      if (!is.null(past_result_df)) {
        combined_results <- rbind(past_result_df, results_df)
      } else {
        combined_results <- results_df
      }
      
      return(list(current_result_list = combined_results, all_empty_df = empty_df))
    }
  }
  
  # Save results to file and return combined results
  saveRDS(results_list, file = "notcombined_results.rds")
  message("Results saved to 'notcombined_results.rds'")
  
  # Combine all the results into a single data frame
  combined_results <- do.call(rbind, results_list)
  
  # Return the combined results and any empty data frames
  return(list(antibody_results = combined_results, accession_research = empty_df))
}

