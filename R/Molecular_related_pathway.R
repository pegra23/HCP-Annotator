#' @title Molecular related pathway
#'
#' @description
#' This script processes information from Uniprot, Brenda, Expansy and KEGG pathway
#' by relying in API calls (Uniprot), install packages (Brenda, KEGG) and manually curated 
#' information obtained by downloading the overall enzyme class from Expasy
#' 
#'
#' @details
#' The overall aim of this script is to retrieve information from multiple Accessions by relying in different curated DBs.
#' For you to use the "full power" of the script it will be necessary to have use the 
#' enzyme class info from Expansy which can be found here: Z:/Students at Alphalyse/Pedro Granjo, Master student 2024/Enzyme classes.xlsx
#'
#' @author Pedro Granjo
#' @date 11-FEB-2025
#'
#'
NULL  

############################# Install Packages if necessary################################################

libraries <- c(
  "httr", "jsonlite", "tidyr","KEGGREST", "readr", "stringr", "dplyr","brendaDb"
)

# Packages to install (those not already installed)
to_install <- setdiff(libraries, rownames(installed.packages()))

# Install missing packages
if(any(c("KEGGREST","brendaDb") %in% to_install)){
  
  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  if("KEGGREST" %in% to_install){
    
    BiocManager::install("KEGGREST")
    
  }
  if("brendaDb" %in% to_install) {
    
    BiocManager::install("brendaDb")
    
  }
  to_install <- setdiff(libraries, rownames(installed.packages()))
}



if (length(to_install) > 0) {
  install.packages(to_install, repos = "http://cran.us.r-project.org")
}

# Load required libraries
lapply(libraries, library, character.only = TRUE)



#source("problematic_hcp_functions.R")


######################################## Uniprot ################################################


#' @title Process Protein Accession Information from UniProt
#' 
#' @description
#' Retrieves information such as EC numbers, sequences, and organisms for proteins from UniProt 
#' based on accession numbers. Organizes this information into two tables for further analysis.
#' 
#' @param protein_info A data frame containing `Protein_name` and `Accession` columns.
#' @return A list containing two data frames: `main_table` (general information) and 
#' `evidence_table` (evidence-related details).
accession_processor <- function(protein_info){
  acession_list <- as.vector(protein_info$Accession)
  base_url <- "https://rest.uniprot.org/uniprotkb/"
  
  ec_info <- data.frame(Accession = character(),protein_name = character(), organism = character(), ec_numbers = character(), linear_sequence = character(), linear_sequence_length = character(), stringsAsFactors = FALSE)
  evidence_info <- data.frame(Accession = character(), enzymatic_evidenceCode = character(), source = character(), id = character(), stringsAsFactors = FALSE)
  
  protein_n <- 1 
  for(acession in acession_list){
    tryCatch({# API Request
      res <- httr::GET(paste0(base_url, acession, ".json"))
      data <- jsonlite::fromJSON(rawToChar(res$content))
      
      # Extract organism, EC numbers and some info of the sequence
      organism <- data$organism$commonName
      ec_numbers <- get_ec_numbers_fast(acession)
      entry_type_uniprot <- data[["entryType"]]
      linear_sequence <- data$sequence$value
      linear_sequence_length <- data$sequence$length
      protein_names <- paste(data[["proteinDescription"]][["recommendedName"]][["fullName"]][["value"]], collapse = "; ")
      protein_existence <- data[["proteinExistence"]]
      genes <- if (!is.null(data$genes$geneName$value)) {
        paste(data$genes$geneName$value, collapse = ",")
      } else {
        ""
      }
      protein_names <- if (!is.null(protein_names) && protein_names != "") {
        protein_names
      } else if (!is.null(data[["proteinDescription"]][["alternativeNames"]][["fullName"]][["value"]]) && 
                 length(data[["proteinDescription"]][["alternativeNames"]][["fullName"]][["value"]]) > 0) {
        
        paste(data[["proteinDescription"]][["alternativeNames"]][["fullName"]][["value"]], collapse = "; Alternative names =")
        
      } else if (!is.null(data[["proteinDescription"]][["submissionNames"]][["fullName"]][["value"]]) && 
                 length(data[["proteinDescription"]][["submissionNames"]][["fullName"]][["value"]]) > 0) {
        
        paste(data[["proteinDescription"]][["submissionNames"]][["fullName"]][["value"]], collapse = "; ")
        
      } else {
        "Not Available"
      }
      
      # Add the row to ec_info
      ec_info <- rbind(ec_info, data.frame(
        Accession = acession,
        protein_name = protein_names,
        organism = ifelse(!is.null(organism),organism, ifelse(!is.null(data$organism$scientificName),data$organism$scientificName,"Not Available")),
        entry_type_uniprot = ifelse(is.null(entry_type_uniprot), "Not Available", entry_type_uniprot), 
        inference_level = ifelse(is.null(protein_existence),"Not Available",  protein_existence),
        gene = genes,
        ec_numbers = ifelse(is.null(ec_numbers),"Not Available", paste(ec_numbers, collapse = ", ")),
        linear_sequence = ifelse(is.null(linear_sequence), "Not Available", paste(linear_sequence, collapse =", ")), 
        linear_sequence_length = ifelse(is.null(linear_sequence_length), "Not Available", paste(linear_sequence_length, collapse =", ")), 
        stringsAsFactors = FALSE
      ))
      
      protein_n <- 1 + protein_n
      # Extract evidence and store in the evidence table
      evidence <- data$proteinDescription$recommendedName$ecNumbers$evidences
      if (!is.null(evidence)) {
        count <- 1 
        for (ev in evidence) {
          id <- paste(acession, count, sep ="_")
          enzymatic_evidenceCode <- if ("enzymatic_evidenceCode" %in% names(ev)) ev$enzymatic_evidenceCode else "Not Available"
          source <- if ("source" %in% names(ev)) ev$source else "Not Available"
          info_ids <- if ("id" %in% names(ev)) ev$id else "Not Available"
          
          # Append to the evidence_info table
          evidence_info <- rbind(evidence_info, data.frame(
            ids = id,
            Accession = acession,
            enzymatic_evidenceCode = enzymatic_evidenceCode,
            enzymatic_source = source,
            info_id = info_ids,
            stringsAsFactors = FALSE
          ))
          count <- count + 1
        }
      } else {
        # If no evidence, we could also add "Not Available" if desired
        evidence_info <- rbind(evidence_info, data.frame(
          ids = paste0(acession,"_1"),
          Accession = acession,
          enzymatic_evidenceCode = "Not Available",
          enzymatic_source = "Not Available",
          info_id = "Not Available",
          stringsAsFactors = FALSE
        ))
      }
      
      if (protein_n %% 1500 == 0) {
        Sys.sleep(30)  # Pause execution to avoid overloading reenzymatic_sources
      }
    },error = function(e) {
      message(e)
      return(list(main_table = clean_and_extract_host_cell(ec_info), evidence_table = evidence_info))
    }
    )
  }
  return(list(main_table = clean_and_extract_host_cell(ec_info), evidence_table = evidence_info))
}


#' @title Retrieve EC Numbers from UniProt for a Given Accession
#'
#' @description This function fetches Enzyme Commission (EC) numbers for a specified
#'   UniProt accession ID by querying the UniProt REST API. It checks multiple
#'   locations within the UniProt entry in a defined hierarchical order to
#'   find the most relevant EC number(s).
#'
#' @param accession A character string representing a UniProt accession ID (e.g., "P00533").
#'
#' @return A character string containing the EC number(s) found, separated by ", ".
#'   If no EC numbers are found after checking all specified locations, it returns
#'   "Not Available". If the API request fails or returns invalid data, it returns
#'   \code{NA_character_}.
#'
#' @details The function prioritizes EC numbers from the following sources in order:
#'   \enumerate{
#'     \item \code{proteinDescription$recommendedName$ecNumbers$value}
#'     \item \code{uniProtKBCrossReferences} (specifically from the "BRENDA" database)
#'     \item \code{comments} (specifically from "reaction" type comments)
#'     \item \code{proteinDescription$submissionNames$ecNumbers$value}
#'   }
#'   It includes error handling for API requests and JSON parsing.
#'
#' @importFrom httr GET http_error
#' @importFrom jsonlite fromJSON
#' @importFrom utils head
#' @importFrom stats na.omit
#'
#' @examples
#' \dontrun{
#' # Requires internet connection
#' get_ec_numbers_fast("P00533")
#' get_ec_numbers_fast("Q9Y223")
#' get_ec_numbers_fast("P12345") # An accession that might not have easily found ECs
#' }
#' @export
get_ec_numbers_fast <- function(accession) {
  base_url <- "https://rest.uniprot.org/uniprotkb/"
  
  url <- paste0(base_url, accession, ".json")
  res <- tryCatch(httr::GET(url, timeout = 10), error = function(e) NULL)
  
  if (is.null(res) || httr::http_error(res)) {
    return(list(accession = accession, ec_numbers = NA_character_))
  }
  
  data <- tryCatch(jsonlite::fromJSON(rawToChar(res$content)), error = function(e) NULL)
  
  if (is.null(data)) {
    return(list(accession = accession, ec_numbers = NA_character_))
  }
  
  # Initialize ec_numbers to "Not Available"
  # This makes the subsequent 'if' conditions cleaner
  current_ec_numbers <- "Not Available"
  
  # 1. Try to get from proteinDescription$recommendedName$ecNumbers$value first (your initial preference)
  if (!is.null(data$proteinDescription) &&
      !is.null(data$proteinDescription$recommendedName) &&
      !is.null(data$proteinDescription$recommendedName$ecNumbers) &&
      length(data$proteinDescription$recommendedName$ecNumbers$value) > 0) {
    current_ec_numbers <- paste(data$proteinDescription$recommendedName$ecNumbers$value, collapse = ", ")
  }
  
  # 2. If still "Not Available", check BRENDA cross-references
  if (current_ec_numbers == "Not Available" &&
      !is.null(data[["uniProtKBCrossReferences"]]) &&
      nrow(data[["uniProtKBCrossReferences"]]) > 0) { # Check for rows before accessing 'database'
    brenda_refs <- data[["uniProtKBCrossReferences"]][data[["uniProtKBCrossReferences"]]$database == "BRENDA", "id"]
    if (length(brenda_refs) > 0) {
      current_ec_numbers <- paste(brenda_refs, collapse = ", ")
    }
  }
  
  # 3. If still "Not Available", check 'reaction' comments
  if (current_ec_numbers == "Not Available" && !is.null(data[["comments"]])) {
    # Filter for 'reaction' type comments
    reaction_comments_list <- data[["comments"]][sapply(data[["comments"]], function(x) !is.null(x$type) && x$type == "reaction")]
    
    if (length(reaction_comments_list) > 0 &&
        !is.null(reaction_comments_list[[1]][["reaction"]]) &&
        !is.null(reaction_comments_list[[1]][["reaction"]][["ecNumber"]])) {
      # Use na.omit just in case, though usually not needed if ecNumber is direct
      current_ec_numbers <- paste(na.omit(reaction_comments_list[[1]][["reaction"]][["ecNumber"]]), collapse = ", ")
    }
  }
  
  # 4. If still "Not Available", check proteinDescription$submissionNames$ecNumbers$value
  if (current_ec_numbers == "Not Available" &&
      !is.null(data[["proteinDescription"]]) &&
      !is.null(data[["proteinDescription"]][["submissionNames"]]) &&
      length(data[["proteinDescription"]][["submissionNames"]]) > 0 &&
      !is.null(data[["proteinDescription"]][["submissionNames"]][[1]][["ecNumbers"]]) &&
      length(data[["proteinDescription"]][["submissionNames"]][[1]][["ecNumbers"]]) > 0 &&
      !is.null(data[["proteinDescription"]][["submissionNames"]][[1]][["ecNumbers"]][[1]][["value"]])) {
    current_ec_numbers <- data[["proteinDescription"]][["submissionNames"]][[1]][["ecNumbers"]][[1]][["value"]]
  }
  
  # Final check and return
  return(current_ec_numbers)
}



#' @title Clean Organism Names and Extract Host Cell Information
#'
#' @description This function processes a dataframe containing organism names,
#'   standardizing common variations of organism names and extracting
#'   host cell information when present within parentheses in the organism string.
#'
#' @param df A data frame that must contain an 'organism' column.
#'
#' @return A data frame with an additional 'host_cell' column. The 'organism'
#'   column will have standardized names, and 'host_cell' will contain
#'   the extracted host cell information or NA if not found.
#'
#' @details The function performs two main operations:
#'   \enumerate{
#'     \item Extracts content within parentheses from the 'organism' column
#'           and assigns it to a new 'host_cell' column. If no parentheses
#'           are found, 'host_cell' will be NA.
#'     \item Standardizes specific common organism names (e.g., "Baker's yeast"
#'           to "Saccharomyces cerevisiae", "Escherichia coli (strain B / BL21-DE3)"
#'           to "Escherichia coli"). This part uses a pre-defined mapping.
#'   }
#'   Note: The organism cleaning logic in the example provided is based on a
#'   hardcoded list of replacements. For larger-scale cleaning, consider
#'   a more dynamic mapping approach (e.g., a lookup table).
#'
#' @importFrom dplyr mutate
#' @importFrom stringr str_extract
#' @importFrom stats na.omit # Used implicitly by some string ops if regex produces NA
#'
#' @examples
#' \dontrun{
#' # Example data frame
#' test_df <- data.frame(
#'   organism = c(
#'     "Homo sapiens (Human)",
#'     "Escherichia coli (strain K12)",
#'     "Baker's yeast",
#'     "Pseudomonas aeruginosa (strain ATCC 15692)",
#'     "Mus musculus",
#'     "Spreading-leaved earth moss",
#'     "Unknown organism"
#'   ),
#'   value = 1:7,
#'   stringsAsFactors = FALSE
#' )
#'
#' cleaned_df <- clean_and_extract_host_cell(test_df)
#' print(cleaned_df)
#' }
#' @export
clean_and_extract_host_cell <- function(df) {
  
  # Ensure the input is a data frame and has an 'organism' column
  if (!is.data.frame(df) || !("organism" %in% names(df))) {
    stop("Input must be a data frame with an 'organism' column.")
  }
  
  # 1. Extract host cell information
  # Use stringr::str_extract for a more robust way to get content in parentheses
  df <- df %>%
    mutate(
      host_cell = dplyr::case_when(
        stringr::str_detect(organism, "\\([^)]+\\)") ~
          stringr::str_extract(organism, "(?<=\\()[^)]+(?=\\))"),
        TRUE ~ NA_character_
      )
    )
  
  # 2. Clean organism names
  # Use a named vector for mapping for better readability and extensibility
  organism_mapping <- c(
    "Spreading-leaved earth moss" = "Physcomitrium patens",
    "Green monkey" = "Chlorocebus sabaeus",
    "Baker's yeast" = "Saccharomyces cerevisiae",
    "Pseudomonas aeruginosa (strain ATCC 15692 / DSM 22644 / CIP 104116 / JCM 14847 / LMG 12228 / 1C / PRS 101 / PAO1)" = "Pseudomonas aeruginosa",
    "Fall armyworm" = "Spodoptera frugiperda",
    "Fruit fly" = "Drosophila melanogaster",
    "Chicken" = "Gallus gallus",
    "Japanese quail" = "Coturnix japonica",
    "Common tobacco" = "Nicotiana tabacum",
    "Chinese hamster" = "CHO",
    "Golden hamster" = "Mesocricetus auratus",
    "Long-tailed dwarf hamster" = "Cricetulus longicaudatus",
    "Bacteriophage T7" = "Escherichia phage T7",
    "Mouse" = "Mus musculus"
  )
  
  # Apply general replacements using dplyr::recode for a cleaner approach
  # This handles the specific strings that need exact replacement.
  df <- df %>%
    mutate(organism = dplyr::recode(organism, !!!organism_mapping))
  
  # Handle Escherichia coli variants specifically as it's a common case
  # This pattern matches any string starting with "Escherichia coli"
  # followed by anything in parentheses, or just "Escherichia coli".
  df$organism[grepl("^Escherichia coli(\\s*\\([^)]*\\))?$", df$organism)] <- "Escherichia coli"
  
  
  return(df)
}


#' @title Merge and Clean Protein Evidence and Main Tables
#' 
#' @description
#' This function merges `main_table` and `evidence_table` from a list, consolidating enzymatic_sources 
#' and IDs for each unique accession number and removing unnecessary columns for final output.
#' 
#' @param info_list A list containing `main_table` and `evidence_table` data frames.
#' @return A data frame with cleaned and merged data for each unique accession, combining 
#'         information from `main_table` and consolidated `evidence_table`.
#' 
#' @examples
#' # Assuming `info_list` is a list with `main_table` and `evidence_table`
#' combined_data <- merging_tables(info_list)
#' print(combined_data)
#' 
#' @export
merging_tables <- function(info_list) {
  
  # Initialize an empty evidence table
  info_list[["evidence_table"]] <- info_list[["evidence_table"]][,-1] #Removing ids column which was the unique identiifier from this table
  updated_evidence_table <- info_list[["evidence_table"]][0,]
  
  # Iterate over unique accession numbers
  for (accession in unique(info_list[["evidence_table"]]$Accession)) {
    
    # Filter for the current accession number
    truncated_info <- info_list[["evidence_table"]][info_list[["evidence_table"]]$Accession == accession, , drop = FALSE]
    
    # Combine enzymatic_source and info_id
    new_enzymatic_sources <- paste(truncated_info$enzymatic_source, truncated_info$info_id, sep = "-")
    nyenzymatic_source <- paste(new_enzymatic_sources, collapse = ",")
    
    # Update the enzymatic_source column
    truncated_info$enzymatic_source <- nyenzymatic_source
    
    # Add the first row of the truncated info to the updated table
    updated_evidence_table <- rbind(updated_evidence_table, truncated_info[1, ])
  }
  
  # Replace specific value in the enzymatic_source column
  updated_evidence_table$enzymatic_source[updated_evidence_table$enzymatic_source == "Not Available-Not Available"] <- "Not Available"
  
  # Remove potential duplicate rows from Uniprot
  updated_evidence_table <- updated_evidence_table[!duplicated(updated_evidence_table), ]
  
  # Replace the old evidence table with the updated one
  info_list[["evidence_table"]] <- updated_evidence_table
  
  # Ensure main_table matches the updated evidence table
  info_list[["main_table"]] <- info_list[["main_table"]][match(info_list[["main_table"]]$Accession, info_list[["evidence_table"]]$Accession),]
  
  # Merge main_table and evidence_table for final use
  #The -4 is to remove the info_ids column which was merge in the enzymatic_source to reduce the overall number of row and allow the possibility of merging both tables
  combo <- merge(info_list[["main_table"]], info_list[["evidence_table"]][,-4], by = "Accession", all = T)
  
  return(combo = combo)
}


################################### Enzymatic Information ################################################

#' @title Compile Enzymatic Information for EC Numbers
#' 
#' @description
#' Processes each EC number, retrieving and formatting reaction information for a detailed analysis from Brenda Database.
#' The dependency of this function is the need to download locally the Brenda Database which will be the enzymatic_source of all info.
#' 
#' @param combo A data frame with EC numbers in the `ec_numbers` column.
#' @param df A data frame with enzymatic details for each EC number.
#' @return A data frame with processed `reaction_type` and `reaction` for each EC number.
#' 
#' @examples
#' enzymatic_info <- process_enzymatic_info(combo, enzymatic_list_df)
#' print(enzymatic_info)
#' 
#' @export
process_ec_numbers <- function(combo, df) {
  
  # Step 1: Clean and extract unique EC numbers from combo
  ec_numbers <- combo$ec_numbers[!(combo$ec_numbers == "Not Available")]
  ec_numbers <- as.vector(trimws(unlist(strsplit(ec_numbers, ", "))))
  ec_numbers <- unique(ec_numbers)
  
  # Step 2: Query Brenda database for information related to these EC numbers
  res <- brendaDb::QueryBrenda(df, EC = unique(ec_numbers), fields = c("REACTION", "REACTION_TYPE", "LOCALIZATION", "SPECIFIC ACTIVITY", "REFERENCE"))
  
  # Step 3: Process the result from Brenda to retrieve reaction_type and reaction
  results <- lapply(res, function(x) {
    enzymatic_reaction_type <- NA
    enzymatic_reaction <- NA
    
    # Process reaction_type
    if (!is.null(x[["nomenclature"]][["reaction.type"]]) && "description" %in% colnames(x[["nomenclature"]][["reaction.type"]])) {
      descriptions <- na.omit(x[["nomenclature"]][["reaction.type"]]$description)
      if (length(descriptions) > 0) {
        enzymatic_reaction_type <- paste("reaction type", paste(unique(na.omit(descriptions)), collapse = ", reaction type - "), sep = " - ")
      }
    }
    
    # Process reaction
    if (!is.null(x[["nomenclature"]][["reaction"]]$description)) {
      enzymatic_reaction <- paste("reaction", x[["nomenclature"]][["reaction"]]$description, sep = " - ", collapse = "; ")
    }
    
    # Return a data frame with the extracted information
    data.frame(enzymatic_reaction_type = enzymatic_reaction_type, enzymatic_reaction = enzymatic_reaction, stringsAsFactors = FALSE)
  })
  
  # Step 4: Combine all processed results into a single data frame, adding the respective EC numbers
  enzymatic_comission_combined_info <- do.call(rbind, lapply(names(results), function(ec_number) {
    df <- cbind(ec_number, results[[ec_number]])
  }))
  
  return(enzymatic_comission_combined_info)
}


#' @title Compile Enzymatic Information for EC Numbers
#' 
#' @description
#' Processes each EC number, retrieving and formatting reaction information for a detailed analysis.
#' This function was developed to curate potential missing information of different EC yourself. For example 
#' A few EC numbers are preliminary, they include an 'n' as part of the fourth (serial) digit (e.g. EC 3.5.1.n3). Thus their info is not in Brenda Database
#' The same can be said for enzyms with EC 1.1.-.-  but they will have their enzyme classification, which can be considered enough
#' 
#' @param combo A data frame with EC numbers in the `ec_numbers` column.
#' @param enzymatic_list_info A data frame with enzymatic details for each EC number.
#' @return A data frame with processed `reaction_type` and `reaction` for each EC number combined with the info from `combo` df
#' 
#' @export
process_enzymatic_info <- function(combo, enzymatic_list_info) {
  # Add a unique identifier column for the combo (to retain traceability)
  
  # Identify unique EC number combinations
  unique_combos <- unique(combo[, c("ec_numbers"), drop = F])
  unique_combos$ec_numbers_split <- lapply(
    strsplit(as.character(unique_combos$ec_numbers), ","),
    function(x) trimws(unique(x))
  )
  
  # Process unique combinations
  processed_combos <- lapply(seq_len(nrow(unique_combos)), function(i) {
    
    ec_numbers <- unique_combos$ec_numbers_split[[i]]

    
    if (length(ec_numbers) == 1 && ec_numbers[1] == "Not Available") {
      # Handle "Not Available" case
      return(data.frame(
        ec_numbers = unique_combos$ec_numbers[i],
        enzymatic_reaction_type = NA,
        enzymatic_reaction = NA,
        stringsAsFactors = FALSE
      ))
    }
    
    classes <- lapply(unique_combos$ec_numbers_split[[i]], function(ec_splits){
      
    
      # Filter enzymatic list info based on EC numbers
      row_info <- enzymatic_list_info[enzymatic_list_info$ec_number %in% ec_splits, ]
    
      if (nrow(row_info) == 0) {
        # No matches found
        return(data.frame(
          ec_numbers = unique_combos$ec_numbers[i],
          enzymatic_reaction_type = NA,
          enzymatic_reaction = NA,
          stringsAsFactors = FALSE
        ))
      }
      
    # Aggregate reaction_type and reaction into a single entry
    aggregated_info <- data.frame(
      ec_numbers = ec_splits,
      enzymatic_reaction_type = paste(unique(na.omit(row_info$enzymatic_reaction_type)), collapse = ", "),
      enzymatic_reaction = paste(unique(na.omit(row_info$enzymatic_reaction)), collapse = ", "),
      stringsAsFactors = FALSE
    )
    
    return(aggregated_info)
  })
    
    res <- as.data.frame(do.call(rbind,classes))
    
    collapsed_df <- res %>% 
      dplyr::summarise(across(everything(), ~ paste(unique(na.omit(gsub("Not Available",NA, .))), collapse = ", ")))
    
  })
    
  # Combine processed results into a single data frame
  processed_combos_df <- do.call(rbind, processed_combos)
  
  # Merge processed results back to the original combo based on ec_numbers
  final_result <- merge(combo, processed_combos_df, by = "ec_numbers", all.x = TRUE)
  
  return(final_result)
}


#' @title Compile Enzymatic Information for EC Numbers
#' 
#' @description
#' Processes each EC number, retrieving and formatting reaction information for a detailed analysis.
#' 
#' @param combo A data frame with EC numbers in the `ec_numbers` column.
#' @param enzymatic_class_df A data frame with enzymatic details for each EC number (their class and subclass)
#' @return A data frame with processed `enzyme_class` and `enzyme_subclass` for each EC number.
#' 
#' @export
process_enzymatic_class <- function(combo, enzymatic_class_df) {
  # Add a unique identifier column for the combo (to retain traceability)
  
  
  # Identify unique EC number combinations
  unique_combos <- unique(combo[, c("ec_numbers"), drop = F])
  unique_combos$ec_numbers_split <- lapply(
    strsplit(as.character(unique_combos$ec_numbers), ","),
    function(x) trimws(unique(x))
  )
  
  # Process unique combinations
  processed_combos <- lapply(seq_len(nrow(unique_combos)), function(i) {
    ec_numbers <- unique_combos$ec_numbers_split[[i]]
    
    if (length(ec_numbers) == 1 && ec_numbers[1] == "Not Available") {
      # Handle "Not Available" case
      return(data.frame(
        og_ec_number = unique_combos$ec_numbers[i],
        enzyme_class = "Not Available",
        enzyme_subclass = "Not Available",
        stringsAsFactors = FALSE
      ))
    }
    
    classes <- lapply(unique_combos$ec_numbers_split[[i]], function(ec_splits){
      
      if(grepl("-",ec_splits)){
        
        ec_class_info <- enzymatic_class_df[enzymatic_class_df$ec_number ==  trimws(ec_splits),]
        
      } else{
        
        unprocessed <- str_split(string = ec_splits, pattern = "\\.")
        ec_class <- paste0(paste(unlist(unprocessed)[1:(length(unlist(unprocessed))-1)], collapse = "."), ".-")
        ec_class_info <- enzymatic_class_df[enzymatic_class_df$ec_number == trimws(ec_class),]
      }
      
      return(ec_class_info)
    })
    
    if(length(ec_numbers) == 1){
      res <- as.data.frame(t(unlist(classes)))
      res_f <- data.frame(og_ec_number = unique_combos$ec_numbers[i],res, stringsAsFactors = FALSE)
      return(res_f)
    }
    else{
      res <- as.data.frame(do.call(rbind,classes))
      
      collapsed_df <- res %>% 
        dplyr::summarise(across(everything(), ~ paste(unique(.), collapse = ", ")))

      res_f <- data.frame(og_ec_number = unique_combos$ec_numbers[i],collapsed_df, stringsAsFactors = FALSE)
    }
    return(res_f)
  })
  
  # Combine processed results into a single data frame
  processed_combos_df <- do.call(bind_rows, processed_combos)
  
  # Merge processed results back to the original combo based on ec_numbers
  final_result <- merge(combo, processed_combos_df %>% select(-ec_number), by.x = "ec_numbers",by.y = "og_ec_number", all.x = TRUE)
  
  return(final_result)
}


################################### KEGG Pathway ################################################

#' @title Retrieve KEGG Pathways for UniProt Accessions
#' 
#' @description
#' Queries KEGG to retrieve pathway information linked to UniProt accessions.
#' 
#' @param uniprot_ids A vector of UniProt IDs.
#' @return A list of pathway information for each UniProt ID.
#' 
#' 
#' @export
retrieve_kegg_pathways_for_hcp <- function(uniprot_ids) {
  library(KEGGREST)
  
  kegg_pathway_info <- list()  # Initialize an empty list to store results
  counter <- 0  # Initialize a counter to track processed IDs
  
  for (uniprot_id in uniprot_ids) {
    counter <- counter + 1  # Increment counter for each processed ID
    
    # Step 1: Convert UniProt ID to kegg gene ID using keggREST
    tryCatch({
      kegg_conversion <- keggConv("genes", paste0("uniprot:", uniprot_id))
      
      if (length(kegg_conversion) > 0) {
        kegg_id <- kegg_conversion[1]  # Extract kegg gene ID
        
        # Step 2: Link kegg gene ID to pathways
        kegg_pathways <- keggLink("pathway", as.character(kegg_id))
        
        if (length(kegg_pathways) > 0) {
          kegg_pathway_info[[uniprot_id]] <- list()  # Create sublist for each UniProt ID
          
          for (pathway_id in names(kegg_pathways)) {
            Sys.sleep(1)
            # Custom request using `httr` for additional control
            response <- GET(paste0("https://rest.kegg.jp/get/", pathway_id),
                            add_headers(`User-Agent` = "R"))
            
            # Check if the request is successful
            if (status_code(response) == 200) {
              pathway_info <- content(response, "text")
              kegg_pathway_info[[uniprot_id]][[pathway_id]] <- pathway_info
            }
          }
        }
      }
    }, error = function(e) {
      message("Error with UniProt ID ", uniprot_id, ": ", e$message)
    })
    
    # Pause for 30 seconds every 1000 processed IDs
    if (counter %% 5000 == 0) {
      message("Processed ", counter, " IDs. Pausing for 30 seconds.")
      Sys.sleep(30)
    }
  }
  
  return(kegg_pathway_info)  # Return the nested list
}


#' @title Process KEGG Pathway Data for UniProt Accession
#' 
#' @description
#' Processes detailed pathway information from KEGG for a given UniProt accession.
#' 
#' @param kegg_data A list of pathway data retrieved for each UniProt accession.
#' @param uniprot_accession_id A single UniProt accession ID.
#' @return A data frame with KEGG ID, hyperlinks, and pathway details.
#' 
process_kegg_info <- function(kegg_data, uniprot_accession_id) {
  
  tryCatch({
    # Primary Extraction Logic
    kegg_id_hyperlink_info <- paste0("https://www.genome.jp/entry/", names(kegg_data[[uniprot_accession_id]]))
    kegg_id <- names(kegg_data[[uniprot_accession_id]])
    pathway_entries <- paste(names(kegg_data[[uniprot_accession_id]][[1]][[1]]$PATHWAY), as.character(kegg_data[[uniprot_accession_id]][[1]][[1]]$PATHWAY), sep = ":")
    pathway_entries_row <- paste(pathway_entries, collapse = ", ")
    pathway_hyperlink <- paste("https://www.kegg.jp/entry/pathway", names(kegg_data[[uniprot_accession_id]][[1]][[1]]$PATHWAY), sep = "+")
    pathway_hyperlink_rows <- paste(pathway_hyperlink, collapse = ", ")
    
    # Compile all info into a data frame
    data.frame(
      Accession = uniprot_accession_id,
      kegg_id = kegg_id,
      kegg_id_hyperlink = kegg_id_hyperlink_info,
      kegg_pathway = pathway_entries_row,
      kegg_pathway_hyperlink = pathway_hyperlink_rows,
      stringsAsFactors = FALSE
    )
    
  }, error = function(e) {
    # Alternative Extraction Logic (from previous improvement)
    tryCatch({
    # Extract kegg ID and hyperlink info
    # Split pathway information by newline
    try_out <- str_split(kegg_data[[uniprot_accession_id]][[1]][[1]], pattern = "\n")[[1]]
    entry_code <- str_split(trimws(gsub("^ENTRY", "", try_out[grepl("^ENTRY", unlist(try_out))])), pattern = " ")[[1]][[1]]
    
    kegg_id <- paste("hsa", entry_code, sep = ":")
    kegg_id_hyperlink_info <- paste0("https://www.genome.jp/entry/", kegg_id)
    
    pathway_lines_idx <- grep("^PATHWAY", try_out)
    if (length(pathway_lines_idx) > 0) {
      # Capture lines between PATHWAY and next non-indented line
      pathway_lines <- gsub("^PATHWAY", "       ", try_out[grep("^PATHWAY", try_out):length(try_out)])
      pathway_end <- which(substr(pathway_lines, 1, 5) != "     ")[1] - 1
      pathway_lines <- pathway_lines[1:pathway_end]
      # Process pathway lines
      cpathway_lines <- c()  # To store cleaned pathway lines
      pathway_codes <- sapply(pathway_lines, function(x) {
      empty_space_removed <- trimws(x)
      cpathway_lines <<- c(cpathway_lines, empty_space_removed)
      str_split(empty_space_removed, pattern = " ")[[1]][1]
    })
      # Final collapsed outputs
      fpathway_lines <- paste(cpathway_lines, collapse = ", ")
      fpathway_hyperlinks <- paste(paste0("https://www.kegg.jp/entry/pathway/", pathway_codes), collapse = ", ")
    } else {
      fpathway_lines <- "Unclassified"
      fpathway_hyperlinks <- "Unclassified"
    }
    
    # Compile alternative info into a data frame
    data.frame(
      Accession = uniprot_accession_id,
      KEGG_id = kegg_id,
      KEGG_id_hyperlink = kegg_id_hyperlink_info,
      KEGG_pathway = fpathway_lines,
      KEGG_pathway_hyperlink = fpathway_hyperlinks,
      stringsAsFactors = FALSE
    )}, error = function(inner_e) {
      data.frame(
        Accession = uniprot_accession_id,
        KEGG_id = "out_of_bounds",
        KEGG_id_hyperlink = "out_of_bounds",
        KEGG_pathway = "out_of_bounds",
        KEGG_pathway_hyperlink = "out_of_bounds",
        stringsAsFactors = FALSE
      )
    })
  })
}


#' @title Merge KEGG Pathway Data with Main Protein Data
#' 
#' @description
#' Retrieves and processes KEGG pathway information for UniProt accessions and merges it 
#' with an existing data frame of protein data.
#' 
#' @param all_collected_info A data frame containing UniProt accession data.
#' @return A data frame with additional KEGG pathway information.
process_kegg_and_merge <- function(all_collected_info) {
  
  # Step 1: Retrieve kegg pathway info for the provided accession numbers
  kegg_pathway_info <- retrieve_kegg_pathways_for_hcp(all_collected_info$Accession)
  
  # Save the kegg_pathway_info object as an RDS file to avoid re-fetching
  saveRDS(kegg_pathway_info, file = "kegg_pathway_info.rds")
  message("kegg_pathway_info saved to 'kegg_pathway_info.rds' for future use.")
  
  # If already saved, you can load with:
  # kegg_pathway_info <- readRDS("kegg_pathway_info.rds")
  
  # Step 2: Process the kegg pathway information
  kegg_compiled_infos <- lapply(names(kegg_pathway_info), function(id) {
    process_kegg_info(kegg_pathway_info, id)
  })
  
  saveRDS(kegg_compiled_infos, file = "kegg_pathway_info_new.rds")
  
  # Step 3: Compile the processed kegg pathway info into a data frame
  kegg_compiled_infos_df <- do.call(rbind, kegg_compiled_infos)
  
  # Step 4: Merge kegg info with the provided all_collected_info data frame
  moreinofcollected <- merge(all_collected_info, kegg_compiled_infos_df, by = "Accession", all = TRUE)
  
  return(moreinofcollected)
}

