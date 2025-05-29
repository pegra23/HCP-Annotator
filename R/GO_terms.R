#' @title GO_terms
#'
#' @description
#' This script processes information from QuickGO website
#' 
#'
#' @details
#' You can retrieve info from GO terms based on query based approach, Uniprot Accession, and
#' essential info from GO terms GO terms group E.G(Biological Process, and Text Definition of the GO terms)
#' 
#'
#' @author Pedro Granjo
#' @date 11-FEB-2025
#'
#'
NULL  

############################# Install Packages if necessary################################################

libraries <- c("httr", "jsonlite", "dplyr", "pryr")

to_install <- setdiff(libraries, rownames(installed.packages()))

# Install missing packages
if (length(to_install) > 0) {
  install.packages(to_install, repos = "http://cran.us.r-project.org")
}

# Load required libraries
lapply(libraries, library, character.only = TRUE)



##################################### GO API Set Up #############################################

#' Retrives detailed info of GO terms per Accession 
#'
#'
#'@param request_url Request url used to retrieve info from QuickGO API
#'
#'@return A processed long formated df with GO terms per Acession
fetch_and_parse_data <- function(request_url) {
  r <- try(GET(request_url, accept("text/gpad")), silent = TRUE)
  
  # Check for HTTP errors
  if (inherits(r, "try-error") || status_code(r) != 200) {
    warning(paste("Failed to fetch data for URL:", request_url))
    return(NULL)
  }
  
  # Process the content
  content <- content(r, as = "text")
  content_lines <- strsplit(content, "\n")[[1]]
  
  # Check for enough lines to process
  if (length(content_lines) <= 9) { # Metadata occupies the first 9 lines
    return(NULL)
  }
  
  # Extract rows of data and convert to a data frame
  data_rows <- lapply(content_lines[-(1:9)], function(x) strsplit(x, "\t")[[1]])
  result_df <- as.data.frame(do.call(rbind, data_rows), stringsAsFactors = FALSE)
  
  # Filter and rename columns
  result_df <- result_df[, -c(1, 8, 9, 11), drop = FALSE]
  colnames(result_df) <- c("Accession", "Qualifier", "GO_term", "Reference", "Evidence", "From", "Assigned By", "Annotation Extension")
  
  return(result_df)
}

#' Retrives detailed info of GO terms based on a list of Accession of interest
#'
#'
#'@param term Uniprot Accession ID used to pull a request to QuickGO API about which GO terms this Uniprot Accession has
#'
#'@return A processed long formated df with GO terms per Acession
accession_go_annotation <- function(term) {
  base_url <- "https://www.ebi.ac.uk/QuickGO/services/annotation/downloadSearch"
  
  # Construct the request URL for Accession
  request_url <- paste0(
    base_url,
    "?geneProductType=protein&geneProductId=", term
  )
  
  # Fetch data
  result_df <- fetch_and_parse_data(request_url)
  
  # Return the result
  return(result_df)
}


#' Retrives detailed info of GO terms based on a list of Accession of interest
#'
#'Constructs the request to Quick GO API based on the parameters of interste
#'Warning: This API retrieves always as much info as it can with a maximum of 10 000 rows per result.
#'So if your results are more 10 000, we advice to do multiple calls where the results can be subsetted based on taxon ids, references, evidence code left over ids might be able to be retrieved
#'
#'Personal advice try to do it manually the number of results that can be retrieved manually iin the Quick GO website. Just be aware of the default settings used by the app
#'
#'@param taxa taxon_id of interest
#'@param refs references information (df formated)
#'@param evidence_code evidence code information (df formated)
#'
#'@return A processed long formated df with GO terms per Accession
construct_request <- function(taxa, refs = NULL, evidence_code = NULL) {
  
  request <- paste0(
    base_url,
    "?goId=", go_term_clean,
    "&geneProductType=protein",
    "&goUsage=exact",
    "&taxonId=", paste(taxa, collapse = "%2C")
  )
  
  
  if (!is.null(refs)) {
    request <- paste0(request, "&reference=", paste(refs, collapse = "%2C"))
  }
  
  
  if (!is.null(evidence_code)) {
    request <- paste0(request, "&evidenceCode=", paste(evidence_code, collapse = "%2C"), "&evidenceCodeUsage=exact")
  }
  
  return(request)
}


#' Retrives detailed info of GO terms based on a list of Accession of interest
#'
#'Warning: This API retrieves always as much info as it can with a maximum of 10 000 rows per result.
#'So if your results are more 10 000, the results are subsetted based on taxon ids, references, evidence code left over ids might be able to be retrieved
#'
#'Personal advice try to do it manually the number of results that can be retrieved manually iin the Quick GO website. Just be aware of the default settings used by the app
#'
#'
#'@param terms GO terms
#'@param references references information (df of references like GO_REF, PMID, for more info explore QuickGO UI)
#'@param max_objects threshold for the size of the df that u want to collected. 10000 is the current limit for the API
#'@param max_refs_per_request maximum of references per request (e.g GO_Ref, PMID, DOI)
#'@param evidence_codes evidence code information (df of evidence codes that are used to annotate GO terms)
#'@param max_ecos_per_request evidence code information (df formated)
#'
#'@return A processed long formated df with GO terms per Accession
go_term_annotation <- function(term, references = NULL, max_objects = 10000, max_refs_per_request = 10, evidence_codes = NULL, max_ecos_per_request = 25) {
  base_url <- "https://www.ebi.ac.uk/QuickGO/services/annotation/downloadSearch"
  all_results <- list()
  overflow_log <- list()
  
  go_term_clean <- gsub(":", "%3A", term)
  taxon_ids <- c("10029", "562", "9606", "4932", "3218", "7108", "1280", "287", 
                 "4097", "60711", "93934", "35621", "9031", "7227", "470", 
                 "469008", "316407", "511693", "83333")
  
  # Helper to construct request URL
  construct_request <- function(taxa, refs = NULL, ecos = NULL) {
    request <- paste0(
      base_url,
      "?goId=", go_term_clean,
      "&geneProductType=protein",
      "&goUsage=exact",
      "&taxonId=", paste(taxa, collapse = "%2C"),
      #"&taxonUsage=exact",
      if (!is.null(refs)) paste0("&reference=", paste(refs, collapse = "%2C")),
      if (!is.null(ecos)) paste0("&evidenceCode=", paste(ecos, collapse = "%2C"), "&evidenceCodeUsage=exact")
    )
    return(request)
  }
  
  # Recursive fetch function with fractionation
  #Created to fragment refs, eco evidences to pull smaller request at the time
  fetch_results <- function(taxa, refs = NULL, ecos = NULL) {
    local_overflow_log <- list()
    if (!is.null(refs)) {
      ref_chunks <- split(refs, ceiling(seq_along(refs) / max_refs_per_request))
      
      # Apply over reference chunks
      results <- lapply(ref_chunks, function(ref_subset) {
        request_url <- construct_request(taxa, ref_subset)
        result_df <- fetch_and_parse_data(request_url)
        
        # Check if result exceeds max_objects; if not, reset ecos to NULL
        if (is.data.frame(result_df) && nrow(result_df) < max_objects) {
          return(result_df)  # No need to split further, return result directly
        }
        
        # If result exceeds max_objects, split further by evidence codes
        if (!is.null(ecos)) {
          
          if ("ECO:0000256" %in% ecos) {
            ecos <- ecos[!(ecos == "ECO:0000256")]
            eco_chunks <- split(ecos, ceiling(seq_along(ecos) / max_ecos_per_request))
            eco_chunks[[length(eco_chunks) + 1]] <- "ECO:0000256"
          }
          
          # Apply over evidence code chunks
          eco_results <- lapply(eco_chunks, function(eco_subset) {
            eco_subset <- gsub(":", "%3A", eco_subset)
            request_url <- construct_request(taxa, ref_subset, eco_subset)
            result_df <- fetch_and_parse_data(request_url)
            
            
            if (is.data.frame(result_df) && nrow(result_df) >=  max_objects) {
              message(sprintf("Reference chunk: %s, Evidence Code chunk: %s exceeded limit!", 
                              paste(ref_subset, collapse = ", "), 
                              paste(eco_subset, collapse = ", ")))
              
              local_overflow_log[[length(local_overflow_log) + 1]] <<- data.frame(GO_term = term,
                                                                                  Reference_Chunk = gsub("%3A", ":", paste(ref_subset, collapse = ", ")),
                                                                                  Evidence_Code_Chunk = gsub("%3A", ":", paste(eco_subset, collapse = ", ")),
                                                                                  stringsAsFactors = FALSE
              )
              return(NULL)
              
            }
            return(result_df)
          })
          
          # Combine results from all evidence code chunks
          return(do.call(rbind, eco_results))
        }
        
        return(result_df)
      })
      
      overflow_log <<- c(overflow_log, do.call(rbind, local_overflow_log))
      
      # Combine results from all reference chunks
      return(do.call(rbind, results))
      
    } else if (!is.null(ecos)) {
      # If no references but evidence codes are provided, split only by evidence codes
      if ("ECO:0000256" %in% ecos) {
        ecos <- ecos[!(ecos == "ECO:0000256")]
        eco_chunks <- split(ecos, ceiling(seq_along(ecos) / max_ecos_per_request))
        eco_chunks[[length(eco_chunks) + 1]] <- "ECO:0000256"
      }
      
      results <- lapply(eco_chunks, function(eco_subset) {
        eco_subset <- gsub(":", "%3A", eco_subset)
        request_url <- construct_request(taxa, refs = NULL, eco_subset)
        result_df <- fetch_and_parse_data(request_url)
        return(result_df)
      })
      
      # Combine results from all evidence code chunks
      return(do.call(rbind, results))
      
    } else {
      # If neither references nor evidence codes are provided, fetch directly by taxa
      request_url <- construct_request(taxa)
      result_df <- fetch_and_parse_data(request_url)
      return(result_df)
    }
  }
  
  
  # Attempt initial fetch with all taxa
  result_df <- fetch_results(taxon_ids)
  
  if (is.data.frame(result_df) && nrow(result_df) >= max_objects) {
    
    taxon_results <- lapply(taxon_ids, function(t) {
      taxon_result <- fetch_results(c(t))
      
      if (is.data.frame(taxon_result) && nrow(taxon_result) >= max_objects && !is.null(references)) {
        ref_results <- fetch_results(c(t), references)
        
        if (is.data.frame(ref_results) && nrow(ref_results) >= max_objects && !is.null(evidence_codes)) {
          eco_results <- fetch_results(c(t), references, evidence_codes)
          
          if (is.null(eco_results) || nrow(eco_results) == 0) {
            return(NULL)
          }
          
          return(eco_results)
        }
        
        return(ref_results)
      }
      
      return(taxon_result)
    })
    
    # Combine non-empty results from all taxa
    result_df <- do.call(rbind, taxon_results[!sapply(taxon_results, is.null)])
  }
  
  # Check if result_df is empty after all processing
  if (is.null(result_df) || nrow(result_df) == 0) {
    warning("No results were obtained for the given GO term, references, and evidence codes.")
    return(NULL)
  }
  
  # Return the main result
  return(list(main_result = result_df, overflow_log = as.data.frame(do.call(rbind, overflow_log))))
}


#' Multiple GO terms Annotation 
#'
#'Warning: Max retrieval from the API is 10 000 rows per results. If subseted with taxon ids, references, evidence code left over ids might be able to be retrieved
#'
#'@param go_terms GO terms 
#'@param references references information (df of references like GO_REF, PMID, for more info explore QuickGO UI)
#'@param evidence_codes evidence code information (df of evidence codes that are used to annotate GO terms)
#'@param batch_size To check how many go terms have been processed already
#'@param sleep_time Request stop time after the divisor of the predefined batch size
#'
#'@return A processed long formated df with GO terms per Accession
multiple_go_terms_annotation <- function(go_terms, references = NULL, evidence_codes = NULL, batch_size = 200, sleep_time = 30) {
  main_results <- list()
  overflow_logs <- list()
  
  # Iterate through the GO terms and annotate
  for (i in seq_along(go_terms)) {
    term <- go_terms[i]
    
    # Call go_term_annotation for the current term
    result <- go_term_annotation(term, references = references, evidence_codes = evidence_codes)
    
    if (!is.null(result)) {
      # Store main result and overflow log
      main_results[[term]] <- result$main_result
      overflow_logs[[term]] <- result$overflow_log
    }
    
    if (i %% 100 == 0) {
      number <- 100
      message("%s HCPs were already processed",number)
      number <- number + 100
    }
    
    # Log progress at every batch size interval
    if (i %% batch_size == 0) {
      Sys.sleep(sleep_time) 
    }
  }
  
  # Return combined results as a list
  return(list(main_results = main_results, overflow_logs = overflow_logs))
}

#' GO terms Fundamental information: Type (MF, CC, BP; Definition)
#'
#'@param go_term Go terms information
#'
#'@return A processed long formated df with GO terms per Accession
go_term_search <- function(go_term) {
  
  # Format the GO term to ensure it's URL-safe
  go_term_clean <- gsub(":", "%3A", go_term)
  requestURL <- paste0("https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/", go_term_clean)
  
  # API Request
  r <- tryCatch({
    
    GET(requestURL, accept("application/json"))
  }, error = function(e) {
    
    return(NULL) # Return NULL on error
  })
  
  # Check if the request was successful
  if (!is.null(r) && status_code(r) == 200) {
    
    # Parse JSON response
    content <- content(r)
    json <- toJSON(content)
    results <- fromJSON(json)$results
    
    return(as.data.frame(results)) 
    
  } else {
    
    warning(paste("Failed to fetch data for", go_term))
    return(go_term)
  }
}


##################################### Query Search Set up ##########################################


#' Keyword Query
#'
#'@param keyword query over each keyword to retrieve information about the GO terms
#'
#'@return Processed df with information about the GO terms
quickgo_search_query <- function(keyword){
  
  keyword_mod <- gsub(" ","%20",keyword)
  
  requestURL <- paste0("https://www.ebi.ac.uk/QuickGO/services/ontology/go/search?query=", keyword_mod)
  
  r <- GET(requestURL, accept("application/json"))
  
  stop_for_status(r)
  json <- toJSON(content(r))
  
  results <- fromJSON(json)$results
  results$keyword <- keyword
  
  
  return(as.data.frame(results))
}

#' Combination of Queries
#'
#'@param terms Apply multiple queries at once
#'
#'@return Processed df with all info in regards to each query
combine_quickgo_queries <- function(terms) {
  
  # Apply the quickgo_search_query function to each term
  df_list <- lapply(terms, quickgo_search_query)
  
  # Filter out any empty results
  df_list <- Filter(function(x) length(x) > 0 , df_list)
  
  # Combine all data frames into one
  df <- do.call(dplyr::bind_rows, df_list)
  
  if ("id" %in% colnames(df)) {
    df$id <- as.character(df$id)
  } else {
    stop("The combined data frame is missing an 'id' column.")
  }

  # Aggregate keywords by 'id'
  df_keywords <- aggregate(
    keyword ~ id, 
    data = df, 
    FUN = function(x) paste(unique(x), collapse = ", ")
  )
  
  df <- df %>% select(-keyword)
  
  df <- merge(df, 
              df_keywords, by = "id")
  non_keyword_columns <- setdiff(colnames(df), "keyword")
  
  df <- df[!(df$id == "NULL"), ]
  
  # Remove duplicate rows
  df <- distinct(df)
  
  df$isObsolete <- unlist(df$isObsolete)
  df$name <- unlist(df$name)
  df$aspect <- unlist(df$aspect)
  colnames(df)[colnames(df) == "definition"] <- "def_text"
  df$def_text <- unlist(df$def_text)
  
  return(df)
}


#' GO_terms_annotator_assembly
#'
#'It gives a final df of GO_terms that can be defined by a specific labels (can be a combination of two or more if they are found in two set of lists)
#'as well as the queries that were used in the first place to retrieve the information 
#'
#'@param dataframes List of dfs retrived for each queried set of keywords
#'@param labels Label for each set of keyword. It can be related to a specific topic
#'
#'@return create a compacted format df with information about the keyword, per GO terms and labels (category)
GO_terms_annotator_assembly <- function(dataframes, labels) {
  # Ensure the inputs are valid
  if (length(dataframes) != length(labels)) {
    stop("The number of data frames must match the number of labels.")
  }
  
  # Add the `Concern` column to each data frame
  for (i in seq_along(dataframes)) {
    if (!is.null(dataframes[[i]]) && nrow(dataframes[[i]]) > 0) {
      dataframes[[i]]$Concern <- labels[i]
    }
  }
  
  # Combine all data frames into a single data frame
  combined_df <- do.call(rbind, dataframes)
  
  # Remove duplicate rows and aggregate the `Concern` and `keyword` columns
  aggregated_df <- combined_df %>%
    group_by(id, name, aspect, def_text) %>%
    summarise(
      keyword = paste(unique(keyword), collapse = ", "),
      Concern = paste(unique(Concern), collapse = ", "),
      .groups = "drop"
    )
  
  return(aggregated_df)
}


#' Retrieve descendants GO terms for a given set of GO terms
#'
#'@param go_terms set of go terms
#'
#'@return a list of children in one subset and descendants on the other
go_children_descendants_query <- function(go_term){
  go_term_clean <- gsub(":", "%3A", go_term)
  requestURL <- paste0("https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/", go_term_clean,"/descendants?relations=is_a%2Cpart_of%2Coccurs_in%2Cregulates")
  r <- GET(requestURL, accept("application/json"))
  stop_for_status(r)
  json <- toJSON(content(r))
  results <- fromJSON(json)$results
  artf <- list(children = results$children[[1]], descendants = results$descendants[[1]])
  return(artf)
}



##################################### GO terms Definition ##########################################


#'Cleans pre-settings of each column of the retrieve df
#'
#'@param go_terms set of go terms
#'
#'@return Df with unlisted values per column
combine_go_terms <- function(go_terms) {
  
  matrix <- lapply(go_terms, function(x){
    go_term_search(x)
  })
  
  matrix <- do.call(bind_rows, matrix)
  matrix_def_text <- do.call(rbind,matrix[["definition"]]$text)
  matrix <- cbind(matrix[,1:9, drop =F], text_def = matrix_def_text)
  matrix <- matrix[,!(colnames(matrix) %in% c("definition","synonyms","children", "comment")),drop = F]
  matrix$id <- unlist(matrix$id)
  matrix$isObsolete <- unlist(matrix$isObsolete)
  matrix$name <- unlist(matrix$name)
  matrix$aspect <- unlist(matrix$aspect)

  return(matrix)
}




##################################### Generate GO columns for database ########################################


#'Formated GO columns stats and compacted into a df
#'
#'@param df_GO_terms_annotator A processed long formated df with GO terms per Accession
#'@param GO_info_new Fundamental info of each GO term
#'
#'@return df formated info of each GO term present in each Accession (compacted) with number per category
generate_GO_columns <- function(df_GO_terms_annotator, GO_info_new) {
  # Process each unique Accession
  GO_columns_info <- apply(
    as.data.frame(unique(df_GO_terms_annotator$Accession)), 
    1, 
    function(accession) {
      # Filter rows related to the current accession
      accession_info <- df_GO_terms_annotator[df_GO_terms_annotator$Accession == accession, , drop = FALSE]
      
      # Extract GO terms and count
      GO_ids <- paste(accession_info$GO_term, collapse = "; ")
      GO_n_counts <- nrow(accession_info)
      
      # Initialize placeholders for GO annotations by aspect
      GO_BP <- c()
      GO_MF <- c()
      GO_CC <- c()
      
      # Process each GO term for the accession
      sapply(accession_info$GO_term, function(go_term) {
        GO_annotation <- GO_info_new[GO_info_new$id == go_term, ]
        if (nrow(GO_annotation) > 0) {  # Ensure the GO term exists in GO_info_new
          if (GO_annotation$aspect == "molecular_function") {
            GO_MF <<- c(GO_MF, paste0(GO_annotation$name, " [", go_term, "]"))
          } else if (GO_annotation$aspect == "biological_process") {
            GO_BP <<- c(GO_BP, paste0(GO_annotation$name, " [", go_term, "]"))
          } else {
            GO_CC <<- c(GO_CC, paste0(GO_annotation$name, " [", go_term, "]"))
          }
        }
      })
      
      # Combine results for the current accession
      return(c(
        Accession = accession,
        GO_ids = GO_ids,
        GO_n_counts = GO_n_counts,
        GO_BP = paste(GO_BP, collapse = "; "),
        GO_MF = paste(GO_MF, collapse = "; "),
        GO_CC = paste(GO_CC, collapse = "; ")
      ))
    }
  )
  
  # Convert the result to a data frame
  GO_columns_info_df <- as.data.frame(t(GO_columns_info), stringsAsFactors = FALSE)
  
  # Set appropriate column names
  colnames(GO_columns_info_df) <- c("Accession", "GO_ids", "GO_n_counts", "GO_BP", "GO_MF", "GO_CC")
  
  return(GO_columns_info_df)
}



####################################GO Concern terms ########################################


#'Addition of labels in a concern column
#'
#'@param df 
#'
#'@return df with a Concern Column (Merged pre-defined labels regarding which HCPs)
concern_column_trimming <- function(df){
  elements <- as.vector(sapply(df$Concern, function(elem){
    elements <- unique(trimws(unlist(strsplit(elem, ", "))))
    elements <- paste(elements, collapse = ", ")
    return(elements)
  }))

  df$Concern <- elements
  return(df)
}


#'Transformation of df into a compacted format 
#'
#'@param df_GO_terms_annotator A processed long formated df with GO terms per Accession
#'@param GO_targets List of GO terms to be found in the df_GO_terms annotator
#'
#'@return df with a Concern Column (Merged pre-defined labels regarding each HCPs)
filter_keyword_compacted_format <- function(df_GO_terms_annotator, GO_targets) {
  
  # Filter the data for relevant GO terms
  filtered_data <- df_GO_terms_annotator[df_GO_terms_annotator$GO_term %in% GO_targets$id, ]
  
  # Add Info column to GO_targets for compact annotation
  GO_targets$Info <- paste0("[", GO_targets$id, "-", GO_targets$name, "]")
  
  # Join filtered data with GO_targets for additional metadata
  filtered_data <- filtered_data %>%
    left_join(
      GO_targets %>% select(all_of(c("id", "Info", "Concern","keyword"))),
      by = c("GO_term" = "id")
    )
  
  # Annotate toxicity or concern-related information
  filtered_data$GO_Concern <- paste(filtered_data$Info, filtered_data$Reference, sep = "_")
  
  # Expand concerns into separate columns
  unique_concerns <- unique(unlist(strsplit(filtered_data$Concern, split = ", ")))
  for (concern in unique_concerns) {
    filtered_data[[concern]] <- grepl(concern, filtered_data$Concern)
  }
  
  # Aggregate concern-related information by Accession
  aggregated <- filtered_data %>%
    group_by(Accession) %>%
    summarise(
      GO_Concern_terms = paste(unique(GO_Concern), collapse = ", "),
      From = paste(unique(From), collapse = ", "),
      Concern = paste(unique(Concern), collapse = ", "),
      across(all_of(unique_concerns), ~ any(.), .names = "{.col}")
    )
  
  # Rename columns for better clarity
  colnames(aggregated)[colnames(aggregated) %in% c("GO_Concern_terms")] <- paste("GO_Concern", "terms", sep = "_")
  
  aggregated <- concern_column_trimming(aggregated)
  return(aggregated)
}




#'Transformation of df into a wider, and compacted format
#'
#'@param df_GO_terms_annotator A processed long formated df with GO terms per Accession
#'@param GO_targets List of GO terms to be found in the df_GO_terms annotator
#'
#'@return df with a Concern Column (Merged pre-defined labels regarding which HCPs). Go terms with specific molecular activity are distributed that way  
filter_keyword_wider_format <- function(df_GO_terms_annotator, GO_targets) {
  # Filter the data for relevant GO terms
  filtered_data <- df_GO_terms_annotator[df_GO_terms_annotator$GO_term %in% GO_targets$id, ]
  
  # Add Info column to GO_targets
  GO_targets$Info <- paste0("[", GO_targets$id, "-", GO_targets$name, "]")
  
  # Join GO targets metadata into the filtered data
  filtered_data <- filtered_data %>%
    left_join(GO_targets %>% select(all_of(c("id", "Info", "Concern","keyword"))), by = c("GO_term" = "id"))
  
  # Get unique concerns
  unique_concerns <- unique(unlist(strsplit(filtered_data$Concern, split = ", ")))
  
  # Add empty columns for each concern
  for (concern in unique_concerns) {
    filtered_data[[concern]] <- NA
  }
  
  # Populate concern columns with appropriate values
  filtered_data <- apply(filtered_data, 1, function(row) {
    row <- as.list(row)  # Convert to list for easier assignment
    concerns <- trimws(unlist(strsplit(row[["Concern"]], split = ",")))
    for (con in concerns) {
      if (!is.na(con)) {
        row[[con]] <- paste(row[["Info"]], row[["Reference"]], sep = "_")
      }
    }
    return(row)
  })
  
  # Convert back to a data frame
  filtered_data <- do.call(rbind.data.frame, filtered_data)
  
  # Aggregate results by Accession
  
  if (!"Accession" %in% colnames(filtered_data)) {
    stop("Column `Accession` is missing after processing. Please debug the intermediate steps.")
  }
  
  filtered_data <- filtered_data %>%
    group_by(across(all_of("Accession"))) %>%  # Use `all_of()` to safely handle column names
    summarise(
      across( 
        #.cols = everything(), 
        .fns = ~ paste(unique(.x[!is.na(.x)]), collapse = ", "),  # Summarize non-NA unique values
      ),
      .groups = "drop"
    ) %>%
    select(! all_of(c("Annotation.Extension","Assigned.By","Qualifier", "Reference", "Evidence", "From")))
  result <- filtered_data 
  result <- concern_column_trimming(result)
  
  return(result)
}



#' @title Annotate UniProt Accessions with Gene Ontology (GO) Terms
#'
#' @description This function iterates through a list of UniProt accession IDs,
#'   calls a separate function (`accession_go_annotation`) to retrieve GO terms
#'   for each accession, and then stores these annotations in a consolidated
#'   data frame. It ensures that only unique GO term annotations are retained.
#'
#' @param accession_list A vector of unique UniProt accession IDs (e.g., from
#'   a column like `main_table_new_hcp$Accession`).
#' @param accession_go_annotation_func A function that takes a single
#'   accession ID as input and returns GO term annotations for that ID.
#'   This function should return NULL if no annotations are found or if an
#'   error occurs during its internal processing.
#'
#' @return A data frame containing unique GO term annotations. The structure
#'   of this data frame depends on the output of `accession_go_annotation_func`.
#'   It typically includes columns like 'Accession', 'GO_ID', 'GO_Term', etc.
#'   Returns an empty data frame if no annotations are successfully retrieved.
#'
#' @details
#' The function loops through each accession, calls the provided
#' `accession_go_annotation_func` for GO term retrieval, and stores the results
#' in an internal list. After processing all accessions, it combines these
#' results into a single data frame and removes any duplicate rows to ensure
#' a unique set of GO term annotations.
#'
#' The `accession_go_annotation_func` is expected to handle its own API calls,
#' error handling, and formatting of the GO term output.
#'
#' @importFrom dplyr distinct
#' @importFrom methods is
#'
#' @examples
#' \dontrun{
#' # Assuming you have a dummy or actual accession_go_annotation function
#' # For demonstration, let's create a simple dummy function:
#' dummy_accession_go_annotation <- function(acc) {
#'   if (acc == "P00533") {
#'     return(data.frame(
#'       Accession = "P00533",
#'       GO_ID = c("GO:0005575", "GO:0006979"),
#'       GO_Term = c("cellular_component", "response to DNA damage stimulus"),
#'       stringsAsFactors = FALSE
#'     ))
#'   } else if (acc == "Q9Y223") {
#'     return(data.frame(
#'       Accession = "Q9Y223",
#'       GO_ID = "GO:0008150",
#'       GO_Term = "biological_process",
#'       stringsAsFactors = FALSE
#'     ))
#'   } else {
#'     return(NULL) # No annotation for other accessions
#'   }
#' }
#'
#' # Create a dummy list of accessions
#' dummy_accessions <- unique(data.frame(Accession = c("P00533", "Q9Y223", "A0A009", "P00533"), stringsAsFactors = FALSE)$Accession)
#'
#' # Run the function
#' go_annotations_df <- annotate_accessions_with_go(
#'   accession_list = dummy_accessions,
#'   accession_go_annotation_func = dummy_accession_go_annotation
#' )
#' print(go_annotations_df)
#' }
#' @export
annotate_accessions_with_go <- function(accession_list) {
  if(is.data.frame(accession_list)){
    accession_list <- accession_list$Accession
  }
  
  # Initialize a list to store GO term results
  GO_terms_results_list <- list()
  
  # Ensure accession_list is a character vector for iteration
  accession_list <- as.character(accession_list)
  
  # Loop through each unique accession to get GO annotations
  for (term_accession in accession_list) {
    # Call the provided GO annotation function
    result <- accession_go_annotation(term_accession)
    
    # If the result is not NULL, store it in the list
    if (!is.null(result) && methods::is(result, "data.frame")) {
      GO_terms_results_list[[term_accession]] <- result
    } else {
      # Optional: message if a result is NULL or not a data.frame
      # message(paste("No GO annotations or invalid format for accession:", term_accession))
    }
  }
  
  # Combine all results from the list into a single data frame
  # Using data.table::rbindlist or purrr::list_rbind could be more efficient for very large lists
  if (length(GO_terms_results_list) > 0) {
    GO_terms_annotator <- do.call(rbind, GO_terms_results_list)
    
    # Ensure uniqueness of the GO term annotations
    df_GO_terms_annotator <- dplyr::distinct(GO_terms_annotator)
  } else {
    # Return an empty data frame with appropriate columns if no data was collected
    # This requires knowing the expected column names from accession_go_annotation_func
    # For a generic solution, you might need to infer columns from the first non-NULL result
    # or pass expected columns as an argument.
    # For now, returning an empty data.frame. User should define initial empty DF if specific columns are needed.
    df_GO_terms_annotator <- data.frame() # Returns an empty dataframe
  }
  
  
  return(df_GO_terms_annotator)
}







