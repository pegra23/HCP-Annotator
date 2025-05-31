#' @title DB_query.R
#'
#' @description
#' This script processes retrieves information from European PMC. Currently is incomplete but it's still capable of retrieve information if the article id is provided.
#' To improve the API calls multi core calls would be suitable
#' 
#'
#' @details
#' The overall aim of this script is to retrieve information from different articles that might of interest retrieve metadata such as
#' author, PubYear, CitationNumber, OpenSource, Title, JournalName, SubmissiondDate, among others
#'
#' @author Pedro Granjo
#' @date 16-APR-2025
#'
#'
NULL  

############################# Install Packages if necessary################################################


library <- c("dplyr","httr", "jsonlite")


libraries <- c(
  "dplyr","httr", "jsonlite"
)

# Packages to install (those not already installed)
to_install <- setdiff(libraries, rownames(installed.packages()))

if (length(to_install) > 0) {
  install.packages(to_install, repos = "http://cran.us.r-project.org")
}

# Load required libraries
lapply(libraries, library, character.only = TRUE)


#######################################################################################################



fetch_europe_pmc_data <- function(Organism = NULL, Accession_ID = NULL, Article_ID = NULL, Cursormarker = "*",
                                  pageSize = 1000, save_path = "europe_pmc_results.csv") {
  # Define API endpoint
  base_url <- "https://www.ebi.ac.uk/europepmc/webservices/rest/searchPOST"
  
  # Construct query string dynamically
  query_parts <- c()
  
  #if (!is.null(Accession_ID)) query_parts <- c(query_parts, paste0('ACCESSION_ID:"', Accession_ID, '"'))
  if (!is.null(Article_ID)) query_parts <- c(query_parts, paste0('EXT_ID:"', Article_ID, '"'))
  # If no parameters provided, stop execution
  if (length(query_parts) == 0) {
    stop("At least one parameter (GOterm, Organism, Accession_ID, Gene_Product) must be provided.")
  }
  
  # Construct the full search query
  search_query <- paste(query_parts)
  
  # Define parameters for the POST request
  query_params <- list(
    query = search_query,
    resultType = "lite",
    synonym = "N",
    format = "json",
    pageSize = pageSize,
    cursorMark = Cursormarker
  )
  
  # Send the request
  response <- POST(
    url = base_url,
    body = query_params,
    encode = "form",
    content_type("application/x-www-form-urlencoded")
  )
  
  # Extract and parse the JSON response
  response_content <- content(response, "text", encoding = "UTF-8")
  parsed_data <- fromJSON(response_content)
  
  # Extract the results
  df <- parsed_data[["resultList"]][["result"]]
  
  # If no results, return NULL
  if (is.null(df) || length(df) == 0) {
    return(NULL)
  }
  
  # Add the keyword combination to the data frame for tracking
  df <- df[,!(colnames(df) %in% c("fullTextIdList","dbCrossReferenceList","tmAccessionTypeList")), drop = F]
  #vector <- unlist(rep(Accession_ID, nrow(df)))
  #df$search_query <- vector
  
  #write.csv(combined_data, save_path, row.names = FALSE)
  return(df)
}



process_row <- function(row) {
  # Extract query parameters dynamically
  #GOterm <- if (!is.null(row["name"])) row["name"] else NULL
  Accession_IDs <-  if (!is.na(row["Accession"]) && row["Accession"] != "") row["Accession"] else NULL
  
  result <- tryCatch({
    fetch_europe_pmc_data(Accession_ID = Accession_IDs)
  }, error = function(e) {
    message(sprintf("Error %s",e$message))
    return(NULL)
  })
  
  # Return result if successful
  if (!is.null(result)){
    return(list(success = TRUE, data = result))
  } 
  
  # Log failed queries if all attempts fail
  return(list(
    success = FALSE,
    failed_queries = Accession_IDs))
}



