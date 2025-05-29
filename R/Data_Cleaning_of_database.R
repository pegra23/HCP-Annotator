table(main_table_2611$organism)



new_accessions_merge <- new_accessions_merge[!(new_accessions_merge$accession_number %in% main_table_2611$Accession),, drop = F]

new_accessions_merge$infoplace <- "GO terms"



names <- do.call(rbind, names)
matching_indices <- which(new_accessions_merge$accession_number %in% names$Accession)
new_accessions_merge[matching_indices, "protein_name"] <- as.vector(names$Protein_name)


host_cell <- lapply(main_table_new_hcp$organism, function(row){
  organism <- row
  
  if (grepl("\\(", organism) && grepl("\\)", organism)) {
    # Extract content between parentheses
    host_cell <-  gsub(".*\\((.*)\\).*", "\\1", organism) 
    #host_cell <- sub("\\).*", "", sub(".*\\(", "", string))
  } else {
    host_cell <- NA
  }
  return(host_cell)
  
})


host_cell <- do.call(rbind, host_cell)

main_table_new_hcp$host_cell <- as.vector(host_cell)
table(main_table_new_hcp$organism)

annotation_table[annotation_table$organism == "Spreading-leaved earth moss",]$organism <- "Physcomitrium patens"
annotation_table[annotation_table$organism == "Green monkey",]$organism <- "Chlorocebus sabaeus"
main_table_new_hcp[main_table_new_hcp$organism == "Baker's yeast",]$organism <- "Saccharomyces cerevisiae"
main_table_new_hcp[main_table_new_hcp$organism == "Pseudomonas aeruginosa (strain ATCC 15692 / DSM 22644 / CIP 104116 / JCM 14847 / LMG 12228 / 1C / PRS 101 / PAO1)",]$organism <- "Pseudomonas aeruginosa" 


annotation_table[annotation_table$organism == "Fall armyworm",]$organism <- "Spodoptera frugiperda"
annotation_table[annotation_table$organism == "Fruit fly",]$organism <- "Drosophila melanogaster"
annotation_table[annotation_table$organism == "Chicken",]$organism <- "Gallus gallus"
annotation_table[annotation_table$organism == "Japanese quail",]$organism <-  "Coturnix japonica"
annotation_table[annotation_table$organism == "Common tobacco",]$organism <- "Nicotiana tabacum"
annotation_table[annotation_table$organism == "Chinese hamster",]$organism <- "CHO"
annotation_table[annotation_table$organism == "Golden hamster",]$organism <- "Mesocricetus auratus"
annotation_table[annotation_table$organism == "Long-tailed dwarf hamster",]$organism <- "Cricetulus longicaudatus"
annotation_table[annotation_table$organism == "Bacteriophage T7",]$organism <- "Escherichia phage T7"
annotation_table[annotation_table$organism == "Mouse",]$organism <- "Mus musculus"
  
HCPAnnotatorDB[["annotation_table"]] <- annotation_table

main_table_new_hcp[main_table_new_hcp$organism %in% c("Escherichia coli", "Escherichia coli (strain B / BL21-DE3)" ,"Escherichia coli (strain B / BL21)",
                                                              "Escherichia coli (strain K12 / W3110 / ATCC 27325 / DSM 5911)",
                                                              "Escherichia coli (strain K12)"),]$organism <- "Escherichia coli"


annot <- HCPAnnotatorDB2404final_2[["annotation_table"]]


miss_col <- setdiff(colnames(annot),colnames(main_table_new_hcp))

main_table_new_hcp <- main_table_new_hcp %>% select(-combo_id)
       

for(m_col in miss_col){
  main_table_new_hcp[m_col] <- NA
}

library(httr)
library(jsonlite)


base_url <- "https://rest.uniprot.org/uniprotkb/"

accession_list <- new_accessions_merge$accession_number[new_accessions_merge$protein_name == ""]
protein_n <- 1

organism <- lapply(updated_main_table_2801[updated_main_table_2801$infoplace != "GO terms","Accession"], function(x){
  res <- GET(paste0(base_url, x, ".json"))
  data <- fromJSON(rawToChar(res$content))
  organism <- data$organism$commonName
  organisms <- ifelse(!is.null(organism),organism, ifelse(!is.null(data$organism$scientificName),data$organism$scientificName,"Not Available"))
  
  retrieved_info <- data.frame(Accession =x , Organism =organisms)
  
  protein_n <<- protein_n + 1 
  
  if (protein_n %% 5000 == 0) {
    Sys.sleep(30) 
  }
  
  return(retrieved_info)
})

final_df_name <- do.call(rbind, names)  

new_accessions_merge[which(new_accessions_merge$accession_number %in% final_df_name$Accession),] <- final_df_name$Protein_name


save(final_df_name, file = "New_Accessions_With_names_v2.RData")
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))





##################################nEED TO SEE WHAT IS THIS 


colnames(new_accessions_merge)[1] <- "Accession"


new_accessions_merge$Accession <- as.character(new_accessions_merge$Accession)
kegg_compilation_f$Accession <- as.character(kegg_compilation_f$Accession)


common_accessions <- intersect(new_accessions_merge$Accession, kegg_compilation_f$Accession)


test_merge <- merge(new_accessions_merge, kegg_compilation_f, by = "Accession", all = FALSE)


new_accessions_df <- new_accessions_df[rowSums(is.na(new_accessions_df[-1])) < ncol(new_accessions_df[-1]), ]


new_accessions_df <- merge(new_accessions_merge, kegg_compilation_f, by = "Accession", all.x = TRUE)


kegg_compilation_f <- kegg_compilation_f[kegg_compilation_f$Accession %in% new_accessions_merge$Accession,] 

new_accessions_df <- merge(new_accessions_merge, kegg_compilation_f, by = "Accession", all = TRUE)


############################################################# Merge Old and New info #######################

library(dplyr)

main_table_2801 <- bind_rows(main_table_2611, new_accessions_merge_kegg)

Accession_GO_annotation <- bind_rows(GO_fundamental_information[["annotator_old_term"]],GO_fundamental_information[["annotator_new_term"]])

GO_combo_wider_total <- filter_keyword_wider_format(Accession_GO_annotation,GO_fundamental_information[["elements_go_expansion"]][["GO_Concern"]])

New_table_info <- GO_combo_wider_total %>%
  select(Accession, Concern)

main_table_2801 <- merge(main_table_2801, New_table_info, by= "Accession", all = T)

ec_numbers_info <- process_enzymatic_info(new_accessions_merge, ec_numbers_list_processed)

main_table_2801 %>% 
  -select(all_of(c("enzymatic_reaction_type", "enzymatic_reaction")))

main_table_2801 <- main_table_2801 %>% 
  select(-all_of(c("enzymatic_reaction_type", "enzymatic_reaction")))



main_table_2801 <- process_enzymatic_info(main_table_2801, ec_numbers_list_processed)

paste(colnames(main_table_2801), collapse =", ")


save(main_table_2801, file= "Main_table2801.RData")


setwd(dirname(rstudioapi::getActiveDocumentContext()$path))







# Reorder the columns of main_table_2801 based on the specified order

updated_main_table_2801 <- main_table_2801[, 
                                           c(
                                             "Accession", 
                                             "protein_name", 
                                             "organism",
                                             "host_cell", 
                                             "ec_numbers",
                                             "enzymatic_reaction_type", 
                                             "enzymatic_reaction",
                                             "evidenceCode",
                                             "source",
                                             "Label",
                                             "KEGG_id",
                                             "KEGG_id_hyperlink",
                                             "KEGG_pathway",
                                             "KEGG_pathway_hyperlink",
                                             "GO_ids",
                                             "GO_n_counts",
                                             "GO_BP", 
                                             "GO_MF", 
                                             "GO_CC",
                                             "Concern",
                                             "linear_sequence",
                                             "linear_sequence_length",
                                             "infoplace",
                                             "reference_id",
                                             "Critical_for", 
                                             "old_Functional_Impact", 
                                             "Function", 
                                             "Function_short"
                                           ), drop = FALSE]



enzyme_numberss <- data.frame("ec_numbers" = c(
  "1.7.2.-", "4.2.1.-", "4.3.99.-", "7.5.2.-"), "enzymatic_reaction_type" =
    c("reaction type - Oxidoreductases. Acting on other nitrogenous compounds as donors. With a cytochrome as acceptor.",
  "reaction type - Carbon-oxygen lyases. Hydro-lyases.",
  "reaction type - Carbon-nitrogen lyases. Other carbon-nitrogen lyases.",
  "reaction type - Translocases. Catalysing the translocation carbohydrates and their derivatives. Linked to the hydrolysis of a nucleoside triphosphate."
)) 


updated_main_table_2801[which(updated_main_table_2801$ec_numbers %in% enzyme_numberss$ec_numbers),"enzymatic_reaction", drop = F] <- enzyme_numberss$enzymatic_reaction_type


save(updated_main_table_2801, main_table_2801, file  = "Main_table2801.RData")

final_df_name <- do.call(rbind,names)

#updated_main_table_2801[match(final_df_name$Accession, updated_main_table_2801$Accession), "protein_name"] <- final_df_name$Protein_name


###################

unique(new_accessions_merge$organism)

colnames(new_accessions_merge)[1] <- "Accession"

new_accessions_merge_fuller <- merge(new_accessions_merge, GO_fundamental_information[["new_go_columns_info"]], by = "Accession")
kegg_compilation_2801 <- kegg_compilation_f[! (kegg_compilation_f$Accession %in% main_table_2611$Accession),, drop = F]


new_accessions_merge_kegg <- merge(new_accessions_merge_fuller, kegg_compilation_2801, by = "Accession",all = T)

colnames(new_accessions_merge_kegg)[20] <- "Label"
new_accessions_merge_kegg$Label <- "Predicted"


library(stringr)
colnames(new_accessions_merge_kegg)[16:19] <- gsub("Kegg", "KEGG", colnames(new_accessions_merge_kegg)[16:19])

save(updated_main_table_2801, file = "main_table_2901.RData")









###################################### Accession Integration

info_list <- accession_processor(hcp_concern_accession)
combo_1 <- merging_tables(info_list) 

database_brenda <- brendaDb::ReadBrenda("brenda_2023_1.txt")
combined_df <- process_ec_numbers(new_accessions_merge,database_brenda)

writexl::write_xlsx(list(
  res = as.data.frame(new_accessions_merge$ec_numbers)
), path = "/work/Alphalyse/ec_numbers_listja2025.xlsx")


ec_numbers_list_processed <- readxl::read_excel("file_show (1)(1).xlsx")

integrated_ec_numbers <- process_enzymatic_info(new_accessions_merge, ec_numbers_list_processed)

main_table_new_hcp <- process_keeg_and_merge(final_combo)
new_accessions_merge_kegg$Label <- "Predicted"
new_accessions_merge_kegg$infoplace <- "Predicted"



Accession <- unique(main_table_2611$Accession)

GO_terms_annotator <- list()

try <- Filter(function(x) nrow(x) > 0 ,results$overflow_logs)


nos <- apply(Accession, 1, function(x) {
  term <- unlist(x) 
  result <- accession_go_annotation(term) 
  if (!is.null(result)) {
    GO_terms_annotator[[term]] <<- result 
  }
})

GO_terms_annotator <- do.call(rbind, GO_terms_annotator)
df_GO_terms_annotator <- distinct(do.call(rbind,GO_terms_annotator))


#matrix
accession_terms_wgo <- as.vector(unique(df_GO_terms_annotator$GO_term))
matrix_go_annotator <- combine_go_terms(accession_terms_wgo)


#New columns for main table of HCPs of Concern Table
go_columns_info <- generate_GO_columns(df_GO_terms_annotator, matrix_go_annotator)



updated_main_table_2801 <- main_table_2801[, 
                                           c(
                                             "Accession", 
                                             "protein_name", 
                                             "organism",
                                             "host_cell", 
                                             "ec_numbers",
                                             "enzymatic_reaction_type", 
                                             "enzymatic_reaction",
                                             "evidenceCode",
                                             "source",
                                             "Label",
                                             "KEGG_id",
                                             "KEGG_id_hyperlink",
                                             "KEGG_pathway",
                                             "KEGG_pathway_hyperlink",
                                             "GO_ids",
                                             "GO_n_counts",
                                             "GO_BP", 
                                             "GO_MF", 
                                             "GO_CC",
                                             "Concern",
                                             "linear_sequence",
                                             "linear_sequence_length",
                                             "infoplace",
                                             "reference_id",
                                             "Critical_for", 
                                             "old_Functional_Impact", 
                                             "Function", 
                                             "Function_short"
                                           ), drop = FALSE]

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("brendaDb")



if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("KEGGREST")


# Load required libraries
library(brendaDb)
library(writexl)
library(readxl)
library(dplyr)


getwd()
load("Molecular_Pathway.R")

annotation_table <- HCPAnnotatorDB[["annotation_table"]]



setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

new_additional_accession <- hcps_table$Accession[!(hcps_table$Accession %in% new_annotation_table$Accession)]
hcp_concern_accession <- data.frame(Accession =new_additional_accession)



annotation_table <- HCPAnnotatorDB[["annotation_table"]]

source("Molecular_related_pathway.R")
source("GO_terms.R")

hcp_concern_accession <- data.frame("Accession" = c("A0A140NAY9", "A0A140NHB7"))



    info_list <- accession_processor(hcp_concern_accession)
    
    message("### Step 2: Merging Tables - CHECK THIS STEP MANUALLY ###")
    combo <- merging_tables(info_list)
    
    message("Please verify the merging process and press Enter to continue...")
    readline()
    
    
    sum(all_accession_enz$ec_numbers != "Not Available")
    message("### Step 3: Read BRENDA Database ###")
    database_brenda <- brendaDb::ReadBrenda("C:/Users/granjo/OneDrive - Alphalyse A S/Documents/Master/Tasks/Task_2 - DB Enrichment/brenda_2023_1.txt")
    
    message("### Step 4: Process EC Numbers ###")
    combined_df <- process_ec_numbers(all_accession_enz, database_brenda)
    
    enzymatic_type_df <- process_enzymatic_info(all_accession_enz, combined_df)
    
    enzymatic_class_df <- 
    process_enzymatic_class(enzymatic_type_df, Enzyme_classes[c("ec_number","enzyme_class","enzyme_subclass")])
    
    
    
    enzyme_new_info <- enzymatic_class_df %>%
      group_by(Accession) %>%
      summarise(across(everything(), ~ paste(unique(na.omit(.x)), collapse = ", ")))
    
    all_prot_info_update <- merge(all_prot_info %>% select(-ec_numbers), enzyme_new_info, by = "Accession", all.x = T)
    
    
    
    save(all_prot_info_update, file = "All_new_columns_updated.RData")
    
    message("### Step 8: Process KEGG and Merge ###")
    main_table_new_hcp <- process_kegg_and_merge(enzymatic_class_df)
    
    message("### Step 9: Assign Labels and Infoplace ###")
    
    # Prompt user for input
    user_label <- readline("Enter value for Label (default: Predicted): ")
    if (user_label == "") user_label <- "Predicted"  # Default if no input
    
    user_infoplace <- readline("Enter value for Infoplace (default: Predicted): ")
    if (user_infoplace == "") user_infoplace <- "Predicted"  # Default if no input
    
    # Assign values
    integrated_ec_numbers$Label <- user_label
    main_table_new_hcp$Label <- ""
    main_table_new_hcp$infoplace <- ""
    
    message("Label set to: ", user_label)
    message("Infoplace set to: ", user_infoplace)
    
    message("### Step 10: GO Terms Annotation ###")
    Accession <- unique(main_table_new_hcp$Accession)
    GO_terms_annotator <- list()
  
    
    nos <- lapply(main_table_new_hcp$Accession, function(x) {
      term <- unlist(x) 
      result <- accession_go_annotation(term) 
      if (!is.null(result)) {
        GO_terms_annotator[[term]] <<- result 
      }
    })
    
    
    GO_terms_annotator <- do.call(rbind, GO_terms_annotator)
    df_GO_terms_annotator <- distinct(GO_terms_annotator)
    
    message("### Step 11: Generate GO Annotations ###")
    accession_terms_wgo <- as.vector(unique(df_GO_terms_annotator$GO_term))
    matrix_go_annotator <- combine_go_terms(accession_terms_wgo)
    
    message("### Step 12: Create GO Columns ###")
    go_columns_info <- generate_GO_columns(df_GO_terms_annotator, matrix_go_annotator)
    
    main_table_new_hcp_v2 <- main_table_new_hcp_v2[, !(colnames(main_table_new_hcp_v2) %in% "GO_ids"), drop = F]
    
    main_table_new_hcp_v2 <- merge(main_table_new_hcp,go_columns_info, all.x = T)
    message("### Step 13: Update Main Table ###")
    main_table_new_hcp_fv <- main_table_new_hcp[, 
                                               c(
                                                 "Accession",
                                                 "protein_name",
                                                 "organism",
                                                 "host_cell",
                                                 "inference_level",
                                                 "gene",                       
                                                 "ec_numbers",
                                                 "enzyme_class",
                                                 "enzyme_subclass",
                                                 "enzymatic_reaction_type",
                                                 "enzymatic_reaction",
                                                 "enzymatic_evidenceCode",
                                                 "enzymatic_source",
                                                 "Label",
                                                 "KEGG_id",
                                                 "KEGG_id_hyperlink",
                                                 "KEGG_pathway",
                                                 "KEGG_pathway_hyperlink",
                                                  "GO_ids",
                                                 "GO_n_counts",
                                                 "GO_BP",
                                                 "GO_MF",
                                                 "GO_CC",
                                                 "Potential_functional_impact",
                                                 "linear_sequence",
                                                 "linear_sequence_length",
                                                 "infoplace",
                                                 "reference_id",
                                                 "Critical_for",
                                                 "old_Functional_Impact",
                                                 "Function",
                                                 "Function_short"  
                                               ), drop = FALSE]
    

colnames(main_table_new_hcp)[10] <- "enzymatic_source"
    
setdiff(colnames(HCPAnnotatorDB2404final_2[["annotation_table"]]) ,colnames(main_table_new_hcp_fv))

main_table_new_hcp_vf <-  merge(main_table_new_hcp, go_columns_info, by = "Accession", all.x = T)


# Run the pipeline
run_pipeline(interval_minutes = 30)


organism <- do.call(rbind, organism)
colnames(organism)

host_cell <- lapply(seq_along(organism$Organism), function(row){
  organism_type <- organism[row,"Organism", drop =F]
  accession <- organism[row,"Accession", drop = F]
  
  if (grepl("\\(", organism_type) && grepl("\\)", organism_type)) {
    # Extract content between parentheses
    Accession <- accession
    host_cell <-  gsub(".*\\((.*)\\).*", "\\1", organism_type)
    main_info <- data.frame(Accession, host_cell)
    #host_cell <- sub("\\).*", "", sub(".*\\(", "", string))
  } else {
    return(NULL)
  }
  return(main_info)
  
})

res <- Filter(function(x){length(x) > 1},host_cell)
res <- do.call(rbind, res)

ecoli_strain1 <- data.frame(Accession = organism$Accession[organism$Organism %in% "Escherichia coli O157:H7"], host_cell = "O157:H7")
ecoli_strain2 <- data.frame(Accession = organism$Accession[organism$Organism %in% "Escherichia coli O1:K1 / APEC"], host_cell = "O1:K1 / APEC")

res$host_cell[res$host_cell == "strain E2348/69 / EPEC"] <- "O127:H6 (strain E2348/69 / EPEC)"
res$host_cell[res$host_cell == "strain CFT073 / ATCC 700928 / UPEC"] <- "O6:H1 (strain CFT073 / ATCC 700928 / UPEC)" 

res <- rbind(res,ecoli_strain1,ecoli_strain2)

updated_main_table_2801[match(res$Accession,updated_main_table_2801$Accession),"host_cell"] <-res$host_cell

main_table_3001 <- updated_main_table_2801

saveRDS(HCPAnnotatorDB3001, file = "HCPAnnotatorDB3001.rds")

View(main_table)

main_table <- 

new_concern_column <- c()
for(concern in main_table[main_table$infoplace == "IEDB","Concern",]){
  res <- unlist(str_split(concern,","))
  if(!(any(grepl("Immune",res)))){
    new_concern <- paste(unique(trimws(c(na.omit(res),"Immune"))), collapse = ", ")
    new_concern_column <- c(new_concern_column, new_concern)
  } else{
    new_concern_column <- c(new_concern_column, concern)
  }
}

HCPAnnotatorDB3001[["main_table"]] <- main_table

library(dplyr)

old_columns <- updated_main_table_2801 %>% dplyr::select(all_of(c("Accession", "Concern")))
old_columns <- old_columns[!duplicated(old_columns$Accession),]


prelunde <- distinct(HCPAnnotatorDB[["annotation_table"]] %>% select(-all_of(c("Concern"))))


annotation_table <- merge(prelunde,old_columns, by = "Accession")


annotation_table <- annotation_table[!duplicated(annotation_table$Accession),]




table(main_table$organism)


host_cell <- lapply(main_table$organism, function(row){
  organism <- row
  
  if (grepl("\\(", organism) && grepl("\\)", organism)) {
    # Extract content between parentheses
    host_cell <-  gsub(".*\\((.*)\\).*", "\\1", organism) 
    #host_cell <- sub("\\).*", "", sub(".*\\(", "", string))
  } else {
    host_cell <- NA
  }
  return(host_cell)
  
})

host_cell <- do.call(rbind, host_cell)

main_table[main_table$organism %in% c("Escherichia coli", "Escherichia coli (strain B / BL21-DE3)" ,"Escherichia coli (strain B / BL21)",
                                                          "Escherichia coli (strain K12 / W3110 / ATCC 27325 / DSM 5911)",
                                                          "Escherichia coli (strain K12)"),]$organism <- "Escherichia coli"
main_table$host_cell <- host_cell




base_url <- "https://rest.uniprot.org/uniprotkb/"

accession_list <- main_table$Accession[is.na(main_table$protein_name)]
protein_n <- 1

library(httr)
library(jsonlite)


protein_names <- lapply(accession_list, function(x){
  res <- GET(paste0(base_url, x, ".json"))
  data <- fromJSON(rawToChar(res$content))
  protein_names_og <- paste(data[["proteinDescription"]][["recommendedName"]][["fullName"]][["value"]], collapse = "; ")
  
  protein_names <- if (!is.null(protein_names_og) && protein_names_og != "") {
    protein_names_og
  } else if (!is.null(data[["proteinDescription"]][["alternativeNames"]][["fullName"]][["value"]]) && 
             length(data[["proteinDescription"]][["alternativeNames"]][["fullName"]][["value"]]) > 0) {
    
    paste(data[["proteinDescription"]][["alternativeNames"]][["fullName"]][["value"]], collapse = "; Alternative names =")
    
  } else if (!is.null(data[["proteinDescription"]][["submissionNames"]][["fullName"]][["value"]]) && 
             length(data[["proteinDescription"]][["submissionNames"]][["fullName"]][["value"]]) > 0) {
    
    paste(data[["proteinDescription"]][["submissionNames"]][["fullName"]][["value"]], collapse = "; ")
    
  } else {
    "Not Available"
  }
  
  retrieved_info <- data.frame(Accession =x , Organism =protein_names)
  
  protein_n <<- protein_n + 1 
  
  if (protein_n %% 5000 == 0) {
    Sys.sleep(30) 
  }
  
  return(retrieved_info)
})

final_df_name <- do.call(rbind, protein_names)  

main_table <- main_table[,!(colnames(main_table) %in% c("protein_names"))]

annotation_table_v3 <- annotation_table_v2 %>%
  mutate(Concern = ifelse(
    grepl("IEDB", infoplace) & !grepl("Immune", Concern), 
    paste0(Concern, ", Immune"), 
    Concern
  ))


table(annotation_table_v3$Concern)
table(annotation_table_v3$infoplace)




lines_without_suffix
library("stringr")


lines_without_suffix_v2 <- lapply(seq_len(nrow(new_accession_df)), function(i) {
  # Extract the row as a list
  row <- lines_without_suffix[i, , drop = F]
  
  id <- row$reference_id  # Access the reference_id column
  ids <- trimws(unlist(str_split(id, ",")))  # Split the IDs by ", "
  # Initialize a vector for valid IDs
  valid_ids <- c()
  
  # Process each ID
  sapply(ids, function(id) {
    tryCatch({
      # If conversion to integer succeeds, add it to valid_ids
      o <- as.integer(id)
      
      if(is.na(o)){
        clean_id <- unlist(str_split(id, "_"))[1]
        valid_ids <<- rbind(valid_ids, clean_id)
      }
    }, error = function(e) {
      print(id)
      # If there's an error, add the ID to valid_ids
      clean_id <- unlist(str_split(id, "_"))[1]
      valid_ids <<- rbind(valid_ids, clean_id)
    })
  })
  
  # Replace `reference_id` with collapsed version of valid IDs
  row$reference_id <- paste(valid_ids, collapse = ", ")
  
  return(row)
})


# Convert the list back to a data frame
lines_without_suffix <- do.call(rbind, lines_without_suffix_v2)


annotation_table_v4 <- merge(annotation_table_v3 %>% select(-reference_id), lines_without_suffix, by = "Accession", all = T)

annotation_table_v5 <- annotation_table_v4[annotation_table_v4$organism != "Not Available",]



annotation_table_final  <- annotation_table_final %>% select(-reference_id)

annotation_table_final_2  <-  merge(annotation_table_v5, lines_without_suffix, by = "Accession", all.x = T)

annotation_table_final  <- annotation_table_final_2[,colnames(HCPAnnotatorDB[["annotation_table"]]), drop = F]


colnames(annotation_table_final_2)[28] <- "reference_id"

HCPAnnotatorDB[["annotation_table"]] <- annotation_table_final

GO_Accession_Annotator <- HCPAnnotatorDB[["GO_Accession_Annotator"]]
GO_info <- HCPAnnotatorDB[["GO_info"]]

GO_Accession_Annotator_v2 <- distinct(rbind(GO_Accession_Annotator, GO_terms_annotator))

GO_Accession_Annotator_v2[duplicated(GO_Accession_Annotator_v2 %>% select(Accession, GO_term, Reference)),]

clean_matrix <-  matrix_go_annotator %>% select(-all_of(c("isObsolete","usage")))
GO_info_v2 <- distinct(bind_rows(GO_info,clean_matrix))

HCPAnnotatorDB[["GO_info"]] <- GO_info_v2
HCPAnnotatorDB[["GO_Accession_Annotator"]] <- GO_Accession_Annotator_v2

saveRDS(HCPAnnotatorDB, file = "HCPAnnotatorDB.rds")


new_annotation_table <- HCPAnnotatorDB[["annotation_table"]] %>% select(-infoplace)
new_annotation_table <- merge(annotation_table_final %>% select(-all_of(c("Concern"))), annotation_table_v2 %>% select(Accession, Concern), by = "Accession")

new_annotation_table <- new_annotation_table[, colnames(HCPAnnotatorDB[["annotation_table"]]),drop = F]


new_annotation_table <- new_annotation_table %>% 
  mutate(Concern = ifelse(
    grepl("IEDB", infoplace) & !grepl("\\bImmune\\b", Concern),  # Add "Immune" only if missing
    str_trim(paste(unique(c(unlist(str_split(Concern, ",\\s*")[[1]]),"Immune")), collapse = ", ")), 
    Concern
  ))


new_annotation_table <- new_annotation_table[new_annotation_table$organism != "Not Available",]

table(new_annotation_table$Concern)
    
HCPAnnotatorDB[["annotation_table"]] <- new_annotation_table

saveRDS(HCPAnnotatorDB, file = "HCPAnnotatorDB.rds")

View(new_annotation_table[grepl("Alphalyse Manual Curation",new_annotation_table$infoplace),])

new_annotation_table$infoplace[new_annotation_table$infoplace == "Alphalyse Manual Curation,"] <- "Alphalyse Manual Curation"

sum(table(new_annotation_table$organism))

sum(table(new_annotation_table$infoplace[new_annotation_table$Accession %in% lines_without_suffix$Accession]))



#' Check Maximum Character Length per Column in Each Table
#'
#' This function iterates through each table in the provided database list and determines the maximum character length of each column.
#'
#' @param db_list A named list of data frames where each data frame represents a table in the database.
#'
#' @return A data frame with table names, column names, and their maximum character lengths.
#' @export
check_max_char_length <- function(db_list) {
  max_length_info <- list()
  
  for (table_name in names(db_list)) {
    table_data <- db_list[[table_name]]
    
    if (is.data.frame(table_data)) {
      for (col_name in names(table_data)) {
        if (is.character(table_data[[col_name]])) {
          max_length <- max(nchar(table_data[[col_name]]), na.rm = TRUE)
          max_length_info <- append(max_length_info, list(data.frame(
            Table = table_name,
            Column = col_name,
            Max_Length = max_length
          )))
        }
      }
    }
  }
  
  result_df <- do.call(rbind, max_length_info)
  return(result_df)
}




charct <- check_max_char_length(HCPAnnotatorDB)

index <- 1
library(stringr)

r <- lapply(seq_along(new_annotation_table$GO_ids), function(index){
  number <- unique(trimws(unlist(str_split(new_annotation_table$GO_ids[[index]], pattern = ";"))))
  new_annotation_table$GO_ids[[index]] <- paste(number, collapse = "; ")
  
  new_annotation_table$GO_BP[[index]]<- paste(unique(trimws(unlist(str_split(new_annotation_table$GO_BP[[index]], pattern = ";")))), collapse = "; ")
  new_annotation_table$GO_MF[[index]] <- paste(unique(trimws(unlist(str_split(new_annotation_table$GO_MF[[index]], pattern = ";")))), collapse = "; ")
  new_annotation_table$GO_CC[[index]] <- paste(unique(trimws(unlist(str_split(new_annotation_table$GO_CC[[index]], pattern = ";")))), collapse = "; ")
  new_annotation_table$GO_n_counts[[index]] <- length(number)
})



main_table_new_hcp_v2 <- main_table_new_hcp



GO_Accession_Annotator_cc <-HCPAnnotatorDB[["GO_Accession_Annotator"]]
GO_info_cc <- HCPAnnotatorDB[["GO_info"]]

colnames(df_GO_terms_annotator) <- gsub(" ", "_",colnames(df_GO_terms_annotator))

GO_Accession_Annotator_cc <- distinct(rbind(GO_Accession_Annotator_cc, df_GO_terms_annotator))
GO_info_cc <- distinct(bind_rows(GO_info_cc, matrix_go_annotator[,colnames(matrix_go_annotator) %in% colnames(GO_info_cc),drop=F]))

HCPAnnotatorDB[["GO_Accession_Annotator"]] <- GO_Accession_Annotator_cc
HCPAnnotatorDB[["GO_info"]] <- GO_info_cc


saveRDS(HCPAnnotatorDB, file = "HCPAnnotatorDB2502.rds") 

library(writexl)
write_xlsx(HCPAnnotatorDB[["GO_info"]], path = "C:/Users/granjo/OneDrive - Alphalyse A S/Documents/Master/Tasks/Scripts_RObjects_January/AllGO.xlsx")


length(unique(c(test_1$Accession)))

length(na.omit(HCPAnnotatorDB$annotation_table$reference_id[HCPAnnotatorDB$annotation_table$Accession %in% unique(c(test_1$Accession, processed_2802$processed$search_query)) ]))


test_1$Reference <- gsub("PMID:", "", test_1$Reference)


columns <- c(
  "Accession", 
  "protein_name", 
  "organism",
  "host_cell",
  "inference_level",
  "gene",
  "ec_numbers",
  "enzyme_class",
  "enzyme_subclass",
  "enzymatic_reaction_type", 
  "enzymatic_reaction",
  "enzymatic_evidenceCode",
  "enzymatic_source",
  "Label",
  "KEGG_id",
  "KEGG_id_hyperlink",
  "KEGG_pathway",
  "KEGG_pathway_hyperlink",
  "Potential_functional_impact",
  "GO_ids",
  "GO_n_counts",
  "GO_BP", 
  "GO_MF", 
  "GO_CC",
  "linear_sequence",
  "linear_sequence_length",
  "immuno_pos_assay_freq",      
  "immuno_pos_neg_ratio", 
  "infoplace",
  "reference_id",
  "Critical_for", 
  "old_Functional_Impact", 
  "Function", 
  "Function_short")

miss_colum <- setdiff(colnames(HCPAnnotatorDB2404final_2[["annotation_table"]]), colnames(main_table_new_hcp))

for(i in miss_colum){
  main_table_new_hcp[i] <- NA
}
