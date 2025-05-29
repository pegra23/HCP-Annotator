# Set working directory to the directory of the current script
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Source necessary R scripts
source("Molecular_related_pathway.R")
source("GO_terms.R")

#Enzyme Class Information
enzyme_classes <- read_rds("enzyme_classes.rds")

#Define accessions which the information should be upploaded
hcp_concern_accession <- data.frame("Accession" = c("A0A140NAY9", "A0A140NHB7"))


#Retrieve Information from Uniprot
info_list <- accession_processor(hcp_concern_accession)
combo <- merging_tables(info_list)

#Retrieve information from Brenda and Expansy(Enzyme classes object)
database_brenda <- brendaDb::ReadBrenda("C:/Users/granjo/OneDrive - Alphalyse A S/Documents/Master/Tasks/Task_2 - DB Enrichment/brenda_2023_1.txt")
combined_df <- process_ec_numbers(combo, database_brenda)
enzymatic_type_df <- process_enzymatic_info(combo, combined_df)
enzymatic_class_df <- 
  process_enzymatic_class(enzymatic_type_df, enzyme_classes)

#Process KEGG information
main_table_new_hcp <- process_kegg_and_merge(enzymatic_class_df)

#Retrieve GO term information from QuickGOAPI
df_GO_terms_annotator <- annotate_accessions_with_go(main_table_new_hcp)

#Make table of GO_info type, etc
accession_terms_wgo <- as.vector(unique(df_GO_terms_annotator$GO_term))
matrix_go_annotator <- combine_go_terms(accession_terms_wgo)

#Generate GO columns for the final database
go_columns_info <- generate_GO_columns(df_GO_terms_annotator, matrix_go_annotator)

#Merge information
main_table_new_hcp <- merge(main_table_new_hcp,go_columns_info, all.x = T)


#Final columns. There are some that were brought from the old Alphalyse internal database annotations~
#They are not maintaned but we keep the information that it was stored

columns <- c(
  "Accession",
  "protein_name",
  "organism",
  "entry_type_uniprot",
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
)

#Add missing columns and order them 
miss_col <- setdiff(columns,colnames(main_table_new_hcp))

if (length(miss_col) > 0) {
  main_table_new_hcp[miss_col] <- NA
}

main_table_new_hcp <- main_table_new_hcp[,columns,drop = F]


#################### Data Retrieval from IEDB ###############


load("IEDB_tools_db_query.RData")

tables_metadata <- list(
  antigen_search = list(endpoint = "antigen_search",
                        valid_parameters = c(
                          "parent_source_antigen_id", "parent_source_antigen_iri", "parent_source_antigen_names", 
                          "parent_source_antigen_source_org_iri", "parent_source_antigen_source_org_name", "e_modifications", 
                          "linear_sequence_lengths", "qualitative_measures", "mhc_allele_evidences", 
                          "antibody_isotypes", "direct_ex_vivo_bool", "iedb_assay_ids", 
                          "iedb_assay_iris", "reference_ids", "reference_iris", 
                          "submission_ids", "pdb_ids", "chebi_ids", 
                          "structure_ids", "structure_iris", "submission_iris", 
                          "structure_types", "structure_descriptions", "curated_source_antigens", 
                          "linear_sequences", "receptor_ids", "receptor_group_ids", 
                          "tcr_receptor_group_ids", "bcr_receptor_group_ids", "receptor_group_iris", 
                          "tcr_receptor_group_iris", "bcr_receptor_group_iris", "receptor_types", 
                          "receptor_names", "receptor_chain1_types", "receptor_chain2_types", 
                          "receptor_chain1_full_seqs", "receptor_chain2_full_seqs", "receptor_chain1_cdr1_seqs", 
                          "receptor_chain2_cdr1_seqs", "receptor_chain1_cdr2_seqs", "receptor_chain2_cdr2_seqs", 
                          "receptor_chain1_cdr3_seqs", "receptor_chain2_cdr3_seqs", "host_organism_iri_search", 
                          "host_organism_iris", "host_organism_names", "source_organism_iri_search", 
                          "source_organism_iris", "source_organism_names", "mhc_allele_iri_search", 
                          "mhc_allele_iris", "mhc_allele_names", "parent_source_antigen_iri_search", 
                          "disease_iri_search", "disease_iris", "disease_names", 
                          "assay_iri_search", "assay_iris", "assay_names", 
                          "epitope_structures_defined", "non_peptidic_molecule_iri_search", "non_peptidic_molecule_iris", 
                          "non_peptidic_molecule_names", "r_object_source_molecule_iri_search", "r_object_source_molecule_iris", 
                          "r_object_source_molecule_names", "r_object_source_organism_iri_search", "r_object_source_organism_iris", 
                          "r_object_source_organism_names", "mhc_classes", "mhc_allele_resolutions", 
                          "e_related_object_types", "tcell_ids", "tcell_iris", 
                          "bcell_ids", "bcell_iris", "elution_ids", 
                          "elution_iris", "journal_names", "reference_types", 
                          "pubmed_ids", "reference_titles", "reference_authors", 
                          "reference_dates"),
                        default_values = list(offset = 0, order = "parent_source_antigen_id"),
                        primary_key = "parent_source_antigen_id"),
  bcr_search = list(endpoint = "bcr_search",
                    valid_parameters = c(
                      "receptor_group_id", "tcr_receptor_group_id", "bcr_receptor_group_id", 
                      "receptor_group_iri", "tcr_receptor_group_iri", "bcr_receptor_group_iri", 
                      "receptor_ids", "receptor_type", "receptor_species_names", 
                      "chain1_cdr3_seq", "chain2_cdr3_seq", "e_modifications", 
                      "linear_sequence_lengths", "reference_ids", "reference_iris", 
                      "submission_ids", "submission_iris", "pdb_ids", 
                      "chebi_ids", "qualitative_measures", "mhc_allele_evidences", 
                      "antibody_isotypes", "direct_ex_vivo_bool", "iedb_assay_ids", 
                      "iedb_assay_iris", "structure_ids", "structure_iris", 
                      "structure_types", "structure_descriptions", "curated_source_antigens", 
                      "linear_sequences", "receptor_names", "receptor_chain1_types", 
                      "receptor_chain2_types", "receptor_chain1_full_seqs", "receptor_chain2_full_seqs", 
                      "receptor_chain1_cdr1_seqs", "receptor_chain2_cdr1_seqs", "receptor_chain1_cdr2_seqs", 
                      "receptor_chain2_cdr2_seqs", "receptor_chain1_cdr3_seqs", "receptor_chain2_cdr3_seqs", 
                      "host_organism_iri_search", "host_organism_iris", "host_organism_names", 
                      "source_organism_iri_search", "source_organism_iris", "source_organism_names", 
                      "mhc_allele_iri_search", "mhc_allele_iris", "mhc_allele_names", 
                      "parent_source_antigen_iri_search", "parent_source_antigen_iris", "parent_source_antigen_names", 
                      "parent_source_antigen_source_org_iris", "parent_source_antigen_source_org_names", "disease_iri_search", 
                      "disease_iris", "disease_names", "assay_iri_search", 
                      "assay_iris", "assay_names", "epitope_structures_defined", 
                      "non_peptidic_molecule_iri_search", "non_peptidic_molecule_iris", "non_peptidic_molecule_names", 
                      "r_object_source_molecule_iri_search", "r_object_source_molecule_iris", "r_object_source_molecule_names", 
                      "r_object_source_organism_iri_search", "r_object_source_organism_iris", "r_object_source_organism_names", 
                      "mhc_classes", "mhc_allele_resolutions", "e_related_object_types", 
                      "tcell_ids", "tcell_iris", "bcell_ids", 
                      "bcell_iris", "elution_ids", "elution_iris", 
                      "journal_names", "reference_types", "pubmed_ids", 
                      "reference_titles", "reference_authors", "reference_dates"
                    ),
                    default_values = list(offset = 0, order = "receptor_group_id"),
                    primary_key = "receptor_group_id"),
  bcell_search = list(endpoint = "bcell_search",
                      valid_parameters = c(
                        "bcell_id", "bcell_iri", "structure_id", 
                        "structure_iri", "linear_sequence", "structure_type", 
                        "structure_description", "curated_source_antigen", "reference_id", 
                        "reference_iri", "reference_type", "pubmed_id", 
                        "reference_authors", "reference_titles", "reference_dates", 
                        "journal_name", "epitope_summary", "reference_summary", 
                        "immunization_description", "antigen_description", "antigen_er", 
                        "assay_description", "e_modification", "linear_sequence_length", 
                        "qualitative_measure", "mhc_class", "mhc_allele_resolution", 
                        "submission_id", "submission_iri", "pdb_id", 
                        "chebi_ids", "mhc_allele_evidence", "antibody_isotype", 
                        "direct_ex_vivo_bool", "receptor_ids", "receptor_group_ids", 
                        "tcr_receptor_group_ids", "bcr_receptor_group_ids", "receptor_group_iris", 
                        "tcr_receptor_group_iris", "bcr_receptor_group_iris", "receptor_types", 
                        "receptor_names", "receptor_chain1_types", "receptor_chain2_types", 
                        "receptor_chain1_full_seqs", "receptor_chain2_full_seqs", "receptor_chain1_cdr1_seqs", 
                        "receptor_chain2_cdr1_seqs", "receptor_chain1_cdr2_seqs", "receptor_chain2_cdr2_seqs", 
                        "receptor_chain1_cdr3_seqs", "receptor_chain2_cdr3_seqs", "host_organism_iri_search", 
                        "host_organism_iri", "host_organism_name", "source_organism_iri_search", 
                        "source_organism_iri", "source_organism_name", "mhc_allele_iri_search", 
                        "mhc_allele_iri", "mhc_allele_name", "parent_source_antigen_iri_search", 
                        "parent_source_antigen_iri", "parent_source_antigen_name", "parent_source_antigen_source_org_iri", 
                        "parent_source_antigen_source_org_name", "disease_iri_search", "disease_iris", 
                        "disease_names", "assay_iri_search", "assay_iris", 
                        "assay_names", "epitope_structure_defined", "non_peptidic_molecule_iri_search", 
                        "non_peptidic_molecule_iri", "non_peptidic_molecule_name", "r_object_source_molecule_iri_search", 
                        "r_object_source_molecule_iri", "r_object_source_molecule_name", "r_object_source_organism_iri_search", 
                        "r_object_source_organism_iri", "r_object_source_organism_name", "e_related_object_type", 
                        "host_mhc_types_present"
                      ),
                      default_values = list(offset = 0, order ='bcell_id'),
                      primary_key = "bcell_id"),
  
  epitope_search = list(endpoint = "epitope_search",
                        valid_parameters = c(
                          "structure_id", "structure_iri", "structure_descriptions", 
                          "curated_source_antigens", "structure_type", "linear_sequence", 
                          "e_modification", "linear_sequence_length", "iedb_assay_ids", 
                          "iedb_assay_iris", "reference_ids", "reference_iris", 
                          "submission_ids", "submission_iris", "pdb_ids", 
                          "chebi_ids", "qualitative_measures", "mhc_allele_evidences", 
                          "antibody_isotypes", "direct_ex_vivo_bool", "receptor_ids", 
                          "receptor_group_ids", "tcr_receptor_group_ids", "bcr_receptor_group_ids", 
                          "receptor_group_iris", "tcr_receptor_group_iris", "bcr_receptor_group_iris", 
                          "receptor_types", "receptor_names", "receptor_chain1_types", 
                          "receptor_chain2_types", "receptor_chain1_full_seqs", "receptor_chain2_full_seqs", 
                          "receptor_chain1_cdr1_seqs", "receptor_chain2_cdr1_seqs", "receptor_chain1_cdr2_seqs", 
                          "receptor_chain2_cdr2_seqs", "receptor_chain1_cdr3_seqs", "receptor_chain2_cdr3_seqs", 
                          "host_organism_iri_search", "host_organism_iris", "host_organism_names", 
                          "source_organism_iri_search", "source_organism_iris", "source_organism_names", 
                          "mhc_allele_iri_search", "mhc_allele_iris", "mhc_allele_names", 
                          "parent_source_antigen_iri_search", "parent_source_antigen_iris", "parent_source_antigen_names", 
                          "parent_source_antigen_source_org_iris", "parent_source_antigen_source_org_names", "disease_iri_search", 
                          "disease_iris", "disease_names", "assay_iri_search", 
                          "assay_iris", "assay_names", "epitope_structures_defined", 
                          "non_peptidic_molecule_iri_search", "non_peptidic_molecule_iris", "non_peptidic_molecule_names", 
                          "r_object_source_molecule_iri_search", "r_object_source_molecule_iris", "r_object_source_molecule_names", 
                          "r_object_source_organism_iri_search", "r_object_source_organism_iris", "r_object_source_organism_names", 
                          "mhc_classes", "mhc_allele_resolutions", "e_related_object_types", 
                          "tcell_ids", "tcell_iris", "bcell_ids", 
                          "bcell_iris", "elution_ids", "elution_iris", 
                          "journal_names", "reference_types", "pubmed_ids", 
                          "reference_titles", "reference_authors", "reference_dates"
                        ),
                        default_values = list(offset = 0, order ='structure_id'),
                        primary_key = "structure_id"),
  mhc_search = list(endpoint = "mhc_search",
                    valid_parameters = c(
                      "elution_id", "elution_iri", "structure_id", 
                      "structure_iri", "linear_sequence", "structure_type", 
                      "structure_description", "curated_source_antigen", "reference_id", 
                      "reference_iri", "reference_type", "pubmed_id", 
                      "reference_authors", "reference_titles", "reference_dates", 
                      "journal_name", "epitope_summary", "reference_summary", 
                      "merged_host_imm_desc", "mhc_restriction", "quantitative_measure", 
                      "assay_description", "e_modification", "linear_sequence_length", 
                      "qualitative_measure", "mhc_class", "mhc_allele_resolution", 
                      "submission_id", "submission_iri", "pdb_id", 
                      "chebi_ids", "mhc_allele_evidence", "antibody_isotype", 
                      "direct_ex_vivo_bool", "receptor_ids", "receptor_group_ids", 
                      "tcr_receptor_group_ids", "bcr_receptor_group_ids", "receptor_group_iris", 
                      "tcr_receptor_group_iris", "bcr_receptor_group_iris", "receptor_types", 
                      "receptor_names", "receptor_chain1_types", "receptor_chain2_types", 
                      "receptor_chain1_full_seqs", "receptor_chain2_full_seqs", "receptor_chain1_cdr1_seqs", 
                      "receptor_chain2_cdr1_seqs", "receptor_chain1_cdr2_seqs", "receptor_chain2_cdr2_seqs", 
                      "receptor_chain1_cdr3_seqs", "receptor_chain2_cdr3_seqs", "host_organism_iri_search", 
                      "host_organism_iri", "host_organism_name", "source_organism_iri_search", 
                      "source_organism_iri", "source_organism_name", "mhc_allele_iri_search", 
                      "mhc_allele_iri", "mhc_allele_name", "parent_source_antigen_iri_search", 
                      "parent_source_antigen_iri", "parent_source_antigen_name", "parent_source_antigen_source_org_iri", 
                      "parent_source_antigen_source_org_name", "disease_iri_search", "disease_iris", 
                      "disease_names", "assay_iri_search", "assay_iris", 
                      "assay_names", "epitope_structure_defined", "non_peptidic_molecule_iri_search", 
                      "non_peptidic_molecule_iri", "non_peptidic_molecule_name", "r_object_source_molecule_iri_search", 
                      "r_object_source_molecule_iri", "r_object_source_molecule_name", "r_object_source_organism_iri_search", 
                      "r_object_source_organism_iri", "r_object_source_organism_name", "e_related_object_type", 
                      "host_mhc_types_present"
                    ),
                    default_values = list(offset = 0, order ='elution_id'),
                    primary_key = "elution_id"),
  reference_search = list(endpoint = "reference_search",
                          valid_parameters = c(
                            "reference_id", "reference_iri", "e_modifications", 
                            "linear_sequence_lengths", "qualitative_measures", "journal_name", 
                            "antibody_isotypes", "direct_ex_vivo_bool", "submission_ids", 
                            "submission_iris", "pdb_ids", "iedb_assay_ids", 
                            "iedb_assay_iris", "structure_ids", "structure_iris", 
                            "structure_types", "linear_sequences", "mhc_classes", 
                            "mhc_allele_resolutions", "parent_source_antigen_iris", "parent_source_antigen_names", 
                            "parent_source_antigen_source_org_iris", "parent_source_antigen_source_org_names", "e_related_object_types", 
                            "tcell_ids", "tcell_iris", "bcell_ids", 
                            "bcell_iris", "elution_ids", "elution_iris", 
                            "receptor_ids", "receptor_group_ids", "tcr_receptor_group_ids", 
                            "bcr_receptor_group_ids", "receptor_group_iris", "tcr_receptor_group_iris", 
                            "bcr_receptor_group_iris", "receptor_types", "receptor_names", 
                            "receptor_chain1_types", "receptor_chain2_types", "receptor_chain1_cdr1_seqs", 
                            "receptor_chain1_cdr2_seqs", "receptor_chain1_cdr3_seqs", "receptor_chain1_full_seqs", 
                            "receptor_chain2_cdr1_seqs", "receptor_chain2_cdr2_seqs", "receptor_chain2_cdr3_seqs", 
                            "receptor_chain2_full_seqs", "reference_title", "reference_title2", 
                            "reference_author", "reference_author2", "reference_date", 
                            "reference_date2", "pubmed_id", "mhc_allele_iri_search", 
                            "mhc_allele_iris", "mhc_allele_names", "assay_iri_search", 
                            "assay_iris", "assay_names", "epitope_structures_defined", 
                            "host_organism_iri_search", "host_organism_iris", "host_organism_names", 
                            "source_organism_iri_search", "source_organism_iris", "source_organism_names", 
                            "disease_iri_search", "disease_iris", "disease_names", 
                            "parent_source_antigen_iri_search", "non_peptidic_molecule_iri_search", "non_peptidic_molecule_iris", 
                            "non_peptidic_molecule_names", "r_object_source_molecule_iri_search", "r_object_source_molecule_iris", 
                            "r_object_source_molecule_names", "r_object_source_organism_iri_search", "r_object_source_organism_iris", 
                            "r_object_source_organism_names", "mhc_allele_evidences", "chebi_ids", 
                            "structure_descriptions"
                          ),
                          default_values = list(offset = 0, order ="reference_id"),
                          primary_key = "reference_id"),
  tcell_search = list(endpoint = 'tcell_search',
                      valid_parameters = c("tcell_id", "tcell_iri", "structure_id",
                                           "structure_iri", "linear_sequence", "structure_type",
                                           "structure_description", "curated_source_antigen", "reference_id",
                                           "reference_iri", "reference_type", "pubmed_id",
                                           "reference_authors", "reference_titles", "reference_dates",
                                           "journal_name", "epitope_summary", "reference_summary",
                                           "immunization_description", "antigen_description", "antigen_er",
                                           "mhc_restriction", "assay_description", "e_modification",
                                           "linear_sequence_length", "qualitative_measure", "mhc_class",
                                           "mhc_allele_resolution", "submission_id", "submission_iri",
                                           "pdb_id", "chebi_ids", "mhc_allele_evidence",
                                           "antibody_isotype", "direct_ex_vivo_bool", "receptor_ids",
                                           "receptor_group_ids", "tcr_receptor_group_ids", "bcr_receptor_group_ids",
                                           "receptor_group_iris", "tcr_receptor_group_iris", "bcr_receptor_group_iris",
                                           "receptor_types", "receptor_names", "receptor_chain1_types",
                                           "receptor_chain2_types", "receptor_chain1_full_seqs", "receptor_chain2_full_seqs",
                                           "receptor_chain1_cdr1_seqs", "receptor_chain2_cdr1_seqs",
                                           "receptor_chain1_cdr2_seqs", "receptor_chain2_cdr2_seqs",
                                           "receptor_chain1_cdr3_seqs", "receptor_chain2_cdr3_seqs",
                                           "host_organism_iri_search", "host_organism_iri", "host_organism_name",
                                           "source_organism_iri_search", "source_organism_iri", "source_organism_name",
                                           "mhc_allele_iri_search", "mhc_allele_iri", "mhc_allele_name",
                                           "parent_source_antigen_iri_search", "parent_source_antigen_iri", "parent_source_antigen_name",
                                           "parent_source_antigen_source_org_iri", "parent_source_antigen_source_org_name", "disease_iri_search",
                                           "disease_iris", "disease_names", "assay_iri_search",
                                           "assay_iris", "assay_names", "epitope_structure_defined",
                                           "non_peptidic_molecule_iri_search", "non_peptidic_molecule_iri", "non_peptidic_molecule_name",
                                           "r_object_source_molecule_iri_search", "r_object_source_molecule_iri", "r_object_source_molecule_name",
                                           "r_object_source_organism_iri_search", "r_object_source_organism_iri", "r_object_source_organism_name",
                                           "e_related_object_type", "host_mhc_types_present"),
                      default_values = list(offset = 0, order ='tcell_id'),
                      primary_key = 'tcell_id'),
  tcr_search = list(endpoint = 'tcr_search',
                    valid_parameters = c(
                      "receptor_group_id", "tcr_receptor_group_id", "bcr_receptor_group_id",
                      "receptor_group_iri", "tcr_receptor_group_iri", "bcr_receptor_group_iri",
                      "receptor_ids", "receptor_type", "receptor_species_names",
                      "chain1_cdr3_seq", "chain2_cdr3_seq", "e_modifications",
                      "linear_sequence_lengths", "reference_ids", "reference_iris",
                      "submission_ids", "submission_iris", "pdb_ids",
                      "chebi_ids", "qualitative_measures", "mhc_allele_evidences",
                      "antibody_isotypes", "direct_ex_vivo_bool", "iedb_assay_ids",
                      "iedb_assay_iris", "structure_ids", "structure_iris",
                      "structure_types", "structure_descriptions", "curated_source_antigens",
                      "linear_sequences", "receptor_names", "receptor_chain1_types",
                      "receptor_chain2_types", "receptor_chain1_full_seqs", "receptor_chain2_full_seqs",
                      "receptor_chain1_cdr1_seqs", "receptor_chain2_cdr1_seqs", "receptor_chain1_cdr2_seqs",
                      "receptor_chain2_cdr2_seqs", "receptor_chain1_cdr3_seqs", "receptor_chain2_cdr3_seqs",
                      "host_organism_iri_search", "host_organism_iris", "host_organism_names",
                      "source_organism_iri_search", "source_organism_iris", "source_organism_names",
                      "mhc_allele_iri_search", "mhc_allele_iris", "mhc_allele_names",
                      "parent_source_antigen_iri_search", "parent_source_antigen_iris", "parent_source_antigen_names",
                      "parent_source_antigen_source_org_iri", "parent_source_antigen_source_org_names", "disease_iri_search",
                      "disease_iris", "disease_names", "assay_iri_search",
                      "assay_iris", "assay_names", "epitope_structures_defined",
                      "non_peptidic_molecule_iri_search", "non_peptidic_molecule_iris", "non_peptidic_molecule_names",
                      "r_object_source_molecule_iri_search", "r_object_source_molecule_iris", "r_object_source_molecule_names",
                      "r_object_source_organism_iri_search", "r_object_source_organism_iris", "r_object_source_organism_names",
                      "mhc_classes", "mhc_allele_resolutions", "e_related_object_types",
                      "tcell_ids", "tcell_iris", "bcell_ids",
                      "bcell_iris", "elution_ids", "elution_iris",
                      "journal_names", "reference_types", "pubmed_ids",
                      "reference_titles", "reference_authors", "reference_dates"),
                    default_values = list(offset = 0, order ='receptor_group_id'),
                    primary_key = 'receptor_group_id')
  
  
)



species_name <- c('Escherichia coli','Homo sapiens','Saccharomyces cerevisiae', 'Spodoptera frugiperda', 'Staphylococcus aureus','Pseudomonas aeruginosa', '
Nicotiana tabacum', 'Gallus gallus','Drosophila melanogaster','Acinetobacter baumannii')

iedb_epi <- process_multiple_species(species_name, tables_metadata)




