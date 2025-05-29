CREATE DATABASE HCP_Annotator_DB;
GO

USE HCP_Annotator_DB;
GO

CREATE TABLE[annotation_table](
 [Accession] varchar(255) UNIQUE NOT NULL,
 [protein_name] NVARCHAR(1025),
 [organism] varchar(255),
 [host_cell] varchar(255),
 [inference_level] varchar(255),
 [gene] varchar(450),
 [ec_numbers] varchar(255), 
 [enzyme_class] varchar(255),
 [enzyme_subclass] varchar(600), 
 [enzymatic_reaction_type] varchar(1000),
 [enzymatic_reaction] varchar(2000),
 [enzymatic_evidenceCode] varchar(255),
 [enzymatic_source] varchar(300),
 [Label] varchar(255),             
 [KEGG_id] varchar(255),
 [KEGG_id_hyperlink] nvarchar(255),
 [KEGG_pathway] varchar(7000),
 [KEGG_pathway_hyperlink] varchar(7000),
 [GO_ids] varchar(5000),    
 [GO_n_counts] INT,
 [GO_BP] varchar(max),
 [GO_MF] varchar(max),
 [GO_CC] varchar(max),
 [Potential_functional_impact] varchar(500),       
 [linear_sequence] varchar(max),
 [linear_sequence_length] varchar(255),
 [infoplace] varchar(255),
 [reference_id] varchar(max),
 [Critical_for] varchar(1500),      
 [old_Functional_Impact] varchar(255),
 [Function] varchar(1000),
 [Function_short] varchar(500),

CONSTRAINT prim_key_annotation PRIMARY KEY(Accession)
)
GO


CREATE TABLE[iedb_table](
 [Accession] varchar(255) UNIQUE NOT NULL,
 [iedb_assay_id] varchar(max),                       
 [reference_id] varchar(4000),                         
 [qualitative_measure] varchar(max),
 [Assay_number] INT,
 [Positive_High_assay_number] INT,
 [Positive_Intermediate_assay_number] INT,
 [Positive_assay_number] INT,
 [Positive_Low_assay_number] INT,            
 [Negative_assay_number] INT,
 [assay_names] varchar(max),                      
 [mhc_allele_name] varchar(max),
 [linear_sequence] varchar(max),
 [linear_sequence_length] varchar(max),   
 [parent_source_antigen_id] varchar(255),
 [parent_source_antigen_source_org_iri] varchar(255),
 [parent_source_antigen_names] varchar(2550),
 [parent_source_antigen_source_org_name] varchar(255),


---Primary Key
CONSTRAINT prim_key_iedb PRIMARY KEY(Accession)
)
GO


CREATE TABLE[GO_info](
 [id] varchar(350) UNIQUE NOT NULL,
 [name] varchar(600),
 [aspect] varchar(800),
 [text_def] varchar(3000),
 [keyword] varchar(255),
 [Potential_functional_impact] varchar(255),
---Primary Key
CONSTRAINT prim_key_go_info PRIMARY KEY(id)
)
GO

CREATE TABLE[GO_Accession_Annotator](
 [Accession] varchar(255) NOT NULL,
 [Qualifier] varchar(255),
 [GO_term] varchar(350) NOT NULL,
 [Reference] varchar(800),
 [Evidence] varchar(255),
 [From] varchar(3000), 
 [Assigned_By] varchar(255),
 [Annotation_Extension] varchar(255), 

 ---Foreign Keys
 CONSTRAINT forn_key_GO_Accession_Annotator_Accession FOREIGN KEY ([Accession]) REFERENCES [annotation_table] ([Accession]),
 CONSTRAINT forn_key_GO_Accession_Annotator_GO FOREIGN KEY ([GO_term]) REFERENCES [GO_info] ([id])
)
GO

CREATE TABLE references_table (
    [id] INT,
    [pmid] INT,
    [source_id] VARCHAR(255),
    [source] VARCHAR(255),
    [pmcid] INT,
    [doi] VARCHAR(255),
    [title] VARCHAR(1500),
    [authors] VARCHAR(4000),
    [journalname] VARCHAR(255),
    [issue] VARCHAR(255),
    [journalVolume] VARCHAR(255),
    [pubYear] INT,
    [pageInfo] VARCHAR(255),
    [pubType] VARCHAR(600),
    [isOpenAccess] VARCHAR(20),
    [inEPMC] VARCHAR(20),
    [inPMC] VARCHAR(20),
    [firstIndexDate] VARCHAR(255),
    [firstPublicationDate] VARCHAR(255),
	[id_iedb] INT
 ---Primary Key
CONSTRAINT prim_key_references PRIMARY KEY(id)
)
GO





