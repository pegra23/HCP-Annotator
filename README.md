# HCP Annotator Database

This repository contains the code  from the master's thesis project "Enhancing the Identification of Host Cell Proteins of Concern in Biopharmaceuticals: Implementing an Integrative Approach." The project's goal is to improve the current centralized information on each Host Cell Protein (HCP), as well as to develop a risk assessment classification system for HCPs found in client samples.

## Table of Contents

* [About the Project](#about-the-project)
* [Project Structure](#project-structure)
* [Installation](#installation)
    * [R Packages](#r-packages)
* [Database](#database)
* [Author](#author)

---

## [About the Project](#about-the-project)

This master's thesis focuses primarily on the characterization of lesser-known organisms, extracting features from each protein that may be relevant when a protein is present in a client sample.  As a result, a database titled "HCP Annotator" was created to serve as the centralized data repository for the information gathered throughout this project.  This repository contains the various R scripts created for gathering information via different APIs.


The project includes:
* **R scripts** for data handling, analysis, and potentially database interaction.
* **SQL scripts** for database creation.

---

## Project Structure

An overview of the main folders and files in this repository:
* **`R/`**: Contains all R source code for the project.
* **`data/`**: A duplicate data file from the result SQL database.
* **`sql/`**: Contains SQL scripts related to your database schema and operations.

## Installation


### R Packages

Open RStudio or your R console and run the following commands to install the necessary packages.

```R
libraries <- c(
  "httr", "jsonlite", "tidyr","KEGGREST", "readr","tibble, "stringr", "dplyr","brendaDb"
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
```

---

### Database

* Ensure you have Microsoft SQL Server installed and running (e.g., SQL Server Express or Developer Edition).
* Make sure you have SQL Server Management Studio (SSMS) or Azure Data Studio installed to easily interact with your SQL Server instance.
* Create the database: Connect to your SQL Server instance and execute the following SQL command to create your database (if it doesn't already exist):

```sql
CREATE DATABASE [HCPAnnotatorDB];
GO
USE [HCPAnnotatorDB];
GO
```
A diagram of the database schema can be found in docs/database_design.pdf (if you have an graphic ER diagram).


---
### Author
Pedro Mauritti Granjo - granjo@alphalyse.com


