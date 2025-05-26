# HCP Annotator Database

This repository contains the code and database for the master's thesis project "Enhancing the Identification of Host Cell Proteins of Concern in Biopharmaceuticals: Implementing an Integrative Approach." The project's goal is to improve the current centralized information on each HCP, as well as to develop a risk assessment classification system for HCP found in client samples.

## Table of Contents

* [About the Project](#about-the-project)
* [Project Structure](#project-structure)
* [Installation](#installation)
    * [R Packages](#r-packages)
* [Usage](#usage)
* [Database](#database)
* [Author](#author)
* [License](#license)

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

To get the project up and running locally, follow these steps:

1.  **Clone the repository:**

    ```bash
    git clone [https://github.com/pegra23/HCP-Annotator-Database.git](https://github.com/YOUR_GITHUB_USERNAME/HCP-Annotator-Database.git)
    cd HCP-Annotator-Database
    ```

### R Packages

Open RStudio or your R console and run the following commands to install the necessary packages.

```R
ibraries <- c(
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

