# Extracting the data from the Harmonized database
# Looking for expression data, mutational and clinical data for TCGA-LGG

datatable(readr::read_csv("https://docs.google.com/spreadsheets/d/1f98kFdj9mxVDc1dv4xTZdx8iWgUiDYO-qiFJINvmTZs/export?format=csv&gid=2046985454",col_types = readr::cols()),
          filter = 'top',
          options = list(scrollX = TRUE, keys = TRUE, pageLength = 40), 
          rownames = FALSE)
# The below reads in the relevant data file and gives the relevant categories needed for your queries

list.data.categories <- readr::read_csv("https://docs.google.com/spreadsheets/d/1f98kFdj9mxVDc1dv4xTZdx8iWgUiDYO-qiFJINvmTZs/export?format=csv&gid=2046985454",col_types = readr::cols())
list.data.categories

library(SummarizedExperiment)
library (TCGAbiolinks)

# Extracting transcriptome/expression data

query.exp1 <- GDCquery(
  project = "TCGA-LGG",
  data.category = "Transcriptome Profiling",
  legacy = FALSE,
  data.type = "Gene Expression Quantification", 
  workflow.type = "HTSeq - FPKM-UQ"
)

# The following lines allow you to look at the information you have about the data you have pulled with your query

query.exp1.results <- getResults(query.exp1)
head(query.exp1.results)

# Pulls/downloads the actual data from the above query

GDCdownload(query.exp1)

# GDCprepare takes the data downloads and prepares it into an R object
# Turns the downloaded files into a summerizedExperiment object (data.exp1)
# Containing both phenotypic/sample data and the actual data

data.exp1 <- GDCprepare(query.exp1)

# Extracts just the phenotype data

pData.exp1 <- as.data.frame(colData(data.exp1))
head(pData.exp1)

# Extracts just the expression data

expr.data.exp1 <- assay(data.exp1)
head(expr.data.exp1)

# Note: rownames of pData.exp1 match the colnames of expr.data.exp1






# Extracting mutational data

query.mut <- GDCquery(
  project = "TCGA-LGG",
  data.category = "Simple Nucleotide Variation",
  legacy = FALSE,
  data.type = "Raw Simple Somatic Mutation"
  
)

query.mut.results <- getResults(query.mut)
head(query.mut.results)

# Gets the results table from the above query
getResults(query.mut)






# Extracting clinical data
query.clin <- GDCquery(
  project = "TCGA-LGG",
  data.category = "Clinical",
  legacy = FALSE,
  data.type = "Clinical Supplement"
  
)

datatable(getResults(query.clin), 
          filter = 'top',
          options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
          rownames = FALSE)

# OR EXTRACT CLINICAL DATA USING THE CODE BELOW I FOUND IN CASE STUDY #3:

dataClin <- GDCquery_clinic(project = "TCGA-LGG")


# Get all patients that have gene expression, mutational and clinical data
# The code below does not work - trying to find the intersection between the patients for which all data is available

common.patients <- intersect(
  substr(getResults(query.mut, cols = "cases"), 1, 12),
  substr(getResults(query.exp1, cols = "cases"), 1, 12)
)

query.exp2 <- GDCquery(
  project = "TCGA-LGG",
  data.category = "Transcriptome Profiling",
  legacy = FALSE,
  data.type = "Gene Expression Quantification", 
  workflow.type = "HTSeq - FPKM-UQ",
  barcode = common.patients
)


head(dataClin)
subsetDataClin <- dataClin[dataClin$age_at_index <30,]
patientsunder30LGG <- subsetDataClin$bcr_patient_barcode






# Questions
# Is this the correct way to go about downloading the data
  # i.e. 3 seperate GDCquery functions
  # Or is there something similar to the GDCquery_clinc function I found in Case Study 1 that I should be using
  # What is the difference between the GDCquery_clinic and GDCquery functions

# Also I do not understand how to apply age filters - and at what step/function should I apply these filters?
  # Should I be applying these filters in each of the GDCqueries?

# How do I make the intersection directly above these questions work?
 
 # Am I even on the right lines with the intersection/do I need this for what you would like me to do


