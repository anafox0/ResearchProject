# Extracting the data from the Harmonized database
# Looking for expression data, mutational and clinical data for TCGA-LGG

datatable(readr::read_csv("https://docs.google.com/spreadsheets/d/1f98kFdj9mxVDc1dv4xTZdx8iWgUiDYO-qiFJINvmTZs/export?format=csv&gid=2046985454",col_types = readr::cols()),
          filter = 'top',
          options = list(scrollX = TRUE, keys = TRUE, pageLength = 40), 
          rownames = FALSE)

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

datatable(getResults(query.exp1), 
          filter = 'top',
          options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
          rownames = FALSE)


# Extracting mutational data

query.mut <- GDCquery(
  project = "TCGA-LGG",
  data.category = "Simple Nucleotide Variation",
  legacy = FALSE,
  data.type = "Raw Simple Somatic Mutation"
  
)

datatable(getResults(query.mut), 
          filter = 'top',
          options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
          rownames = FALSE)

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
  substr(getResults(query.exp, cols = "cases")),
  substr(getResults(query.clin, cols = "cases"))
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
