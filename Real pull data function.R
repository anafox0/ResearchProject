# Creating a singular function for data extraction of:
# Clinical data
# Transcriptomic/expression data
# Mutational data

real.pull.data <- function(
  project = "",
  tumor = "" ,
  workflow.type = "HTSeq - FPKM-UQ",
  save.filename = ""
){
  
  # Downloading expression data
  
  query.exp <- GDCquery(
    project = project,
    data.category = "Transcriptome Profiling",
    legacy = FALSE,
    data.type = "Gene Expression Quantification",
    workflow.type = workflow.type
  )
  
  GDCdownload(query.exp)
  expression <- GDCprepare(query = query.exp,
                           save = TRUE,
                           save.filename = save.filename)
  
  
  # Downloading clinical data
  
  query.clin <- GDCquery(
    project = project,
    file.type = "xml",
    data.category = "Clinical",
    legacy = FALSE,
    data.type = "Clinical Supplement"
  )
  
  GDCdownload(query.clin)
  clinical <- GDCprepare_clinic(query = query.clin, clinical.info = "patient")
  

  # Downloading mutational data
  
  library(maftools)
  
  mutational <- GDCquery_Maf(tumor = tumor, pipelines = "mutect2")

  return(list(expression = expression,
              clinical = clinical,
              mutational = mutational))
  
  
}


lgg.data <- real.pull.data(project = "TCGA-LGG",
                           tumor = "LGG",
                           save.filename = "LGG.Exp.rda"
)







