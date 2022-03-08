# End of day: Mon 7/02/2022 - Code extracts expression and clinical data for those under 30
                            # GDCquery_Maf doesn't allow barcode filter - currently extracting all mutational data
                            # Barcodes are not universal across the platforms from what I can see


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
  
  # From the clinical data, subsetting only for the patients we are interested in
  subsetclinical <- clinical[clinical$age_at_initial_pathologic_diagnosis <30,]
  clinpatientsunder30 <- subsetclinical$bcr_patient_barcode
  
  # Downloading expression data only for those subsetted clinical patients
  
  query.exp <- GDCquery(
    project = project,
    data.category = "Transcriptome Profiling",
    legacy = FALSE,
    data.type = "Gene Expression Quantification",
    workflow.type = "HTSeq - FPKM-UQ",
    barcode = clinpatientsunder30
  )
  
  GDCdownload(query.exp)
  
  expression <- GDCprepare(query = query.exp,
                           save = TRUE,
                           save.filename = save.filename)

  # Downloading mutational data
  
  library(maftools)
  
  mutational.pre <- GDCquery_Maf(tumor = tumor, pipelines = "mutect2")
  
  mutational <- intersect(substr(mutational.pre$Tumor_Sample_Barcode, 1, 12),
                          substr(subsetclinical$bcr_patient_barcode, 1, 12)
  )



  return(list(expression = expression,
              clinical = clinical,
              mutational = mutational))
  
  
}

)

lgg.data <- real.pull.data(project = "TCGA-LGG", 
                           tumor = "LGG", 
                           save.filename = "lgg.Exp.rda")


gbm.data <- real.pull.data(project = "TCGA-GBM",
                           tumor = "GBM",
                           save.filename = "gbm.Exp.rda")

