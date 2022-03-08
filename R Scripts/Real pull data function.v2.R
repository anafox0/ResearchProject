# End of day: Mon 7/02/2022 - Code extracts expression and clinical data for those under 30
                            # GDCquery_Maf doesn't allow barcode filter - currently extracting all mutational data
                            # Barcodes are not universal across the platforms from what I can see


# Creating a singular function for data extraction of:
# Clinical data
# Transcriptomic/expression data
# Mutational data

real.pull.all.data <- function(
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

  # Downloading expression data only for those subsetted clinical patients
  
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

  # Downloading mutational data
  
  library(maftools)
  
  mutational <- GDCquery_Maf(tumor = tumor, pipelines = "mutect2")
  
  return(list(expression = expression,
              clinical = clinical,
              mutational = mutational))
  
}


lgg.all.data <- real.pull.all.data(project = "TCGA-LGG", 
                           tumor = "LGG", 
                           save.filename = "lgg.Exp.rda")

gbm.all.data <- real.pull.all.data(project = "TCGA-GBM",
                           tumor = "GBM",
                           save.filename = "gbm.Exp.rda")

str(lgg.all.data, max.level = 1)

### create a vector of patients we want to keep 
selected.patients <- lgg.all.data$clinical[lgg.all.data$clinical$age_at_initial_pathologic_diagnosis <30,"bcr_patient_barcode"]
### create a temporary copy of all data
lgg.all.data.filt.match <- lgg.all.data
### stepwise replace the indiviudal slots with filtered/matched data
### create an index for the expression slot where the numbers are columns of patients who are also in selected patients
which(lgg.all.data$expression$patient%in%selected.patients) -> exp.idx
lgg.all.data$expression[,exp.idx] -> lgg.all.data.filt.match$expression
### remove duplicated patients
lgg.all.data.filt.match$expression <- lgg.all.data.filt.match$expression[,which(!duplicated(lgg.all.data.filt.match$expression$patient))]

### create an match for the clinical data 
### remove exact duplicate lines
lgg.all.data$clinical <- lgg.all.data$clinical[which(!duplicated(lgg.all.data$clinical$bcr_patient_barcode)),]
match(lgg.all.data.filt.match$expression$patient, lgg.all.data$clinical$bcr_patient_barcode) -> clin.idx
clin.idx <- clin.idx[!is.na(clin.idx)]
lgg.all.data$clinical[clin.idx,] -> lgg.all.data.filt.match$clinical

lgg.all.data.filt.match$mutational <- lgg.all.data$mutational[substr(lgg.all.data$mutational$Tumor_Sample_Barcode,1,12)%in%selected.patients,]

### show that the two tables match
lgg.all.data.filt.match$clinical$bcr_patient_barcode
colnames(lgg.all.data.filt.match$expression)
lgg.all.data.filt.match$expression$patient
identical(lgg.all.data.filt.match$clinical$bcr_patient_barcode,
          lgg.all.data.filt.match$expression$patient)




#### spare code

# From the clinical data, subsetting only for the patients we are interested in
# subsetclinical <- clinical[clinical$age_at_initial_pathologic_diagnosis <30,]
# clinpatientsunder30 <- subsetclinical$bcr_patient_barcode


#rowSums(table(lgg.all.data.filt.match$expression$patient, lgg.all.data.filt.match$expression$barcode))