# 1 - Introduction --------------------------------------------------------


#TCGAbiolinks is able to access The National Cancer Institute (NCI) Genomic Data Commons (GDC) thorough its
#GDC Application Programming Interface (API) to search, download and prepare relevant data for analysis in R.

#Installation of the stable version

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("TCGAbiolinks")

a
#Libraries required for this tutorial:

library(TCGAbiolinks)
library(dplyr)
library(DT)


# Searching GDC Database --------------------------------------------------

# GDC consists of two sources: Harmonized and Legacy Archive
# GDC Legacy Archive:
#     Provides access to unmodified copy of data that was previously stored in CGHub and in the TCGA Data Portal
#     hosted by the TCGA Data Coordinating Center, in which uses as references 
#     GRCh37 (hg19) and GRCh36 (hg18)

# GDC Harmonized Database
#     Data available was harmonized against GRCh38 (hg38) using GDC Bioinformatics Pipelines
#     providing methods to the standardisation of biospecimen and clinical data

# Understanding the barcode
# TCGA barcode = composed of collectio nof identifiers
# Each specifically identifies a TCGA data element
# Aliquot Barcode - TCGA-G4-6317-02A-11D-2064-05
# Participant Barcode: TCGA-G4-6317
# Sample Barcode: TCGA-G4-6317-02



# Searching Arguments
# Easily search GDC data using the GDCquery function
# Function uses a summary of filters as used in the TCGA portal to use the following arguments

# Data.category - project (TCGA-BRCA, GENIE-MSK, GENIE-VICC, GENIE-UHN, CPTAC-2, CMI-ASC etc.)
# data.type
# workflow.type
# legacy
# access
# platform
# file.type
# barcode
# experimental.strategy
# sample.type - TP, TR, TB, TRBM etc.

# Harmonized data options (legacy = FALSE)

datatable(readr::read_csv("https://docs.google.com/spreadsheets/d/1f98kFdj9mxVDc1dv4xTZdx8iWgUiDYO-qiFJINvmTZs/export?format=csv&gid=2046985454",col_types = readr::cols()),
          filter = 'top',
          options = list(scrollX = TRUE, keys = TRUE, pageLength = 40), 
          rownames = FALSE)

# Harmonized Database Examples

# DNA methylation data: Recurrent tumor samples
# In this example we will access the harmonized database (legacy = FALSE) and search for all DNA methylation data for recurrent glioblastoma multiform (GBM) and low grade gliomas (LGG) samples.

query <- GDCquery(
  project = c("TCGA-GBM", "TCGA-LGG"),
  data.category = "DNA Methylation",
  legacy = FALSE,
  platform = c("Illumina Human Methylation 450"),
  sample.type = "Recurrent Tumor"
)
datatable(getResults(query), 
          filter = 'top',
          options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
          rownames = FALSE)

# Samples with DNA methylation and gene expression data
# In this example we will access the harmonized database (legacy = FALSE) and search for all patients with DNA methylation (platform HumanMethylation450k) and gene expression data for Colon Adenocarcinoma tumor (TCGA-COAD).

query.met <- GDCquery(
  project = "TCGA-COAD",
  data.category = "DNA Methylation",
  legacy = FALSE,
  platform = c("Illumina Human Methylation 450")
)
query.exp <- GDCquery(
  project = "TCGA-COAD",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification", 
  workflow.type = "HTSeq - FPKM-UQ"
)

# Get all patients that have DNA methylation and gene expression.
common.patients <- intersect(
  substr(getResults(query.met, cols = "cases"), 1, 12),
  substr(getResults(query.exp, cols = "cases"), 1, 12)
)

# Only select the first 5 patients
query.met <- GDCquery(
  project = "TCGA-COAD",
  data.category = "DNA Methylation",
  legacy = FALSE,
  platform = c("Illumina Human Methylation 450"),
  barcode = common.patients[1:5]
)
query.exp <- GDCquery(
  project = "TCGA-COAD",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification", 
  workflow.type = "HTSeq - FPKM-UQ",
  barcode = common.patients[1:5]
)


datatable(
  getResults(query.met, cols = c("data_type","cases")),
  filter = 'top',
  options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
  rownames = FALSE
)
datatable(
  getResults(query.exp, cols = c("data_type","cases")), 
  filter = 'top',
  options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
  rownames = FALSE
)

# Raw Sequencing Data: Finding the match between file names and barcode for Controlled data.
# This example shows how the user can search for breast cancer Raw Sequencing Data (“Controlled”) and verify the name of the files and the barcodes associated with it.

query <- GDCquery(
  project = c("TCGA-BRCA"),
  data.category = "Sequencing Reads",  
  sample.type = "Primary Tumor"
)
# Only first 100 to make render faster
datatable(
  getResults(query, rows = 1:100,cols = c("file_name","cases")), 
  filter = 'top',
  options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
  rownames = FALSE
)




# Legacy Archive Data Options (legacy = TRUE)

datatable(readr::read_csv("https://docs.google.com/spreadsheets/d/1f98kFdj9mxVDc1dv4xTZdx8iWgUiDYO-qiFJINvmTZs/export?format=csv&gid=1817673686",col_types = readr::cols()),
          filter = 'top',
          options = list(scrollX = TRUE, keys = TRUE, pageLength = 40), 
          rownames = FALSE)

# Legacy archive examples
# DNA methylation
    #Array-based assays
    #This example shows how the user can search for glioblastoma multiform (GBM) and DNA methylation data for platform Illumina Human Methylation 450 and Illumina Human Methylation 27.

query <- GDCquery(
  project = c("TCGA-GBM"),
  legacy = TRUE,
  data.category = "DNA methylation",
  platform = c("Illumina Human Methylation 450", "Illumina Human Methylation 27")
)
datatable(
  getResults(query, rows = 1:100), 
  filter = 'top',
  options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
  rownames = FALSE
)

    # Whole-genome bisulfite sequencing (WGBS)

query <- GDCquery(
  project = c("TCGA-LUAD"),
  legacy = TRUE,
  data.category = "DNA methylation",
  data.type = "Methylation percentage",
  experimental.strategy = "Bisulfite-Seq"
)

# VCF - controlled data
query <- GDCquery(
  project = c("TCGA-LUAD"),
  legacy = TRUE,
  data.category = "DNA methylation",
  data.type = "Bisulfite sequence alignment",
  experimental.strategy = "Bisulfite-Seq"
)


# WGBS BAM files - controlled data
query <- GDCquery(
  project = c("TCGA-LUAD"),
  legacy = TRUE,
  data.type = "Aligned reads",
  data.category = "Raw sequencing data",
  experimental.strategy = "Bisulfite-Seq"
)

# Gene expression
# This example shows how the user can search for glioblastoma multiform (GBM) gene expression data with the normalized results for expression of a gene. For more information about file.types check GDC TCGA file types

# Gene expression aligned against hg19.
query.exp.hg19 <- GDCquery(
  project = "TCGA-GBM",
  data.category = "Gene expression",
  data.type = "Gene expression quantification",
  platform = "Illumina HiSeq", 
  file.type  = "normalized_results",
  experimental.strategy = "RNA-Seq",
  barcode = c("TCGA-14-0736-02A-01R-2005-01", "TCGA-06-0211-02A-02R-2005-01"),
  legacy = TRUE
)

datatable(
  getResults(query.exp.hg19), 
  filter = 'top',
  options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
  rownames = FALSE
)


# Get Manifest file
# If you want to get the manifest file from the query object you can use the function getManifest. If you set save to TRUEm a txt file that can be used with GDC-client Data transfer tool (DTT) or with its GUI version ddt-ui will be created.

getManifest(query.exp.hg19,save = FALSE) 

# ATAC-seq data

datatable(
  getResults(TCGAbiolinks:::GDCquery_ATAC_seq())[,c("file_name","file_size")], 
  filter = 'top',
  options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
  rownames = FALSE
)

# You can use the function GDCquery_ATAC_seq filter the manifest table and use GDCdownload to save the data locally.

query <- TCGAbiolinks:::GDCquery_ATAC_seq(file.type = "rds") 
GDCdownload(query,method = "client")

query <- TCGAbiolinks:::GDCquery_ATAC_seq(file.type = "bigWigs") 
GDCdownload(query,method = "client")


# Summary of available files per patient
# Retrieve the numner of files under each data_category + data_type + experimental_strategy + platform. Almost like https://portal.gdc.cancer.gov/exploration

tab <-  getSampleFilesSummary(project = "TCGA-ACC")
datatable(
  head(tab),
  filter = 'top',
  options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
  rownames = FALSE
)



# 3 - Downloading and preparing files for analysis ------------------------

# 2 methods to download GDC data using TCGAbiolinks
# Client - this method creates a MANIFEST file and downloads the data using the GDC Data Transfer Tool
    # More reliable method but slightly slower vs api method

# api method - Uses GDC Application Programming Interface (API) to download the data
    # Creates a MANIFEST file and downloaded data is compressed into tar.gz file
    # If the size/number of files are too big this tar.gz will be too big
    # --> high probability of download failure
    # Solve this via files.per.chunk argument - splits files into small chunks
    # If chunks.per.download = 10, we will download 10 files inside each tar.gz

# DATA PREPARED - SUMMARIZEDEXPERIMENT OBJECT

# A SummarizedExperiment object has three main matrices that can be accessed using the SummarizedExperiment package):

# 1 : Sample matrix information is accessed via colData(data): stores sample information. TCGAbiolinks will add indexed clinical data and subtype information from marker TCGA papers.
# 2 : Assay matrix information is accessed via assay(data): stores molecular data
# 3 : Feature matrix information (gene information) is accessed via rowRanges(data): stores metadata about the features, including their genomic ranges

# SUMMARIZED EXPERIMENT: ANNOTATION INFORMATION

# When using the function GDCprepare there is an argument called SummarizedExperiment which defines the output type of a Summarized Experiment (default option) or a data frame. 

# To create a summarized Experiment object we annotate the data with genomic positions with last patch release version of the genome available. 
# For legacy data (data aligned to hg19) TCGAbiolinks is using GRCh37.p13 
# For harmonized data (data aligned to hg38) now it is using GRCh38.p7 (May 2017).

# Unfortunately, some of the updates changes/remove gene symbols, change coordinates, etc. Which might introduce some loss of data. 
# For example, if the gene was removed we cannot map it anymore and that information will be lost in the SummarizedExperiment.

# If you set SummarizedExperiment to FALSE, you will get the data unmodified just as they are in the files and ad your own annotation.

# Also, there are no updated for DNA methylation data. But the last metadata available can be found here: http://zwdzwd.github.io/InfiniumAnnotation

# GDCdownload Arguments
# query, token.file, method, directory, files.per.chunk

# GDCprepare Arguments
# query, save, save.filename, directory, summarizedExperiment, remove.files.prepared, add.gistic2.mut, mut.pipeline, mutant_variant_classification

# Search and download data from legacy databses using GDC api method

# In this example we will download gene expression data from legacy database (data aligned against genome of reference hg19) using GDC api method and we will show object data and metadata.

query <- GDCquery(
  project = "TCGA-GBM",
  data.category = "Gene expression",
  data.type = "Gene expression quantification",
  platform = "Illumina HiSeq", 
  file.type  = "normalized_results",
  experimental.strategy = "RNA-Seq",
  barcode = c("TCGA-14-0736-02A-01R-2005-01", "TCGA-06-0211-02A-02R-2005-01"),
  legacy = TRUE
)
GDCdownload(query, method = "api", files.per.chunk = 10)
data <- GDCprepare(query)

# Gene expression aligned against hg19.
datatable(as.data.frame(colData(data)), 
          options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
          rownames = FALSE)


# Only first 100 to make render faster
datatable(assay(data)[1:100,], 
          options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
          rownames = TRUE)

rowRanges(data)

# Search and download data for two samples from database
# In this example we will download gene expression quantification from harmonized database (data aligned against genome of reference hg38). Also, it shows the object data and metadata.

# Gene expression aligned against hg38
query <- GDCquery(
  project = "TCGA-GBM",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification", 
  workflow.type = "HTSeq - FPKM-UQ",
  barcode = c("TCGA-14-0736-02A-01R-2005-01", "TCGA-06-0211-02A-02R-2005-01")
)
GDCdownload(query)
data <- GDCprepare(query)
datatable(as.data.frame(colData(data)), 
          options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
          rownames = FALSE)


datatable(assay(data)[1:100,], 
          options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
          rownames = TRUE)


# GDCprepare: Outputs
# This function is still under development, it is not working for all cases. See the tables below with the status. Examples of query, download, prepare can be found in this gist.
# See this page to view which cases are working within Harmonized and/or Legacy data

# Examples

# Legacy archive
# DNA methylation: Get all TCGA IDAT files

# Example to idat files from TCGA projects
projects <- TCGAbiolinks:::getGDCprojects()$project_id
projects <- projects[grepl('^TCGA',projects,perl=T)]
match.file.cases.all <- NULL
for(proj in projects){
  print(proj)
  query <- GDCquery(
    project = proj,
    data.category = "Raw microarray data",
    data.type = "Raw intensities", 
    experimental.strategy = "Methylation array", 
    legacy = TRUE,
    file.type = ".idat",
    platform = "Illumina Human Methylation 450"
  )
  match.file.cases <- getResults(query,cols=c("cases","file_name"))
  match.file.cases$project <- proj
  match.file.cases.all <- rbind(match.file.cases.all,match.file.cases)
  tryCatch(
    GDCdownload(query, method = "api", files.per.chunk = 20),
    error = function(e) GDCdownload(query, method = "client")
  )
}

# This will create a map between idat file name, cases (barcode) and project
readr::write_tsv(match.file.cases.all, path =  "idat_filename_case.txt")

# code to move all files to local folder
for(file in dir(".",pattern = ".idat", recursive = T)){
  TCGAbiolinks::move(file,basename(file))
}




# DNA methylation: aligned against hg19
query_meth.hg19 <- GDCquery(
  project= "TCGA-LGG", 
  data.category = "DNA methylation", 
  platform = "Illumina Human Methylation 450", 
  barcode = c("TCGA-HT-8111-01A-11D-2399-05","TCGA-HT-A5R5-01A-11D-A28N-05"), 
  legacy = TRUE
)
GDCdownload(query_meth.hg19)
data.hg19 <- GDCprepare(query_meth.hg19)

# Protein expression
query <- GDCquery(
  project = "TCGA-GBM",
  data.category = "Protein expression",
  legacy = TRUE, 
  barcode = c("TCGA-OX-A56R-01A-21-A44T-20","TCGA-08-0357-01A-21-1898-20")
)
GDCdownload(query)
data <- GDCprepare(
  query, save = TRUE, 
  save.filename = "gbmProteinExpression.rda",
  remove.files.prepared = TRUE
)

# Gene expression: aligned against hg19
# Aligned against Hg19
query.exp.hg19 <- GDCquery(
  project = "TCGA-GBM",
  data.category = "Gene expression",
  data.type = "Gene expression quantification",
  platform = "Illumina HiSeq", 
  file.type  = "normalized_results",
  experimental.strategy = "RNA-Seq",
  barcode = c("TCGA-14-0736-02A-01R-2005-01", "TCGA-06-0211-02A-02R-2005-01"),
  legacy = TRUE
)
GDCdownload(query.exp.hg19)
data <- GDCprepare(query.exp.hg19)

# Harmonized database examples
# Copy Number
query <- GDCquery(
  project = "TCGA-ACC", 
  data.category = "Copy Number Variation",
  data.type = "Copy Number Segment",
  barcode = c( "TCGA-OR-A5KU-01A-11D-A29H-01", "TCGA-OR-A5JK-01A-11D-A29H-01")
)
GDCdownload(query)
data <- GDCprepare(query)

# GISTIC2
query <- GDCquery(
  project = "TCGA-ACC",
  data.category = "Copy Number Variation",
  data.type = "Gene Level Copy Number Scores",              
  access = "open"
)
GDCdownload(query)
data <- GDCprepare(query)

# Gene expression: aligned against hg38
# For more examples, please check: http://rpubs.com/tiagochst/TCGAbiolinks_RNA-seq_new_projects

# mRNA pipeline: https://gdc-docs.nci.nih.gov/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline/
query.exp.hg38 <- GDCquery(
  project = "TCGA-GBM", 
  data.category = "Transcriptome Profiling", 
  data.type = "Gene Expression Quantification", 
  workflow.type = "HTSeq - FPKM-UQ",
  barcode =  c("TCGA-14-0736-02A-01R-2005-01", "TCGA-06-0211-02A-02R-2005-01")
)
GDCdownload(query.exp.hg38)
expdat <- GDCprepare(
  query = query.exp.hg38,
  save = TRUE, 
  save.filename = "exp.rda"
)

# miRNA
library(TCGAbiolinks)
query.mirna <- GDCquery(
  project = "TARGET-AML", 
  experimental.strategy = "miRNA-Seq",
  data.category = "Transcriptome Profiling", 
  barcode = c("TARGET-20-PATDNN","TARGET-20-PAPUNR"),
  data.type = "miRNA Expression Quantification"
)
GDCdownload(query.mirna)
mirna <- GDCprepare(
  query = query.mirna,
  save = TRUE, 
  save.filename = "mirna.rda"
)


query.isoform <- GDCquery(
  project = "TARGET-AML", 
  experimental.strategy = "miRNA-Seq",
  data.category = "Transcriptome Profiling", 
  barcode = c("TARGET-20-PATDNN","TARGET-20-PAPUNR"),
  data.type = "Isoform Expression Quantification"
)
GDCdownload(query.isoform)

isoform <- GDCprepare(
  query = query.isoform,
  save = TRUE, 
  save.filename = "mirna-isoform.rda"
)
# DNA methylation: aligned against hg38

# DNA methylation data

# DNA methylation aligned to hg38
query_met.hg38 <- GDCquery(
  project= "TCGA-LGG", 
  data.category = "DNA Methylation", 
  platform = "Illumina Human Methylation 450", 
  barcode = c("TCGA-HT-8111-01A-11D-2399-05","TCGA-HT-A5R5-01A-11D-A28N-05")
)
GDCdownload(query_met.hg38)
data.hg38 <- GDCprepare(query_met.hg38)

# DNA methylation IDAT
# Using sesame  http://bioconductor.org/packages/sesame/
# Please cite 10.1093/nar/gky691 and doi: 10.1093/nar/gkt090.
library(TCGAbiolinks)
proj <- "TCGA-ACC"
query <- GDCquery(
  project = proj,
  data.category = "Raw microarray data",
  data.type = "Raw intensities", 
  experimental.strategy = "Methylation array", 
  legacy = TRUE,
  barcode = c("TCGA-OR-A5JT","CGA-OR-A5LG","TCGA-OR-A5JX"),
  file.type = ".idat",
  platform = "Illumina Human Methylation 450"
)
tryCatch(
  GDCdownload(query, method = "api", files.per.chunk = 20),
  error = function(e) GDCdownload(query, method = "client")
)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("sesame")

a
betas <- GDCprepare(query)


# 4 - Clinical Data -------------------------------------------------------

# Clinical Data can be retireved from the GDC database from different sources:
# Indexed Clinical: Refined clinical data created using the XML files
# XML files: Original source of the data
# BCR Biotab: tsv files parsed from XML files

# 2 main differences between indexed clinical and XML files:
# XML - more information: radiation, drugs information, follow-ups, biospecimen etc.
        # Indexed = subset of XML files
# Indexed - contains updated data with the follow-up information 


# BCR Biotab

# Clinical
# In this example, we will fetch clinical data from BCR Biotab files
query <- GDCquery(project = "TCGA-ACC", 
                  data.category = "Clinical",
                  data.type = "Clinical Supplement", 
                  data.format = "BCR Biotab")
GDCdownload(query)
clinical.BCRtab.all <- GDCprepare(query)
names(clinical.BCRtab.all)

query <- GDCquery(project = "TCGA-ACC", 
                  data.category = "Clinical",
                  data.type = "Clinical Supplement", 
                  data.format = "BCR Biotab",
                  file.type = "radiation")
GDCdownload(query)
clinical.BCRtab.radiation <- GDCprepare(query)

clinical.BCRtab.all$clinical_drug_acc  %>% 
  head  %>% 
  DT::datatable(options = list(scrollX = TRUE, keys = TRUE))

# In this example we will fetch all BRCA BCR Biotab files, and look for the ER status.

library(TCGAbiolinks)
query <- GDCquery(project = "TCGA-BRCA", 
                  data.category = "Clinical",
                  data.type = "Clinical Supplement", 
                  data.format = "BCR Biotab")
GDCdownload(query)
clinical.BCRtab.all <- GDCprepare(query)

# All available tables
names(clinical.BCRtab.all)

## [1] "clinical_follow_up_v4.0_brca"     "clinical_follow_up_v1.5_brca"    
## [3] "clinical_follow_up_v4.0_nte_brca" "clinical_omf_v4.0_brca"          
## [5] "clinical_nte_brca"                "clinical_radiation_brca"         
## [7] "clinical_follow_up_v2.1_brca"     "clinical_patient_brca"           
## [9] "clinical_drug_brca"

# colnames from clinical_patient_brca
tibble::tibble(sort(colnames(clinical.BCRtab.all$clinical_patient_brca)))

# ER status count
plyr::count(clinical.BCRtab.all$clinical_patient_brca$er_status_by_ihc)

# ER content 
er.cols <- grep("^er",colnames(clinical.BCRtab.all$clinical_patient_brca))
clinical.BCRtab.all$clinical_patient_brca[,c(2,er.cols)] %>% 
  DT::datatable(options = list(scrollX = TRUE))

# All columns content first rows
clinical.BCRtab.all$clinical_patient_brca %>% 
  head  %>% 
  DT::datatable(options = list(scrollX = TRUE, keys = TRUE))


# Biospecimen
# Biospecimen BCR Biotab
query.biospecimen <- GDCquery(project = "TCGA-BRCA", 
                              data.category = "Biospecimen",
                              data.type = "Biospecimen Supplement", 
                              data.format = "BCR Biotab")
GDCdownload(query.biospecimen)
biospecimen.BCRtab.all <- GDCprepare(query.biospecimen)

# All available tables
names(biospecimen.BCRtab.all)
##  [1] "ssf_tumor_samples_brca"             "biospecimen_sample_brca"           
##  [3] "biospecimen_shipment_portion_brca"  "biospecimen_slide_brca"            
##  [5] "biospecimen_analyte_brca"           "ssf_normal_controls_brca"          
##  [7] "biospecimen_diagnostic_slides_brca" "biospecimen_portion_brca"          
##  [9] "biospecimen_aliquot_brca"           "biospecimen_protocol_brca"
biospecimen.BCRtab.all$ssf_normal_controls_ov  %>% 
  head  %>% 
  DT::datatable(options = list(scrollX = TRUE, keys = TRUE))

# Clinical indexed data
# In this example we will fetch clinical indexed data (same as showed in the data portal).

clinical <- GDCquery_clinic(project = "TCGA-LUAD", type = "clinical")
clinical %>%
  head %>% 
  DT::datatable(filter = 'top', 
                options = list(scrollX = TRUE, keys = TRUE, pageLength = 5),  
                rownames = FALSE)

clinical <- GDCquery_clinic(project = "BEATAML1.0-COHORT", type = "clinical")
clinical %>% 
  head %>% 
  DT::datatable(filter = 'top', 
                options = list(scrollX = TRUE, keys = TRUE, pageLength = 5),  
                rownames = FALSE)

clinical <- GDCquery_clinic(project = "CPTAC-2", type = "clinical")
clinical %>% 
  head %>% 
  DT::datatable(filter = 'top', 
                options = list(scrollX = TRUE, keys = TRUE, pageLength = 5),  
                rownames = FALSE)

clinical <- GDCquery_clinic(project = "GENIE-MSK", type = "clinical")
clinical %>% 
  head %>% 
  DT::datatable(filter = 'top', 
                options = list(scrollX = TRUE, keys = TRUE, pageLength = 5),  
                rownames = FALSE)


# XML clinical data
# The process to get data directly from the XML are: 
# 1. Use GDCquery and GDCDownload functions to search/download either biospecimen or clinical XML files 
# 2. Use GDCprepare_clinic function to parse the XML files.

# The relation between one patient and other clinical information are 1:n, one patient could have several radiation treatments. 
# For that reason, we only give the option to parse individual tables (only drug information, only radiation information,…) 
# The selection of the table is done by the argument clinical.info.

# Review the page to look at the clinical.info options to parse information for each data category

# Below are several examples fetching clinical data directly from the clinical XML files.

query <- GDCquery(project = "TCGA-COAD", 
                  data.category = "Clinical", 
                  file.type = "xml", 
                  barcode = c("TCGA-RU-A8FL","TCGA-AA-3972"))
GDCdownload(query)
clinical <- GDCprepare_clinic(query, clinical.info = "patient")
clinical %>% 
  datatable(filter = 'top', 
            options = list(scrollX = TRUE, keys = TRUE, pageLength = 5),  
            rownames = FALSE)

clinical.drug <- GDCprepare_clinic(query, clinical.info = "drug")
clinical.drug %>% 
  datatable(filter = 'top', 
            options = list(scrollX = TRUE, keys = TRUE, pageLength = 5),  
            rownames = FALSE)

clinical.radiation <- GDCprepare_clinic(query, clinical.info = "radiation")
clinical.radiation %>% 
  datatable(filter = 'top', 
            options = list(scrollX = TRUE, keys = TRUE, pageLength = 5),  
            rownames = FALSE)

clinical.admin <- GDCprepare_clinic(query, clinical.info = "admin")
clinical.admin %>% 
  datatable(filter = 'top', 
            options = list(scrollX = TRUE, keys = TRUE, pageLength = 5),  
            rownames = FALSE)

# Microsatellite data
# MSI-Mono-Dinucleotide Assay is performed to test a panel of four mononucleotide repeat loci (polyadenine tracts BAT25, BAT26, BAT40, and transforming growth factor receptor type II) and three dinucleotide repeat loci (CA repeats in D2S123, D5S346, and D17S250). 
# Two additional pentanucleotide loci (Penta D and Penta E) are included in this assay to evaluate sample identity. 
# Multiplex fluorescent-labeled PCR and capillary electrophoresis were used to identify MSI if a variation in the number of microsatellite repeats was detected between tumor and matched non-neoplastic tissue or mononuclear blood cells. 
# Equivocal or failed markers were re-evaluated by singleplex PCR.

# Classifications: microsatellite-stable (MSS), low level MSI (MSI-L) if less than 40% of markers were altered and high level MSI (MSI-H) if greater than 40% of markers were altered.

# Reference: TCGA wiki

# Level 3 data is included in BCR clinical-based submissions and can be downloaded as follows:

query <- GDCquery(project = "TCGA-COAD", 
                  data.category = "Other",
                  legacy = TRUE,
                  access = "open",
                  data.type = "Auxiliary test",
                  barcode = c("TCGA-AD-A5EJ","TCGA-DM-A0X9"))  
GDCdownload(query)
msi_results <- GDCprepare_clinic(query, "msi")
msi_results %>% DT::datatable(options = list(scrollX = TRUE, keys = TRUE))


# Tissue Slide Image (SVS format)

# Tissue slide image files from legacy database
query.legacy <- GDCquery(project = "TCGA-COAD", 
                         data.category = "Clinical", 
                         data.type = "Tissue slide image",
                         legacy = TRUE,
                         barcode = c("TCGA-RU-A8FL","TCGA-AA-3972")) 

# Tissue slide image files from harmonized database
query.harmonized <- GDCquery(project = "TCGA-OV",
                             data.category = "Biospecimen",
                             data.type = 'Slide Image')
query.legacy %>% 
  getResults %>% 
  DT::datatable(options = list(scrollX = TRUE, keys = TRUE))


query.harmonized  %>% 
  getResults %>% 
  head  %>% 
  DT::datatable(options = list(scrollX = TRUE, keys = TRUE))

# Diagnostic Slide (SVS format)

# Pathology report from harmonized portal 
query.harmonized <- GDCquery(project = "TCGA-COAD", 
                             data.category = "Biospecimen", 
                             data.type = "Slide Image",
                             experimental.strategy = "Diagnostic Slide",
                             barcode = c("TCGA-RU-A8FL","TCGA-AA-3972"))  
query.harmonized  %>% 
  getResults %>% 
  head  %>% 
  DT::datatable(options = list(scrollX = TRUE, keys = TRUE))

# Legacy archive files
# The clinical data types available in legacy database are:
  
# 1 - Biospecimen data (Biotab format)
# 2 - Tissue slide image (SVS format)
# 3 - Clinical Supplement (XML format)
# 4 - Pathology report (PDF)
# 5 - Clinical data (Biotab format)

# Pathology report (PDF)
# Pathology report from legacy portal 
query.legacy <- GDCquery(project = "TCGA-COAD", 
                         data.category = "Clinical", 
                         data.type = "Pathology report",
                         legacy = TRUE,
                         barcode = c("TCGA-RU-A8FL","TCGA-AA-3972"))  
query.legacy %>% 
  getResults %>% 
  DT::datatable(options = list(scrollX = TRUE, keys = TRUE))

# Tissue slide image (SVS format)
# Tissue slide image
query <- GDCquery(project = "TCGA-COAD", 
                  data.category = "Clinical", 
                  data.type = "Tissue slide image",
                  legacy = TRUE,
                  barcode = c("TCGA-RU-A8FL","TCGA-AA-3972")) 
query %>% 
  getResults %>% 
  DT::datatable(options = list(scrollX = TRUE, keys = TRUE))

# Clinical Supplement (XML format)
# Clinical Supplement
query <- GDCquery(project = "TCGA-COAD", 
                  data.category = "Clinical", 
                  data.type = "Clinical Supplement",
                  legacy = TRUE,
                  barcode = c("TCGA-RU-A8FL","TCGA-AA-3972")) 
query %>% 
  getResults %>% 
  DT::datatable(options = list(scrollX = TRUE, keys = TRUE))

# Clinical data (Biotab format)
# Clinical data
query <- GDCquery(project = "TCGA-COAD", 
                  data.category = "Clinical", 
                  data.type = "Clinical data",
                  legacy = TRUE,
                  file.type = "txt")  
query %>% 
  getResults %>% 
  select(-matches("cases"))%>% 
  DT::datatable(options = list(scrollX = TRUE, keys = TRUE))
  
GDCdownload(query)
clinical.biotab <- GDCprepare(query)
names(clinical.biotab)
## [1] "clinical_radiation_coad"          "clinical_nte_coad"               
## [3] "clinical_patient_coad"            "clinical_drug_coad"              
## [5] "clinical_follow_up_v1.0_nte_coad" "clinical_omf_v4.0_coad"          
## [7] "clinical_follow_up_v1.0_coad"
datatable(clinical.biotab$clinical_radiation_coad, options = list(scrollX = TRUE, keys = TRUE))

# Filter Functions
# Function TCGAquery_SampleTypes will filter barcodes based on the argument typesample

# type sample - character vector indicating the tissue type to query
    # e.g. TP (primary tumour), TR (recurrent tumour)

# TCGAquery_MatchedCoupledSampleTypes - filters the samples that have all the typesample provided as an argument
# e.g. if TP and TR set as typesample, the function will return the barcodes of a patient if it has both types
# So if it has a TP, but not a TR, no barcode will be returned. 
# If it has a TP AND a TR, then BOTH barcodes wil be returned


bar <- c("TCGA-G9-6378-02A-11R-1789-07", "TCGA-CH-5767-04A-11R-1789-07",  
         "TCGA-G9-6332-60A-11R-1789-07", "TCGA-G9-6336-01A-11R-1789-07",
         "TCGA-G9-6336-11A-11R-1789-07", "TCGA-G9-7336-11A-11R-1789-07",
         "TCGA-G9-7336-04A-11R-1789-07", "TCGA-G9-7336-14A-11R-1789-07",
         "TCGA-G9-7036-04A-11R-1789-07", "TCGA-G9-7036-02A-11R-1789-07",
         "TCGA-G9-7036-11A-11R-1789-07", "TCGA-G9-7036-03A-11R-1789-07",
         "TCGA-G9-7036-10A-11R-1789-07", "TCGA-BH-A1ES-10A-11R-1789-07",
         "TCGA-BH-A1F0-10A-11R-1789-07", "TCGA-BH-A0BZ-02A-11R-1789-07",
         "TCGA-B6-A0WY-04A-11R-1789-07", "TCGA-BH-A1FG-04A-11R-1789-08",
         "TCGA-D8-A1JS-04A-11R-2089-08", "TCGA-AN-A0FN-11A-11R-8789-08",
         "TCGA-AR-A2LQ-12A-11R-8799-08", "TCGA-AR-A2LH-03A-11R-1789-07",
         "TCGA-BH-A1F8-04A-11R-5789-07", "TCGA-AR-A24T-04A-55R-1789-07",
         "TCGA-AO-A0J5-05A-11R-1789-07", "TCGA-BH-A0B4-11A-12R-1789-07",
         "TCGA-B6-A1KN-60A-13R-1789-07", "TCGA-AO-A0J5-01A-11R-1789-07",
         "TCGA-AO-A0J5-01A-11R-1789-07", "TCGA-G9-6336-11A-11R-1789-07",
         "TCGA-G9-6380-11A-11R-1789-07", "TCGA-G9-6380-01A-11R-1789-07",
         "TCGA-G9-6340-01A-11R-1789-07", "TCGA-G9-6340-11A-11R-1789-07")

S <- TCGAquery_SampleTypes(bar,"TP")
S2 <- TCGAquery_SampleTypes(bar,"NB")

# Retrieve multiple tissue types  NOT FROM THE SAME PATIENTS
SS <- TCGAquery_SampleTypes(bar,c("TP","NB"))

# Retrieve multiple tissue types  FROM THE SAME PATIENTS
SSS <- TCGAquery_MatchedCoupledSampleTypes(bar,c("NT","TP"))



# To get all the information for TCGA samples you can use the code below:

# This code will get all clinical indexed data from TCGA
library(data.table)
library(dplyr)
library(regexPipes)
clinical <- TCGAbiolinks:::getGDCprojects()$project_id %>% 
  regexPipes::grep("TCGA",value=T) %>% 
  sort %>% 
  plyr::alply(1,GDCquery_clinic, .progress = "text") %>% 
  rbindlist
readr::write_csv(clinical,path = paste0("all_clin_indexed.csv"))

# This code will get all clinical XML data from TCGA
getclinical <- function(proj){
  message(proj)
  while(1){
    result = tryCatch({
      query <- GDCquery(project = proj, data.category = "Clinical",file.type = "xml")
      GDCdownload(query)
      clinical <- GDCprepare_clinic(query, clinical.info = "patient")
      for(i in c("admin","radiation","follow_up","drug","new_tumor_event")){
        message(i)
        aux <- GDCprepare_clinic(query, clinical.info = i)
        if(is.null(aux) || nrow(aux) == 0) next
        # add suffix manually if it already exists
        replicated <- which(grep("bcr_patient_barcode",colnames(aux), value = T,invert = T) %in% colnames(clinical))
        colnames(aux)[replicated] <- paste0(colnames(aux)[replicated],".",i)
        if(!is.null(aux)) clinical <- merge(clinical,aux,by = "bcr_patient_barcode", all = TRUE)
      }
      readr::write_csv(clinical,path = paste0(proj,"_clinical_from_XML.csv")) # Save the clinical data into a csv file
      return(clinical)
    }, error = function(e) {
      message(paste0("Error clinical: ", proj))
    })
  }
}
clinical <- TCGAbiolinks:::getGDCprojects()$project_id %>% 
  regexPipes::grep("TCGA",value=T) %>% sort %>% 
  plyr::alply(1,getclinical, .progress = "text") %>% 
  rbindlist(fill = TRUE) %>% setDF %>% subset(!duplicated(clinical))

readr::write_csv(clinical,path = "all_clin_XML.csv")
# result: https://drive.google.com/open?id=0B0-8N2fjttG-WWxSVE5MSGpva1U
# Obs: this table has multiple lines for each patient, as the patient might have several followups, drug treatments,
# new tumor events etc...



# 5 - Mutation Data -------------------------------------------------------

# Search and Download
# TCGAbiolinks has provided a few functions to download mutation data from GDC. There are two options to download the data:
  
# 1. Use GDCquery_Maf which will download MAF aligned against hg38
# 2. Use GDCquery, GDCdownload and GDCpreprare to download MAF aligned against hg19
# 3. Use getMC3MAF(), to download MC3 MAF from https://gdc.cancer.gov/about-data/publications/mc3-2017

# Mutation data (hg38)
# This example will download MAF (mutation annotation files) for variant calling pipeline muse. Pipelines options are: muse, varscan2, somaticsniper, mutect. For more information please access GDC docs.

maf <- GDCquery_Maf("CHOL", pipelines = "muse")
# Only first 50 to make render faster
datatable(maf[1:20,],
          filter = 'top',
          options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
          rownames = FALSE)

# Mutation data (hg19)
# This example will download MAF (mutation annotation files) aligned against hg19 (Old TCGA maf files)
query.maf.hg19 <- GDCquery(project = "TCGA-CHOL", 
                           data.category = "Simple nucleotide variation", 
                           data.type = "Simple somatic mutation",
                           access = "open", 
                           legacy = TRUE)
# Check maf availables
datatable(dplyr::select(getResults(query.maf.hg19),-contains("cases")),
          filter = 'top',
          options = list(scrollX = TRUE, keys = TRUE, pageLength = 10), 
          rownames = FALSE)

query.maf.hg19 <- GDCquery(project = "TCGA-CHOL", 
                           data.category = "Simple nucleotide variation", 
                           data.type = "Simple somatic mutation",
                           access = "open", 
                           file.type = "bcgsc.ca_CHOL.IlluminaHiSeq_DNASeq.1.somatic.maf",
                           legacy = TRUE)
GDCdownload(query.maf.hg19)
maf <- GDCprepare(query.maf.hg19)

# Only first 50 to make render faster
datatable(maf[1:20,],
          filter = 'top',
          options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
          rownames = FALSE)

# Mutation data MC3 file
# This will download the MC3 MAF file from https://gdc.cancer.gov/about-data/publications/mc3-2017, and add project each sample belongs.

maf <- getMC3MAF()

# Visualize the data
# To visualize the data you can use the Bioconductor package maftools. For more information, please check its vignette.

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("maftools")

library(maftools)
library(dplyr)
maf <- GDCquery_Maf("CHOL", pipelines = "muse") %>% read.maf
datatable(getSampleSummary(maf),
          filter = 'top',
          options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
          rownames = FALSE)
plotmafSummary(maf = maf, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE)

oncoplot(maf = maf, top = 10, removeNonMutated = TRUE)
titv = titv(maf = maf, plot = FALSE, useSyn = TRUE)
#plot titv summary
plotTiTv(res = titv)


# 6 - Compilation of TCGA Molecular Subtypes ------------------------------

# TCGAbiolinks retrieved molecular subtypes information from TCGA samples. 
# The functions PanCancerAtlas_subtypes and TCGAquery_subtype can be used to get the information tables.

# While the PanCancerAtlas_subtypes function gives access to a curated table retrieved from synapse (probably with the most updated molecular subtypes)
# TCGAquery_subtype function has the complete table also with sample information retrieved from the TCGA marker papers.


# PanCancerAtlas_subtypes: Curated molecular subtypes.
# Data and description retrieved from synapse (https://www.synapse.org/#!Synapse:syn8402849)
                                               
    # Synapse has published a single file with all available molecular subtypes that have been described by TCGA (all tumor types and all molecular platforms)
    # Can be accessed using the PanCancerAtlas_subtypes function as below:
subtypes <- PanCancerAtlas_subtypes()
DT::datatable(subtypes,
      filter = 'top',
      options = list(scrollX = TRUE, keys = TRUE, pageLength = 5),
      rownames = FALSE)

# TCGAquery_subtype: Working with Molecular Subtypes Data
#The Cancer Genome Atlas (TCGA) Research Network has reported integrated genome-wide studies of various diseases

#These subtypes will be automatically added in the summarizedExperiment object through GDCprepare. But you can also use the TCGAquery_subtype function to retrieve this information.
lgg.gbm.subtype <- TCGAquery_subtype(tumor = "lgg")

#Session Information

sessionInfo()

# Analysing and Visualising TCGA Data -------------------------------------

# TCGAanalyze: Analyze data from TCGA

# TCGAanalyze_Preprocessing: Preprocessing of Gene Expression data (IlluminaHiSeq_RNASeqV2)
# You can easily search TCGA samples, download and prepare a matrix of gene expression.

# You can define a list of samples to query and download providing relative TCGA barcodes.
listSamples <- c("TCGA-E9-A1NG-11A-52R-A14M-07","TCGA-BH-A1FC-11A-32R-A13Q-07",
                 "TCGA-A7-A13G-11A-51R-A13Q-07","TCGA-BH-A0DK-11A-13R-A089-07",
                 "TCGA-E9-A1RH-11A-34R-A169-07","TCGA-BH-A0AU-01A-11R-A12P-07",
                 "TCGA-C8-A1HJ-01A-11R-A13Q-07","TCGA-A7-A13D-01A-13R-A12P-07",
                 "TCGA-A2-A0CV-01A-31R-A115-07","TCGA-AQ-A0Y5-01A-11R-A14M-07")

# Query platform Illumina HiSeq with a list of barcode 
query <- GDCquery(project = "TCGA-BRCA", 
                  data.category = "Gene expression",
                  data.type = "Gene expression quantification",
                  experimental.strategy = "RNA-Seq",
                  platform = "Illumina HiSeq",
                  file.type = "results",
                  barcode = listSamples, 
                  legacy = TRUE)

# Download a list of barcodes with platform IlluminaHiSeq_RNASeqV2
GDCdownload(query)

# Prepare expression matrix with geneID in the rows and samples (barcode) in the columns
# rsem.genes.results as values
BRCARnaseqSE <- GDCprepare(query)

BRCAMatrix <- assay(BRCARnaseqSE,"raw_count") # or BRCAMatrix <- assay(BRCARnaseqSE,"raw_count")

# For gene expression if you need to see a boxplot correlation and AAIC plot to define outliers you can run
BRCARnaseq_CorOutliers <- TCGAanalyze_Preprocessing(BRCARnaseqSE)

# TCGAanalyze_DEA() & TCGAanalyze_LevelTab(): Differential Expression Analysis

# Perform DEA (Differential expression analysis) to identify differentially expressed genes (DEGs) using the TCGAanalyze_DEA function.

# TCGAanalyze_DEA performs DEA using following functions from R :
  
# edgeR::DGEList - Converts the count matrix into an edgeR object.
# edgeR::estimateCommonDisp - Each gene gets assigned the same dispersion estimate.
# edgeR::exactTest - Performs pair-wise tests for differential expression between two groups.
# edgeR::topTags - Takes the output from exactTest(), adjusts the raw p-values using the False Discovery Rate (FDR) correction, and returns the top differentially expressed genes.

# This function receives as arguments:
  
# mat1 - The matrix of the first group (in the example, group 1 is the normal samples),
# mat2 - The matrix of the second group (in the example, group 2 is tumor samples)
# Cond1type - Label for group 1
# Cond1type - Label for group 2

#Next, we filter the output of dataDEGs by abs(LogFC) >=1
# And use the TCGAanalyze_LevelTab function to create a table with:
# DEGs (differentially expressed genes)
# log Fold Change (FC)
# False discovery rate (FDR) 
# Gene expression level for samples in Cond1type, and Cond2type, 
# Delta value (the difference of gene expression between the two conditions multiplied logFC).

# Downstream analysis using gene expression data  
# TCGA samples from IlluminaHiSeq_RNASeqV2 with type rsem.genes.results
# save(dataBRCA, geneInfo , file = "dataGeneExpression.rda")
library(TCGAbiolinks)

# normalization of genes
dataNorm <- TCGAanalyze_Normalization(tabDF = dataBRCA, geneInfo =  geneInfo)

# quantile filter of genes
dataFilt <- TCGAanalyze_Filtering(tabDF = dataNorm,
                                  method = "quantile", 
                                  qnt.cut =  0.25)

# selection of normal samples "NT"
samplesNT <- TCGAquery_SampleTypes(barcode = colnames(dataFilt),
                                   typesample = c("NT"))

# selection of tumor samples "TP"
samplesTP <- TCGAquery_SampleTypes(barcode = colnames(dataFilt), 
                                   typesample = c("TP"))

# Diff.expr.analysis (DEA)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("edgeR")

dataDEGs <- TCGAanalyze_DEA(mat1 = dataFilt[,samplesNT],
                            mat2 = dataFilt[,samplesTP],
                            Cond1type = "Normal",
                            Cond2type = "Tumor",
                            fdr.cut = 0.01 ,
                            logFC.cut = 1,
                            method = "glmLRT")

# DEGs table with expression values in normal and tumor samples
dataDEGsFiltLevel <- TCGAanalyze_LevelTab(dataDEGs,"Tumor","Normal",
                                          dataFilt[,samplesTP],dataFilt[,samplesNT])


# HTSeq data: Downstream analysis BRCA

CancerProject <- "TCGA-BRCA"
DataDirectory <- paste0("../GDC/",gsub("-","_",CancerProject))
FileNameData <- paste0(DataDirectory, "_","HTSeq_Counts",".rda")

query <- GDCquery(project = CancerProject,
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification", 
                  workflow.type = "HTSeq - Counts")

samplesDown <- getResults(query,cols=c("cases"))

dataSmTP <- TCGAquery_SampleTypes(barcode = samplesDown,
                                  typesample = "TP")



dataSmNT <- TCGAquery_SampleTypes(barcode = samplesDown,
                                  typesample = "NT")
dataSmTP_short <- dataSmTP[1:10]
dataSmNT_short <- dataSmNT[1:10]

queryDown <- GDCquery(project = CancerProject, 
                      data.category = "Transcriptome Profiling",
                      data.type = "Gene Expression Quantification", 
                      workflow.type = "HTSeq - Counts", 
                      barcode = c(dataSmTP_short, dataSmNT_short))

GDCdownload(query = queryDown,
            directory = DataDirectory)

dataPrep <- GDCprepare(query = queryDown, 
                       save = TRUE, 
                       directory =  DataDirectory,
                       save.filename = FileNameData)

dataPrep <- TCGAanalyze_Preprocessing(object = dataPrep, 
                                      cor.cut = 0.6,
                                      datatype = "HTSeq - Counts")                      

dataNorm <- TCGAanalyze_Normalization(tabDF = dataPrep,
                                      geneInfo = geneInfoHT,
                                      method = "gcContent") 

boxplot(dataPrep, outline = FALSE)

boxplot(dataNorm, outline = FALSE)

dataFilt <- TCGAanalyze_Filtering(tabDF = dataNorm,
                                  method = "quantile", 
                                  qnt.cut =  0.25)   

dataDEGs <- TCGAanalyze_DEA(mat1 = dataFilt[,dataSmTP_short],
                            mat2 = dataFilt[,dataSmNT_short],
                            Cond1type = "Normal",
                            Cond2type = "Tumor",
                            fdr.cut = 0.01 ,
                            logFC.cut = 1,
                            method = "glmLRT")  


# miRNA expression data: Downstream analysis BRCA

require(TCGAbiolinks)

CancerProject <- "TCGA-BRCA"
DataDirectory <- paste0("../GDC/",gsub("-","_",CancerProject))
FileNameData <- paste0(DataDirectory, "_","miRNA_gene_quantification",".rda")

query.miR <- GDCquery(project = CancerProject, 
                      data.category = "Gene expression",
                      data.type = "miRNA gene quantification",
                      file.type = "hg19.miRNA",
                      legacy = TRUE)

samplesDown.miR <- getResults(query.miR,cols=c("cases"))

dataSmTP.miR <- TCGAquery_SampleTypes(barcode = samplesDown.miR,
                                      typesample = "TP")

dataSmNT.miR <- TCGAquery_SampleTypes(barcode = samplesDown.miR,
                                      typesample = "NT")
dataSmTP_short.miR <- dataSmTP.miR[1:10]
dataSmNT_short.miR <- dataSmNT.miR[1:10]

queryDown.miR <- GDCquery(project = CancerProject, 
                          data.category = "Gene expression",
                          data.type = "miRNA gene quantification",
                          file.type = "hg19.mirna",
                          legacy = TRUE,
                          barcode = c(dataSmTP_short.miR, dataSmNT_short.miR))

GDCdownload(query = queryDown.miR,
            directory = DataDirectory)

dataAssy.miR <- GDCprepare(query = queryDown.miR, 
                           save = TRUE, 
                           save.filename = FileNameData, 
                           summarizedExperiment = TRUE,
                           directory =DataDirectory )
rownames(dataAssy.miR) <- dataAssy.miR$miRNA_ID

# using read_count's data 
read_countData <-  colnames(dataAssy.miR)[grep("count", colnames(dataAssy.miR))]
dataAssy.miR <- dataAssy.miR[,read_countData]
colnames(dataAssy.miR) <- gsub("read_count_","", colnames(dataAssy.miR))

dataFilt <- TCGAanalyze_Filtering(tabDF = dataAssy.miR,
                                  method = "quantile", 
                                  qnt.cut =  0.25)   

dataDEGs <- TCGAanalyze_DEA(mat1 = dataFilt[,dataSmNT_short.miR],
                            mat2 = dataFilt[,dataSmTP_short.miR],
                            Cond1type = "Normal",
                            Cond2type = "Tumor",
                            fdr.cut = 0.01 ,
                            logFC.cut = 1,
                            method = "glmLRT")  

# TCGAanalyze_EAcomplete() & TCGAvisualize_EAbarplot(): Enrichment Analysis

# Researchers, in order to better understand the underlying biological processes, often want to retrieve a functional profile of a set of genes that might have an important role. This can be done by performing an enrichment analysis.

# We will perform an enrichment analysis on gene sets using the TCGAanalyze_EAcomplete function. Given a set of genes that are up-regulated under certain conditions, an enrichment analysis will identify classes of genes or proteins that are over-represented using annotations for that gene set.

#To view the results you can use the TCGAvisualize_EAbarplot function as shown below.

library(TCGAbiolinks)
# Enrichment Analysis EA
# Gene Ontology (GO) and Pathway enrichment by DEGs list
Genelist <- rownames(dataDEGsFiltLevel)

system.time(ansEA <- TCGAanalyze_EAcomplete(TFname="DEA genes Normal Vs Tumor",Genelist))

# Enrichment Analysis EA (TCGAVisualize)
# Gene Ontology (GO) and Pathway enrichment barPlot

TCGAvisualize_EAbarplot(tf = rownames(ansEA$ResBP), 
                        GOBPTab = ansEA$ResBP,
                        GOCCTab = ansEA$ResCC,
                        GOMFTab = ansEA$ResMF,
                        PathTab = ansEA$ResPat,
                        nRGTab = Genelist, 
                        nBar = 10)


# TCGAanalyze_survival(): Survival Analysis

# When analyzing survival, different problems come up than the ones discussed so far. 
# One question is how do we deal with subjects dropping out of a study.
# For example, assuming that we test a new cancer drug. 
# While some subjects die, others may believe that the new drug is not effective, and decide to drop out of the study before the study is finished. 
# A similar problem would be faced when we investigate how long a machine lasts before it breaks down.

# Using the clinical data, it is possible to create a survival plot with the function TCGAanalyze_survival as follows:
  
  clin.gbm <- GDCquery_clinic("TCGA-GBM", "clinical")
TCGAanalyze_survival(clin.gbm,
                     "gender",
                     main = "TCGA Set\n GBM",height = 10, width=10)

#The arguments of TCGAanalyze_survival are:
  
# clinical_patient - TCGA Clinical patient with the information days_to_death
# clusterCol - Column with groups to plot. This is a mandatory field, the caption will be based in this column
# legend - Legend title of the figure
# xlim - xlim x axis limits e.g. xlim = c(0, 1000). Present narrower X axis, but not affect survival estimates.
# main - main title of the plot
# ylab - y-axis text of the plot
# xlab - x-axis text of the plot
# filename - The name of the pdf file
# color - Define the colors of the lines.
# pvalue - Show pvalue in the plot.
# risk.table - Show or not the risk table
# conf.int - Show confidence intervals for point estimates of survival curves.


# TCGAanalyze_SurvivalKM(): Correlating gene expression and survival analysis

library(TCGAbiolinks)
# Survival Analysis SA

clinical_patient_Cancer <- GDCquery_clinic("TCGA-BRCA","clinical")
dataBRCAcomplete <- log2(BRCA_rnaseqv2)

tokenStop<- 1

tabSurvKMcomplete <- NULL

for( i in 1: round(nrow(dataBRCAcomplete)/100)){
  message( paste( i, "of ", round(nrow(dataBRCAcomplete)/100)))
  tokenStart <- tokenStop
  tokenStop <-100*i
  tabSurvKM<-TCGAanalyze_SurvivalKM(clinical_patient_Cancer,
                                    dataBRCAcomplete,
                                    Genelist = rownames(dataBRCAcomplete)[tokenStart:tokenStop],
                                    Survresult = F,
                                    ThreshTop=0.67,
                                    ThreshDown=0.33)
  
  tabSurvKMcomplete <- rbind(tabSurvKMcomplete,tabSurvKM)
}

tabSurvKMcomplete <- tabSurvKMcomplete[tabSurvKMcomplete$pvalue < 0.01,]
tabSurvKMcomplete <- tabSurvKMcomplete[order(tabSurvKMcomplete$pvalue, decreasing=F),]

tabSurvKMcompleteDEGs <- tabSurvKMcomplete[
  rownames(tabSurvKMcomplete) %in% dataDEGsFiltLevel$mRNA,
]


#TCGAanalyze_DMR: Differentially methylated regions analysis



# We will search for differentially methylated CpG sites using the TCGAanalyze_DMR function. 
#In order to find these regions we use the beta-values (methylation values ranging from 0.0 to 1.0) to compare two groups.

# First, it calculates the difference between the mean DNA methylation of each group for each probe.

# Second, it test for differential expression between two groups using the wilcoxon test adjusting by the Benjamini-Hochberg method.
# The default arguments was set to require a minimum absolute beta-values difference of 0.2 and an adjusted p-value of < 0.01.

# After these tests, we save a volcano plot (x-axis:diff mean methylation, y-axis: statistical significance) that will help the user identify the differentially methylated CpG sites, then the results are saved in a csv file (DMR_results.groupCol.group1.group2.csv) and finally the object is returned with the calculus in the rowRanges.

# Please view page to see the arguments of TCGAanalyze_DMR

data <- TCGAanalyze_DMR(data, groupCol = "methylation_subtype",
                        group1 = "CIMP.H",
                        group2="CIMP.L",
                        p.cut = 10^-5,
                        diffmean.cut = 0.25,
                        legend = "State",
                        plot.filename = "coad_CIMPHvsCIMPL_metvolcano.png")


# Also, the TCGAanalyze_DMR function will save the plot as pdf and return the same SummarizedExperiment that was given as input with the values of p-value, p-value adjusted, diffmean and the group it belongs in the graph (non significant, hypomethylated, hypermethylated) in the rowRanges. 
# The columns will be (where group1 and group2 are the names of the groups):
  
# diffmean.group1.group2 (mean.group2 - mean.group1)
# diffmean.group2.group1 (mean.group1 - mean.group2)
# p.value.group1.group2
# p.value.adj.group1.group2
# status.group1.group2 (Status of probes in group2 in relation to group1)
# status.group2.group1 (Status of probes in group1 in relation to group2)
#This values can be view/acessed using the rowRanges accessesor (rowRanges(data)).

#Observation: Calling the same function again, with the same arguments will only plot the results, as it was already calculated. 
# If you want to have them recalculated, please set overwrite to TRUE or remove the calculated columns.


#TCGAvisualize: Visualize results from analysis functions with TCGA's data

#TCGAvisualize_Heatmap: Create heatmaps with cluster bars
#In order to have a better view of clusters, we normally use heatmaps. 
#TCGAvisualize_Heatmap will plot a heatmap and add to each sample bars representing different features. 
# This function is a wrapper to the package ComplexHeatmap package,

# The arguments of this function are:
  
# data - The object with the heatmap data (expression, methylation)
# col.metadata - Metadata for the columns (patients). It should have the column bcr_patient_barcode or patient or ID with the patients barcodes.
# row.metadata - Metadata for the rows genes (expression) or probes (methylation)
# col.colors - A list of names colors
# row.colors - A list of named colors
# show_column_names - Show column names names? Default: FALSE
# show_row_names - Show row names? Default: FALSE
# cluster_rows - Cluster rows ? Default: FALSE
# cluster_columns - Cluster columns ? Default: FALSE
# sortCol - Name of the column to be used to sort the columns
# title - Title of the plot
# type - Select the colors of the heatmap values. Possible values are “expression” (default), “methylation”
# scale - Use z-score to make the heamat? If we want to show differences between genes, it is good to make Z-score by samples (force each sample to have zero mean and standard deviation=1). If we want to show differences between samples, it is good to make Z-score by genes (force each gene to have zero mean and standard deviation=1). Possibilities: “row”, “col. Default”none"

# TCGAVisualize_volcano(): Create volcano plot
# Creates a volcano plot for DNA methylation or expression
# Please view page for arguments of this function

# TCGAvisualize_PCA: Principal Component Analysis plot for differentially expressed genes
# In order to better understand our genes, we can perform a PCA to reduce the number of dimensions of our gene set. 
# The function TCGAvisualize_PCA will plot the PCA for different groups.

# The arguments of this function are:
  
# dataFilt - The expression matrix after normalization and quantile filter
# dataDEGsFiltLevel - The TCGAanalyze_LevelTab output
# ntopgenes - Number of DEGs genes to plot in PCA
# ** group1 a string containing the barcode list of the samples in control group
# ** group2 a string containing the barcode list of the samples in disease group

# normalization of genes
dataNorm <- TCGAbiolinks::TCGAanalyze_Normalization(dataBRCA, geneInfo)

# quantile filter of genes
dataFilt <- TCGAanalyze_Filtering(tabDF = dataNorm,
                                  method = "quantile", 
                                  qnt.cut =  0.25)

# selection of normal samples "NT" 
group1 <- TCGAquery_SampleTypes(colnames(dataFilt), typesample = c("NT"))
# selection of normal samples "TP" 
group2 <- TCGAquery_SampleTypes(colnames(dataFilt), typesample = c("TP"))

# Principal Component Analysis plot for ntop selected DEGs
pca <- TCGAvisualize_PCA(dataFilt,dataDEGsFiltLevel, ntopgenes = 200, group1, group2)

# TCGAvisualize_meanMethylation: Mean DNA Methylation Analysis
# Using the data and calculating the mean DNA methylation per group, it is possible to create a mean DNA methylation boxplot with the function TCGAvisualize_meanMethylation as follows:
  
  query <- GDCquery(project = "TCGA-GBM",
                    data.category = "DNA methylation",
                    platform = "Illumina Human Methylation 27",
                    legacy = TRUE, 
                    barcode = c("TCGA-02-0058-01A-01D-0186-05", "TCGA-12-1597-01B-01D-0915-05",
                                "TCGA-12-0829-01A-01D-0392-05", "TCGA-06-0155-01B-01D-0521-05",
                                "TCGA-02-0099-01A-01D-0199-05", "TCGA-19-4068-01A-01D-1228-05",
                                "TCGA-19-1788-01A-01D-0595-05", "TCGA-16-0848-01A-01D-0392-05"))
GDCdownload(query, method = "api")
data <- GDCprepare(query)

# TCGAvisualize_starburst: Integration of gene expression and DNA methylation data
# The starburst plot is proposed to combine information from two volcano plots, and is applied for a study of DNA methylation and gene expression. It first introduced in 2010 (Noushmehr et al. 2010). In order to reproduce this plot, we will use the TCGAvisualize_starburst function.

# The function creates Starburst plot for comparison of DNA methylation and gene expression. 
# The log10 (FDR-corrected P value) for DNA methylation is plotted in the x axis, and for gene expression in the y axis, for each gene. 
# The black dashed line shows the FDR-adjusted P value of 0.01.

# Please view page to see arguments for this function

starburst <- TCGAvisualize_starburst(coad.SummarizeExperiment, 
                                     different.experssion.analysis.data,
                                     group1 = "CIMP.H",
                                     group2 = "CIMP.L",
                                     met.platform = "450K",
                                     genome = "hg19",
                                     met.p.cut = 10^-5, 
                                     exp.p.cut = 10^-5,
                                     names = TRUE)




