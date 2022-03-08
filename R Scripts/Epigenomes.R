# TCGA Workflow: Analyze cancer genomics and epigenomics data using Bioconductor packages

if (!"BiocManager" %in% rownames(installed.packages()))
  install.packages("BiocManager")
BiocManager::install("TCGAWorkflow")
BiocManager::install("TCGAWorkflowData")

library(TCGAWorkflowData)
library(DT)
library(TCGAbiolinks)
library(SummarizedExperiment)

# Obs: The data in the legacy database has been aligned to hg19

query.exp1 <- GDCquery(
  project = "TCGA-LGG",
  data.category = "Transcriptome Profiling",
  legacy = FALSE,
  data.type = "Gene Expression Quantification", 
  workflow.type = "HTSeq - FPKM-UQ"
)

GDCdownload(query.exp1)
exp.lgg <- GDCprepare(query = query.exp1, save = TRUE, save.filename = "lggExp.rda")

data(LGGIllumina_HiSeq)
data <- assay(lgg.exp)
datatable(
  data = data[1:10,],
  options = list(scrollX = TRUE, keys = TRUE, pageLength = 5),
  rownames = TRUE
)

genes.info <- rowRanges(lgg.exp)
genes.info

sample.info <- colData(lgg.exp)
datatable(
  data = as.data.frame(sample.info),
  options = list(scrollX = TRUE, keys = TRUE, pageLength = 5),
  rownames = FALSE
)

dataClin1 <- GDCquery_clinic(project = "TCGA-LGG", type = "Clinical")
datatable(dataClin[1:10,], options = list(scrollX = TRUE, keys = TRUE), rownames = FALSE)



LGGmut <- GDCquery_Maf(tumor = "LGG", pipelines = "mutect2")
data(mafMutect2LGG)
datatable(LGGmut, options = list(scrollX = TRUE, keys = TRUE), rownames = FALSE)


lgg.subtypes <- TCGAquery_subtype(tumor = "lgg")
datatable(lgg.subtypes[1:10,],options = list(scrollX = TRUE, keys = TRUE), rownames = FALSE)


library(maftools)

LGGmut <- GDCquery_Maf(tumor = "LGG", pipelines = "mutect2")
data(mafMutect2LGG)

LGGmut <- GDCquery_Maf(tumor = "LGG", pipelines = "mutect2", )
data(mafMutect2LGG)

colnames(dataClin1)[1] <- "Tumor_Sample_Barcode"

plyr::count(dataClin1$vital_status)n1$vital_status)TRUE
dataClin1$Overall_Survival_Status <- 1 #dead
dataClin1$Overall_Survival_Status[which(dataClin1$vital_status != "Dead")] <- 0

dataClin1$time <- dataClin1$days_to_death
dataClin1$time[is.na(dataClin1$days_to_death)] <- dataClin1$days_to_last_follow_up[is.na(dataClin1$days_to_death)]

maf <- read.maf(maf = LGGmut, clinicalData = dataClin1, isTCGA = TRUE)

plotmafSummary(
  maf = maf,
  rmOutlier = TRUE,
  addStat = 'median',
  dashboard = TRUE
)


oncoplot(
  maf = maf,
  top = 20,
  legendFontSize = 8,
  clinicalFeatures = c("tissue_or_organ_of_origin")
)

plot <- mafSurvival(
  maf = maf,
  genes = "TP53",
  time = 'time',
  Status = 'Overall_Survival_Status',
  isTCGA = TRUE
)



# Transcriptomic Analysis

query.exp1$results[[1]] <- query.exp1$results[[1]][1:20,]
GDCdownload(query.exp1)

data("LGGIllumina_HiSeq")

dataPrep_LGG <- TCGAanalyze_Preprocessing(
  object = lgg.exp,
  cor.cut = 0.6,
  datatype = "raw_count",
  filename = "LGG_IlluminaHiSeq_RNASeqV2.png"
)


dataNorm <- TCGAanalyze_Normalization(
  tabDF = dataPrep_LGG, 
  geneInfo = TCGAbiolinks::geneInfo,
  method = "gcContent"
)

dataFilt <- TCGAanalyze_Filtering(
  tabDF = dataNorm,
  method = "quantile",
  qnt.cut = 0.25
)

save(dataFilt, file = paste0("LGG_Norm_IlluminaHISeq.rda"))

dataFiltLGG <- subset(
  dataFilt,
  select = substr(colnames(dataFilt), 1, 12) %in% dataClin1$bcr_patient_barcode
)

samplesNT <- TCGAquery_SampleTypes(barcode = colnames(dataFilt),
                                   typesample = c("NT"))

samplesTP <- TCGAquery_SampleTypes(barcode = colnames(dataFilt),
                                   typesample = c("TP"))


#This line is not working below
dataDEGs <- TCGAanalyze_DEA(mat1 = dataFilt[,samplesNT],
                            mat2 = dataFilt[,samplesTP],
                            Cond1type = "Normal",
                            Cond2type = "Tumor",
                            fdr.cut = 0.01 ,
                            logFC.cut = 1,
                            method = "glmLRT")


dataNorm <- TCGAanalyze_Normalization(tabDF = dataPrep_LGG, geneInfo =  geneInfo)

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



















# Enrichment Analysis
ansEA <- TCGAanalyze_EAcomplete(
  TFname = "DEA genes LGG",
  RegulonList = rownames(dataDEGs)
)

TCGAvisualize_EAbarplot(
  tf = rownames(ansEA$ResBP),
  filename = NULL,
  GOBPTab = ansEA$ResBP,
  nRGTab = rownames(dataDEGs),
  nBar = 20
)


TCGAvisualize_EAbarplot(
  tf = rownames(ansEA$ResBP),
  filename = NULL,
  GOCCTab = ansEA$ResCC,
  nRGTab = rownames(dataDEGs),
  nBar = 20
)

TCGAvisualize_EAbarplot(
  tf = rownames(ansEA$ResBP),
  filename = NULL,
  GOMFTab = ansEA$ResMF,
  nRGTab = rownames(dataDEGs),
  nBar = 20
)

TCGAvisualize_EAbarplot(
  tf = rownames(ansEA$ResBP),
  filename = NULL,
  GOMFTab = ansEA$ResPat,
  nRGTab = rownames(dataDEGs),
  nBar = 20
)

