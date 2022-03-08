## ---------------------------
##
##
##        _ _____  _____  _    _ 
##       | |  __ \|  __ \| |  | |
##       | | |__) | |__) | |__| |
##   _   | |  ___/|  _  /|  __  |
##  | |__| | |    | | \ \| |  | |
##   \____/|_|    |_|  \_\_|  |_|
##
##
## Script name: TCGAbiolinksTEST.R
##
## Purpose of script:
## To practice with biolinks
##
## Author: James Hacking
##
## Date Created: 02-11-2020 
##
## Copyright (c) James Hacking, 2020
## Email: james.hacking@ncl.ac.uk
##
## ---------------------------
##
## Notes: TESTING GIT
##
## Case studies are using TCGAbiolinks_2.18.0
## Compsvr currently running TCGAbiolinks_2.10.5
##
## Add in this: 
##
## ---------------------------

# source("functions/packages.R")       # loads up all the packages we need

## ---------------------------

## load up our functions into memory

# source("functions/summarise_data.R") 

## ---------------------------

### Just some code that shows how to install if not there.
#if(!require(somepackage)){
#  install.packages("somepackage")
#  library(somepackage)
#}



# http://bioconductor.org/packages/release/bioc/vignettes/TCGAbiolinks/inst/doc/index.html

#Introduction

if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("TCGAbiolinks")

if(!require("DT")){
  install.packages("DT")
  library("DT")}
 

library(TCGAbiolinks)
library(dplyr)
library(DT)

version

packageVersion("TCGAbiolinks")

#Case study Pan Cancer downstream analysis LGG

BiocManager::install("EDASeq")
BiocManager::install("genefilter")
BiocManager::install("ConsensusClusterPlus")
BiocManager::install("survminer")
BiocManager::install("ComplexHeatmap")

library(TCGAbiolinks)
library(SummarizedExperiment)


query.exp <- GDCquery(project = "TCGA-LGG", 
                      legacy = TRUE,
                      data.category = "Gene expression",
                      data.type = "Gene expression quantification",
                      platform = "Illumina HiSeq", 
                      file.type = "results",
                      experimental.strategy = "RNA-Seq",
                      sample.type = "Primary Tumor")
GDCdownload(query.exp)
lgg.exp <- GDCprepare(query = query.exp, save = TRUE, save.filename = "lggExp.rda")

library(dplyr)

# get subtype information
dataSubt <- TCGAquery_subtype(tumor = "BRCA")

# Which samples are Primary Tumor
dataSmTP <- TCGAquery_SampleTypes(getResults(query.exp,cols="cases"),"TP") 

# which samples are solid tissue normal
dataSmNT <- TCGAquery_SampleTypes(getResults(query.exp,cols="cases"),"NT")

dataPrep <- TCGAanalyze_Preprocessing(object = lgg.exp, cor.cut = 0.6)
dataNorm <- TCGAanalyze_Normalization(tabDF = dataPrep,
                                      geneInfo = geneInfo,
                                      method = "gcContent")

datFilt <- dataNorm %>% TCGAanalyze_Filtering(method = "varFilter") %>%
  TCGAanalyze_Filtering(method = "filter1") %>%  TCGAanalyze_Filtering(method = "filter2",foldChange = 0.2)

data_Hc2 <- TCGAanalyze_Clustering(tabDF = datFilt,
                                   method = "consensus",
                                   methodHC = "ward.D2") 
# Add  cluster information to Summarized Experiment
colData(lgg.exp)$groupsHC <- paste0("EC",data_Hc2[[4]]$consensusClass)

#The next steps will be to visualize the data. First, we created the survival plot.

TCGAanalyze_survival(data = colData(lgg.exp),
                     clusterCol = "groupsHC",
                     main = "TCGA kaplan meier survival plot from consensus cluster",
                     legend = "RNA Group",height = 10,
                     risk.table = T,conf.int = F,
                     color = c("black","red","blue","green3"),
                     filename = "survival_lgg_expression_subtypes.png")

#We will also, create a heatmap of the expression.

TCGAvisualize_Heatmap(t(datFilt),
                      col.metadata =  colData(lgg.exp)[,c("barcode",
                                                          "groupsHC",
                                                          "subtype_Histology",
                                                          "subtype_IDH.codel.subtype")],
                      col.colors =  list(
                        groupsHC = c("EC1"="black",
                                     "EC2"="red",
                                     "EC3"="blue",
                                     "EC4"="green3")),
                      sortCol = "groupsHC",
                      type = "expression", # sets default color
                      scale = "row", # use z-scores for better visualization. Center gene expression level around 0.
                      title = "Heatmap from concensus cluster", 
                      filename = "case2_Heatmap.png",
                      cluster_rows = TRUE,
                      color.levels = colorRampPalette(c("green", "black", "red"))(n = 11),
                      exatrems =seq(-5,5,1),
                      cluster_columns = FALSE,
                      width = 1000,
                      height = 1000)

sessionInfo()

lgg.exp$