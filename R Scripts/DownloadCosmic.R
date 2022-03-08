### cosmic

### download your file and untar using treminal

### tar -xvzf cosmic.file.tar.gz
# read.delim(file = "./CosmicGSData/CosmicGenomeScreensMutantExport.tsv")

# cosmicGS <- read.delim(file = "./CosmicGSData/CosmicGenomeScreensMutantExport.tsv")
# saveRDS(cosmicGS, file = "./CosmicGSData/CosmicGenomeScreensMutantExport.rds")


#setwd("/home/anafox/ResearchProject/")
#cosmicGS <- readRDS(file = "./CosmicGSData/CosmicGenomeScreensMutantExport.rds")
head(cosmicGS)

library(dplyr)

### use that to explore the data
head(cosmicGS)
table(cosmicGS$Primary.histology)

### create a subset with only the types that you want. Check whether you need top add oligodendroglioma or others add to below?
cosmicGS %>%
filter(Primary.histology=="glioma") %>%
 filter(Histology.subtype.1=="astrocytoma_Grade_II"|Histology.subtype.1=="astrocytoma_Grade_III")   -> cosmicGS.lgg

#saveRDS(cosmicGS.lgg, file = "./CosmicGSData/cosmicGS.lgg.rds")
# cosmicGS.lgg <- readRDS(file = "./CosmicGSData/cosmicGS.lgg.rds")

# ALternative
# cosmicGS[cosmicGS$Primary.histology=="glioma",]

table(cosmicGS.glioma$Histology.subtype.1)
head(cosmicGS.lgg)
dim(cosmicGS.lgg)


table(cosmicGS.lgg$Gene.name)

#### only looking for rows which contain letters AASDH in the Gene.name column
cosmicGS.lgg[grep("TP53",cosmicGS.lgg$Gene.name),]

# count number of genes with unique ID_tumour, and unique  GENOMIC_MUTATION_ID

### trim everything after the underscore
gsub('_.*', '',cosmicGS.lgg$Gene.name) -> gene.names.trimmed
### add an extra column to the subsetted dataframe
cosmicGS.lgg$gene.names.trimmed <- gene.names.trimmed

### creates a concatenated string of gene name, patient name, mutation name
paste(cosmicGS.lgg$gene.names.trimmed, cosmicGS.lgg$ID_tumour, cosmicGS.lgg$GENOMIC_MUTATION_ID) -> barcode
#### add it as a column back to the table
cosmicGS.lgg$barcode <- barcode
### remove any rows which have a duplicate barcode
cosmicGS.lgg[!duplicated(cosmicGS.lgg$barcode),] -> cosmicGS.lgg.dedup

### count the number of times the genes are found in the column gene names trimmed
table(cosmicGS.lgg.dedup$gene.names.trimmed) -> freq.gene.mut
### order them in decreasing order
freq.gene.mut[order(freq.gene.mut, decreasing = T)] -> ordered.freq.gene.mut

plot(ordered.freq.gene.mut)

head(cosmicGS.lgg)

#saveRDS(cosmicGS.lgg, file = "./CosmicGSData/cosmicGS.lgg.dedup.rds")
# cosmicGS.lgg <- readRDS(file = "./CosmicGSData/cosmicGS.lgg.dedup.rds"))


#########################################

# save(cosmicGS, file = "./CosmicGSData/CosmicGenomeScreensMutantExport.Rdata")
# load(file = "./CosmicGSData/CosmicGenomeScreensMutantExport.Rdata")

R.utils::gunzip("CosmicGenomeScreensMutantExport.tsv")

list.files("./CosmicGSData/", full.names = T) -> files

list.files("./cosmic_data/", pattern = ".csv", full.names = T) -> csv.files

read.csv(csv.files[1])

lapply(csv.files, read.csv) -> csv.data
str(csv.data)
names(csv.data) <- basename(csv.files)

list.files("./CosmicGSData/", pattern = ".tsv", full.names = T) -> tsv.files
lapply(tsv.files, read.delim) -> tsv.data
names(tsv.data) <- basename(tsv.files)
list.files()


list.files("./cosmic_data/", pattern = ".vcf", full.names = T) -> vcf.files
library("VariantAnnotation")
readVcf(vcf.files[1])
lapply(vcf.files, readVcf) -> vcf.data
names(vcf.data) <- basename(vcf.files)

all.data <- c(csv.data,tsv.data, vcf.data)


