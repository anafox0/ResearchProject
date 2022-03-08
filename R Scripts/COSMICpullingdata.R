
# BRCA --------------------------------------------------------------------

head(cosmicGS)

library(dplyr)

### use that to explore the data
head(cosmicGS)
table(cosmicGS$Primary.histology)

### create a subset with only the types that you want. Check whether you need top add oligodendroglioma or others add to below?
cosmicGS %>%
   filter(Primary.histology=="carcinoma") %>%
   filter(Histology.subtype.1=="ductal_carcinoma_in_situ"|Histology.subtype.1=="lobular_carcinoma") -> cosmicGS.brca



#saveRDS(cosmicGS.brca, file = "./CosmicGSData/cosmicGS.brca.rds")
# cosmicGS.brca <- readRDS(file = "./CosmicGSData/cosmicGS.brca.rds")

# ALternative
# cosmicGS[cosmicGS$Primary.histology=="glioma",]

table(cosmicGS.brca$Histology.subtype.1)
head(cosmicGS.brca)
dim(cosmicGS.brca)


table(cosmicGS.brca$Gene.name)

#### only looking for rows which contain letters AASDH in the Gene.name column
cosmicGS.brca[grep("TP53",cosmicGS.brca$Gene.name),]

# count number of genes with unique ID_tumour, and unique  GENOMIC_MUTATION_ID

### trim everything after the underscore
gsub('_.*', '',cosmicGS.brca$Gene.name) -> gene.names.trimmed
gsub('.{2}$', '',cosmicGS.brca$Accession.Number) -> AN.trimmed1
gsub('ENST', 'ENSG', AN.trimmed1) -> AN.trimmed

### add an extra column to the subsetted dataframe
cosmicGS.brca$gene.names.trimmed <- gene.names.trimmed
cosmicGS.brca$AN.trimmed <- AN.trimmed

### creates a concatenated string of gene name, patient name, mutation name
paste(cosmicGS.brca$gene.names.trimmed, cosmicGS.brca$ID_tumour, cosmicGS.brca$GENOMIC_MUTATION_ID) -> barcode
#### add it as a column back to the table
cosmicGS.brca$barcode <- barcode
### remove any rows which have a duplicate barcode
cosmicGS.brca[!duplicated(cosmicGS.brca$barcode),] -> cosmicGS.brca.dedup

# Create a subsetted dataframe from this deduped dataframe containing mutated genes
# New dataframe will contain only those genes predicted to be pathogenic (FATHMM.prediction)
# Convert character values in FATHMM.prediction column to factors
library(purrr)
brca.cosmic.dedup <- as.data.frame(map_if(cosmicGS.brca.dedup, is.character, as.factor))
summary(brca.cosmic.dedup$FATHMM.prediction)
# 1891 = PATHOGNEIC, NEUTRAL = 41370, BLANK = 3489
# Subset data to extract only those that are predicted to be pathogenic
brca.cosmic.path <- brca.cosmic.dedup[brca.cosmic.dedup$FATHMM.prediction == "PATHOGENIC", ]
dim(brca.cosmic.path)
#Should be the same dimensions are dim(summary(brca.cosmic.dedup$FATHMM.prediction)

#Now extracting the columns that we are interested in for our annotations

brca.cosmic <- brca.cosmic.path[,c(13:14,20:22,28:31,39,41)]
dim(brca.cosmic)
head(brca.cosmic)
#Should be 10 columns

#From here, we will make
brca.cosmic[brca.cosmic$gene.names.trimmed=="TP53",]
summary(brca.cosmic$gene.names.trimmed)
head(brca.cosmic[brca.cosmic$gene.names.trimmed=="TP53",])
aggregate(FATHMM.score~gene.names.trimmed, data = brca.cosmic, mean) -> Mean.FATHMM.score
aggregate(FATHMM.score~gene.names.trimmed, data = brca.cosmic, length) -> num.muts
##Format Mean.FATHMM.score as a data.frame and change column names to allow input into master annotation file
as.data.frame(Mean.FATHMM.score)
data.table::setnames(Mean.FATHMM.score,'FATHMM.score','mean.FATHMM.score')

### count the number of times the genes are found in the column gene names trimmed
table(brca.cosmic$gene.names.trimmed) -> freq.gene.mut
### order them in decreasing order
freq.gene.mut[order(freq.gene.mut, decreasing = T)] -> ordered.freq.gene.mut
plot(ordered.freq.gene.mut)

#Format freq.gene.mut as a data.frame and change column names to allow input into master annotation file
data.table::as.data.table(freq.gene.mut) -> freq.gene.mut
data.table::setnames(freq.gene.mut,'V1','gene.names.trimmed')
data.table::setnames(freq.gene.mut,'N','Mutation Frequency')


#Use the full_join() to join the three databases together

brca.cosmic %>% full_join(Mean.FATHMM.score, by = "gene.names.trimmed") -> brca.cosmic1
brca.cosmic1 %>% full_join(freq.gene.mut, by = "gene.names.trimmed") -> brca.cosmic

saveRDS(brca.cosmic, file = "./CosmicGSData/brca.cosmic.rds")

head(cosmicGS.brca)


write.table(brca.cosmic, file = "./cosmicdata/cosmic.brca.txt", sep = "\t", quote = F, row.names = F)
write.csv(brca.cosmic, file = "./cosmicdata/cosmic.brca.csv", row.names = F)



# GBM ---------------------------------------------------------------------

head(cosmicGS)

library(dplyr)

### use that to explore the data
head(cosmicGS)
table(cosmicGS$Primary.histology)

### create a subset with only the types that you want. Check whether you need top add oligodendroglioma or others add to below?
#cosmicGS %>%
 # filter(Primary.histology=="glioma") %>%
  #filter(Histology.subtype.1=="astrocytoma_Grade_IV")  -> cosmicGS.gbm

#saveRDS(cosmicGS.gbm, file = "./CosmicGSData/cosmicGS.gbm.rds")
# cosmicGS.gbm <- readRDS(file = "./CosmicGSData/cosmicGS.gbm.rds")

# ALternative
# cosmicGS[cosmicGS$Primary.histology=="glioma",]

table(cosmicGS.gbm$Histology.subtype.1)
head(cosmicGS.gbm)
dim(cosmicGS.gbm)


table(cosmicGS.gbm$Gene.name)

#### only looking for rows which contain letters AASDH in the Gene.name column
cosmicGS.gbm[grep("TP53",cosmicGS.gbm$Gene.name),]

# count number of genes with unique ID_tumour, and unique  GENOMIC_MUTATION_ID

### trim everything after the underscore
gsub('_.*', '',cosmicGS.gbm$Gene.name) -> gene.names.trimmed
gsub('\\..*', '', cosmicGS.gbm$Accession.Number) -> AN.trimmed1
gsub('ENST', 'ENSG', AN.trimmed1) -> AN.trimmed

### add an extra column to the subsetted dataframe
cosmicGS.gbm$gene.names.trimmed <- gene.names.trimmed
cosmicGS.gbm$AN.trimmed <- AN.trimmed
head(cosmicGS.gbm)
### creates a concatenated string of gene name, patient name, mutation name
paste(cosmicGS.gbm$gene.names.trimmed, cosmicGS.gbm$ID_tumour, cosmicGS.gbm$GENOMIC_MUTATION_ID) -> barcode
#### add it as a column back to the table
cosmicGS.gbm$barcode <- barcode
### remove any rows which have a duplicate barcode
cosmicGS.gbm[!duplicated(cosmicGS.gbm$barcode),] -> cosmicGS.gbm.dedup

# Create a subsetted dataframe from this deduped dataframe containing mutated genes
# New dataframe will contain only those genes predicted to be pathogenic (FATHMM.prediction)
# Convert character values in FATHMM.prediction column to factors
library(purrr)
gbm.cosmic.dedup <- as.data.frame(map_if(cosmicGS.gbm.dedup, is.character, as.factor))
summary(gbm.cosmic.dedup$FATHMM.prediction)
# 62986 = PATHOGNEIC, NEUTRAL = 67349, BLANK = 33507
# Subset data to extract only those that are predicted to be pathogenic
gbm.cosmic.path <- gbm.cosmic.dedup[gbm.cosmic.dedup$FATHMM.prediction == "PATHOGENIC", ]
dim(gbm.cosmic.path)
#Should be the same dimensions are dim(summary(gbm.cosmic.dedup$FATHMM.prediction)

#Now extracting the columns that we are interested in for our annotations
head(gbm.cosmic.path)
gbm.cosmic <- gbm.cosmic.path[,c(13:14,20:22,28:31,39,41)]
dim(gbm.cosmic)
head(gbm.cosmic)
#Should be 11 columns

#From here, we will make
gbm.cosmic[gbm.cosmic$gene.names.trimmed=="TP53",]
summary(gbm.cosmic$gene.names.trimmed)
head(gbm.cosmic[gbm.cosmic$gene.names.trimmed=="TP53",])
aggregate(FATHMM.score~gene.names.trimmed, data = gbm.cosmic, mean) -> Mean.FATHMM.score
aggregate(FATHMM.score~gene.names.trimmed, data = gbm.cosmic, length) -> num.muts
##Format Mean.FATHMM.score as a data.frame and change column names to allow input into master annotation file
as.data.frame(Mean.FATHMM.score)
data.table::setnames(Mean.FATHMM.score,'FATHMM.score','mean.FATHMM.score')

### count the number of times the genes are found in the column gene names trimmed
table(gbm.cosmic$gene.names.trimmed) -> freq.gene.mut
### order them in decreasing order
freq.gene.mut[order(freq.gene.mut, decreasing = T)] -> ordered.freq.gene.mut
plot(ordered.freq.gene.mut)

#Format freq.gene.mut as a data.frame and change column names to allow input into master annotation file
data.table::as.data.table(freq.gene.mut) -> freq.gene.mut
data.table::setnames(freq.gene.mut,'V1','gene.names.trimmed')
data.table::setnames(freq.gene.mut,'N','Mutation Frequency') -> freq.gene.mut
#Use the full_join() to join the three databases together

gbm.cosmic %>% full_join(Mean.FATHMM.score, by = "gene.names.trimmed") -> gbm.cosmic1
gbm.cosmic1 %>% full_join(freq.gene.mut, by = "gene.names.trimmed") -> gbm.cosmic





plot(ordered.freq.gene.mut)

head(cosmicGS.gbm)

#saveRDS(cosmicGS.gbm, file = "./CosmicGSData/cosmicGS.gbm.dedup.rds")
# cosmicGS.gbm <- readRDS(file = "./CosmicGSData/cosmicGS.gbm.dedup.rds")
write.table(gbm.cosmic, file = "./cosmicdata/cosmic.gbm.txt", sep = "\t", quote = F, row.names = F)
write.csv(gbm.cosmic, file = "./cosmicdata/cosmic.gbm.csv", row.names = F)
# LAML ---------------------------------------------------------------------

head(cosmicGS)

library(dplyr)

### use that to explore the data
head(cosmicGS)
table(cosmicGS$Primary.histology)

### create a subset with only the types that you want. Check whether you need top add oligodendroglioma or others add to below?
#cosmicGS %>%
 # filter(Primary.histology=="haematopoietic_neoplasm") %>%
  #filter(Histology.subtype.1=="acute_myeloid_leukaemia")  -> cosmicGS.haemoneo

#saveRDS(cosmicGS.haemoneo, file = "./CosmicGSData/cosmicGS.laml.rds")
# cosmicGS.laml <- readRDS(file = "./CosmicGSData/cosmicGS.laml.rds")

# ALternative
# cosmicGS[cosmicGS$Primary.histology=="glioma",]

table(cosmicGS.laml$Histology.subtype.1)
head(cosmicGS.laml)
dim(cosmicGS.laml)


table(cosmicGS.laml$Gene.name)

#### only looking for rows which contain letters AASDH in the Gene.name column
cosmicGS.laml[grep("TP53",cosmicGS.laml$Gene.name),]

# count number of genes with unique ID_tumour, and unique  GENOMIC_MUTATION_ID

### trim everything after the underscore
gsub('_.*', '',cosmicGS.laml$Gene.name) -> gene.names.trimmed
gsub('\\..*', '', cosmicGS.laml$Accession.Number) -> AN.trimmed1
gsub('ENST', 'ENSG', AN.trimmed1) -> AN.trimmed

### add an extra column to the subsetted dataframe
cosmicGS.laml$gene.names.trimmed <- gene.names.trimmed
cosmicGS.laml$AN.trimmed <- AN.trimmed
### creates a concatenated string of gene name, patient name, mutation name
paste(cosmicGS.laml$gene.names.trimmed, cosmicGS.laml$ID_tumour, cosmicGS.laml$GENOMIC_MUTATION_ID) -> barcode
#### add it as a column back to the table
cosmicGS.laml$barcode <- barcode
### remove any rows which have a duplicate barcode
cosmicGS.laml[!duplicated(cosmicGS.laml$barcode),] -> cosmicGS.laml.dedup

# Create a subsetted dataframe from this deduped dataframe containing mutated genes
# New dataframe will contain only those genes predicted to be pathogenic (FATHMM.prediction)
# Convert character values in FATHMM.prediction column to factors
library(purrr)
laml.cosmic.dedup <- as.data.frame(map_if(cosmicGS.laml.dedup, is.character, as.factor))
summary(laml.cosmic.dedup$FATHMM.prediction)
# 23581 = PATHOGNEIC, NEUTRAL = 175704, BLANK = 6938
# Subset data to extract only those that are predicted to be pathogenic
laml.cosmic.path <- laml.cosmic.dedup[laml.cosmic.dedup$FATHMM.prediction == "PATHOGENIC", ]
dim(laml.cosmic.path)
#Should be the same dimensions are dim(summary(laml.cosmic.dedup$FATHMM.prediction)

#Now extracting the columns that we are interested in for our annotations

laml.cosmic <- laml.cosmic.path[,c(13:14,20:22,28:31,39,41)]
dim(laml.cosmic)
head(laml.cosmic)
#Should be 10 columns

#From here, we will make
laml.cosmic[laml.cosmic$gene.names.trimmed=="TP53",]
summary(laml.cosmic$gene.names.trimmed)
head(laml.cosmic[laml.cosmic$gene.names.trimmed=="TP53",])
aggregate(FATHMM.score~gene.names.trimmed, data = laml.cosmic, mean) -> Mean.FATHMM.score
aggregate(FATHMM.score~gene.names.trimmed, data = laml.cosmic, length) -> num.muts
##Format Mean.FATHMM.score as a data.frame and change column names to allow input into master annotation file
as.data.frame(Mean.FATHMM.score)
data.table::setnames(Mean.FATHMM.score,'FATHMM.score','mean.FATHMM.score')

### count the number of times the genes are found in the column gene names trimmed
table(laml.cosmic$gene.names.trimmed) -> freq.gene.mut
### order them in decreasing order
freq.gene.mut[order(freq.gene.mut, decreasing = T)] -> ordered.freq.gene.mut
plot(ordered.freq.gene.mut)

#Format freq.gene.mut as a data.frame and change column names to allow input into master annotation file
data.table::as.data.table(freq.gene.mut) -> freq.gene.mut
data.table::setnames(freq.gene.mut,'V1','gene.names.trimmed')
data.table::setnames(freq.gene.mut,'N','Mutation Frequency')


#Use the full_join() to join the three databases together

laml.cosmic %>% full_join(Mean.FATHMM.score, by = "gene.names.trimmed") -> laml.cosmic1
laml.cosmic1 %>% full_join(freq.gene.mut, by = "gene.names.trimmed") -> laml.cosmic

saveRDS(laml.cosmic, file = "./CosmicGSData/laml.cosmic.rds")
write.csv(laml.cosmic, file = "./cosmicdata/cosmic.laml.csv", row.names = F)



# KIRC ---------------------------------------------------------------------

head(cosmicGS)

library(dplyr)

### use that to explore the data
head(cosmicGS)
table(cosmicGS$Primary.histology)

### create a subset with only the types that you want. Check whether you need top add oligodendroglioma or others add to below?
 cosmicGS %>%
   filter(Primary.histology=="carcinoma") %>%
   filter(Histology.subtype.1=="clear_cell_renal_cell_carcinoma")  -> cosmicGS.kirc

#saveRDS(cosmicGS.kirc, file = "./CosmicGSData/cosmicGS.kirc.rds")
#cosmicGS.kirc <- readRDS(file = "./CosmicGSData/cosmicGS.kirc.rds")

# ALternative
# cosmicGS[cosmicGS$Primary.histology=="glioma",]

table(cosmicGS.kirc$Histology.subtype.1)
head(cosmicGS.kirc)
dim(cosmicGS.kirc)


table(cosmicGS.kirc$Gene.name)

#### only looking for rows which contain letters AASDH in the Gene.name column
cosmicGS.kirc[grep("TP53",cosmicGS.kirc$Gene.name),]

# count number of genes with unique ID_tumour, and unique  GENOMIC_MUTATION_ID

### trim everything after the underscore
gsub('_.*', '',cosmicGS.kirc$Gene.name) -> gene.names.trimmed
gsub('.{2}$', '',cosmicGS.kirc$Accession.Number) -> AN.trimmed1
gsub('ENST', 'ENSG', AN.trimmed1) -> AN.trimmed

### add an extra column to the subsetted dataframe
cosmicGS.kirc$gene.names.trimmed <- gene.names.trimmed
cosmicGS.kirc$AN.trimmed <- AN.trimmed
### creates a concatenated string of gene name, patient name, mutation name
paste(cosmicGS.kirc$gene.names.trimmed, cosmicGS.kirc$ID_tumour, cosmicGS.kirc$GENOMIC_MUTATION_ID) -> barcode
#### add it as a column back to the table
cosmicGS.kirc$barcode <- barcode
### remove any rows which have a duplicate barcode
cosmicGS.kirc[!duplicated(cosmicGS.kirc$barcode),] -> cosmicGS.kirc.dedup

# Create a subsetted dataframe from this deduped dataframe containing mutated genes
# New dataframe will contain only those genes predicted to be pathogenic (FATHMM.prediction)
# Convert character values in FATHMM.prediction column to factors
library(purrr)
kirc.cosmic.dedup <- as.data.frame(map_if(cosmicGS.kirc.dedup, is.character, as.factor))
summary(kirc.cosmic.dedup$FATHMM.prediction)
# 43687 = PATHOGNEIC, NEUTRAL = 49200, BLANK = 11061
# Subset data to extract only those that are predicted to be pathogenic
kirc.cosmic.path <- kirc.cosmic.dedup[kirc.cosmic.dedup$FATHMM.prediction == "PATHOGENIC", ]
dim(kirc.cosmic.path)
#Should be the same dimensions are dim(summary(kirc.cosmic.dedup$FATHMM.prediction)

#Now extracting the columns that we are interested in for our annotations

kirc.cosmic <- kirc.cosmic.path[,c(13:14,20:22,28:31,39,41)]
dim(kirc.cosmic)
head(kirc.cosmic)
#Should be 10 columns

#From here, we will make
kirc.cosmic[kirc.cosmic$gene.names.trimmed=="TP53",]
summary(kirc.cosmic$gene.names.trimmed)
head(kirc.cosmic[kirc.cosmic$gene.names.trimmed=="TP53",])
aggregate(FATHMM.score~gene.names.trimmed, data = kirc.cosmic, mean) -> Mean.FATHMM.score
aggregate(FATHMM.score~gene.names.trimmed, data = kirc.cosmic, length) -> num.muts
##Format Mean.FATHMM.score as a data.frame and change column names to allow input into master annotation file
as.data.frame(Mean.FATHMM.score)
data.table::setnames(Mean.FATHMM.score,'FATHMM.score','mean.FATHMM.score')

### count the number of times the genes are found in the column gene names trimmed
table(kirc.cosmic$gene.names.trimmed) -> freq.gene.mut
### order them in decreasing order
freq.gene.mut[order(freq.gene.mut, decreasing = T)] -> ordered.freq.gene.mut
plot(ordered.freq.gene.mut)

#Format freq.gene.mut as a data.frame and change column names to allow input into master annotation file
data.table::as.data.table(freq.gene.mut) -> freq.gene.mut
data.table::setnames(freq.gene.mut,'V1','gene.names.trimmed')
data.table::setnames(freq.gene.mut,'N','Mutation Frequency')


#Use the full_join() to join the three databases together

kirc.cosmic %>% full_join(Mean.FATHMM.score, by = "gene.names.trimmed") -> kirc.cosmic1
kirc.cosmic1 %>% full_join(freq.gene.mut, by = "gene.names.trimmed") -> kirc.cosmic

saveRDS(kirc.cosmic, file = "./CosmicGSData/kirc.cosmic.rds")
write.table(kirc.cosmic, file = "./cosmicdata/cosmic.kirc.txt", sep = "\t", quote = F, row.names = F)
write.csv(kirc.cosmic, file = "./cosmicdata/cosmic.kirc.csv", row.names = F)

# LGG ---------------------------------------------------------------------

head(cosmicGS)

library(dplyr)

### use that to explore the data
head(cosmicGS)
table(cosmicGS$Primary.histology)

### create a subset with only the types that you want. Check whether you need top add oligodendroglioma or others add to below?
cosmicGS %>%
  filter(Primary.histology=="glioma") %>%
  filter(Histology.subtype.1=="astrocytoma_Grade_II"|Histology.subtype.1=="astrocytoma_Grade_III"|Histology.subtype.1=="oligoastrocytoma_Grade_II"|Histology.subtype.1=="oligoastrocytoma_Grade_III"|Histology.subtype.1=="oligodendroglioma_Grade_II"|Histology.subtype.1=="oligodendroglioma_Grade_III") -> cosmicGS.lgg

#saveRDS(cosmicGS.lgg, file = "./CosmicGSData/cosmicGS.lgg.rds")
cosmicGS.lgg <- readRDS(file = "./CosmicGSData/cosmicGS.lgg.rds")

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
gsub('.{2}$', '',cosmicGS.lgg$Accession.Number) -> AN.trimmed1
gsub('ENST', 'ENSG', AN.trimmed1) -> AN.trimmed

### add an extra column to the subsetted dataframe
cosmicGS.lgg$gene.names.trimmed <- gene.names.trimmed
cosmicGS.lgg$AN.trimmed <- AN.trimmed
### creates a concatenated string of gene name, patient name, mutation name
paste(cosmicGS.lgg$gene.names.trimmed, cosmicGS.lgg$ID_tumour, cosmicGS.lgg$GENOMIC_MUTATION_ID) -> barcode
#### add it as a column back to the table
cosmicGS.lgg$barcode <- barcode
### remove any rows which have a duplicate barcode
cosmicGS.lgg[!duplicated(cosmicGS.lgg$barcode),] -> cosmicGS.lgg.dedup

# Create a subsetted dataframe from this deduped dataframe containing mutated genes
# New dataframe will contain only those genes predicted to be pathogenic (FATHMM.prediction)
# Convert character values in FATHMM.prediction column to factors
library(purrr)
lgg.cosmic.dedup <- as.data.frame(map_if(cosmicGS.lgg.dedup, is.character, as.factor))
summary(lgg.cosmic.dedup$FATHMM.prediction)
# 6628 = PATHOGNEIC, NEUTRAL = 5903, BLANK = 1918
# Subset data to extract only those that are predicted to be pathogenic
lgg.cosmic.path <- lgg.cosmic.dedup[lgg.cosmic.dedup$FATHMM.prediction == "PATHOGENIC", ]
dim(lgg.cosmic.path)
#Should be the same dimensions are dim(summary(lgg.cosmic.dedup$FATHMM.prediction)

#Now extracting the columns that we are interested in for our annotations

lgg.cosmic <- lgg.cosmic.path[,c(13:14,20:22,28:31,39,41)]
dim(lgg.cosmic)
head(lgg.cosmic)
#Should be 10 columns

#From here, we will make
lgg.cosmic[lgg.cosmic$gene.names.trimmed=="TP53",]
summary(lgg.cosmic$gene.names.trimmed)
head(lgg.cosmic[lgg.cosmic$gene.names.trimmed=="TP53",])
aggregate(FATHMM.score~gene.names.trimmed, data = lgg.cosmic, mean) -> Mean.FATHMM.score
aggregate(FATHMM.score~gene.names.trimmed, data = lgg.cosmic, length) -> num.muts
##Format Mean.FATHMM.score as a data.frame and change column names to allow input into master annotation file
as.data.frame(Mean.FATHMM.score)
data.table::setnames(Mean.FATHMM.score,'FATHMM.score','mean.FATHMM.score')

### count the number of times the genes are found in the column gene names trimmed
table(lgg.cosmic$gene.names.trimmed) -> freq.gene.mut
### order them in decreasing order
freq.gene.mut[order(freq.gene.mut, decreasing = T)] -> ordered.freq.gene.mut
plot(ordered.freq.gene.mut)

#Format freq.gene.mut as a data.frame and change column names to allow input into master annotation file
data.table::as.data.table(freq.gene.mut) -> freq.gene.mut
data.table::setnames(freq.gene.mut,'V1','gene.names.trimmed')
data.table::setnames(freq.gene.mut,'N','Mutation Frequency')


#Use the full_join() to join the three databases together

lgg.cosmic %>% full_join(Mean.FATHMM.score, by = "gene.names.trimmed") -> lgg.cosmic1
lgg.cosmic1 %>% full_join(freq.gene.mut, by = "gene.names.trimmed") -> lgg.cosmic

saveRDS(lgg.cosmic, file = "./CosmicGSData/lgg.cosmic.rds")
write.table(lgg.cosmic, file = "./cosmicdata/cosmic.lgg.txt", sep = "\t", quote = F, row.names = F)
write.csv(lgg.cosmic, file = "./cosmicdata/cosmic.lgg.csv", row.names = F)
