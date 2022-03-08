# LAML
# GBM
# KIRC
# LGG
# DLBC
source(file = "./Real pull data function.v2.source.R")

# Acute Myeloid Leukaemia -------------------------------------------------


laml.all.data <- real.pull.all.data(project = "TCGA-LAML", 
                                   tumor = "LAML", 
                                   save.filename = "LAML.Exp.rda")

saveRDS(laml.all.data, file = "./laml.all.data.rds")
laml.all.data <- readRDS(file = "./laml.all.data.rds")

laml.exp.matrix <- assay(laml.all.data$expression)
laml.exp.matrix[1:10,1:10]
dim(laml.exp.matrix)

laml.exp.matrix.filtered.gp <- laml.exp.matrix[apply(laml.exp.matrix,1,gp.style.filter,fold.change=5, delta=500, prop=0.05, base=10000, prop.base=0.05, na.rm = TRUE, neg.rm = TRUE),]
dim(laml.exp.matrix.filtered.gp)

### make the first column gene the same as the row names
laml.exp.matrix.aracne <- data.frame(gene = rownames(laml.exp.matrix.filtered.gp),laml.exp.matrix.filtered.gp)

### write out as a tab delimitted text file
write.table(laml.exp.matrix.aracne, file = "./aracne_input/laml_test/laml.exp.matrix.aracne.txt", sep = "\t", quote = F, row.names = F)

write.table(laml.exp.matrix.aracne$gene[1:1000], file = "./aracne_input/laml_test/laml.exp.matrix.aracne.gene.names.txt", sep = "\t", quote = F, row.names = F, col.names = F)








selected.patients <- laml.all.data$clinical[laml.all.data$clinical$age_at_initial_pathologic_diagnosis <30,"bcr_patient_barcode"]

laml.all.data.filt.match <- laml.all.data

which(laml.all.data$expression$patient%in%selected.patients) -> exp.idx
laml.all.data$expression[,exp.idx] -> laml.all.data.filt.match$expression

laml.all.data.filt.match$expression <- laml.all.data.filt.match$expression[,which(!duplicated(laml.all.data.filt.match$expression$patient))]


laml.all.data$clinical <- laml.all.data$clinical[which(!duplicated(laml.all.data$clinical$bcr_patient_barcode)),]
match(laml.all.data.filt.match$expression$patient, laml.all.data$clinical$bcr_patient_barcode) -> clin.idx
clin.idx <- clin.idx[!is.na(clin.idx)]
laml.all.data$clinical[clin.idx,] -> laml.all.data.filt.match$clinical

laml.all.data.filt.match$mutational <- laml.all.data$mutational[substr(laml.all.data$mutational$Tumor_Sample_Barcode,1,12)%in%selected.patients,]









# Glioblastoma Multiforme -------------------------------------------------


gbm.all.data <- real.pull.all.data(project = "TCGA-GBM", 
                                    tumor = "GBM", 
                                    save.filename = "GBM.Exp.rda")

saveRDS(gbm.all.data, file = "./gbm.all.data.rds")
gbm.all.data <- readRDS(file = "./gbm.all.data.rds")

gbm.exp.matrix <- assay(gbm.all.data$expression)
gbm.exp.matrix[1:10,1:10]
dim(gbm.exp.matrix)

gbm.exp.matrix.filtered.gp <- gbm.exp.matrix[apply(gbm.exp.matrix,1,gp.style.filter,fold.change=5, delta=500, prop=0.05, base=10000, prop.base=0.05, na.rm = TRUE, neg.rm = TRUE),]
dim(gbm.exp.matrix.filtered.gp)

### make the first column gene the same as the row names
gbm.exp.matrix.aracne <- data.frame(gene = rownames(gbm.exp.matrix.filtered.gp),gbm.exp.matrix.filtered.gp)

### write out as a tab delimitted text file
write.table(gbm.exp.matrix.aracne, file = "./aracne_input/gbm_test/gbm.exp.matrix.aracne.txt", sep = "\t", quote = F, row.names = F)

write.table(gbm.exp.matrix.aracne$gene[1:16304], file = "./aracne_input/gbm_test/gbm.exp.matrix.aracne.gene.names.txt", sep = "\t", quote = F, row.names = F, col.names = F)









selected.patients <- gbm.all.data$clinical[gbm.all.data$clinical$age_at_initial_pathologic_diagnosis <30,"bcr_patient_barcode"]

gbm.all.data.filt.match <- gbm.all.data

which(gbm.all.data$expression$patient%in%selected.patients) -> exp.idx
gbm.all.data$expression[,exp.idx] -> gbm.all.data.filt.match$expression

gbm.all.data.filt.match$expression <- gbm.all.data.filt.match$expression[,which(!duplicated(gbm.all.data.filt.match$expression$patient))]


gbm.all.data$clinical <- gbm.all.data$clinical[which(!duplicated(gbm.all.data$clinical$bcr_patient_barcode)),]
match(gbm.all.data.filt.match$expression$patient, gbm.all.data$clinical$bcr_patient_barcode) -> clin.idx
clin.idx <- clin.idx[!is.na(clin.idx)]
gbm.all.data$clinical[clin.idx,] -> gbm.all.data.filt.match$clinical

gbm.all.data.filt.match$mutational <- gbm.all.data$mutational[substr(gbm.all.data$mutational$Tumor_Sample_Barcode,1,12)%in%selected.patients,]









# 	Kidney Renal Clear Cell Carcinoma --------------------------------------

kirc.all.data <- real.pull.all.data(project = "TCGA-KIRC", 
                                    tumor = "KIRC", 
                                    save.filename = "KIRC.Exp.rda")

saveRDS(kirc.all.data, file = "./kirc.all.data.rds")
kirc.all.data <- readRDS(file = "./kirc.all.data.rds")

kirc.exp.matrix <- assay(kirc.all.data$expression)
kirc.exp.matrix[1:10,1:10]
dim(kirc.exp.matrix)

kirc.exp.matrix.filtered.gp <- kirc.exp.matrix[apply(kirc.exp.matrix,1,gp.style.filter,fold.change=5, delta=500, prop=0.05, base=10000, prop.base=0.05, na.rm = TRUE, neg.rm = TRUE),]
dim(kirc.exp.matrix.filtered.gp)

### make the first column gene the same as the row names
kirc.exp.matrix.aracne <- data.frame(gene = rownames(kirc.exp.matrix.filtered.gp),kirc.exp.matrix.filtered.gp)

### write out as a tab delimitted text file
write.table(kirc.exp.matrix.aracne, file = "./aracne_input/kirc_test/kirc.exp.matrix.aracne.txt", sep = "\t", quote = F, row.names = F)

write.table(kirc.exp.matrix.aracne$gene[1:1000], file = "./aracne_input/kirc_test/kirc.exp.matrix.aracne.gene.names.txt", sep = "\t", quote = F, row.names = F, col.names = F)
























selected.patients <- kirc.all.data$clinical[kirc.all.data$clinical$age_at_initial_pathologic_diagnosis <30,"bcr_patient_barcode"]

kirc.all.data.filt.match <- kirc.all.data

which(kirc.all.data$expression$patient%in%selected.patients) -> exp.idx
kirc.all.data$expression[,exp.idx] -> kirc.all.data.filt.match$expression

kirc.all.data.filt.match$expression <- kirc.all.data.filt.match$expression[,which(!duplicated(kirc.all.data.filt.match$expression$patient))]


kirc.all.data$clinical <- kirc.all.data$clinical[which(!duplicated(kirc.all.data$clinical$bcr_patient_barcode)),]
match(kirc.all.data.filt.match$expression$patient, kirc.all.data$clinical$bcr_patient_barcode) -> clin.idx
clin.idx <- clin.idx[!is.na(clin.idx)]
kirc.all.data$clinical[clin.idx,] -> kirc.all.data.filt.match$clinical

kirc.all.data.filt.match$mutational <- kirc.all.data$mutational[substr(kirc.all.data$mutational$Tumor_Sample_Barcode,1,12)%in%selected.patients,]








# Lower Grade Glioma ------------------------------------------------------

lgg.all.data <- real.pull.all.data(project = "TCGA-LGG", 
                                    tumor = "LGG", 
                                    save.filename = "LGG.Exp.rda")

saveRDS(lgg.all.data, file = "./lgg.all.data.rds")
lgg.all.data <- readRDS(file = "./lgg.all.data.rds")

lgg.exp.matrix <- assay(lgg.all.data$expression)
lgg.exp.matrix[1:10,1:10]
dim(lgg.exp.matrix)

##### ARACNE
### remove genes that are invariant or not expressing at a high value
lgg.exp.matrix.filtered.gp <- lgg.exp.matrix[apply(lgg.exp.matrix,1,gp.style.filter,fold.change=5, delta=500, prop=0.05, base=10000, prop.base=0.05, na.rm = TRUE, neg.rm = TRUE),]
dim(lgg.exp.matrix.filtered.gp)

### make the first column gene the same as the row names
lgg.exp.matrix.aracne <- data.frame(gene = rownames(lgg.exp.matrix.filtered.gp),lgg.exp.matrix.filtered.gp)

### write out as a tab delimitted text file
write.table(lgg.exp.matrix.aracne, file = "./aracne_input/lgg_test/lgg.exp.matrix.aracne.txt", sep = "\t", quote = F, row.names = F)

write.table(lgg.exp.matrix.aracne$gene[1:1000], file = "./aracne_input/lgg_test/lgg.exp.matrix.aracne.gene.names.txt", sep = "\t", quote = F, row.names = F, col.names = F)


##################

selected.patients <- lgg.all.data$clinical[lgg.all.data$clinical$age_at_initial_pathologic_diagnosis <30,"bcr_patient_barcode"]

lgg.all.data.filt.match <- lgg.all.data

which(lgg.all.data$expression$patient%in%selected.patients) -> exp.idx
lgg.all.data$expression[,exp.idx] -> lgg.all.data.filt.match$expression

lgg.all.data.filt.match$expression <- lgg.all.data.filt.match$expression[,which(!duplicated(lgg.all.data.filt.match$expression$patient))]

lgg.all.data$clinical <- lgg.all.data$clinical[which(!duplicated(lgg.all.data$clinical$bcr_patient_barcode)),]
match(lgg.all.data.filt.match$expression$patient, lgg.all.data$clinical$bcr_patient_barcode) -> clin.idx
clin.idx <- clin.idx[!is.na(clin.idx)]
lgg.all.data$clinical[clin.idx,] -> lgg.all.data.filt.match$clinical

lgg.all.data.filt.match$mutational <- lgg.all.data$mutational[substr(lgg.all.data$mutational$Tumor_Sample_Barcode,1,12)%in%selected.patients,]






# Lymphoid Neoplasm Diffuse Large B-cell Lymphoma -------------------------


dlbc.all.data <- real.pull.all.data(project = "TCGA-DLBC", 
                                    tumor = "DLBC", 
                                    save.filename = "DLBC.Exp.rda")

saveRDS(dlbc.all.data, file = "./dlbc.all.data.rds")
dlbc.all.data <- readRDS(file = "./dlbc.all.data.rds")

dlbc.exp.matrix <- assay(dlbc.all.data$expression)
dlbc.exp.matrix[1:10,1:10]
dim(dlbc.exp.matrix)

dlbc.exp.matrix.filtered.gp <- dlbc.exp.matrix[apply(dlbc.exp.matrix,1,gp.style.filter,fold.change=5, delta=500, prop=0.05, base=10000, prop.base=0.05, na.rm = TRUE, neg.rm = TRUE),]
dim(dlbc.exp.matrix.filtered.gp)

### make the first column gene the same as the row names
dlbc.exp.matrix.aracne <- data.frame(gene = rownames(dlbc.exp.matrix.filtered.gp),dlbc.exp.matrix.filtered.gp)

### write out as a tab delimitted text file
write.table(dlbc.exp.matrix.aracne, file = "./aracne_input/dlbc_test/dlbc.exp.matrix.aracne.txt", sep = "\t", quote = F, row.names = F)

write.table(dlbc.exp.matrix.aracne$gene[1:1000], file = "./aracne_input/dlbc_test/dlbc.exp.matrix.aracne.gene.names.txt", sep = "\t", quote = F, row.names = F, col.names = F)



















