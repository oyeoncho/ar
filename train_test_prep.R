library(readxl)
library(dplyr)
library(moonBook)
library(lubridate)
library(tidyverse)
library(edgeR)
library (Hmisc)
library (leaps)
library(xlsx)
# load all data 

cx <- readxl::read_excel(path= "dataset/cx.xlsx", sheet="cx_test", col_names=TRUE)


load(file='dataset/smRa.Rdata')
miRa <- smRa[["miRa"]]

miRa1 <- miRa %>% select(Mature_ID, contains(cx$ID)) %>% select(Mature_ID, ends_with("01_Read_Count"), ends_with("Pre_Read_Count"), ends_with("-1_Read_Count"),
                         ends_with("02_Read_Count"), ends_with("After_Read_Count"), ends_with("-2_Read_Count")) 


## TMN normalization
mirna2 <- column_to_rownames(miRa1, var ="Mature_ID") 
y <- DGEList(counts=mirna2)
y <- calcNormFactors(y)
cpms <- cpm(y, log=TRUE)

## H- clustering
f1 <- cpms %>% as.data.frame() %>% t
f <- hclust(dist(f1, method="euclidean"), method = "average")
plot(f)

m <- cx  %>% filter(ID != "25580130") %>% select(ID, ar) 

miRa1 <- miRa %>% select(Mature_ID, contains(m$ID)) %>% 
  select(Mature_ID, ends_with("01_Read_Count"), ends_with("Pre_Read_Count"), ends_with("-1_Read_Count"),
         ends_with("02_Read_Count"), ends_with("After_Read_Count"), ends_with("-2_Read_Count")) 

num <- c() 
for (i in 1:nrow(miRa1)) {
  out <- sum(miRa1[i,]==0)
  num <- rbind(num, out)
} # num - NA number

miRa1_num <- cbind(miRa1, num) %>% 
  as_tibble() 

mirna <- miRa1_num %>% filter(num<nrow(m)) %>% select(-num)

names(mirna)
nrow(m)
mirna %>% filter(Mature_ID=="hsa-miR-150-3p")

for (i in 2:(nrow(m)+1)) {
  mirna1 <- column_to_rownames(mirna, var ="Mature_ID") %>% 
    select(colnames(mirna)[[i]], colnames(mirna)[[i+nrow(m)]])
  dge <- DGEList(counts=mirna1, group=c(1,2))
  dge <- calcNormFactors(dge)
  design <- model.matrix(~group, dge$samples)
  rownames(design) <- rownames(dge$samples)
  dge_dispersion <- estimateGLMCommonDisp(dge)
  dge_tagwise <- estimateGLMTagwiseDisp(dge_dispersion)
  plotBCV(dge_tagwise)
  fit <- glmFit(dge_tagwise, dispersion = dge_tagwise$tagwise.dispersion)
  lrt <- glmLRT(fit)
  norm <- round(cpm(dge_tagwise, normalized.lib.sizes = T))
  normEXP <- merge(norm, topTags(lrt, n=10000), by.x="row.names", by.y="row.names")
  names(normEXP)[1]
  head(normEXP)
  mirna_fc <- data.frame(normEXP$Row.names, normEXP$logFC, normEXP$PValue)
  colnames(mirna_fc)[1] <- "Mature_ID"
  colnames(mirna_fc)[2] <- str_c(str_sub(colnames(mirna)[[i]], 1, 8), "_FC")
  colnames(mirna_fc)[3] <- str_c(str_sub(colnames(mirna)[[i]], 1, 8), "_pvalue")
  mirna <- merge(mirna, mirna_fc, by="Mature_ID")
}

miR <- mirna %>% as_tibble() %>% select(Mature_ID, contains("FC")) %>% column_to_rownames(var="Mature_ID") %>% t %>% 
  as.data.frame() %>% rownames_to_column(var="ID") %>% as_tibble() %>% mutate(ID=substr(ID,1,8)) 

miR <- m %>% inner_join(miR, by="ID")


###
mRa <- smRa[["mRa"]]
mRa1 <- mRa %>% select(Gene_Symbol, contains(m$ID)) %>% select(Gene_Symbol, ends_with("01_Read_Count"), ends_with("Pre_Read_Count"), ends_with("-1_Read_Count"),
                                                               ends_with("02_Read_Count"), ends_with("After_Read_Count"), ends_with("-2_Read_Count")) 



num <- c() 
for (i in 1:nrow(mRa1)) {
  out <- sum(mRa1[i,]==0)
  num <- rbind(num, out)
} # num - NA number

mRa1_num <- cbind(mRa1, num) %>% 
  as_tibble() 

mrna <- mRa1_num %>% filter(num<nrow(m)) %>% select(-num)

for (i in 2:(nrow(m)+1)) {
  mrna1 <- column_to_rownames(mrna, var ="Gene_Symbol") %>% 
    select(colnames(mrna)[[i]], colnames(mrna)[[i+nrow(m)]])
  dge <- DGEList(counts=mrna1, group=c(1,2))
  dge <- calcNormFactors(dge)
  design <- model.matrix(~group, dge$samples)
  rownames(design) <- rownames(dge$samples)
  dge_dispersion <- estimateGLMCommonDisp(dge)
  dge_tagwise <- estimateGLMTagwiseDisp(dge_dispersion)
  plotBCV(dge_tagwise)
  fit <- glmFit(dge_tagwise, dispersion = dge_tagwise$tagwise.dispersion)
  lrt <- glmLRT(fit)
  norm <- round(cpm(dge_tagwise, normalized.lib.sizes = T))
  normEXP <- merge(norm, topTags(lrt, n=30000), by.x="row.names", by.y="row.names")
  names(normEXP)[1]
  head(normEXP)
  mrna_fc <- data.frame(normEXP$Row.names, normEXP$logFC, normEXP$PValue)
  colnames(mrna_fc)[1] <- "Gene_Symbol"
  colnames(mrna_fc)[2] <- str_c(str_sub(colnames(mrna)[[i]], 1, 8), "_FC")
  colnames(mrna_fc)[3] <- str_c(str_sub(colnames(mrna)[[i]], 1, 8), "_pvalue")
  mrna <- merge(mrna, mrna_fc, by="Gene_Symbol")
}

mR <- mrna %>% as_tibble() %>% select(Gene_Symbol, contains("FC")) %>% column_to_rownames(var="Gene_Symbol") %>% t %>% 
  as.data.frame() %>% rownames_to_column(var="ID") %>% as_tibble() %>% mutate(ID=substr(ID,1,8)) 

all_rna <- miR %>% inner_join(mR, by="ID")

a3 <- all_rna[,2:ncol(all_rna)]
x<-select_if (a3, is.numeric)
rx <- rcorr (as.matrix (x))$r
result<-ifelse (rcorr (as.matrix (x))$P<0.05, rx, 0)

result <-  result %>% as_tibble() %>% mutate(ID=colnames(a3))

write.csv(result, file="dataset/cor.csv", row.names = FALSE)
write.csv(mR, file="dataset/fc_mR.csv", row.names = FALSE)
write.xlsx(miR, file="dataset/fc_miR.xlsx", sheetName="ar",append=TRUE)
