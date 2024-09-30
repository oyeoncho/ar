library(readxl)
library(dplyr)
library(moonBook)
library(lubridate)
library(tidyverse)
library(edgeR)
library (Hmisc)
library (leaps)
library(xlsx)
library(corrplot)
library(caret)

## tv idff testset vs val set
cx_val <- readxl::read_excel(path= "dataset/cx.xlsx", sheet="cx_val", col_names=TRUE) %>% filter(ID != "25580130" & tv0 >4) %>% 
  mutate(group="val") %>% select(ID, tv0, tv2, group)
cx_test <- readxl::read_excel(path= "dataset/cx.xlsx", sheet="cx_test", col_names=TRUE) %>% filter(ID != "25580130" & tv0 >4) %>%
  mutate(group="test") %>% select(ID, tv0, tv2, group)

tv_diff <- cx_val %>% bind_rows(cx_test)

mytable(group~., data=tv_diff, method=2, digits=2)

# load all data 

cx1 <- readxl::read_excel(path= "dataset/cx.xlsx", sheet="cx_all", col_names=TRUE) 
cx1 %>% ggplot(aes(tv0)) + geom_density(position = "identity",alpha=0.3, size=1.5) +scale_x_continuous(breaks=seq(0,250,20))


z1 <- c()
for (i in 4:200){
  m <- cx1 %>% filter(tv0>i & ID != "25580130") %>% nrow()
  z <- c(i, m)
  z1 <- rbind(z1,z)
}
colnames(z1) <- c("tv", "number")

a0 <- z1 %>% as_tibble() %>% ggplot(aes(tv, number)) +geom_line()+ labs(
  x = "Initial tumor volume (cm³)", 
  y = "Patients whose initial tumor volume \nexceeds a certain threshold"
) +
  theme(
    axis.title = element_text(size = 15, face = "bold"), 
    title = element_text(size = 15, face = "bold"),
    plot.title = element_text(hjust = 0.5, vjust = -1, size=15), axis.text = element_text(size=13))+ geom_vline(xintercept=40, color="red", lty=2)+
    scale_y_continuous(breaks=seq(0,40,2))+ scale_x_continuous(breaks=seq(0,200,20))
  
  
  
  

tv <- seq(4, 48, 4)

load(file='dataset/smRa.Rdata')
miRa <- smRa[["miRa"]]

miRa1 <- miRa %>%  select(Mature_ID, ends_with("01_Read_Count"), ends_with("Pre_Read_Count"), ends_with("-1_Read_Count"),
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


final <- c()

for (k in tv) {
  cx <- cx1 %>% filter(tv0>k,ID != "25580130") 
  
  miRa1 <- miRa %>% select(Mature_ID, contains(cx$ID)) %>% 
    select(Mature_ID, ends_with("01_Read_Count"), ends_with("Pre_Read_Count"), ends_with("-1_Read_Count"),
           ends_with("02_Read_Count"), ends_with("After_Read_Count"), ends_with("-2_Read_Count")) 
  
  num <- c() 
  for (i in 1:nrow(miRa1)) {
    out <- sum(miRa1[i,]==0)
    num <- rbind(num, out)
  } # num - NA number
  
  miRa1_num <- cbind(miRa1, num) %>% 
    as_tibble() 
  
  mirna <- miRa1_num %>% filter(num<nrow(cx)) %>% select(-num)
  
  for (i in 2:(nrow(cx)+1)) {
    mirna1 <- column_to_rownames(mirna, var ="Mature_ID") %>% 
      select(colnames(mirna)[[i]], colnames(mirna)[[i+nrow(cx)]])
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
  
  
  miR <- cx %>% select(ID, ar) %>% inner_join(miR, by="ID")
  
  mRa <- smRa[["mRa"]]
  mRa1 <- mRa %>% select(Gene_Symbol, contains(cx$ID)) %>% select(Gene_Symbol, ends_with("01_Read_Count"), ends_with("Pre_Read_Count"), ends_with("-1_Read_Count"),
                                                                 ends_with("02_Read_Count"), ends_with("After_Read_Count"), ends_with("-2_Read_Count")) 
  
  
  
  num <- c() 
  for (i in 1:nrow(mRa1)) {
    out <- sum(mRa1[i,]==0)
    num <- rbind(num, out)
  } # num - NA number
  
  mRa1_num <- cbind(mRa1, num) %>% 
    as_tibble() 
  
  mrna <- mRa1_num %>% filter(num<nrow(cx)) %>% select(-num)
  
  
  for (i in 2:(nrow(cx)+1)) {
    mrna1 <- column_to_rownames(mrna, var ="Gene_Symbol") %>% 
      select(colnames(mrna)[[i]], colnames(mrna)[[i+nrow(cx)]])
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
  
  mR <- mrna %>% as_tibble()  %>% select(Gene_Symbol, contains("FC")) %>% column_to_rownames(var="Gene_Symbol") %>% t %>% 
    as.data.frame() %>% rownames_to_column(var="ID") %>% as_tibble() %>% mutate(ID=substr(ID,1,8)) 
  
  all_rna <- miR %>% inner_join(mR, by="ID") %>% select(ar, `hsa-miR-150-3p`, `hsa-miR-424-3p`, PRDM1, NMT2) %>% 
    mutate(var1=`hsa-miR-150-3p`+ NMT2+PRDM1, var2= `hsa-miR-150-3p`+ PRDM1)
  
  ctrl <- trainControl(method = "LOOCV")
  model0 <- train(ar ~ `hsa-miR-150-3p`+ PRDM1+ NMT2, data = all_rna, method = "lm", trControl = ctrl)
  model1 <- train(ar ~ var1, data = all_rna, method = "lm", trControl = ctrl)
  model2 <- train(ar ~ var2, data = all_rna, method = "lm", trControl = ctrl)
  model3 <- train(ar ~ `hsa-miR-150-3p`, data = all_rna, method = "lm", trControl = ctrl)
  model4 <- train(ar ~ PRDM1, data = all_rna, method = "lm", trControl = ctrl)
  model5 <- train(ar ~ NMT2, data = all_rna, method = "lm", trControl = ctrl)
  
  result <- bind_rows(model0$results, model1$results, model2$results, model3$results,model4$results,model5$results) %>% as_tibble %>% 
    mutate(volume=k, var=c("lr", "miR-150-3p+NMT2+PRDM1 (log2FC)", "miR-150-3p+PRDM1 (log2FC)","miR-150-3p (log2FC)", "PRDM1 (log2FC)", "NMT2 (log2FC)"), num=nrow(all_rna))
  
  final <- bind_rows(final, result)
}


write.xlsx(final, file="process/loocv.xlsx", sheetName="final",append=TRUE)

