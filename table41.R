library(readxl)
library(tidyverse)
library(moonBook)
library(xlsx)
library(edgeR)

####
cx <- readxl::read_excel(path= "dataset/cx.xlsx", sheet="cx_all", col_names=TRUE) 

load(file='dataset/smRa.Rdata')
miRa <- smRa[["miRa"]]

miRa1 <- miRa %>% select(Mature_ID, ends_with("01_Read_Count"), ends_with("Pre_Read_Count"), ends_with("-1_Read_Count"),
                         ends_with("02_Read_Count"), ends_with("After_Read_Count"), ends_with("-2_Read_Count")) %>% select(Mature_ID, contains(cx$ID))

## TMN normalization
mirna2 <- column_to_rownames(miRa1, var ="Mature_ID") 
y <- DGEList(counts=mirna2)
y <- calcNormFactors(y)
cpms <- cpm(y, log=TRUE)

## H- clustering
f1 <- cpms %>% as.data.frame() %>% t
f <- hclust(dist(f1, method="euclidean"), method = "average")
plot(f)


###
cx <- readxl::read_excel(path= "dataset/cx.xlsx", sheet="cx_all", col_names=TRUE) %>% filter(ID != "25580130")

load(file='dataset/smRa.Rdata')
miRa <- smRa[["miRa"]]

miRa1 <- miRa %>% select(Mature_ID, ends_with("01_Read_Count"), ends_with("Pre_Read_Count"), ends_with("-1_Read_Count"),
                         ends_with("02_Read_Count"), ends_with("After_Read_Count"), ends_with("-2_Read_Count")) %>% select(Mature_ID, contains(cx$ID))

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

names(mirna)
nrow(cx)
mirna %>% filter(Mature_ID=="hsa-miR-150-3p")

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


miR_data <- mirna  %>% column_to_rownames(var="Mature_ID") %>%  select(contains("FC")) %>% t %>% as.data.frame() %>%  
  rownames_to_column(var="ID") %>% as_tibble() %>% select(ID, 'hsa-miR-150-3p') %>% mutate(ID = substr(ID,1,8))



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

mR_data <- mrna  %>% column_to_rownames(var="Gene_Symbol") %>%  select(contains("FC")) %>% t %>% as.data.frame() %>%  rownames_to_column(var="ID") %>% as_tibble() %>%
  mutate(ID = substr(ID,1,8)) %>% select(ID, PRDM1, NMT2)

test <- cx %>% select(ID, age, pathology, stage, tv0, tv2, ar, dose0:date2, alc0, alc1, alc2, RT_field, ICR, fraction_size, EQD2) %>% 
  inner_join(miR_data, by="ID") %>% 
  inner_join(mR_data, by="ID")  %>% 
  mutate(total = PRDM1+`hsa-miR-150-3p`, total1 = NMT2+PRDM1+`hsa-miR-150-3p`, miR=`hsa-miR-150-3p`,
         EQD2_m = ifelse(EQD2<median(EQD2), 1, 0), ICR_m = ifelse(is.na(ICR)==TRUE, "N/A", 
                                                                  ifelse(fraction_size<4, "EBRT", paste0(ICR,"/",ICR/fraction_size))),
         Age_m=ifelse(age<50,1,0), stage1=ifelse(stage=="3C2"|stage=="4A"|stage=="4B","IIIC2-IVB", 
                                                 ifelse(stage=="1B", "IB", "IIB-IIIC1")))


test %>%
  mutate(`Initial tumor volume` = tv0) %>%
  ggplot(aes(ar, total1, color = `Initial tumor volume`)) +
  geom_point(size = 3) +
  scale_color_gradient(low = "blue", high = "red") +  # 연속형 색상 그라데이션 설정
  theme_minimal() +
  geom_smooth(method = "lm", se = TRUE) +
  labs(x = "Acute tumor response", y = "miR-150-3p+PRDM1+NMT2(log2FC)") +
  theme(
    axis.title = element_text(size = 13, face = "bold"),
    title = element_text(size = 13, face = "bold"),
    plot.title = element_text(hjust = 0.5, vjust = -2.5),
    axis.text = element_text(size = 13, face = "bold"),
    legend.title = element_text(size = 13, face = "bold"), # Adjusts legend title text
    legend.position = "top" # Adjusts legend position, can be "top", "bottom", "left", or "right"
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray", size = 1) +
  ylim(-5, 3) +
  scale_y_continuous(limits = c(-6, 6))

  


summary(test)
out = mytable( ~ tv0+tv2+ar+Age_m+pathology+stage1+RT_field+EQD2_m+ICR_m+alc0+alc1+alc2+miR+PRDM1+NMT2+total+total1, data=test, method=2, digits=2)
mycsv(out, file="table/table_41.csv")
