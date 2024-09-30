library(readxl)
library(dplyr)
library(moonBook)
library(lubridate)
library(tidyverse)
library(edgeR)
library(ggpubr)

# load all data 

cx <- readxl::read_excel(path= "dataset/cx.xlsx", sheet="cx_val", col_names=TRUE)


load(file='dataset/smRa.Rdata')
miRa <- smRa[["miRa"]]


miRa1 <- miRa %>% select(Mature_ID, ends_with("01_Read_Count"), ends_with("Pre_Read_Count"), ends_with("-1_Read_Count"),
                         ends_with("02_Read_Count"), ends_with("After_Read_Count"), ends_with("-2_Read_Count")) %>% select(Mature_ID, contains(cx$ID))


m <- cx %>%  filter(tv0 > 4) 

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


miR_data <- mirna  %>% column_to_rownames(var="Mature_ID") %>%  select(contains("FC")) %>% t %>% as.data.frame() %>%  
  rownames_to_column(var="ID") %>% as_tibble() %>% select(ID, 'hsa-miR-150-3p') %>% mutate(ID = substr(ID,1,8))
   
  

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

mrna %>% filter(Gene_Symbol=="PRDM1")

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

mR_data <- mrna  %>% column_to_rownames(var="Gene_Symbol") %>%  select(contains("FC")) %>% t %>% as.data.frame() %>%  rownames_to_column(var="ID") %>% as_tibble() %>%
  mutate(ID = substr(ID,1,8)) %>% select(ID, PRDM1, NMT2)

test <- m %>% mutate(`tumor response`=ifelse(ar<0.3,"good","poor"), ar_m=ifelse(ar<0.3,0,1)) %>% 
  select(ID, age, pathology, stage, tv0, tv2, ar, `tumor response`, ar_m, dose0:date2, alc0, alc1, alc2, RT_field, ICR, fraction_size, EQD2) %>% inner_join(miR_data, by="ID") %>% 
  inner_join(mR_data, by="ID") 

library(moonBook)
library(ggpubr)

test$mir <- test$`hsa-miR-150-3p`
test$total <- test$`hsa-miR-150-3p`+test$PRDM1
test$total1 <- test$`hsa-miR-150-3p`+test$PRDM1+test$NMT2

summary(lm((ar~mir), data=test))
summary(lm((ar~PRDM1), data=test))
summary(lm((ar~total), data=test))
summary(lm((ar~total1), data=test))


a1 <-   test %>%  ggplot(aes(mir, ar, color=` tumor response`))+geom_point(size=3)+
  geom_abline(intercept= 0.21822, slope=-0.02944 , color='black', linewidth = 1.5)+
  annotate("text",-1,0.5,label=expression(paste(R^2,"=0.138")), size=10)+
  theme(axis.title =element_text(face="bold", size=15), axis.text = element_text(size=15, face="bold"), legend.position = "top")+
  ylab("Acute tumor response")+xlab("miR-150-3p (log2FC)")+
  scale_y_continuous(breaks=seq(0,1,0.1))+
  scale_x_continuous(breaks=seq(-6,6,1))+
  ggtitle("Correlation between tumor response and miR-150-3p (validation)")

a2 <- test %>%  ggplot(aes(PRDM1, ar, color=`tumor response`))+geom_point(size=3)+
  geom_abline(intercept= 0.22986, slope=-0.05197 , color='black', linewidth = 1.5)+
  annotate("text",-1,0.5,label=expression(paste(R^2,"=0.3753")), size=10)+
  theme(axis.title =element_text(face="bold", size=15), axis.text = element_text(size=15, face="bold"), legend.position = "top")+
  ylab("Acute tumor response")+xlab("PRDM1 (log2FC)")+
  scale_y_continuous(breaks=seq(0,1,0.1))+
  scale_x_continuous(breaks=seq(-6,6,1))+
  ggtitle("Correlation between tumor response and PRDM1 (validation)")

a3 <- test %>%  ggplot(aes(total, ar, color=`tumor response`))+geom_point(size=3)+
  geom_abline(intercept= 0.215346, slope=-0.036698  , color='black', linewidth = 1.5)+
  annotate("text",-1,0.5,label=expression(paste(R^2,"=0.437")), size=10)+
  theme(axis.title =element_text(face="bold", size=15), axis.text = element_text(size=15, face="bold"), legend.position = "top")+
  ylab("Acute tumor response")+xlab("miR-150-3p + PRDM1 (log2FC)")+
  scale_y_continuous(breaks=seq(0,1,0.1))+
  scale_x_continuous(breaks=seq(-6,6,1))+
  ggtitle("Correlation between tumor response and miR-150-3p + PRDM1 (validation)")

a4 <- test %>%  ggplot(aes(total1, ar, color=`tumor response`))+geom_point(size=3)+
  geom_abline(intercept= 0.223246, slope=-0.03427  , color='black', linewidth = 1.5)+
  annotate("text",-1,0.5,label=expression(paste(R^2,"=0.4957")), size=10)+
  theme(axis.title =element_text(face="bold", size=15), axis.text = element_text(size=15, face="bold"), legend.position = "top")+
  ylab("Acute tumor response")+xlab("miR-150-3p + PRDM1 + NMT2 (log2FC)")+
  scale_y_continuous(breaks=seq(0,1,0.1))+
  scale_x_continuous(breaks=seq(-6,6,1))+
  ggtitle("Correlation between tumor response and miR-150-3p + PRDM1 + NMT2 (validation)")


a5 <- test %>%  ggplot(aes(total, alc1, color=`tumor response`))+geom_point(size=3)+
  theme(axis.title =element_text(face="bold", size=15), axis.text = element_text(size=15, face="bold"), legend.position = "top")+
  xlab("miR-150-3p + PRDM1 (log2FC)")+ylab("Abolute lymphocyte counts 1 week (ceels/uL)")+
  scale_y_continuous(breaks=seq(0,7000,1000))+
  scale_x_continuous(breaks=seq(-6,6,1))+
  ggtitle("Correlation between abolute lymphocyte counts and miR-150-3p + PRDM1 (validation)")


a6 <- test %>%  ggplot(aes(total1, alc1, color=`tumor response`))+geom_point(size=3)+
  theme(axis.title =element_text(face="bold", size=15), axis.text = element_text(size=15, face="bold"), legend.position = "top")+
  xlab("miR-150-3p + PRDM1 +NMT2 (log2FC)")+ylab("Abolute lymphocyte counts 1 week (ceels/uL)")+
  scale_y_continuous(breaks=seq(0,7000,1000))+
  scale_x_continuous(breaks=seq(-6,6,1))+
  ggtitle("Correlation between abolute lymphocyte counts and miR-150-3p + PRDM1 +NMT2 (validation)")
  

ggarrange(a3,a4,a5,a6, labels=c("A","B","C","D"), nrow=2, ncol=2, font.label = list(size = 30, color = "black"))


b3 <- ggplot(data = test, aes(x = total1, y =ar)) +
  geom_point(size=3)+
  geom_hline(yintercept=0, linetype='dashed', color='red', size=1)+
  geom_smooth(method="lm")+
  labs(x="miR-150-3p+PRDM1+NMT2(log2FC)", y="Acute tumor response")+
  theme(axis.title = element_text(size=15, face="bold"), title=element_text(size=15, face="bold"),
        plot.title=element_text(hjust = 0.5, vjust=-1), axis.text=element_text(size=15, face="bold"))+ 
  annotate("text", x=0, y=0.5, label=expression(R^2==0.4957)) 


p <- ggboxplot(test, x = "tumor response", y = "total1",
               color = "tumor response", palette = c("#0000FF", "#FF0000"),
               add = "jitter", shape = "tumor response")
b4 <- p+ stat_compare_means(label.x = 1.5, label.y = 4, method="wilcox.test", label = "p.signif", size=7) + 
  labs(x="Acute tumor response", y="miR-150-3p+PRDM1+NMT2(log2FC)")+theme(axis.title = element_text(size=13, face="bold"),
                                                                          title=element_text(size=13, face="bold"), plot.title=element_text(hjust = 0.5, vjust=-2.5), 
                                                                          axis.text=element_text(size=13, face="bold"))+geom_hline(yintercept = 0, linetype="dashed", color="gray", size=1)+ylim(-5,3)+
                                                                          scale_y_continuous(limits = c(-6, 6))



ggarrange(b3, b4, labels=c("E","F"), nrow=1, ncol=2,font.label = list(size = 16, color = "black", face = "bold"))


summary(test)

test <- test %>% mutate(stage1=ifelse(stage=="3C2","IIIC2", 
                              ifelse(stage=="1B", "IB", "IIB-IIIC1")), EQD2_m = ifelse(EQD2<median(EQD2), 1, 0), ICR_m = paste0(ICR,"/",ICR/fraction_size),
                        Age_m=ifelse(age<50,1,0)) 
  

summary(test)


out = mytable(ar_m ~ tv0+tv2+ar+Age_m+pathology+stage1+RT_field+EQD2_m+ICR_m+alc0+alc1+alc2+mir+PRDM1+NMT2+total+total1, data=test, method=2, digits=2)
mycsv(out, file="table/table_val.csv")

out1 = mytable( ~ tv0+tv2+ar+Age_m+pathology+stage1+RT_field+EQD2_m+ICR_m+alc0+alc1+alc2+mir+PRDM1+NMT2+total+total1, data=test, method=2, digits=2)
mycsv(out1, file="table/table_val_all.csv")

out2 = mytable( ~ tv0+tv2+ar+age+EQD2+alc0+alc1+alc2+mir+PRDM1+NMT2+total+total1, data=test, method=2, digits=2)
mycsv(out2, file="table/median_val.csv")


