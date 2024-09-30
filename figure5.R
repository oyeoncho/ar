library(readxl)
library(tidyverse)
library(lubridate)
library(moonBook)
library(survival)
library(survminer)
library(xlsx)
library(edgeR)
library(ggpubr)
library(gridExtra)          

cx1 <- readxl::read_excel(path= "dataset/cx.xlsx", sheet="cx_all", col_names=TRUE) 

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


final <- read_excel(path= "process/loocv.xlsx", sheet="final", col_names=TRUE) %>% select(-`...1`)

a1 <- final %>% 
  filter(var != "lr") %>% 
  ggplot(aes(volume, Rsquared, color = var)) +
  geom_line(size = 2) +
  scale_y_continuous(limits=c(0,1), breaks=seq(0,1,0.1)) +
  theme(legend.position = c(0.95, 0.95), legend.justification = c("right", "top")) +
  scale_color_discrete(breaks = c("miR-150-3p+NMT2+PRDM1 (log2FC)", "miR-150-3p+PRDM1 (log2FC)", "PRDM1 (log2FC)", "miR-150-3p (log2FC)", "NMT2 (log2FC)")) +
  theme(legend.title = element_blank()) +
  labs(
    x = "Patients with initial tumor volume \nabove a certain threshold in cm³", 
    y = "R square of LOOCV"
  ) +
  theme(
    axis.title = element_text(size = 15, face = "bold"), 
    title = element_text(size = 15, face = "bold"),
    plot.title = element_text(hjust = 0.5, vjust = -1, size=15), 
    axis.text = element_text(size = 15, face = "bold")
  )

a2 <- final %>% 
  filter(var != "lr") %>% 
  mutate(var = ifelse(var == "var2", "miR-150-3p+NMT2+PRDM1", 
                      ifelse(var == "var1", "miR-150-3p+PRDM1", var))) %>%
  ggplot(aes(volume, RMSE, color = var)) +
  geom_line(size = 2) +
  scale_y_continuous(limits=c(0,0.25), breaks=seq(0,0.25,0.05)) +
  theme(legend.position = c(0.95, 0.95), legend.justification = c("right", "top")) +
  scale_color_discrete(breaks = c("miR-150-3p+NMT2+PRDM1 (log2FC)", "miR-150-3p+PRDM1 (log2FC)", "PRDM1 (log2FC)", "miR-150-3p (log2FC)", "NMT2 (log2FC)")) +
  theme(legend.title = element_blank()) +
  labs(
    x = "Patients with initial tumor volume \nabove a certain threshold in cm³", 
    y = "RMSE of LOOCV")+
  theme(
    axis.title = element_text(size = 15, face = "bold"), 
    title = element_text(size = 15, face = "bold"),
    plot.title = element_text(hjust = 0.5, vjust = -1, size=15), 
    axis.text = element_text(size = 15, face = "bold")
  )



cx <- readxl::read_excel(path= "dataset/cx.xlsx", sheet="cx_all", col_names=TRUE) %>% filter(ID != "25580130" & EQD2 >= 50) 

load(file='dataset/smRa.Rdata')
miRa <- smRa[["miRa"]]

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

all_rna <- miR %>% inner_join(mR, by="ID") %>% select(ID, `hsa-miR-150-3p`,  PRDM1, NMT2) %>% mutate(miR=`hsa-miR-150-3p`,var1=`hsa-miR-150-3p`+ PRDM1, var2= `hsa-miR-150-3p`+ PRDM1 +NMT2)

test <- cx %>% inner_join(all_rna, by="ID") %>% mutate(var1_m = ifelse(var1 < median(var1), 1, 0), var2_m = ifelse(var2 < median(var2), 1, 0)
                                                       , PRDM1_m = ifelse(PRDM1 < median(PRDM1), 1, 0), recur_m = ifelse(recur1>0,1,0),
                                                        miR_m = ifelse(miR < median(miR), 1, 0), NMT2_m = ifelse(NMT2 < median(NMT2), 1, 0)) %>%
  mutate(age_m=ifelse(age<median(age),1,0), stage_m = ifelse(stage=="1B"|stage=="2B"|stage=="3C1",0,1), td_m = ifelse(EQD2<median(EQD2),1,0), path_m = ifelse(pathology=="sqcc",0,1) )

summary(test)

library(ggpubr)
library(gridExtra)
library(survminer)
library(survival)
library(dplyr)

# 예제 데이터 준비
fit0 <- survfit(Surv(recur_date, recur_m == 1) ~ miR_m, data = test) # 3PFS 80% vs. 69.6 p=0.5
fit1 <- survfit(Surv(recur_date, recur_m == 1) ~ NMT2_m, data = test) # 3PFS 80% vs. 69.6 p=0.54
fit2 <- survfit(Surv(recur_date, recur_m == 1) ~ PRDM1_m, data = test) # 3PFS 70% vs. 79.3 p=0.42


# 각 생존 곡선 플롯 및 위험표 생성
p0 <- ggsurvplot(fit0, xlab = "Months", ylab = "Progression free survival rate", pval = TRUE, pval.size = 6, fun = "pct", conf.int = FALSE, risk.table = TRUE, size = 1, linetype = c("dashed","solid"), palette = c("blue", "red"), legend.labs = c("miR-150-3p ≥ -0.67", "miR-150-3p < -0.67"), break.time.by = 12, xlim = c(0, 48), risk.table.height = 0.3, ylim = c(0, 100), ggtheme = theme_classic2(base_size = 10, base_family = "Arial"), font.family = "Arial")
p1 <- ggsurvplot(fit1, xlab = "Months", ylab = "Progression free survival rate", pval = TRUE, pval.size = 6, fun = "pct", conf.int = FALSE, risk.table = TRUE, size = 1, linetype = c("dashed","solid"), palette = c("blue", "red"), legend.labs = c("NMT2 ≥ -0.17", "NMT2 < -0.17"), break.time.by = 12, xlim = c(0, 48), risk.table.height = 0.3, ylim = c(0, 100), ggtheme = theme_classic2(base_size = 10, base_family = "Arial"), font.family = "Arial")
p2 <- ggsurvplot(fit2, xlab = "Months", ylab = "Progression free survival rate", pval = TRUE, pval.size = 6, fun = "pct", conf.int = FALSE, risk.table = TRUE, size = 1, linetype = c("dashed","solid"), palette = c("blue", "red"), legend.labs = c("PRDM1 ≥ 0.13", "PRDM1 < 0.13"), break.time.by = 12, xlim = c(0, 48), risk.table.height = 0.3, ylim = c(0, 100), ggtheme = theme_classic2(base_size = 10, base_family = "Arial"), font.family = "Arial")

# 위험표를 포함한 ggsurvplot 결합
g0 <- arrangeGrob(ggplotGrob(p0$plot), ggplotGrob(p0$table), ncol = 1, heights = c(1.5, 0.5)) # plot과 table 결합
g1 <- arrangeGrob(ggplotGrob(p1$plot), ggplotGrob(p1$table), ncol = 1, heights = c(1.5, 0.5))
g2 <- arrangeGrob(ggplotGrob(p2$plot), ggplotGrob(p2$table), ncol = 1, heights = c(1.5, 0.5))

# ggforest 생성
test2 <- test %>% mutate(TS = Surv(recur_date, recur_m == 1)) %>% 
  rename("Age < 50" = age_m, "EQD2 < 76.25" = td_m, "Non-SqCC" = path_m, "PRDM1 < 0.07" = PRDM1_m, "stage IIIC2-IVB" = stage_m, "miR-150-3p < -0.64" = miR_m, "NMT2 < -0.27" = NMT2_m) 
test2_m <- test2 %>% select(TS, `Age < 50`, `EQD2 < 76.25`, `stage IIIC2-IVB`, `Non-SqCC`, `PRDM1 < 0.07`, `NMT2 < -0.27`, `miR-150-3p < -0.64`) %>% 
  as.data.frame()
out <- mycph(TS ~ ., data = test2_m)
result <- coxph(TS ~ ., data = test2_m)
finalmodel <- step(result, direction = "backward")
summary(finalmodel)
a3 <- ggforest(finalmodel, data = test2_m, fontsize = 0.8, refLabel = "Reference", noDigits = 2) + 
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) # 폰트 크기와 여백 조정


# ggarrange를 사용하여 결합
ggarrange(
  a0, a1, a2, g0, g1, g2, 
  labels = c("A", "B", "C", "D", "E", "F"), 
  nrow = 2, 
  ncol = 3, 
  font.label = list(size = 18, color = "black", face = "bold")
)


mytable(PRDM1_m~path_m+age_m+stage_m+RT_field+ td_m+ alc1+recur1+survival, data=test)
mytable(miR_m~path_m+age_m+stage_m+RT_field+ td_m+ alc1+recur1+survival, data=test)
mytable(var1_m~path_m+age_m+stage_m+RT_field+ td_m+ alc1+recur1+survival, data=test)


