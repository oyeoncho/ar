library(tidyverse)
library(ggpubr)
go <- read.table("process/GO_Biological_Process_2023_table.txt", sep="\t", header = T) %>% as_tibble() %>% select(Term, Overlap, P.value)
go150 <- read.table("process/GO_Biological_Process_2023_table_150.txt", sep="\t", header = T) %>% as_tibble() %>% select(Term, Overlap, P.value)
go_prdm1 <- read.table("process/GO_Biological_Process_2023_prdm1.txt", sep="\t", header = T) %>% as_tibble() %>% select(Term, Overlap, P.value)
go_nmt2 <- read.table("process/GO_Biological_Process_2023_nmt2.txt", sep="\t", header = T) %>% as_tibble() %>% select(Term, Overlap, P.value)
go_all <- read.table("process/GO_Biological_Process_2023_table_all.txt", sep="\t", header = T) %>% as_tibble() %>% select(Term, Overlap, P.value)

cell <- read.table("process/HuBMAP_ASCTplusB_augmented_2022_table.txt", sep="\t", header = T) %>% as_tibble() %>% select(Term, Overlap, P.value)
cell150 <- read.table("process/HuBMAP_ASCTplusB_augmented_2022_table_150.txt", sep="\t", header = T) %>% as_tibble() %>% select(Term, Overlap, P.value)
cell_prdm1 <- read.table("process/HuBMAP_ASCTplusB_augmented_2022_prdm1.txt", sep="\t", header = T) %>% as_tibble() %>% select(Term, Overlap, P.value)
cell_nmt2 <- read.table("process/HuBMAP_ASCTplusB_augmented_2022_nmt2.txt", sep="\t", header = T) %>% as_tibble() %>% select(Term, Overlap, P.value)
cellall <- read.table("process/HuBMAP_ASCTplusB_augmented_2022_table_all.txt", sep="\t", header = T) %>% as_tibble() %>% select(Term, Overlap, P.value)


a1 <- go150 %>% slice(1:5) %>% mutate(logp=round(-log(P.value, 10),2)) %>% ggplot(aes(x=reorder(Term, logp), y=logp, fill=Term)) + 
  geom_bar(stat='identity',width=0.8)+coord_flip()+geom_text(aes(label=Term), hjust=1,color="black", size=5)+
  theme(axis.text.y=element_blank(), legend.position = "none",axis.text.x = element_text(size=15), axis.title.x = element_text(size=17,color="red"), axis.title.y = element_text(size=17)
        ,plot.title = element_text(size = 17, face = "bold"))+
  geom_text(aes(label=logp), hjust=-0.1,color="red", size=5, fontface="italic")+
  xlab("Term")+ylab(expression(-log[10](P-value)))+
  geom_hline(yintercept = -log(0.05,10),color="gray", lty=2, lwd=1)+
  ggtitle("Ontologies (GO Biological Process 2023) miR-150-3p-mR network")

a2 <- cell150 %>% slice(1:5) %>% mutate(logp=round(-log(P.value, 10),2)) %>% ggplot(aes(x=reorder(Term, logp), y=logp, fill=Term)) + 
  geom_bar(stat='identity',width=0.8)+coord_flip()+geom_text(aes(label=Term), hjust=1,color="black", size=5)+
  theme(axis.text.y=element_blank(), legend.position = "none",axis.text.x = element_text(size=15), axis.title.x = element_text(size=17,color="red"), axis.title.y = element_text(size=17)
        ,plot.title = element_text(size = 17, face = "bold"))+
  geom_text(aes(label=logp), hjust=-0.3,color="red", size=5, fontface="italic")+
  xlab("Term")+ylab(expression(-log[10](P-value)))+
  geom_hline(yintercept = -log(0.05,10),color="gray", lty=2, lwd=1)+
  ggtitle("Cell types (HuBMAP 2022) miR-150-3p-mR network")

ggarrange(a1,a2, nrow=2, ncol=1)

a1 <- go %>% slice(1:5) %>% mutate(logp=round(-log(P.value, 10),2)) %>% ggplot(aes(x=reorder(Term, logp), y=logp, fill=Term)) + 
  geom_bar(stat='identity',width=0.8)+coord_flip()+geom_text(aes(label=Term), hjust=1,color="black", size=5)+
  theme(axis.text.y=element_blank(), legend.position = "none",axis.text.x = element_text(size=15), axis.title.x = element_text(size=17,color="red"), axis.title.y = element_text(size=17)
        ,plot.title = element_text(size = 17, face = "bold"))+
  geom_text(aes(label=logp), hjust=-0.1,color="red", size=5, fontface="italic")+
  xlab("Term")+ylab(expression(-log[10](P-value)))+
  geom_hline(yintercept = -log(0.05,10),color="gray", lty=2, lwd=1)+
  ggtitle("Ontologies (GO Biological Process 2023) miR-mR network")

a2 <-cell %>% slice(1:5) %>% mutate(logp=round(-log(P.value, 10),2)) %>% ggplot(aes(x=reorder(Term, logp), y=logp, fill=Term)) + 
  geom_bar(stat='identity',width=0.8)+coord_flip()+geom_text(aes(label=Term), hjust=1,color="black", size=5)+
  theme(axis.text.y=element_blank(), legend.position = "none",axis.text.x = element_text(size=15), axis.title.x = element_text(size=17,color="red"), axis.title.y = element_text(size=17)
        ,plot.title = element_text(size = 17, face = "bold"))+
  geom_text(aes(label=logp), hjust=-0.3,color="red", size=5, fontface="italic")+
  xlab("Term")+ylab(expression(-log[10](P-value)))+
  geom_hline(yintercept = -log(0.05,10),color="gray", lty=2, lwd=1)+
  ggtitle("Cell types (HuBMAP 2022) miR-mR network")

ggarrange(a1,a2, nrow=2, ncol=1)

a1 <- go_prdm1 %>% slice(1:5) %>% mutate(logp=round(-log(P.value, 10),2)) %>% ggplot(aes(x=reorder(Term, logp), y=logp, fill=Term)) + 
  geom_bar(stat='identity',width=0.8)+coord_flip()+geom_text(aes(label=Term), hjust=1,color="black", size=5)+
  theme(axis.text.y=element_blank(), legend.position = "none",axis.text.x = element_text(size=15), axis.title.x = element_text(size=17,color="red"), axis.title.y = element_text(size=17)
        ,plot.title = element_text(size = 17, face = "bold"))+
  geom_text(aes(label=logp), hjust=-0.1,color="red", size=5, fontface="italic")+
  xlab("Term")+ylab(expression(-log[10](P-value)))+
  geom_hline(yintercept = -log(0.05,10),color="gray", lty=2, lwd=1)+
  ggtitle("Ontologies (GO Biological Process 2023) PRDM1-mR network")

a2 <- cell_prdm1 %>% slice(1:5) %>% mutate(logp=round(-log(P.value, 10),2)) %>% ggplot(aes(x=reorder(Term, logp), y=logp, fill=Term)) + 
  geom_bar(stat='identity',width=0.8)+coord_flip()+geom_text(aes(label=Term), hjust=1,color="black", size=5)+
  theme(axis.text.y=element_blank(), legend.position = "none",axis.text.x = element_text(size=15), axis.title.x = element_text(size=17,color="red"), axis.title.y = element_text(size=17)
        ,plot.title = element_text(size = 17, face = "bold"))+
  geom_text(aes(label=logp), hjust=-0.3,color="red", size=5, fontface="italic")+
  xlab("Term")+ylab(expression(-log[10](P-value)))+
  geom_hline(yintercept = -log(0.05,10),color="gray", lty=2, lwd=1)+
  ggtitle("Cell types (HuBMAP 2022) PRDM1-mR network")

ggarrange(a1,a2, nrow=2, ncol=1)

a1 <- go_nmt2 %>% slice(1:5) %>% mutate(logp=round(-log(P.value, 10),2)) %>% ggplot(aes(x=reorder(Term, logp), y=logp, fill=Term)) + 
  geom_bar(stat='identity',width=0.8)+coord_flip()+geom_text(aes(label=Term), hjust=1,color="black", size=5)+
  theme(axis.text.y=element_blank(), legend.position = "none",axis.text.x = element_text(size=15), axis.title.x = element_text(size=17,color="red"), axis.title.y = element_text(size=17)
        ,plot.title = element_text(size = 17, face = "bold"))+
  geom_text(aes(label=logp), hjust=-0.1,color="red", size=5, fontface="italic")+
  xlab("Term")+ylab(expression(-log[10](P-value)))+
  geom_hline(yintercept = -log(0.05,10),color="gray", lty=2, lwd=1)+
  ggtitle("Ontologies (GO Biological Process 2023) NMT2-mR network")

a2 <- cell_nmt2 %>% slice(1:5) %>% mutate(logp=round(-log(P.value, 10),2)) %>% ggplot(aes(x=reorder(Term, logp), y=logp, fill=Term)) + 
  geom_bar(stat='identity',width=0.8)+coord_flip()+geom_text(aes(label=Term), hjust=1,color="black", size=5)+
  theme(axis.text.y=element_blank(), legend.position = "none",axis.text.x = element_text(size=15), axis.title.x = element_text(size=17,color="red"), axis.title.y = element_text(size=17)
        ,plot.title = element_text(size = 17, face = "bold"))+
  geom_text(aes(label=logp), hjust=-0.3,color="red", size=5, fontface="italic")+
  xlab("Term")+ylab(expression(-log[10](P-value)))+
  geom_hline(yintercept = -log(0.05,10),color="gray", lty=2, lwd=1)+
  ggtitle("Cell types (HuBMAP 2022) NMT2-mR network")

ggarrange(a1,a2, nrow=2, ncol=1)

a1 <- go_all %>% slice(1:5) %>% mutate(logp=round(-log(P.value, 10),2)) %>% ggplot(aes(x=reorder(Term, logp), y=logp, fill=Term)) + 
  geom_bar(stat='identity',width=0.8)+coord_flip()+geom_text(aes(label=Term), hjust=1,color="black", size=5)+
  theme(axis.text.y=element_blank(), legend.position = "none",axis.text.x = element_text(size=15), axis.title.x = element_text(size=17,color="red"), axis.title.y = element_text(size=17)
        ,plot.title = element_text(size = 17, face = "bold"))+
  geom_text(aes(label=logp), hjust=-0.1,color="red", size=5, fontface="italic")+
  xlab("Term")+ylab(expression(-log[10](P-value)))+
  geom_hline(yintercept = -log(0.05,10),color="gray", lty=2, lwd=1)+
  ggtitle("Ontologies (GO Biological Process 2023) of mRNAs relevant to all of miR-150-3p, NMT2, and PRDM1")

a2 <- cellall %>% slice(1:5) %>% mutate(logp=round(-log(P.value, 10),2)) %>% ggplot(aes(x=reorder(Term, logp), y=logp, fill=Term)) + 
  geom_bar(stat='identity',width=0.8)+coord_flip()+geom_text(aes(label=Term), hjust=1,color="black", size=5)+
  theme(axis.text.y=element_blank(), legend.position = "none",axis.text.x = element_text(size=15), axis.title.x = element_text(size=17,color="red"), axis.title.y = element_text(size=17)
        ,plot.title = element_text(size = 17, face = "bold"))+
  geom_text(aes(label=logp), hjust=-0.3,color="red", size=5, fontface="italic")+
  xlab("Term")+ylab(expression(-log[10](P-value)))+
  geom_hline(yintercept = -log(0.05,10),color="gray", lty=2, lwd=1)+
  ggtitle("Cell types (HuBMAP 2022) of mRNAs relevant to all of miR-150-3p, NMT2, and PRDM1")

ggarrange(a1,a2, nrow=2, ncol=1)
