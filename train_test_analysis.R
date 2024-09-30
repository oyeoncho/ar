library(readxl)
library(dplyr)
library(moonBook)
library(lubridate)
library(tidyverse)
library(edgeR)
library (Hmisc)
library (leaps)
library(xlsx)
library(ggpubr)
library(igraph)
library(ggraph)
library(tidygraph)
library(ggtext)

# analysis

result <- read.csv("dataset/cor.csv", header=T, na.strings=c("-999", "NA", "")) %>% as_tibble()
colnames(result)[2:ncol(result)] <- result$ID
miR_rna <- read_excel(path= "dataset/fc_miR.xlsx", sheet="ar", col_names=TRUE) %>% select(-`...1`)
mR_rna <- read.csv("dataset/fc_mR.csv", header=T, na.strings=c("-999", "NA", "")) %>% as_tibble() %>% 
  mutate(ID=as.character(ID))
colnames(mR_rna) <- gsub("[.]","-", colnames(mR_rna))


### ar ~ miRNA
sel <- result %>% as_tibble() %>% select(ID, ar) %>% filter(abs(ar) > 0.6 & substr(ID,1,4)=="hsa-")
c <- miR_rna %>% select(ar, sel$ID)
colnames(c) <- gsub("hsa-","", colnames(c))

g <-  regsubsets(ar ~ .,
                 data = c,
                 nbest = 1,       # 1 best model for each number of predictors
                 nvmax = NULL,    # NULL for no limit on number of variables
                 force.in = NULL, force.out = NULL,
                 method = "exhaustive")
out <- summary(g)

plot(g, scale = "adjr2", main = "Subgroup selection for acute tumor response")

variables <- apply(out$which[, -1], 1, function(row) paste(names(c)[-1][row], collapse = " + "))
df <- data.frame(
  Model = variables,      # 변수명 조합으로 된 모델
  adjr2= out$adjr2      # 각 모델의 조정된 R-제곱 값
)


a0 <- ggplot(df, aes(x = reorder(Model, adjr2), y = adjr2, fill=Model)) + # 조정된 R-제곱에 따라 모델명 정렬
  geom_bar(stat='identity',width=0.8)+coord_flip()+geom_text(aes(label=Model), hjust=1.3,color="black", size=4) + # 막대 그래프 생성
  geom_text(aes(label = round(adjr2, 3)), hjust = 0.5, size = 4, color="red",fontface="italic") + # 각 막대에 조정된 R-제곱 값 표시
  theme(axis.text.y=element_blank(), legend.position = "none",axis.text.x = element_text(size=13), axis.title.x = element_text(size=13,color="red"), axis.title.y = element_text(size=13)
        ,plot.title = element_text(size = 13, face = "bold"))+
  labs(title = "Subgroup selection for acute tumor response",
       x = "Variable Combination",
       y = "Adjusted R-Squared") + scale_y_continuous(breaks = seq(0, 1, by = 0.1))+
  coord_flip() # x축과 y축을 뒤집어 가독성을 높임

c <- miR_rna %>% select(ar, sel$ID) %>% mutate(`tumor response` = ifelse(ar<0.2 , "good", "poor"))
c$miR <- c$`hsa-miR-150-3p`-c$`hsa-miR-424-3p`
summary(lm(ar~miR, data=c))

a1 <- ggplot(data = c, aes(x = miR, y =ar)) +
  geom_point(size=3)+
  geom_hline(yintercept=0, linetype='dashed', color='red', size=1)+
  geom_smooth(method="lm")+
  labs(x="miR-150-3p-miR-424-3p (log2FC)", y="Acute tumor response")+
  theme(axis.title = element_text(size=15, face="bold"), title=element_text(size=15, face="bold"),
        plot.title=element_text(hjust = 0.5, vjust=-1), axis.text=element_text(size=15, face="bold"))+ 
  annotate("text", x=-1.5, y=0.4, label=expression(R^2==0.7101)) 


p <- ggboxplot(c, x = "tumor response", y = "miR",
               color = "tumor response", palette = c("#0000FF", "#FF0000"),
               add = "jitter", shape = "tumor response")
a2 <- p+ stat_compare_means(label.x = 1.5, label.y = 2, method="wilcox.test", label = "p.signif", size=7)+
    labs(x="Acute tumor response", y="miR-150-3p-miR-424-3p (log2FC)")+theme(axis.title = element_text(size=13, face="bold"),
                                                                           title=element_text(size=13, face="bold"), plot.title=element_text(hjust = 0.5, vjust=-2.5), 
                                                                           axis.text=element_text(size=13, face="bold"))+geom_hline(yintercept = 0, linetype="dashed", color="gray", size=1)+ylim(-5,3)

p <- ggboxplot(c, x = "tumor response", y = "hsa-miR-150-3p",
               color = "tumor response", palette = c("#0000FF", "#FF0000"),
               add = "jitter", shape = "tumor response")
a3 <- p+ stat_compare_means(label.x = 1.5, label.y = 2, method="wilcox.test", label = "p.signif", size=7)+
    labs(x="Acute tumor response", y="miR-150-3p(log2FC)")+theme(axis.title = element_text(size=13, face="bold"),
                                                               title=element_text(size=13, face="bold"), plot.title=element_text(hjust = 0.5, vjust=-2.5), 
                                                               axis.text=element_text(size=13, face="bold"))+geom_hline(yintercept = 0, linetype="dashed", color="red", size=1)+ylim(-5,3)

p <- ggboxplot(c, x = "tumor response", y = "hsa-miR-424-3p",
               color = "tumor response", palette = c("#0000FF", "#FF0000"),
               add = "jitter", shape = "tumor response")
a4 <- p+ stat_compare_means(label.x = 1.5, label.y = 2, method="wilcox.test", label = "p.signif", size=7)+
  labs(x="Acute tumor response", y="miR-424-3p(log2FC)")+theme(axis.title = element_text(size=13, face="bold"),
                                                               title=element_text(size=13, face="bold"), plot.title=element_text(hjust = 0.5, vjust=-2.5), 
                                                               axis.text=element_text(size=13, face="bold"))+geom_hline(yintercept = 0, linetype="dashed", color="red", size=1)+ylim(-5,3)




## miRNA associated miRNA
sel1 <- result %>% as_tibble() %>% select(ID, `hsa-miR-150-3p`, `hsa-miR-424-3p`) %>% 
  filter(abs(`hsa-miR-150-3p`) >0.4 & abs(`hsa-miR-424-3p`) >0.4 & ID !="ar" & substr(ID,1,4)=="hsa-")


# 원본 데이터 전처리
c <- miR_rna %>% select(sel$ID, sel1$ID)
colnames(c) <- gsub("hsa-", "", colnames(c))
c <- data.frame(t(c))

# 그래프 생성
g <- graph.adjacency(
  as.matrix(as.dist(cor(t(c), method = "pearson"))),
  mode = "undirected",
  weighted = TRUE,
  diag = FALSE
)
g <- simplify(g, remove.multiple = TRUE, remove.loops = TRUE)

# 간선 및 정점 속성 설정
E(g)[which(E(g)$weight < 0)]$color <- "darkblue"  # 음의 상관관계: 푸른색
E(g)[which(E(g)$weight > 0)]$color <- "darkred"   # 양의 상관관계: 붉은색
E(g)$weight <- abs(E(g)$weight)
g <- delete_edges(g, E(g)[which(E(g)$weight < 0.4)])
g <- delete.vertices(g, degree(g) == 0)

# 최소 신장 트리 생성 및 군집 탐지
mst <- mst(g, algorithm = "prim")
mst.communities <- edge.betweenness.community(mst, weights = NULL, directed = FALSE)

# tidygraph 객체로 변환하고 군집 멤버십을 factor로 변환
tg <- as_tbl_graph(mst) %>%
  mutate(community = factor(mst.communities$membership))

# 강조할 노드 목록 정의
highlight_nodes <- c("miR-150-3p", "miR-424-3p")

# 정점의 크기와 색상 설정
tg <- tg %>%
  mutate(
    size = ifelse(name %in% highlight_nodes, 10, 5),  # 강조된 노드의 크기 설정
    color = ifelse(name %in% highlight_nodes, "red", as.character(community))  # 강조된 노드 색상 설정
  )

# ggraph를 사용하여 그래프 그리기 (Fruchterman-Reingold 레이아웃 사용)
a5 <- ggraph(tg, layout = "fr") +  # 'fr' 레이아웃 사용
  geom_edge_link(aes(width = E(mst)$weight * 4, color = color),  # 간선 설정
                 alpha = 0.6) + 
  geom_node_point(aes(size = size, fill = color), shape = 21, color = "white") +  # 정점 설정
  geom_node_text(aes(label = name), repel = TRUE, size = 3) +  # 노드 라벨 설정
  scale_edge_width(range = c(0.5, 2)) +  # 간선 너비 범위 조정
  scale_fill_identity() +  # 사용자 지정 색상을 사용
  scale_edge_color_manual(values = c("darkblue", "darkred")) +  # 간선 색상 설정
  theme_void() +  # 배경 및 축 제거
  guides(size = "none", edge_width = "none", edge_color = "none") +
  theme(plot.title = element_text(hjust = 0.5, vjust = -1, size = 13, face = "bold")) +
  ggtitle("miR-miR Network")  # 그래프 제목

central <- data.frame(degree(mst))
central$ID <- rownames(central)
miR_network <- central

ggarrange(a0,a1,a2,a3, a4, a5, labels=c("A","B","C","D","E","F"), nrow=2, ncol=3,font.label = list(size = 16, color = "black", face = "bold"))

miR_network$ID <- paste0("hsa-", miR_network$ID)

########### miR -150 -3p realted mR
mir150 <- result %>% filter(abs(`hsa-miR-150-3p`) > 0.6) %>%  filter(ID !="ar" & substr(ID, 1,4) != "hsa-") %>% select(ID)

c <- miR_rna %>% inner_join(mR_rna, by="ID") %>% select(`hsa-miR-150-3p`, mir150$ID)
colnames(c) <- gsub("hsa-","", colnames(c))
c <- data.frame(t(c))

g <- graph.adjacency(
  as.matrix(as.dist(cor(t(c), method = "pearson"))),
  mode = "undirected",
  weighted = TRUE,
  diag = FALSE
)

g <- simplify(g, remove.multiple = TRUE, remove.loops = TRUE)

# 간선 및 정점 속성 설정
E(g)[which(E(g)$weight < 0)]$color <- "darkblue"
E(g)[which(E(g)$weight > 0)]$color <- "darkred"
E(g)$weight <- abs(E(g)$weight)
g <- delete_edges(g, E(g)[which(E(g)$weight < 0.6)])
g <- delete.vertices(g, degree(g) == 0)

# 최소 신장 트리 생성 및 군집 탐지
mst <- mst(g, algorithm = "prim")
mst.communities <- edge.betweenness.community(mst, weights = NULL, directed = FALSE)

# tidygraph 객체로 변환하고 군집 멤버십을 factor로 변환
tg <- as_tbl_graph(mst) %>%
  mutate(community = factor(mst.communities$membership))

# 강조할 노드 목록 정의
highlight_nodes <- c("miR-150-3p")

# 정점의 크기와 색상 설정
tg <- tg %>%
  mutate(
    size = ifelse(name %in% highlight_nodes, 10, 5),  # 강조된 노드의 크기 설정
    color = ifelse(name %in% highlight_nodes, "red", as.character(community)),  # 강조된 노드 색상 설정
    label_size = ifelse(name %in% highlight_nodes, 10, 5),
    label_bold = ifelse(name %in% highlight_nodes, "bold", "plain")  # 강조된 노드의 라벨 글꼴 설정
  )

# ggraph를 사용하여 그래프 그리기 (Fruchterman-Reingold 레이아웃 사용)
 ggraph(tg, layout = "fr") +  # 'fr' 레이아웃 사용
   geom_edge_link(aes(width = abs(weight) * 4, color = I(E(tg)$color)),  # 간선 설정
                  alpha = 0.6)+ # 간선 설정
  geom_node_point(aes(size = size, fill = color), shape = 21, color = "white") +  # 정점 설정
  geom_node_text(aes(label = name, fontface = label_bold), 
                 repel = TRUE, color = "black", size=3) +  # 라벨 글꼴만 설정
  scale_edge_width(range = c(0.5, 2)) +  # 간선 너비 범위 조정
  scale_fill_identity() +  # 사용자 지정 색상을 사용
  scale_edge_color_manual(values = c("darkblue", "darkred")) +  # 간선 색상 설정
  theme_void() +  # 배경 및 축 제거
  guides(size = "none", edge_width = "none", edge_color = "none") +
  theme(plot.title = element_markdown(hjust = 0.5, vjust=-1, size = 15, face = "bold"))+
   ggtitle("<span style='color:red;'>miR-150-3p</span>-mR Network")  # 그래프 제목

central <- data.frame(degree(mst))
central$ID <- rownames(central)
central$closeness_centrality <- closeness(mst)
central$betweeness_centrality <- betweenness(mst)
central$eigen_centrality <-  eigen_centrality(mst)$vector
central$community <- mst.communities$membership
mir150_network <- central %>% arrange(community) %>% as_tibble()

write.csv(mir150_network, "process/mir150.csv")


## ar~mR +miR-150-3p associated mR + miR
ar_mR <- result %>% filter(abs(ar) > 0.6) %>% select(ID, ar) %>%
  filter(ID !="ar" & substr(ID, 1,4) != "hsa-")

primary_mR <- result %>% filter(abs(`hsa-miR-150-3p`) > 0.6 & abs(ar) > 0.4) %>% select(ID, `hsa-miR-150-3p`, ar) %>% 
  filter(ID !="ar" & substr(ID, 1,4) != "hsa-")


c <- miR_rna %>% inner_join(mR_rna, by="ID") %>% select(miR_network$ID, primary_mR$ID, ar_mR$ID)
colnames(c) <- gsub("hsa-","", colnames(c))
c <- data.frame(t(c))

g <- graph.adjacency(
  as.matrix(as.dist(cor(t(c), method = "pearson"))),
  mode = "undirected",
  weighted = TRUE,
  diag = FALSE
)
g <- simplify(g, remove.multiple = TRUE, remove.loops = TRUE)

# 간선 및 정점 속성 설정
E(g)[which(E(g)$weight < 0)]$color <- "darkblue"
E(g)[which(E(g)$weight > 0)]$color <- "darkred"
E(g)$weight <- abs(E(g)$weight)
g <- delete_edges(g, E(g)[which(E(g)$weight < 0.6)])
g <- delete.vertices(g, degree(g) == 0)

# 최소 신장 트리 생성 및 군집 탐지
mst <- mst(g, algorithm = "prim")
mst.communities <- edge.betweenness.community(mst, weights = NULL, directed = FALSE)

# tidygraph 객체로 변환하고 군집 멤버십을 factor로 변환
tg <- as_tbl_graph(mst) %>%
  mutate(community = factor(mst.communities$membership))

# 강조할 노드 목록 정의
highlight_nodes <- c("FCER2", "PRDM1", "NMT2", "miR-150-3p")

# 정점의 크기와 색상 설정
tg <- tg %>%
  mutate(
    size = ifelse(name %in% highlight_nodes, 12, 5),  # 강조된 노드의 크기 설정
    color = ifelse(name %in% highlight_nodes, "red", as.character(community)),  # 강조된 노드 색상 설정
    label_bold = ifelse(name %in% highlight_nodes, "bold", "plain")  # 강조된 노드의 라벨 글꼴 설정
  )
# ggraph를 사용하여 그래프 그리기 (Fruchterman-Reingold 레이아웃 사용)
ggraph(tg, layout = "fr") +  # 'fr' 레이아웃 사용
  geom_edge_link(aes(width = abs(weight) * 4, color = I(E(tg)$color)),  # 간선 설정
                 alpha = 0.6)+ # 간선 설정
  geom_node_point(aes(size = size, fill = color), shape = 21, color = "white") +  # 정점 설정
  geom_node_text(aes(label = name, fontface = label_bold), 
                 repel = TRUE, color = "black", size = 3) +  # 라벨 글꼴만 설정
  scale_edge_width(range = c(0.5, 2)) +  # 간선 너비 범위 조정
  scale_fill_identity() +  # 사용자 지정 색상을 사용
  scale_edge_color_manual(values = c("darkblue", "darkred")) +  # 간선 색상 설정
  theme_void() +  # 배경 및 축 제거
  guides(size = "none", edge_width = "none", edge_color = "none") +
  theme(plot.title = element_markdown(hjust = 0.5, vjust=-1, size = 15, face = "bold"))+
  ggtitle("miR-mR Network <span style='color:red;'> (lymphocyte </span> <span style='color:red;'> enriched </span> <span style='color:red;'> RNAs) </span>")  # 그래프 제목


central <- data.frame(degree(mst))
central$ID <- rownames(central)
central$closeness_centrality <- closeness(mst)
central$betweeness_centrality <- betweenness(mst)
central$eigen_centrality <-  eigen_centrality(mst)$vector
central$community <- mst.communities$membership
mR_network <- central %>% arrange(community) %>% as_tibble()
mR_network %>% filter(substr(ID,1,3)=="miR") %>% arrange(desc(closeness_centrality))
write.csv(mR_network, "process/mR_network.csv")


mR_can1 <- result %>% filter(abs(ar) > 0.7) %>% select(ID, ar) %>%
  filter(ID !="ar" & substr(ID, 1,4) != "hsa-")
mR_can2 <-result %>% filter(abs(`hsa-miR-150-3p`) > 0.7 & abs(ar)>0.4) %>% select(ID, `hsa-miR-150-3p`, ar) %>% 
  filter(ID !="ar" & substr(ID, 1,4) != "hsa-")
write.xlsx(mR_can1, file="process/candidate.xlsx", sheetName="ar",append=TRUE)
write.xlsx(mR_can2, file="process/candidate.xlsx", sheetName="150",append=TRUE)

mR_can1


# 데이터 전처리
c <- miR_rna %>% inner_join(mR_rna, by = "ID") %>% 
  select(FCER2, PRDM1, NMT2, `hsa-miR-150-3p`, `hsa-miR-424-3p`)
colnames(c) <- gsub("hsa-", "", colnames(c))
c <- data.frame(t(c))

# 그래프 생성
g <- graph.adjacency(
  as.matrix(as.dist(cor(t(c), method = "pearson"))),
  mode = "undirected",
  weighted = TRUE,
  diag = FALSE
)
g <- simplify(g, remove.multiple = TRUE, remove.loops = TRUE)

# 간선 및 정점 속성 설정
E(g)$color <- ifelse(E(g)$weight < 0, "darkblue", "darkred")  # 음의 상관관계: 푸른색, 양의 상관관계: 붉은색
E(g)$weight <- abs(E(g)$weight)
g <- delete_edges(g, E(g)[which(E(g)$weight < 0.6)])
g <- delete.vertices(g, degree(g) == 0)

# 최소 신장 트리 생성 및 군집 탐지
mst <- mst(g, algorithm = "prim")
mst.communities <- edge.betweenness.community(mst, weights = NULL, directed = FALSE)

# tidygraph 객체로 변환하고 군집 멤버십을 factor로 변환
tg <- as_tbl_graph(mst) %>%
  mutate(community = factor(mst.communities$membership))

# 강조할 노드 목록 정의
highlight_nodes <- c("miR-150-3p", "miR-424-3p")

# 정점의 크기와 색상 설정
tg <- tg %>%
  mutate(
    size = ifelse(name %in% highlight_nodes, 10, 5),  # 강조된 노드의 크기 설정
    color = ifelse(name %in% highlight_nodes, "red", as.character(community))  # 강조된 노드 색상 설정
  )

# ggraph를 사용하여 그래프 그리기 (Fruchterman-Reingold 레이아웃 사용)
b0 <- ggraph(tg, layout = "fr") +  # 'fr' 레이아웃 사용
  geom_edge_link(aes(width = abs(weight) * 4, color = I(E(tg)$color)),  # 간선 설정
                 alpha = 0.6) + 
  geom_node_point(aes(size = size, fill = color), shape = 21, color = "white") +  # 정점 설정
  geom_node_text(aes(label = name), repel = TRUE, size = 5) +  # 노드 라벨 설정
  scale_edge_width(range = c(0.5, 2)) +  # 간선 너비 범위 조정
  scale_fill_identity() +  # 사용자 지정 색상을 사용
  theme_void() +  # 배경 및 축 제거
  guides(size = "none", edge_width = "none", edge_color = "none") +
  theme(plot.title = element_text(hjust = 0.5, vjust = -1, size = 17, face = "bold")) +
  ggtitle("Simplified Network")  # 그래프 제목

central <- data.frame(degree(mst))
central$ID <- rownames(central)
network <- central %>% arrange(ID)
write.csv(network, "process/network.csv")

network$ID <- ifelse(substr(network$ID,1,4) =="miR-", paste0("hsa-", network$ID), network$ID)

##
c <- miR_rna %>% inner_join(mR_rna, by="ID") %>% select(ar, network$ID)
colnames(c) <- gsub("hsa-","", colnames(c))

g <-  regsubsets(ar ~.,
                 data = c,
                 nbest = 1,       # 1 best model for each number of predictors
                 nvmax = NULL,    # NULL for no limit on number of variables
                 force.in = NULL, force.out = NULL,
                 method = "exhaustive")
out <- summary(g)

plot(g, scale = "adjr2", main = "Subgroup selection for acute tumor response")


variables <- apply(out$which[, -1], 1, function(row) paste(names(c)[-1][row], collapse = " + "))
df <- data.frame(
  Model = variables,      # 변수명 조합으로 된 모델
  adjr2= out$adjr2      # 각 모델의 조정된 R-제곱 값
)


b1 <- ggplot(df, aes(x = reorder(Model, adjr2), y = adjr2, fill=Model)) + # 조정된 R-제곱에 따라 모델명 정렬
  geom_bar(stat='identity',width=0.8)+coord_flip()+geom_text(aes(label=Model), hjust=1.5,color="black", size=4) + # 막대 그래프 생성
  geom_text(aes(label = round(adjr2, 3)), hjust = 0.5, size = 4, color="red",fontface="italic") + # 각 막대에 조정된 R-제곱 값 표시
  theme(axis.text.y=element_blank(), legend.position = "none",axis.text.x = element_text(size=13), axis.title.x = element_text(size=13,color="red"), axis.title.y = element_text(size=13)
        ,plot.title = element_text(size = 13, face = "bold"))+
  labs(title = "Subgroup selection for acute tumor response",
       x = "Variable Combination",
       y = "Adjusted R-Squared") + scale_y_continuous(breaks = seq(0, 1, by = 0.1))+
  coord_flip() # x축과 y축을 뒤집어 가독성을 높임

c <- miR_rna %>% inner_join(mR_rna, by="ID") %>% select(ar, network$ID) %>% mutate(`tumor response` = ifelse(ar<0.2 , "good", "poor"))
c$miR <- c$`hsa-miR-150-3p`+c$NMT2+c$PRDM1
summary(lm(ar~miR, data=c))

b2 <- ggplot(data = c, aes(x = miR, y =ar)) +
  geom_point(size=3)+
  geom_hline(yintercept=0, linetype='dashed', color='red', size=1)+
  geom_smooth(method="lm")+
  labs(x="miR-150-3p+PRDM1+NMT2(log2FC)", y="Acute tumor response")+
  theme(axis.title = element_text(size=15, face="bold"), title=element_text(size=15, face="bold"),
        plot.title=element_text(hjust = 0.5, vjust=-1), axis.text=element_text(size=15, face="bold"))+ 
  annotate("text", x=-1.5, y=0.4, label=expression(R^2==0.8305)) +
  scale_x_continuous(limits = c(-7, 3))


p <- ggboxplot(c, x = "tumor response", y = "miR",
               color = "tumor response", palette = c("#0000FF", "#FF0000"),
               add = "jitter", shape = "tumor response")
b3 <- p+ stat_compare_means(label.x = 1.5, label.y = 4, method="wilcox.test", label = "p.signif", size=7) + 
  labs(x="Acute tumor response", y="miR-150-3p+PRDM1+NMT2(log2FC)")+
  theme(axis.title = element_text(size=13, face="bold"), title=element_text(size=13, face="bold"), plot.title=element_text(hjust = 0.5, vjust=-2.5),
        axis.text=element_text(size=13, face="bold"))+geom_hline(yintercept = 0, linetype="dashed", color="gray", size=1)+
  scale_y_continuous(limits = c(-7, 6))
  

ggarrange(b0, b1, b2, b3, labels=c("A","B","C","D"), nrow=2, ncol=2,font.label = list(size = 16, color = "black", face = "bold"))

####

########### prdm1 realted mR
prdm1 <- result %>% filter(abs(PRDM1) > 0.6) %>%  filter(ID !="ar" & substr(ID, 1,4) != "hsa-") %>% select(ID)

c <- miR_rna %>% inner_join(mR_rna, by="ID") %>% select(PRDM1, prdm1$ID)
colnames(c) <- gsub("hsa-","", colnames(c))
c <- data.frame(t(c))

g <- graph.adjacency(
  as.matrix(as.dist(cor(t(c), method = "pearson"))),
  mode = "undirected",
  weighted = TRUE,
  diag = FALSE
)

g <- simplify(g, remove.multiple = TRUE, remove.loops = TRUE)

# 간선 및 정점 속성 설정
E(g)[which(E(g)$weight < 0)]$color <- "darkblue"
E(g)[which(E(g)$weight > 0)]$color <- "darkred"
E(g)$weight <- abs(E(g)$weight)
g <- delete_edges(g, E(g)[which(E(g)$weight < 0.6)])
g <- delete.vertices(g, degree(g) == 0)

# 최소 신장 트리 생성 및 군집 탐지
mst <- mst(g, algorithm = "prim")
mst.communities <- edge.betweenness.community(mst, weights = NULL, directed = FALSE)

# tidygraph 객체로 변환하고 군집 멤버십을 factor로 변환
tg <- as_tbl_graph(mst) %>%
  mutate(community = factor(mst.communities$membership))

# 강조할 노드 목록 정의
highlight_nodes <- c("PRDM1")

# 정점의 크기와 색상 설정
tg <- tg %>%
  mutate(
    size = ifelse(name %in% highlight_nodes, 10, 5),  # 강조된 노드의 크기 설정
    color = ifelse(name %in% highlight_nodes, "red", as.character(community)),  # 강조된 노드 색상 설정
    label_size = ifelse(name %in% highlight_nodes, 10, 5),
    label_bold = ifelse(name %in% highlight_nodes, "bold", "plain")  # 강조된 노드의 라벨 글꼴 설정
  )

# ggraph를 사용하여 그래프 그리기 (Fruchterman-Reingold 레이아웃 사용)
ggraph(tg, layout = "fr") +  # 'fr' 레이아웃 사용
  geom_edge_link(aes(width = abs(weight) * 4, color = I(E(tg)$color)),  # 간선 설정
                 alpha = 0.6)+  # 간선 설정
  geom_node_point(aes(size = size, fill = color), shape = 21, color = "white") +  # 정점 설정
  geom_node_text(aes(label = name, fontface = label_bold), 
                 repel = TRUE, color = "black", size=3) +  # 라벨 글꼴만 설정
  scale_edge_width(range = c(0.5, 2)) +  # 간선 너비 범위 조정
  scale_fill_identity() +  # 사용자 지정 색상을 사용
  scale_edge_color_manual(values = c("darkblue", "darkred")) +  # 간선 색상 설정
  theme_void() +  # 배경 및 축 제거
  guides(size = "none", edge_width = "none", edge_color = "none") +
  theme(plot.title = element_markdown(hjust = 0.5, vjust=-1, size = 15, face = "bold"))+
  ggtitle("<span style='color:red;'>PRDM1</span>-mR Network")  # 그래프 제목

central <- data.frame(degree(mst))
central$ID <- rownames(central)
central$closeness_centrality <- closeness(mst)
central$betweeness_centrality <- betweenness(mst)
central$eigen_centrality <-  eigen_centrality(mst)$vector
central$community <- mst.communities$membership
prdm1_network <- central %>% arrange(community) %>% as_tibble()

write.csv(prdm1_network, "process/prdm1.csv")

########### nmt2 realted mR
nmt2 <- result %>% filter(abs(NMT2) > 0.6) %>%  filter(ID !="ar" & substr(ID, 1,4) != "hsa-") %>% select(ID)

c <- miR_rna %>% inner_join(mR_rna, by="ID") %>% select(NMT2, nmt2$ID)
colnames(c) <- gsub("hsa-","", colnames(c))
c <- data.frame(t(c))

g <- graph.adjacency(
  as.matrix(as.dist(cor(t(c), method = "pearson"))),
  mode = "undirected",
  weighted = TRUE,
  diag = FALSE
)

g <- simplify(g, remove.multiple = TRUE, remove.loops = TRUE)

# 간선 및 정점 속성 설정
E(g)[which(E(g)$weight < 0)]$color <- "darkblue"
E(g)[which(E(g)$weight > 0)]$color <- "darkred"
E(g)$weight <- abs(E(g)$weight)
g <- delete_edges(g, E(g)[which(E(g)$weight < 0.6)])
g <- delete.vertices(g, degree(g) == 0)

# 최소 신장 트리 생성 및 군집 탐지
mst <- mst(g, algorithm = "prim")
mst.communities <- edge.betweenness.community(mst, weights = NULL, directed = FALSE)

# tidygraph 객체로 변환하고 군집 멤버십을 factor로 변환
tg <- as_tbl_graph(mst) %>%
  mutate(community = factor(mst.communities$membership))

# 강조할 노드 목록 정의
highlight_nodes <- c("NMT2")

# 정점의 크기와 색상 설정
tg <- tg %>%
  mutate(
    size = ifelse(name %in% highlight_nodes, 10, 5),  # 강조된 노드의 크기 설정
    color = ifelse(name %in% highlight_nodes, "red", as.character(community)),  # 강조된 노드 색상 설정
    label_size = ifelse(name %in% highlight_nodes, 10, 5),
    label_bold = ifelse(name %in% highlight_nodes, "bold", "plain")  # 강조된 노드의 라벨 글꼴 설정
  )

# ggraph를 사용하여 그래프 그리기 (Fruchterman-Reingold 레이아웃 사용)
ggraph(tg, layout = "fr") +  # 'fr' 레이아웃 사용
  geom_edge_link(aes(width = abs(weight) * 4, color = I(E(tg)$color)),  # 간선 설정
                 alpha = 0.6)+  # 간선 설정
  geom_node_point(aes(size = size, fill = color), shape = 21, color = "white") +  # 정점 설정
  geom_node_text(aes(label = name, fontface = label_bold), 
                 repel = TRUE, color = "black", size=3) +  # 라벨 글꼴만 설정
  scale_edge_width(range = c(0.5, 2)) +  # 간선 너비 범위 조정
  scale_fill_identity() +  # 사용자 지정 색상을 사용
  scale_edge_color_manual(values = c("darkblue", "darkred")) +  # 간선 색상 설정
  theme_void() +  # 배경 및 축 제거
  guides(size = "none", edge_width = "none", edge_color = "none") +
  theme(plot.title = element_markdown(hjust = 0.5, vjust=-1, size = 15, face = "bold"))+
  ggtitle("<span style='color:red;'>NMT2</span>-mR Network")  # 그래프 제목

central <- data.frame(degree(mst))
central$ID <- rownames(central)
central$closeness_centrality <- closeness(mst)
central$betweeness_centrality <- betweenness(mst)
central$eigen_centrality <-  eigen_centrality(mst)$vector
central$community <- mst.communities$membership
nmt2_network <- central %>% arrange(community) %>% as_tibble()

write.csv(nmt2_network, "process/nmt2.csv")


#########
all100 <- result %>% filter(abs(`hsa-miR-150-3p`) > 0 & abs(PRDM1) > 0 & abs(NMT2) > 0) %>%  filter(ID !="ar" & substr(ID, 1,4) != "hsa-") %>% select(ID, `hsa-miR-150-3p`, PRDM1, NMT2)
write.xlsx(all100, file="process/sel_network.xlsx", sheetName="all",append=TRUE)


####

cx <- readxl::read_excel(path= "dataset/cx.xlsx", sheet="cx_test", col_names=TRUE)

cx <- cx %>% mutate( ar_m=ifelse(ar<0.2,0,1)) %>% select(ID, age, pathology, stage, ar,tv0, tv2,  stage, scc0:alc2, RT_field, ICR, fraction_size, EQD2, ar_m)

test <- miR_rna %>% inner_join(mR_rna, by="ID") %>% select(ID, `hsa-miR-150-3p`, NMT2, PRDM1) 


test1 <-  cx %>% inner_join(test, by="ID") %>% 
  mutate(total = PRDM1+`hsa-miR-150-3p`, total1 = NMT2+PRDM1+`hsa-miR-150-3p`, miR=`hsa-miR-150-3p`,
         EQD2_m = ifelse(EQD2<median(EQD2), 1, 0), ICR_m = ifelse(is.na(ICR)==TRUE, "N/A", 
                                                                  ifelse(fraction_size<4, "EBRT", paste0(ICR,"/",ICR/fraction_size))),
         Age_m=ifelse(age<50,1,0), stage1=ifelse(stage=="3C2"|stage=="4A"|stage=="4B","IIIC2-IVB", 
                                                          ifelse(stage=="1B", "IB", "IIIC1")))

out = mytable(ar_m ~ tv0+tv2+ar+Age_m+pathology+stage1+RT_field+EQD2_m+ICR_m+alc0+alc1+alc2+miR+PRDM1+NMT2+total+total1, data=test1, method=2, digits=2)

mycsv(out, file="table/table_test.csv")

out1 = mytable( ~ tv0+tv2+ar+Age_m+pathology+stage1+RT_field+EQD2_m+ICR_m+alc0+alc1+alc2+miR+PRDM1+NMT2+total+total1, data=test1, method=2, digits=2)

mycsv(out1, file="table/table_all.csv")

out2 = mytable( ~ tv0+tv2+ar+age+EQD2+alc0+alc1+alc2+miR+PRDM1+NMT2+total+total1, data=test1, method=2, digits=2)

mycsv(out2, file="table/median_test.csv")
