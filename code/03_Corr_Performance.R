library(plotly)
library(UpSetR)
require(tidyverse)
require(data.table)
require(scales)
require(dplyr)
require(ggplot2)
require(ggrepel)
require(magrittr)
require(tidyr)
require(igraph)
require(patchwork)
library(showtext)
require(Cairo)
require(progress)
require(ROCR)
`%ni%` <- Negate(`%in%`)

load("../data/output/graphs.RData")
load("../data/output/Correlation.RData")

MM2$bip = ifelse(is.na(MM2$bip), 0, MM2$bip)
genes_cor = unique(unique(MM2$from), unique(MM2$to))

protein_coding = Annot_Dic %>%
  filter(Category_cl %in% c("Protein Coding", "TF"))


TFs = Annot_Dic %>%
  filter(Category_cl %in% c("TF"))


PC_interactions = MM2 %>%
  filter( from %in% protein_coding$Symbol & 
            to %in% protein_coding$Symbol)


TF_interactions = MM2 %>%
  filter( from %in% TFs$Symbol & 
            to %in% TFs$Symbol)

rm(MM2)
PC_interactions$bip = ifelse(PC_interactions$bip > 0, 1, 0)
PC_interactions$bipppi = ifelse(PC_interactions$bip > 0 | PC_interactions$PCPC > 0, 1, 0)


TF_interactions$bip = ifelse(TF_interactions$bip > 0, 1, 0)
TF_interactions$bipppi = ifelse(TF_interactions$bip > 0 | TF_interactions$PCPC > 0, 1, 0)


binds = list()

binds[[1]] = PC_interactions %>%
  filter(PCPC > 0) %>%
  dplyr::select(from,to,wto,cor) %>%
  mutate(type = "Direct") %>%
  pivot_longer(cols = c( 'wto','cor' ))

binds[[2]] = PC_interactions %>%
  filter(bip > 0) %>%
  dplyr::select(from,to,wto,cor) %>%
  mutate(type = "Indirect") %>%
  pivot_longer(cols = c( 'wto','cor' ))

binds[[3]] = PC_interactions %>%
  filter(bip > 0 | PCPC > 0) %>%
  dplyr::select(from,to,wto,cor) %>%
  mutate(type = "Indirect + Direct") %>%
  pivot_longer(cols = c( 'wto','cor' ))

binds[[4]] = PC_interactions %>%
  filter(bip == 0 | PCPC == 0) %>%
  dplyr::select(from,to,wto,cor) %>%
  mutate(type = "None") %>%
  pivot_longer(cols = c( 'wto','cor' ))

binds %<>% bind_rows()
binds$value = round(binds$value, 3)

KW_wto = binds %>% 
  mutate(type = factor(type, levels = c("None","Direct", "Indirect","Indirect + Direct")))%>% 
  filter(name == "wto") %>% 
  kruskal.test(abs(value) ~ type, .)

Dunn_wto = binds %>% 
  mutate(type = factor(type, levels = c("None","Direct", "Indirect","Indirect + Direct")))%>% 
  filter(name == "wto") %>% 
  FSA::dunnTest(abs(value) ~ type, .)

Dunn_wto$res %>%
  as.data.frame() %>% 
  # mutate(padj = p.adjust(None, method = "fdr")) %>% 
  mutate(KW = KW_wto$p.value) %>% 
  fwrite("../data/output/Dunn_test_wto.csv")

rm(KW_wto)
rm(Dunn_wto)

KW_cor = binds %>% 
  mutate(type = factor(type, levels = c("None","Direct", "Indirect","Indirect + Direct")))%>% 
  filter(name == "cor") %>% 
  kruskal.test(abs(value) ~ type, .)

Dunn_cor = binds %>% 
  mutate(type = factor(type, levels = c("None","Direct", "Indirect","Indirect + Direct")))%>% 
  filter(name == "cor") %>% 
  FSA::dunnTest(abs(value) ~ type, .)

Dunn_cor$res %>%
  as.data.frame() %>% 
  # mutate(padj = p.adjust(None, method = "fdr")) %>% 
  mutate(KW = KW_cor$p.value) %>% 
  fwrite("../data/output//Dunn_test_cor.csv")

rm(KW_cor)
rm(Dunn_cor)

############# TF TF analysis
############# 
############# 
binds_TF = list()
binds_TF[[1]] = TF_interactions %>%
  filter(PCPC > 0) %>%
  dplyr::select(from,to,wto,cor) %>%
  mutate(type = "Direct") %>%
  pivot_longer(cols = c( 'wto','cor' ))

binds_TF[[2]] = TF_interactions %>%
  filter(bip > 0) %>%
  dplyr::select(from,to,wto,cor) %>%
  mutate(type = "Indirect") %>%
  pivot_longer(cols = c( 'wto','cor' ))

binds_TF[[3]] = TF_interactions %>%
  filter(bip > 0 | PCPC > 0) %>%
  dplyr::select(from,to,wto,cor) %>%
  mutate(type = "Indirect + Direct") %>%
  pivot_longer(cols = c( 'wto','cor' ))

binds_TF[[4]] = TF_interactions %>%
  filter(bip == 0 | PCPC == 0) %>%
  dplyr::select(from,to,wto,cor) %>%
  mutate(type = "None") %>%
  pivot_longer(cols = c( 'wto','cor' ))

binds_TF %<>% bind_rows()
binds_TF$value = round(binds_TF$value, 3)

KW_wto_TF = binds_TF %>% 
  mutate(type = factor(type, levels = c("None","Direct", "Indirect","Indirect + Direct")))%>% 
  filter(name == "wto") %>% 
  kruskal.test(abs(value) ~ type, .)

Dunn_wto_TF = binds_TF %>% 
  mutate(type = factor(type, levels = c("None","Direct", "Indirect","Indirect + Direct")))%>% 
  filter(name == "wto") %>% 
  FSA::dunnTest(abs(value) ~ type, .)

Dunn_wto_TF$res %>%
  as.data.frame() %>% 
  # mutate(padj = p.adjust(None, method = "fdr")) %>% 
  mutate(KW = KW_wto_TF$p.value) %>% 
  fwrite("../data/output/Dunn_test_wto_TF.csv")

rm(KW_wto_TF)
rm(Dunn_wto_TF)

KW_cor_TF = binds_TF %>% 
  mutate(type = factor(type, levels = c("None","Direct", "Indirect","Indirect + Direct")))%>% 
  filter(name == "cor") %>% 
  kruskal.test(abs(value) ~ type, .)

Dunn_cor_TF = binds_TF %>% 
  mutate(type = factor(type, levels = c("None","Direct", "Indirect","Indirect + Direct")))%>% 
  filter(name == "cor") %>% 
  FSA::dunnTest(abs(value) ~ type, .)

Dunn_cor_TF$res %>%
  as.data.frame() %>% 
  # mutate(padj = p.adjust(None, method = "fdr")) %>% 
  mutate(KW = KW_cor_TF$p.value) %>% 
  fwrite("../data/output//Dunn_test_cor_TF.csv")

rm(KW_cor_TF)
rm(Dunn_cor_TF)

colors_defined = c("Direct" = "#DBA9CA",
                   
                   "Indirect" = "#BCF1E9",
                   
                   "Indirect + Direct" = "#FAF69E",
                   
                   "None" = "#B8BECC")

Cairo::CairoPDF("../figs/Fig04E_FigS06B.pdf",
                width = 10, height = 8)
p0 = binds %>%
  ggplot() +
  aes(x = abs(value),
      y = type,
      fill = type,
      colour = type) +
  geom_boxplot(alpha = 0.8,
               outlier.size = 0,
               outlier.alpha = 0) +
  scale_fill_manual(values = colors_defined) +
  scale_color_manual(values = colors_defined) +
  coord_flip() +
  theme_minimal() +
  theme(legend.position = "bottom",
        text = element_text(size = 18),
        axis.title  = element_text(face = "bold",
                                   hjust = 0.5),
        axis.title.x.top = element_text(face = "bold",
                                        hjust = 0.5),
        strip.text = element_text(face = "bold",
                                  hjust = 0.5)) +
  facet_wrap(vars(name)) +
  labs (x = "Absolute co-expression weight",
        y = "Interaction type",
        color = "",
        fill = "")
p0
dev.off()

Cairo::CairoPDF("../figs/Fig04E_FigS06B_TF.pdf",
                width = 10, height = 8)
p0_TF = binds_TF %>%
  ggplot() +
  aes(x = abs(value),
      y = type,
      fill = type,
      colour = type) +
  geom_boxplot(alpha = 0.8,
               outlier.size = 0,
               outlier.alpha = 0) +
  scale_fill_manual(values = colors_defined) +
  scale_color_manual(values = colors_defined) +
  coord_flip() +
  theme_minimal() +
  theme(legend.position = "bottom",
        text = element_text(size = 18),
        axis.title  = element_text(face = "bold",
                                   hjust = 0.5),
        axis.title.x.top = element_text(face = "bold",
                                        hjust = 0.5),
        strip.text = element_text(face = "bold",
                                  hjust = 0.5)) +
  facet_wrap(vars(name)) +
  labs (x = "Absolute co-expression weight",
        y = "Interaction type",
        color = "",
        fill = "")
p0_TF
dev.off()


####
pred_PPI_cor = ROCR::prediction(predictions = abs(PC_interactions$cor),
                                labels = PC_interactions$PCPC)

pred_PPI_wto = ROCR::prediction(predictions = abs(PC_interactions$wto),
                                labels = PC_interactions$PCPC)

pred_BIP_cor = ROCR::prediction(predictions = abs(PC_interactions$cor),
                                labels = PC_interactions$bip)

pred_BIP_wto = ROCR::prediction(predictions = abs(PC_interactions$wto),
                                labels = PC_interactions$bip)

pred_BIPPPI_cor = ROCR::prediction(predictions = abs(PC_interactions$cor),
                                   labels = PC_interactions$bipppi)

pred_BIPPPI_wto = ROCR::prediction(predictions = abs(PC_interactions$wto),
                                   labels = PC_interactions$bipppi)

require(pROC)
PC_interactions$abscor = abs(PC_interactions$cor)
PC_interactions$abswto = abs(PC_interactions$wto)

roc_ppi_wto = roc(PC_interactions, PCPC, abswto, ci = T); roc_ppi_wto
roc_ppi_cor = roc(PC_interactions, PCPC, abscor, ci = T); roc_ppi_cor

roc_bip_wto = roc(PC_interactions, bip, abswto, ci = T); roc_bip_wto
roc_bip_cor = roc(PC_interactions, bip, abscor, ci = T); roc_bip_cor

roc_bipppi_wto = roc(PC_interactions, bipppi, abswto, ci = T); roc_bipppi_wto
roc_bipppi_cor = roc(PC_interactions, bipppi, abscor, ci = T); roc_bipppi_cor

get_auc = function(rocs, base = "xxx", measure = 'wto'){
  o = rocs$ci %>% 
    as.vector() %>% 
    t %>%
    as.data.frame()  %>%
    mutate(base = base, 
           measure = measure) 
  names(o)[1:3] = c("CI_l", "AUC", "CI_h")
  
  return(o)
}

perf_sum_rocs = list()
perf_sum_rocs[[1]] = roc_ppi_wto %>%
  get_auc(., base = "Direct", 
          measure = "wTO")

perf_sum_rocs[[2]] = roc_ppi_cor %>%
  get_auc(., base = "Direct", 
          measure = "cor")

perf_sum_rocs[[3]] = roc_bip_wto %>%
  get_auc(., base = "Indirect", 
          measure = "wTO")

perf_sum_rocs[[4]] = roc_bip_cor %>%
  get_auc(., base = "Indirect", 
          measure = "wTO")

perf_sum_rocs[[5]] = roc_bipppi_wto %>%
  get_auc(., base = "Direct + Indirect", 
          measure = "wTO")

perf_sum_rocs[[6]] = roc_bipppi_cor %>%
  get_auc(., base = "Direct + Indirect", 
          measure = "cor")

perf_sum_rocs %<>% bind_rows()



performances = function(data, label){
  roc <- performance(data,"tpr","fpr")
  roc = data.frame(x = as.numeric(roc@x.values[[1]]),
                   y = as.numeric(roc@y.values[[1]]),
                   type = "ROC")

  auc <- performance(data,"auc")
  auc = auc@y.values

  prec_rec <- performance(data, "prec", "rec")
  prec_rec = data.frame(x = as.numeric(prec_rec@x.values[[1]]),
                        y = as.numeric(prec_rec@y.values[[1]]),
                        type = "ROCPR")


  aucpr <- performance(data, "aucpr")
  aucpr = aucpr@y.values

  out =  bind_rows(roc, prec_rec) %>%
    mutate(base = label) %>%
    mutate(AUC = auc) %>%
    mutate(AUCPR = aucpr)

  return(out)
}
perfs = list()

perfs[[1]] = pred_PPI_cor  %>%
  performances(., label = "PPI cor")

perfs[[2]] = pred_PPI_wto  %>%
  performances(., label = "PPI wto")

perfs[[3]] = pred_BIP_cor  %>%
  performances(., label = "BIP cor")

perfs[[4]] = pred_BIP_wto  %>%
  performances(., label = "BIP wto")

perfs[[5]] = pred_BIPPPI_cor  %>%
  performances(., label = "BIPPPI cor")

perfs[[6]] = pred_BIPPPI_wto  %>%
  performances(., label = "BIPPPI wto")

rm(pred_BIP_cor)
rm(pred_BIP_wto)

rm(pred_BIPPPI_cor)
rm(pred_BIPPPI_wto)

rm(pred_PPI_cor)
rm(pred_PPI_wto)

perfs %<>% bind_rows()

fwrite(perfs, "../data/output/Performance.csv")

perfs = fread("../data/output/Performance.csv")

perf_sum = perfs %>%
  select(base, AUC,AUCPR) %>%
  unique() %>% 
  pivot_longer(cols = c("AUC", "AUCPR")) 

tmp = perf_sum$base %>% 
  stringr::str_split(., " ", simplify = T) %>%
  as.data.frame()
perf_sum$base = tmp$V1
perf_sum$measure = tmp$V2


p1_1 = perf_sum %>% 
  filter(name == "AUC") %>% 
  mutate(base = factor(base, 
                       levels = c("PPI", 
                                  "BIP", 
                                  "BIPPPI"), 
                       labels = c("Direct",
                                  "Indirect",
                                  "Indirect + Direct"))) %>% 
  ggplot() +
  aes(x = base, 
      fill = base, 
      colour = base, 
      weight = value) +
  geom_bar() +
  scale_fill_manual(values = colors_defined) +
  scale_color_manual(values = colors_defined) +
  theme_minimal() +
  # ylim(c(0.5,1)) +
  facet_grid(vars(name), vars(measure))+ 
  theme(legend.position = "bottom", 
        text = element_text(size = 18),
        axis.title  = element_text(face = "bold",
                                   hjust = 0.5), 
        axis.title.x.top = element_text(face = "bold", 
                                        hjust = 0.5), 
        strip.text = element_text(face = "bold", 
                                  hjust = 0.5)) +
  coord_cartesian(ylim = c(0.5, 0.7))+
  labs(y = NULL, #"AUC / AUCPR", 
       x = "Base", 
       fill = "Base", 
       color = "Base")



Cairo::CairoPDF("../figs/Fig_4F_Fig_S6B.pdf", 
                width = 10, height = 8)
p1_1
dev.off()

########
######## Performance TFs
######## 


pred_PPI_cor_TF = ROCR::prediction(predictions = abs(TF_interactions$cor),
                                labels = TF_interactions$PCPC)

pred_PPI_wto_TF = ROCR::prediction(predictions = abs(TF_interactions$wto),
                                labels = TF_interactions$PCPC)

pred_BIP_cor_TF = ROCR::prediction(predictions = abs(TF_interactions$cor),
                                labels = TF_interactions$bip)

pred_BIP_wto_TF = ROCR::prediction(predictions = abs(TF_interactions$wto),
                                labels = TF_interactions$bip)

pred_BIPPPI_cor_TF = ROCR::prediction(predictions = abs(TF_interactions$cor),
                                   labels = TF_interactions$bipppi)

pred_BIPPPI_wto_TF = ROCR::prediction(predictions = abs(TF_interactions$wto),
                                   labels = TF_interactions$bipppi)

require(pROC)
TF_interactions$abscor = abs(TF_interactions$cor)
TF_interactions$abswto = abs(TF_interactions$wto)

roc_ppi_wto_TF = roc(TF_interactions, PCPC, abswto, ci = T); roc_ppi_wto_TF
roc_ppi_cor_TF = roc(TF_interactions, PCPC, abscor, ci = T); roc_ppi_cor_TF

roc_bip_wto_TF = roc(TF_interactions, bip, abswto, ci = T); roc_bip_wto_TF
roc_bip_cor_TF = roc(TF_interactions, bip, abscor, ci = T); roc_bip_cor_TF

roc_bipppi_wto_TF = roc(TF_interactions, bipppi, abswto, ci = T); roc_bipppi_wto_TF
roc_bipppi_cor_TF = roc(TF_interactions, bipppi, abscor, ci = T); roc_bipppi_cor_TF

get_auc = function(rocs, base = "xxx", measure = 'wto'){
  o = rocs$ci %>% 
    as.vector() %>% 
    t %>%
    as.data.frame()  %>%
    mutate(base = base, 
           measure = measure) 
  names(o)[1:3] = c("CI_l", "AUC", "CI_h")
  
  return(o)
}

perf_sum_rocs_TF = list()
perf_sum_rocs_TF[[1]] = roc_ppi_wto_TF %>%
  get_auc(., base = "Direct", 
          measure = "wTO")

perf_sum_rocs_TF[[2]] = roc_ppi_cor_TF %>%
  get_auc(., base = "Direct", 
          measure = "cor")

perf_sum_rocs_TF[[3]] = roc_bip_wto_TF %>%
  get_auc(., base = "Indirect", 
          measure = "wTO")

perf_sum_rocs_TF[[4]] = roc_bip_cor_TF %>%
  get_auc(., base = "Indirect", 
          measure = "wTO")

perf_sum_rocs_TF[[5]] = roc_bipppi_wto_TF %>%
  get_auc(., base = "Direct + Indirect", 
          measure = "wTO")

perf_sum_rocs_TF[[6]] = roc_bipppi_cor_TF %>%
  get_auc(., base = "Direct + Indirect", 
          measure = "cor")

perf_sum_rocs_TF %<>% bind_rows()



performances = function(data, label){
  roc <- performance(data,"tpr","fpr")
  roc = data.frame(x = as.numeric(roc@x.values[[1]]),
                   y = as.numeric(roc@y.values[[1]]),
                   type = "ROC")
  
  auc <- performance(data,"auc")
  auc = auc@y.values
  
  prec_rec <- performance(data, "prec", "rec")
  prec_rec = data.frame(x = as.numeric(prec_rec@x.values[[1]]),
                        y = as.numeric(prec_rec@y.values[[1]]),
                        type = "ROCPR")
  
  
  aucpr <- performance(data, "aucpr")
  aucpr = aucpr@y.values
  
  out =  bind_rows(roc, prec_rec) %>%
    mutate(base = label) %>%
    mutate(AUC = auc) %>%
    mutate(AUCPR = aucpr)
  
  return(out)
}
perfs_TF = list()

perfs_TF[[1]] = pred_PPI_cor_TF  %>%
  performances(., label = "PPI cor")

perfs_TF[[2]] = pred_PPI_wto_TF  %>%
  performances(., label = "PPI wto")

perfs_TF[[3]] = pred_BIP_cor_TF  %>%
  performances(., label = "BIP cor")

perfs_TF[[4]] = pred_BIP_wto_TF  %>%
  performances(., label = "BIP wto")

perfs_TF[[5]] = pred_BIPPPI_cor_TF  %>%
  performances(., label = "BIPPPI cor")

perfs_TF[[6]] = pred_BIPPPI_wto_TF  %>%
  performances(., label = "BIPPPI wto")

rm(pred_BIP_cor_TF)
rm(pred_BIP_wto_TF)

rm(pred_BIPPPI_cor_TF)
rm(pred_BIPPPI_wto_TF)

rm(pred_PPI_cor_TF)
rm(pred_PPI_wto_TF)

perfs_TF %<>% bind_rows()

fwrite(perfs_TF, "../data/output/Performance_TFs.csv")

perfs_TF = fread("../data/output/Performance_TFs.csv")

perf_sum_TF = perfs_TF %>%
  select(base, AUC,AUCPR) %>%
  unique() %>% 
  pivot_longer(cols = c("AUC", "AUCPR"), 
               values_transform = list(value = as.numeric)) 

tmp = perf_sum_TF$base %>% 
  stringr::str_split(., " ", simplify = T) %>%
  as.data.frame()
perf_sum_TF$base = tmp$V1
perf_sum_TF$measure = tmp$V2


p1_1_TF = perf_sum_TF %>% 
  filter(name == "AUC") %>% 
  mutate(base = factor(base, 
                       levels = c("PPI", 
                                  "BIP", 
                                  "BIPPPI"), 
                       labels = c("Direct",
                                  "Indirect",
                                  "Indirect + Direct"))) %>% 
  ggplot() +
  aes(x = base, 
      fill = base, 
      colour = base, 
      weight = value) +
  geom_bar() +
  scale_fill_manual(values = colors_defined) +
  scale_color_manual(values = colors_defined) +
  theme_minimal() +
  # ylim(c(0.5,1)) +
  facet_grid(vars(name), vars(measure))+ 
  theme(legend.position = "bottom", 
        text = element_text(size = 18),
        axis.title  = element_text(face = "bold",
                                   hjust = 0.5), 
        axis.title.x.top = element_text(face = "bold", 
                                        hjust = 0.5), 
        strip.text = element_text(face = "bold", 
                                  hjust = 0.5)) +
  coord_cartesian(ylim = c(0.5, 0.7))+
  labs(y = NULL, #"AUC / AUCPR", 
       x = "Base", 
       fill = "Base", 
       color = "Base")



Cairo::CairoPDF("../figs/Fig_4F_Fig_S6B_TF.pdf", 
                width = 10, height = 8)
p1_1_TF
dev.off()
