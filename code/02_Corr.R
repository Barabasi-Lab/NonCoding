require(dplyr)
require(data.table)
require(magrittr)
require(tidyr)
require(progress)
require(igraph)
require(ggplot2)

#####################################
### Clean GTEx to contain only blood samples
### Group transcripts that map to the same gene
### Save it
#####################################

load("../data/output/graphs.RData")
n = Annot_Dic
rm(Annot_Dic)
gtex = fread("../data/input/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz")

gtex_samples = fread("../data/input/GTEx_Annotations_SamplesSubjectsMerged.txt")
Annot_Dic = fread("~/Dropbox (CCNR)/Biology/99_Toolbox/data/gene_annotations/data/out/multiple_hgnc_dict_set_2022-04-04.txt")

ids_blood = gtex_samples %>% 
  select(SAMPID, SMTSD) %>% 
  filter(SMTSD %in% "Whole Blood") %>%
  pull(SAMPID)

n %<>% filter( Symbol %in% V(gPPInc)$name)

xx = rowSums(gtex[,-c(1:2)])
names(xx) = gtex$Name
keep = subset(xx, xx > quantile(xx, 0.1) & xx < quantile(xx, .9))

gtex_blood = gtex %>% 
  filter(Name %in% names(keep)) %>% 
  select(Description, all_of(ids_blood)) %>%
  dplyr::inner_join(Annot_Dic,., by = c("alias" = "Description")) %>%
  dplyr::inner_join(n, ., by = c("Symbol"="HGNC_Symbol")) %>% 
  select(-c(color, alias)) %>% 
  rename(HGNC_Symbol = Symbol) %>% 
  pivot_longer( - c(Category_cl, HGNC_Symbol)) %>%
  filter(HGNC_Symbol %in% V(gPPInc)$name) %>%
  group_by(name, HGNC_Symbol, Category_cl) %>%
  summarise(value = sum(value)) %>%
  ungroup() %>%
  group_by(HGNC_Symbol) %>%
  mutate(sd = sd(value), 
         sum = sum(value), 
         zeros = sum(value == 0)) 

gtex_blood %>% 
  filter(zeros < length(ids_blood)*.80 & sd > 0.01) %>% 
  select(HGNC_Symbol, Category_cl) %>%
  unique() %>% 
  group_by(Category_cl) %>%  
  summarise(n = n())

gtex_blood %>% 
  filter(HGNC_Symbol %in% GDA$hgnc_symbol) %>% 
  filter(zeros < length(ids_blood)*.80 & sd > 0.001) %>% 
  select(HGNC_Symbol, Category_cl) %>%
  unique() %>% 
  group_by(Category_cl) %>%  
  summarise(n = n())

rm(gtex)

gtex_bloodw = gtex_blood %>%
  filter(zeros < length(ids_blood)*.80 & sd > 0.001) %>% 
  ungroup() %>% 
  pivot_wider(names_from = name, 
              values_from = value)

gtex_bloodw %<>% as.data.frame()
row.names(gtex_bloodw) = gtex_bloodw$HGNC_Symbol
gtex_bloodw = gtex_bloodw[,-c(1:5)]

dim(gtex_bloodw)
rm(gtex_blood)

gtex_bloodw %>% 
  fwrite("../data/input/gtex_bloodw.csv.gz", row.names = T)


#####################################
### Calculate Correlation and wTO
#####################################

rm(list = ls())
gtex_bloodw = fread("../data/input/gtex_bloodw.csv.gz") %>% 
  as.data.frame
row.names(gtex_bloodw) = gtex_bloodw$V1
gtex_bloodw = gtex_bloodw  %>% 
  select(-V1)

cor = gtex_bloodw %>% 
  as.matrix() %>% 
  wTO::CorrelationOverlap(., Overlap = row.names(gtex_bloodw), 
                          method = "p")
fwrite(cor, "../data/output/cor_blood.csv")

wto = cor %>% 
  wTO::wTO(., sign = 'sign')
fwrite(wto, "../data/output/WTO_blood.csv")

rm(list = ls())

require(dplyr)
require(data.table)
require(magrittr)
require(tidyr)
require(progress)
require(igraph)
`%ni%` <- Negate(`%in%`)

OrderNames = function(M){
  n1 = apply(M[, 1:2], 1, min)
  n2 = apply(M[, 1:2], 1, max)
  M[,1]<- n1
  M[,2]<- n2
  
  return(M)
}

OrderNames_chop = function(d, y = 2563){
  lines = nrow(d)
  i = lines/y
  
  O = list()
  
  pb <- progress_bar$new(
    format = "(:spin) Computed [:bar] :percent eta: :eta",
    total = y, 
    clear = FALSE,
    width= 60)
  
  for (k in 1:y){
    pb$tick()
    ini = (i*(k-1) + 1)
    end = i*k
    O[[k]] = d[ini:end, ]%>%
      OrderNames()
  }
  O %<>% bind_rows()
  return(O )
}


wto = fread("../data/output/WTO_blood.csv") %>% 
  as.data.frame()
row.names(wto) = names(wto)
M = as.matrix(wto)
M[upper.tri(M, diag = T)] <- NA

M %<>% reshape2::melt()
M %<>% na.exclude()
M %<>% OrderNames_chop()
rm(wto)

cor = fread("../data/output/cor_blood.csv") %>% 
  as.data.frame()
row.names(cor) = names(cor)
C = as.matrix(cor)
C[upper.tri(C, diag = T)] <- NA

C %<>% reshape2::melt()
C %<>% na.exclude()
C %<>% OrderNames_chop()

rm(cor)

load("../data/output/graphs.RData")
diseases = GDA$DiseaseName %>% unique()

PPINC = gPPInc %>% 
  as_data_frame(., "edges") %>%
  select(from, to) %>%
  filter(from != to) %>% 
  unique() %>%
  OrderNames() %>%
  unique() %>%
  mutate(PPINC = 1)

PT = Annot_Dic %>%
  filter(Category_cl %in% c("TF", "Protein Coding" )) %>%
  pull(Symbol)

C %<>% 
  select(from = Var1, 
         to = Var2, 
         cor = value)

MM = M %>%
  select(from = Var1, 
         to = Var2, 
         wto = value) %>% 
  dplyr::left_join(C) %>% 
  dplyr::filter(from %in% V(gPPInc)$name & to %in% V(gPPInc)$name) %>% 
  dplyr::left_join(PPINC) %>%
  mutate(PPINC = ifelse(!is.na(PPINC), 1, 0)) %>% 
  mutate(PCPC = ifelse(from %in% PT & to %in% PT & PPINC == 1, 1, 0))


rm(C)
rm(M)

##################
#### Create the projections
################## 
bip = PPINC %>%
  mutate(type1 = ifelse(from %in% PT, 1, 0)) %>% 
  mutate(type2 = ifelse(to %in% PT, 1, 0)) %>% 
  mutate(type = type1 + type2) %>%
  filter(type == 1)
Annot_Dic$type = ifelse(Annot_Dic$Symbol %in% PT, 1, 0)

gbip = bip %>%
  graph_from_data_frame(., directed = F, vertices = Annot_Dic) %>%
  delete.vertices(., degree(.) == 0)
B = gbip %>% as_adjacency_matrix()
require(Matrix)
xx = B %*% t(B)

PT_B = PT[PT %in% V(gbip)$name]
xx_pt = xx[PT_B, PT_B]
xx_pt[upper.tri(xx_pt, diag = T)]<-NA

A = xx_pt %>% 
  as.matrix() %>%  
  reshape2::melt() %>% 
  na.exclude() 

for(i in 1:10000){
  if(nrow(A)%%i == 0){
    print(i)
  }
}
# y = 2563

A %<>% 
  OrderNames_chop(., y = 1100)

names(A) = c("from", "to", "bip")
MM = full_join(A, MM)
save(MM, 
    file = "../data/output/Correlation.RData", 
    compress = "xz")


