RequiredPackages <- c(
  
  "data.table",
  "dplyr",
  "ggplot2",
  "magrittr",
  "tidyr",
  "tibble",
  "purrr",
  "igraph",
  "NetSci",
  "parallel"
)


for (i in RequiredPackages) { 
  if (!require(i, character.only = TRUE)) install.packages(i)
}

`%ni%` <- Negate(`%in%`)

OrderNames = function(M){
  n1 = apply(M[, 1:2], 1, min)
  n2 = apply(M[, 1:2], 1, max)
  M[,1]<- n1
  M[,2]<- n2
  
  return(M)
}

load("../data/output/graphs.RData")
set.seed(123)

#### 
#### LCC & Sab
#### 

Hyper = function(j){
  require(magrittr)
  require(dplyr)
  require(NetSci)
  n_B = B %>% 
    filter(Var1 %in% diseases[j] & Var2 %in% diseases[j]) %>%
    pull(value)
  
  n_C = GDA_tmp$hgnc_symbol %>% 
    unique() %>% 
    length()
  hytmp = list()
  for(k in j:length(diseases)){
    n_AB = B %>% 
      filter((Var1 %in% diseases[j] & Var2 %in% diseases[k])|
               (Var1 %in% diseases[k] & Var2 %in% diseases[j])) %>%
      pull(value)
    
    n_A = B %>% 
      filter(Var1 %in% diseases[k] & Var2 %in% diseases[k]) %>%
      pull(value)
    
    htmp = 1 - phyper(n_AB, 
                      n_B, 
                      n_C-n_B, 
                      n_A, 
                      lower.tail = FALSE)
    hytmp[[k]] = data.frame(x = diseases[j], 
                            y = diseases[k], 
                            p = htmp)
  }
  hytmp %<>% bind_rows()
  hytmp$graph = PPIs[i]
  hytmp %<>% OrderNames()
  
  return(Hyper_tmp = hytmp)
}

LCCs_Calc = function(j){
  require(magrittr)
  require(dplyr)
  require(NetSci)
  genes_dis = GDA_tmp %>% 
    filter(DiseaseName == diseases[j]) %>%
    pull(hgnc_symbol)
  
  LCC_aux = LCC_Significance(N = N,
                             Targets = genes_dis,
                             G = x,
                             bins = Bins,
                             min_per_bin = min_per_bin)
  
  LCC_tmp = data.frame(Targets = length(genes_dis),
                       LCC = LCC_aux$LCC,
                       p = LCC_aux$emp_p,
                       Z = LCC_aux$Z,
                       graph = PPIs[i],
                       disease = diseases[j])
  return(LCC = LCC_tmp)
}

Gene_Overlap = list()
LCC = list()

PPIs = c("gPPI",
         "gPPInc" )

N = 1000
Bins = 100
min_per_bin = 1

SAB_JAC = LCCs_Output = SAB_JAC_HYP = list()

for(i in 1:length(PPIs)){
  x = get(PPIs[i])
  message(paste("\n\nStarting Network:", PPIs[i]))
  
  GDA_tmp = GDA %>%
    filter(hgnc_symbol %in% V(x)$name) %>%
    group_by(DiseaseName) %>%
    mutate(n = n()) %>%
    filter(n > 5) 
  
  message("\tDiseases found: ", length(unique(GDA_tmp$DiseaseName)))
  Threads = round(igraph::vcount(x)/ igraph::ecount(x) * length(unique(GDA_tmp$DiseaseName))) + 30
  
  
  diseases = GDA_tmp$DiseaseName %>% unique()
  
  g = GDA_tmp %>% 
    select(DiseaseName , hgnc_symbol) %>%
    graph_from_data_frame(., directed = F) 
  
  igraph::V(g)$type <- igraph::bipartite_mapping(g)$type
  message("\tCalculating Adj.\n")
  A = g %>%
    as_incidence_matrix()
  
  B = A %*% t(A) %>%
    as.matrix() %>% 
    reshape2::melt() %>% 
    OrderNames() %>%
    unique()
  
  LCC_tmp = list()
  
  
  cl = parallel::makeCluster(10)
  clusterExport(cl, "diseases")
  clusterExport(cl, "B")
  clusterExport(cl, "i")
  clusterExport(cl, "GDA_tmp")
  clusterExport(cl, "N")
  clusterExport(cl, "Bins")
  clusterExport(cl, "min_per_bin")
  clusterExport(cl, "PPIs")
  clusterExport(cl, "OrderNames")
  clusterExport(cl, "x")
  
  message("\tCalculating LCC.\n")
  LCC_tmp = clusterApplyLB(cl, 1:length(diseases), LCCs_Calc)
  stopCluster(cl)
  
  LCC_tmp %<>% bind_rows()
  LCCs_Output[[i]] = LCC_tmp 
}

LCCs_Output %<>% bind_rows()

fwrite(LCCs_Output, "../data/output/LCC_100bins.csv")
