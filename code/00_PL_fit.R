require(poweRlaw)
require(dplyr)
require(tidyr)
require(igraph)
require(ggplot2)
require(magrittr)
category = "Protein Coding"
load("../data/output/graphs.RData")


pk <- function(kcut, kprime, ksat, gamma){
  p_k = 1 / (sum((kprime + ksat)^-gamma * exp (-kprime/kcut))) * (kprime + ksat)^-gamma * exp (-kprime/kcut) 
  return(p_k)
}

values =  c("Protein Coding" = "#ef476f", 
            "TF" = "#f78c6b" ,
            "Pseudogene" = "#06d6a0",
            "miRNA" = "#83d483",
            "other" = "#ffd166",
            "lncRNA" = "#0cb0a9", 
            "ncRNA" = "#118ab2")

PPI_colors = c(PPI = "#8D3B72",
               `PPI & NCI` = "#72E1D1")


PLD = function(g = gPPI, 
               category = "Protein Coding", 
               col = NULL){
  threads = parallel::detectCores()
  message("We detected ", threads, " threads. All will be used.")
  
  if(is.null(col)){
    col = values[names(values) %in% category] 
  }
  
  out = list()
  
  data.dist = igraph::as_data_frame(g, what = "vertices") %>% 
    filter(Category_cl %in% category) %>% 
    group_by(degree) %>% 
    summarise(f = n()) %>% 
    ungroup() %>%
    mutate(p_k = f/sum(f)) %>% 
    mutate(cump_k = cumsum(p_k)) %>% 
    select(k = degree, p_k, cump_k) %>% 
    mutate(S_k = 1 - cump_k + p_k)
  
  out$data = data.dist
  
  data = igraph::as_data_frame(g, what = "vertices") %>% 
    filter(Category_cl %in% category)
  
  out$degree = data
  
  data.s <- data.dist$k
  m_pl <- displ$new(data$degree)
  est_1 = estimate_xmin(m_pl)
  
  message("Defining k_cut and k_sat.")
  #####
  # #generate kmin & kmax pairs
  k_mins = data.s[data.s <= 100]
  k_max = c(data.s[(length(data.s)-2) : length(data.s)])
  out$pairs <- expand.grid(k_mins, k_max)
  message("Grid: ", nrow(out$pairs))
  
  find_ks = function(i){
    require(poweRlaw)
    m_pl$setXmin(out$pairs[i,1])
    est = estimate_xmin(m_pl, 
                        xmins = out$pairs[i,1],
                        xmax = out$pairs[i,2],
                        distance = "ks")
    
    tmp = data.frame(xmin = out$pairs[i,1],
                     xmax = out$pairs[i,2],
                     D = est$gof,
                     gamma = est$pars)
    return(tmp)
  }
  
  k_scan = new.env()
  assign("m_pl", m_pl, envir = k_scan)
  assign("out", out, envir = k_scan)
  
  require(parallel)
  
  cl = parallel::makePSOCKcluster(threads)
  clusterExport(cl, "out", envir = k_scan)
  clusterExport(cl, "m_pl", envir = k_scan)
  out$pairs = clusterApplyLB(cl, 1:nrow(out$pairs), find_ks) %>%
    bind_rows() %>%
    mutate(range = xmax - xmin) %>%
    arrange(D)
  
  stopCluster(cl)
  
  return(out)
}


pars = list()

### PPI 
PLD_PPI = PLD(category = c("Protein Coding", "TF"))
i = 1
PLD_PPI$data$p_k_fit = pk(kcut = PLD_PPI$pairs$xmax[i], 
                             ksat = PLD_PPI$pairs$xmin[i], 
                             gamma = PLD_PPI$pairs$gamma[i],
                             kprime = PLD_PPI$data$k)

pars$PLD_PPI = data.frame(kcut = PLD_PPI$pairs$xmax[i], 
                         ksat = PLD_PPI$pairs$xmin[i], 
                         gamma = PLD_PPI$pairs$gamma[i]) %>% 
  mutate(Network = "PPI", 
         type = "PPI")
PLD_PPI$data %<>%   mutate(Network = "PPI", 
                           type = "PPI")
  
p_PPI = PLD_PPI$data %>%
  ggplot() +
  aes(x = k, y = p_k) +
  geom_point(shape = "circle", size = 1.4, aes(color = type), alpha = 0.3) +
  geom_line(aes(y = p_k_fit), col = "red", linewidth = 1) + 
  scale_x_continuous(trans = "log10") +
  scale_y_continuous(trans = "log10", limits = c(0.00005,0.05)) +
  labs(x = "k", y = "p(k)", title = "PPI") +
  theme_minimal() +
  scale_color_manual(values = PPI_colors) + 
  theme(legend.position = "none") 

PLD_PPI_PC = PLD()
i = 1
PLD_PPI_PC$data$p_k_fit = pk(kcut = PLD_PPI_PC$pairs$xmax[i], 
   ksat = PLD_PPI_PC$pairs$xmin[i], 
   gamma = PLD_PPI_PC$pairs$gamma[i],
   kprime = PLD_PPI_PC$data$k)


pars$PLD_PPI_PC = data.frame(kcut = PLD_PPI_PC$pairs$xmax[i], 
                          ksat = PLD_PPI_PC$pairs$xmin[i], 
                          gamma = PLD_PPI_PC$pairs$gamma[i]) %>% 
  mutate(Network = "PPI", 
         type = "Protein Coding")
PLD_PPI_PC$data %<>%   mutate(Network = "PPI", 
                           type = "Protein Coding")

p_PPI_PC = PLD_PPI_PC$data %>%
  ggplot() +
  aes(x = k, y = p_k) +
  geom_point(shape = "circle", size = 1.4, aes(color = type), alpha = 0.3) +
  geom_line(aes(y = p_k_fit), col = "red", linewidth = 1) + 
  scale_x_continuous(trans = "log10") +
  scale_y_continuous(trans = "log10", limits = c(0.00005,0.05)) +
  labs(x = "k", y = "p(k)", title = "PPI - Protein Coding") +
  theme_minimal() +
  scale_color_manual(values = values) + 
  theme(legend.position = "none") ; p_PPI_PC



PLD_PPI_TF = PLD(g = gPPI, category = "TF")
i = 1
PLD_PPI_TF$data$p_k_fit = pk(kcut = PLD_PPI_TF$pairs$xmax[i], 
                             ksat = PLD_PPI_TF$pairs$xmin[i], 
                             gamma = PLD_PPI_TF$pairs$gamma[i],
                             kprime = PLD_PPI_TF$data$k)

pars$PLD_PPI_TF = data.frame(kcut = PLD_PPI_TF$pairs$xmax[i], 
                             ksat = PLD_PPI_TF$pairs$xmin[i], 
                             gamma = PLD_PPI_TF$pairs$gamma[i]) %>% 
  mutate(Network = "PPI", 
         type = "TF")
PLD_PPI_TF$data %<>%   mutate(Network = "PPI", 
                              type = "TF")

p_PPI_TF = PLD_PPI_TF$data %>%
  ggplot() +
  aes(x = k, y = p_k) +
  geom_point(shape = "circle", size = 1.4, aes(color = type), alpha = 0.3) +
  geom_line(aes(y = p_k_fit), col = "red", linewidth = 1) + 
  scale_x_continuous(trans = "log10") +
  scale_y_continuous(trans = "log10", limits = c(0.0003,0.05)) +
  labs(x = "k", y = "p(k)", title = "PPI - TF") +
  theme_minimal() +
  scale_color_manual(values = values) + 
  theme(legend.position = "none") ; p_PPI_TF

############################
# PPI & NCI 

PLD_PPINCI = PLD(category = unique(V(gPPInc)$Category_cl), g = gPPInc)
i = 1
PLD_PPINCI$data$p_k_fit = pk(kcut = PLD_PPINCI$pairs$xmax[i], 
                          ksat = PLD_PPINCI$pairs$xmin[i], 
                          gamma = PLD_PPINCI$pairs$gamma[i],
                          kprime = PLD_PPINCI$data$k)

pars$PLD_PPINCI = data.frame(kcut = PLD_PPINCI$pairs$xmax[i], 
                          ksat = PLD_PPINCI$pairs$xmin[i], 
                          gamma = PLD_PPINCI$pairs$gamma[i]) %>% 
  mutate(Network = "PPI & NCI", 
         type = "PPI & NCI")
PLD_PPINCI$data %<>%   mutate(Network = "PPI & NCI", 
                           type = "PPI & NCI")

p_PPINCI = PLD_PPINCI$data %>%
  ggplot() +
  aes(x = k, y = p_k) +
  geom_point(shape = "circle", size = 1.4, aes(color = type), alpha = 0.3) +
  geom_line(aes(y = p_k_fit), col = "red", linewidth = 1) + 
  scale_x_continuous(trans = "log10") +
  scale_y_continuous(trans = "log10", limits = c(0.00005,0.05)) +
  labs(x = "k", y = "p(k)", title = "PPI & NCI") +
  theme_minimal() +
  scale_color_manual(values = PPI_colors) + 
  theme(legend.position = "none") ; p_PPINCI


PLD_PPINCI_PC = PLD(g = gPPInc, category = "Protein Coding")
i = 1
PLD_PPINCI_PC$data$p_k_fit = pk(kcut = PLD_PPINCI_PC$pairs$xmax[i], 
                             ksat = PLD_PPINCI_PC$pairs$xmin[i], 
                             gamma = PLD_PPINCI_PC$pairs$gamma[i],
                             kprime = PLD_PPINCI_PC$data$k)

pars$PLD_PPINCI_PC = data.frame(kcut = PLD_PPINCI_PC$pairs$xmax[i], 
                             ksat = PLD_PPINCI_PC$pairs$xmin[i], 
                             gamma = PLD_PPINCI_PC$pairs$gamma[i]) %>% 
  mutate(Network = "PPI & NCI", 
         type = "Protein Coding")
PLD_PPINCI_PC$data %<>% mutate(Network = "PPI & NCI", 
                              type = "Protein Coding")

p_PPINCI_PC = PLD_PPINCI_PC$data %>%
  ggplot() +
  aes(x = k, y = p_k) +
  geom_point(shape = "circle", size = 1.4, aes(color = type), alpha = 0.3) +
  geom_line(aes(y = p_k_fit), col = "red", linewidth = 1) + 
  scale_x_continuous(trans = "log10") +
  scale_y_continuous(trans = "log10", limits = c(0.00005,0.05)) +
  labs(x = "k", y = "p(k)", title = "PPI & NCI - Protein Coding") +
  theme_minimal() +
  scale_color_manual(values = values) + 
  theme(legend.position = "none") ; p_PPINCI_PC



PLD_PPINCI_TF = PLD(g = gPPInc, category = "TF")
i = 3
PLD_PPINCI_TF$data$p_k_fit = pk(kcut = PLD_PPINCI_TF$pairs$xmax[i], 
                                ksat = PLD_PPINCI_TF$pairs$xmin[i], 
                                gamma = PLD_PPINCI_TF$pairs$gamma[i],
                                kprime = PLD_PPINCI_TF$data$k)


pars$PLD_PPINCI_TF = data.frame(kcut = PLD_PPINCI_TF$pairs$xmax[i], 
                                ksat = PLD_PPINCI_TF$pairs$xmin[i], 
                                gamma = PLD_PPINCI_TF$pairs$gamma[i]) %>% 
  mutate(Network = "PPI & NCI", 
         type = "TF")
PLD_PPINCI_TF$data %<>% mutate(Network = "PPI & NCI", 
                               type = "TF")

p_PPINCI_TF = PLD_PPINCI_TF$data %>%
  ggplot() +
  aes(x = k, y = p_k) +
  geom_point(shape = "circle", size = 1.4, aes(color = type), alpha = 0.3) +
  geom_line(aes(y = p_k_fit), col = "red", linewidth = 1) + 
  scale_x_continuous(trans = "log10") +
  scale_y_continuous(trans = "log10", limits = c(0.0005,0.05)) +
  labs(x = "k", y = "p(k)", title = "PPI & NCI - TF") +
  theme_minimal() +
  scale_color_manual(values = values) + 
  theme(legend.position = "none") ; p_PPINCI_TF


PLD_PPINCI_lncRNA = PLD(g = gPPInc, category = "lncRNA")
i = 221

PLD_PPINCI_lncRNA$data$p_k_fit = pk(kcut = PLD_PPINCI_lncRNA$pairs$xmax[i], 
                                ksat = PLD_PPINCI_lncRNA$pairs$xmin[i], 
                                gamma = PLD_PPINCI_lncRNA$pairs$gamma[i],
                                kprime = PLD_PPINCI_lncRNA$data$k)


pars$PLD_PPINCI_lncRNA = data.frame(kcut = PLD_PPINCI_lncRNA$pairs$xmax[i], 
                                ksat = PLD_PPINCI_lncRNA$pairs$xmin[i], 
                                gamma = PLD_PPINCI_lncRNA$pairs$gamma[i]) %>% 
  mutate(Network = "PPI & NCI", 
         type = "lncRNA")
PLD_PPINCI_lncRNA$data %<>% mutate(Network = "PPI & NCI", 
                               type = "lncRNA")

p_PPINCI_lncRNA = PLD_PPINCI_lncRNA$data %>%
  ggplot() +
  aes(x = k, y = p_k) +
  geom_point(shape = "circle", size = 1.4, aes(color = type), alpha = 0.3) +
  geom_line(aes(y = p_k_fit), col = "red", linewidth = 1) + 
  scale_x_continuous(trans = "log10") +
  scale_y_continuous(trans = "log10", limits = c(0.0003,0.1)) +
  labs(x = "k", y = "p(k)", title = "PPI & NCI - lncRNA") +
  theme_minimal() +
  scale_color_manual(values = values) + 
  theme(legend.position = "none") ; p_PPINCI_lncRNA



PLD_PPINCI_miRNA = PLD(g = gPPInc, category = "miRNA")
i = 1
PLD_PPINCI_miRNA$data$p_k_fit = pk(kcut = PLD_PPINCI_miRNA$pairs$xmax[i], 
                                ksat = PLD_PPINCI_miRNA$pairs$xmin[i], 
                                gamma = PLD_PPINCI_miRNA$pairs$gamma[i],
                                kprime = PLD_PPINCI_miRNA$data$k)


pars$PLD_PPINCI_miRNA = data.frame(kcut = PLD_PPINCI_miRNA$pairs$xmax[i], 
                                ksat = PLD_PPINCI_miRNA$pairs$xmin[i], 
                                gamma = PLD_PPINCI_miRNA$pairs$gamma[i]) %>% 
  mutate(Network = "PPI & NCI", 
         type = "miRNA")
PLD_PPINCI_miRNA$data %<>% mutate(Network = "PPI & NCI", 
                               type = "miRNA")

p_PPINCI_miRNA = PLD_PPINCI_miRNA$data %>%
  ggplot() +
  aes(x = k, y = p_k) +
  geom_point(shape = "circle", size = 1.4, aes(color = type), alpha = 0.3) +
  geom_line(aes(y = p_k_fit), col = "red", linewidth = 1) + 
  scale_x_continuous(trans = "log10") +
  scale_y_continuous(trans = "log10", limits = c(0.0005,0.05)) +
  labs(x = "k", y = "p(k)", title = "PPI & NCI - miRNA") +
  theme_minimal() +
  scale_color_manual(values = values) + 
  theme(legend.position = "none") ; p_PPINCI_miRNA


PLD_PPINCI_Pseudogenes = PLD(g = gPPInc, category = "Pseudogene")
i = 1
PLD_PPINCI_Pseudogenes$data$p_k_fit = pk(kcut = PLD_PPINCI_Pseudogenes$pairs$xmax[i], 
                                ksat = PLD_PPINCI_Pseudogenes$pairs$xmin[i], 
                                gamma = PLD_PPINCI_Pseudogenes$pairs$gamma[i],
                                kprime = PLD_PPINCI_Pseudogenes$data$k)


pars$PLD_PPINCI_Pseudogenes = data.frame(kcut = PLD_PPINCI_Pseudogenes$pairs$xmax[i], 
                                ksat = PLD_PPINCI_Pseudogenes$pairs$xmin[i], 
                                gamma = PLD_PPINCI_Pseudogenes$pairs$gamma[i]) %>% 
  mutate(Network = "PPI & NCI", 
         type = "Pseudogenes")
PLD_PPINCI_Pseudogenes$data %<>% mutate(Network = "PPI & NCI", 
                               type = "Pseudogenes")

p_PPINCI_Pseudogenes = PLD_PPINCI_Pseudogenes$data %>%
  ggplot() +
  aes(x = k, y = p_k) +
  geom_point(shape = "circle", size = 1.4, aes(color = type), alpha = 0.3) +
  geom_line(aes(y = p_k_fit), col = "red", linewidth = 1) + 
  scale_x_continuous(trans = "log10") +
  scale_y_continuous(trans = "log10", limits = c(0.0001,0.4)) +
  labs(x = "k", y = "p(k)", title = "PPI & NCI - Pseudogenes") +
  theme_minimal() +
  scale_color_manual(values = values) + 
  theme(legend.position = "none") ; p_PPINCI_Pseudogenes

require(patchwork)
Cairo::CairoPDF("../figs/FigS4.pdf", width = 15, height = 10)
(p_PPI + p_PPI_PC + p_PPI_TF) / (p_PPINCI + p_PPINCI_PC + p_PPINCI_TF + p_PPINCI_lncRNA +  p_PPINCI_miRNA + p_PPINCI_Pseudogenes)
dev.off()

pars %>% 
  bind_rows()
