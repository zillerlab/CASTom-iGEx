options(stringsAsFactors=F)
options(max.print=1000)
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(qvalue))
suppressPackageStartupMessages(library(pROC))
suppressPackageStartupMessages(library(pryr))
suppressPackageStartupMessages(library(umap))
suppressPackageStartupMessages(library(igraph))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(SparseM))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(gep2pep))
options(bitmapType = 'cairo', device = 'png')

setwd('OUTPUT_GTEx/predict_CAD/Liver/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/CAD_HARD_clustering/update_corrPCs/')
download.file("http://dsea.tigem.it/data/Cmap_MSigDB_v6.1_PEPs.tar.gz","Cmap_MSigDB_v6.1_PEPs.tar.gz")
untar("Cmap_MSigDB_v6.1_PEPs.tar.gz")
rpBig <- openRepository("Cmap_MSigDB_v6.1_PEPs")

Reactome_coll <- loadCollection(rpBig,"C2_CP:REACTOME")
GO_MF_coll <- loadCollection(rpBig,"C5_MF")
GO_BP_coll <- loadCollection(rpBig,"C5_BP")
GO_CC_coll <- loadCollection(rpBig,"C5_CC")
drug_atc <- read.csv("/psycl/g/mpsziller/lucia/drug_targeting/WHO ATC-DDD 2021-12-03.csv")
# load CAD results:
path_CAD <- get(load('pathOriginal_filtJS0.2_corrPCs_tscoreClusterCases_featAssociation.RData'))
# liver_res <- path_CAD$test_feat[[10]] %>% filter(pval_corr <= 0.01)
tot_res <- do.call(rbind, path_CAD$test_feat) %>% filter(pval_corr <= 0.01)
# gr_res <- lapply(1:5, function(x) liver_res[liver_res$comp == sprintf('gr%i_vs_all', x),])
gr_res <- lapply(1:5, function(x) tot_res[tot_res$comp == sprintf('gr%i_vs_all', x),])
# remove discordant results in sign
for(i in 1:5){
  tmp <- gr_res[[i]]
  dup_path <- names(which(table(tmp$feat) > 1))
  if(length(dup_path)>0){
    rm_path <- c()
    for(j in 1:length(dup_path)){
      tmp_path <- tmp %>% filter(feat == dup_path[j])
      if(!(all(tmp_path$estimates > 0) | all(tmp_path$estimates < 0))){
        rm_path <- c(rm_path, dup_path[j])
      }
    }
    gr_res[[i]] <- gr_res[[i]][!gr_res[[i]]$feat %in% rm_path,]
  }
}


# Reactome_names <- path_CAD$res_pval[[10]]$path[path_CAD$res_pval[[10]]$name == 'Reactome']
Reactome_names <- do.call(rbind, path_CAD$res_pval) %>% filter(name == 'Reactome')
Reactome_names <- Reactome_names$path 
Reactome_names <- Reactome_names$path[!duplicated(Reactome_names$path)]

df_names_R <- data.frame(original = Reactome_names)
# filter
new_names <- Reactome_names %>% toupper() %>% 
  str_replace_all(pattern =  "-",  replacement = "_") %>%
  str_replace_all(pattern =  "/",  replacement = "_") %>%
  str_replace_all(pattern =  " ",  replacement = "_") %>% 
  str_replace_all(pattern =  "[(]",  replacement = "") %>% 
  str_replace_all(pattern =  "[)]",  replacement = "") %>% 
  str_replace_all(pattern =  "___",  replacement = "_")
new_names <- str_c("REACTOME_", new_names)
df_names_R$updated <- new_names

##
GO_names <-  do.call(rbind, path_CAD$res_pval) %>% filter(name == 'GO')
GO_names<- GO_names$path[!duplicated(GO_names$path)]

df_names_GO <- data.frame(original = GO_names)
# filter
new_names <- GO_names %>% toupper() %>% 
  str_replace_all(pattern =  "-",  replacement = "_") %>%
  str_replace_all(pattern =  "/",  replacement = "_") %>%
  str_replace_all(pattern =  " ",  replacement = "_") %>% 
  str_replace_all(pattern =  "[(]",  replacement = "") %>% 
  str_replace_all(pattern =  "[)]",  replacement = "") %>% 
  str_replace_all(pattern =  "___",  replacement = "_")
new_names <- str_c("GO_", new_names)
df_names_GO$updated <- new_names

##################
#### function ####
##################

psea_group <- function(df_names, collection, coll_name){
  
  psea_out <- list()
  for(i in 1:length(gr_res)){
    print(i)
    tmp <- df_names[df_names$original %in% gr_res[[i]]$feat[gr_res[[i]]$estimates < 0],]
    id <- which(sapply(collection, function(x) x@setName) %in% tmp$updated)
    db_filt <- collection[id]
    if(length(id)>0){
      psea <- PathSEA(rpBig, db_filt, collections=c(coll_name))
      psea_out[[i]] <- psea$PathSEA[[coll_name]]
      psea_out[[i]]$gr <- i
      psea_out[[i]] <- psea_out[[i]] %>% 
        mutate(drug = rownames(psea_out[[i]]), sign_p = -log10(PV) * sign(ES)) %>% 
        arrange(desc(sign_p)) %>% filter(sign_p > 0) %>%  filter(PV <= 0.01)
      psea_out[[i]]$atc_code <- drug_atc$atc_code[match(psea_out[[i]]$drug, drug_atc$atc_name)]
    }
  }
  psea_neg <- do.call(rbind, psea_out) 
  
  psea_out <- list()
  for(i in 1:length(gr_res)){
    print(i)
    tmp <- df_names[df_names$original %in% gr_res[[i]]$feat[gr_res[[i]]$estimates > 0],]
    id <- which(sapply(collection, function(x) x@setName) %in% tmp$updated)
    db_filt <- collection[id]
    if(length(id)>0){
      psea <- PathSEA(rpBig, db_filt, collections=c(coll_name))
      psea_out[[i]] <- psea$PathSEA[[coll_name]]
      psea_out[[i]]$gr <- i
      psea_out[[i]] <- psea_out[[i]] %>% 
        mutate(drug = rownames(psea_out[[i]]), sign_p = -log10(PV) * sign(ES)) %>% 
        arrange(sign_p)  %>% filter(sign_p < 0)  %>%  filter(PV <= 0.01) 
      psea_out[[i]]$atc_code <- drug_atc$atc_code[match(psea_out[[i]]$drug, drug_atc$atc_name)]
    }
  }
  psea_pos <- do.call(rbind, psea_out)
  
  return(list(neg = psea_neg, pos = psea_pos))
}

###################################################################

Reactome_out <- psea_group(df_names_R, collection = Reactome_coll, coll_name = "C2_CP:REACTOME")
Reactome_out$neg %>% filter(!is.na(atc_code)) %>% filter(str_starts(atc_code, "C")) %>% filter(PV <= 10^-3)
Reactome_out$pos %>% filter(!is.na(atc_code)) %>% filter(str_starts(atc_code, "C")) %>% filter(PV <= 10^-3)

GO_MF_out <- psea_group(df_names_GO, collection = GO_MF_coll, coll_name = "C5_MF")
GO_MF_out$neg %>% filter(!is.na(atc_code)) %>% filter(str_starts(atc_code, "C")) %>% filter(PV <= 10^-3)
GO_MF_out$pos %>% filter(!is.na(atc_code)) %>% filter(str_starts(atc_code, "C")) %>% filter(PV <= 10^-3)

GO_BP_out <- psea_group(df_names_GO, collection = GO_BP_coll, coll_name = "C5_BP")
GO_BP_out$neg %>% filter(!is.na(atc_code)) %>% filter(str_starts(atc_code, "C")) %>% filter(PV <= 10^-3)
GO_BP_out$pos %>% filter(!is.na(atc_code)) %>% filter(str_starts(atc_code, "C")) %>% filter(PV <= 10^-3)

GO_CC_out <- psea_group(df_names_GO, collection = GO_CC_coll, coll_name = "C5_CC")
GO_CC_out$neg %>% filter(!is.na(atc_code)) %>% 
  filter(str_starts(atc_code, c("C")) | str_starts(atc_code, c( "B"))) %>% 
  filter(PV <= 10^-3)
GO_CC_out$pos %>% filter(!is.na(atc_code)) %>% 
  filter(str_starts(atc_code, c("C")) | str_starts(atc_code, c( "B"))) %>% 
  filter(PV <= 10^-3)


