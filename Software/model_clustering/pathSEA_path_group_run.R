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

parser <- ArgumentParser(description="path SEA analysis (cmap data)")
parser$add_argument("--pathCluster_file", type = "character", help = "file to be loaded")
parser$add_argument("--atc_file", type = "character", help = "")
parser$add_argument("--cmap_fold", type = "character", help = "")
parser$add_argument("--type_cluster", type = "character", default = 'All', help = "All, Cases, Controls")
parser$add_argument("--outFold", type="character", help = "Output file [basename only]")

args <- parser$parse_args()
pathCluster_file <- args$pathCluster_file
atc_file <- args$atc_file
cmap_fold <- args$cmap_fold
type_cluster <- args$type_cluster
outFold <- args$outFold

####################################################################################################################
# pathCluster_file <- 'OUTPUT_GTEx/predict_CAD/Liver/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/CAD_HARD_clustering/update_corrPCs/pathOriginal_filtJS0.2_corrPCs_tscoreClusterCases_featAssociation.RData'
# cmap_fold <- '/psycl/g/mpsziller/lucia/refData/Cmap_MSigDB_v6.1_PEPs'
# atc_file <- '/psycl/g/mpsziller/lucia/drug_targeting/WHO ATC-DDD 2021-12-03.csv'
####################################################################################################################

rpBig <- openRepository(cmap_fold)
Reactome_coll <- loadCollection(rpBig,"C2_CP:REACTOME")
GO_MF_coll <- loadCollection(rpBig,"C5_MF")
GO_BP_coll <- loadCollection(rpBig,"C5_BP")
GO_CC_coll <- loadCollection(rpBig,"C5_CC")
drug_atc <- read.csv(atc_file)

# load results:
path_feat <- get(load(pathCluster_file))
tot_res <- do.call(rbind, path_feat$test_feat) %>% filter(pval_corr <= 0.01)
gr <- sort(unique(tot_res$comp))
n_gr <- length(gr)

gr_res <- lapply(1:n_gr, function(x) tot_res[tot_res$comp == sprintf('gr%i_vs_all', x),])
# remove discordant results in sign
for(i in 1:n_gr){
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

Reactome_names <- do.call(rbind, path_feat$res_pval) %>% filter(name == 'Reactome')
Reactome_names <- Reactome_names$path 
Reactome_names <- Reactome_names[!duplicated(Reactome_names)]

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
GO_names <-  do.call(rbind, path_feat$res_pval) %>% filter(name == 'GO')
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
        arrange(desc(sign_p)) %>% filter(sign_p > 0) %>%  
        mutate(FDR = p.adjust(PV,method = 'BH'))
        #filter(FDR <= 0.05) 
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
        arrange(sign_p)  %>% filter(sign_p < 0)  %>%  
        mutate(FDR = p.adjust(PV,method = 'BH')) 
        #filter(FDR <= 0.05) 
      psea_out[[i]]$atc_code <- drug_atc$atc_code[match(psea_out[[i]]$drug, drug_atc$atc_name)]
    }
  }
  psea_pos <- do.call(rbind, psea_out)
  
  return(list(neg = psea_neg, pos = psea_pos))
}

###################################################################

Reactome_out <- psea_group(df_names_R, collection = Reactome_coll, coll_name = "C2_CP:REACTOME")
GO_MF_out <- psea_group(df_names_GO, collection = GO_MF_coll, coll_name = "C5_MF")
GO_BP_out <- psea_group(df_names_GO, collection = GO_BP_coll, coll_name = "C5_BP")
GO_CC_out <- psea_group(df_names_GO, collection = GO_CC_coll, coll_name = "C5_CC")

# combine 
pos <- rbind(Reactome_out$pos %>% mutate(db = 'Reactome'), 
             GO_MF_out$pos %>% mutate(db = 'GO_MF'), 
             GO_BP_out$pos %>% mutate(db = 'GO_BP'), 
             GO_CC_out$pos %>% mutate(db = 'GO_CC')) %>% 
  filter(FDR <= 0.05)

neg <- rbind(Reactome_out$neg %>% mutate(db = 'Reactome'), 
             GO_MF_out$neg %>% mutate(db = 'GO_MF'), 
             GO_BP_out$neg %>% mutate(db = 'GO_BP'), 
             GO_CC_out$neg %>% mutate(db = 'GO_CC')) %>% 
  filter(FDR <= 0.05)

tot <- rbind(pos %>% mutate(type = 'up-reg pathways'), 
             neg %>% mutate(type = 'down-reg pathways'))

tot$atc_meaning3 <- drug_atc$atc_name[match(substr(tot$atc_code, 1, 3), drug_atc$atc_code)]
tot$atc_meaning1 <- drug_atc$atc_name[match(substr(tot$atc_code, 1, 1), drug_atc$atc_code)]

rownames(tot) <- 1:nrow(tot)



write.table(tot, col.names = T, row.names = F, sep = '\t', quote = F, 
            file = sprintf('%spathSEA_corrPCs_tscoreCluster%s_featAssociation.txt', outFold, type_cluster))
