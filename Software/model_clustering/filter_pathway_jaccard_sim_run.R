#!/usr/bin/env Rscript

# Consider all pathways, filter to obtain a restricted set of 
# not-overlapping genes that retain the pathway with highest coverage
# for each tissue separately since genes are not exactly the same across pathways

options(stringsAsFactors=F)
options(max.print=1000)
suppressPackageStartupMessages(library(argparse))
options(bitmapType = 'cairo', device = 'png')

parser <- ArgumentParser(description="select pathways")
parser$add_argument("--pvalresFile", type = "character", help = "file with pvalue results")
parser$add_argument("--thr_js", type = "double", default = 0.3, help = "")
parser$add_argument("--outFold", type="character", help = "Output file [basename only]")

args <- parser$parse_args()
pvalresFile <- args$pvalresFile
thr_js <- args$thr_js
outFold <- args$outFold

########################################################################################################################
# outFold <- '/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/OUTPUT_GTEx/predict_CAD/Whole_Blood/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/'
# pvalresFile <- sprintf('%spval_CAD_pheno_covCorr.RData', outFold)
# thr_js <- 0.2
# #########################################################################################################################


# get pathway structure
get_path_list <- function(path_structure, path_structure_name){
  df <- data.frame(path = sapply(path_structure[[1]], function(x) x$path$path), 
                   ngenes_tscore = sapply(path_structure[[1]], function(x) x$path$ngenes_tscore), 
                   ngenes_path = sapply(path_structure[[1]], function(x) x$path$ngenes_path))
  df$coverage <- df$ngenes_tscore/df$ngenes_path
  
  df$gene_name <- sapply(path_structure[[1]], function(x) paste0(x$tscore$external_gene_name, collapse = ','))
  df$name <- path_structure_name
  return(df)  
}

compute_js <- function(p1, p2){
  
  n_int <- length(intersect(p1, p2))
  n_union <- length(union(p1, p2))
  
  return(n_int/n_union)
}

## clumping pathways ##
clumping_pathways <- function(path, thr_js, min_ngenes = 3, max_genes = 200){
  
  path <- path[path$ngenes_tscore >= min_ngenes, ]
  path <- path[path$ngenes_tscore <= max_genes & path$ngenes_path <= max_genes, ]
  
  path <- path[order(path$coverage, path$ngenes_tscore, decreasing = T), ]
  path$id <- paste('path', 1:nrow(path), sep = '_')
  
  # compute js total
  path_gene_list <- lapply(path$gene_name, function(x) strsplit(x, split = '[,]')[[1]])
  
  js_mat <- sapply(path_gene_list, function(y) 
    sapply(path_gene_list, function(x) compute_js(x,y)))
  rownames(js_mat) <- colnames(js_mat) <- path$id
  print('total Jaccard similarity computed')
  
  element_rm <- c()
  stop_cond <- F
  list_feat <- path$id
  
  while(!stop_cond){
    
    id <- which(js_mat[,list_feat[1]]>thr_js)
    if(length(id)>1){
      element_rm <- c(element_rm, names(id)[-which.max(path$coverage[match(names(id),path$id)])])
    }
    list_feat <- list_feat[!list_feat %in% names(id)]
    # print(length(list_feat))
    stop_cond <- length(list_feat) == 0
    
  }
  
  element_rm <- unique(element_rm)
  path <- path[!path$id %in% element_rm, -ncol(path)]
  
  return(path)
}

########################################
res_pval <- get(load(pvalresFile))

path_reactome <- get_path_list(res_pval$info_pathScore_reactome, 'Reactome')
path_GO <- get_path_list(res_pval$info_pathScore_GO, 'GO')

tot_path <- rbind(path_reactome, path_GO)

path_filt <- clumping_pathways(path = tot_path, thr_js = thr_js)

# save
write.table(path_filt, file = sprintf('%sselected_pathways_JSthr%s.txt', outFold, 
                                      as.character(thr_js)), 
            quote = T, sep = '\t', col.names = T, row.names = F)




