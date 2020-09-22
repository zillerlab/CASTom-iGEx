# find relevant genes/pathways to cluster structure 
# use GLM

options(stringsAsFactors=F)
options(max.print=1000)
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(qvalue))
suppressPackageStartupMessages(library(igraph))
suppressPackageStartupMessages(library(MASS))
suppressPackageStartupMessages(library(lmtest))
suppressPackageStartupMessages(library(pryr))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(rstatix))
suppressPackageStartupMessages(library(umap))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(RColorBrewer))
options(bitmapType = 'cairo', device = 'png')

parser <- ArgumentParser(description="Find relevant genes for a cluster comparison")
parser$add_argument("--sampleAnnFile", type = "character", help = "file with samples to be used")
parser$add_argument("--clusterFile", type = "character", help = "file with clustering structure")
parser$add_argument("--split_tot", type = "integer", default = 0, help = "if 0 then inpuntFile load alone, otherwise splitted version")
parser$add_argument("--inputFile", type = "character", default = 'NA', help = "file to be loaded (predicted tscore or pathScore)")
parser$add_argument("--type_cluster", type = "character", default = 'All', help = "All, Cases, Controls")
parser$add_argument("--functR", type = "character", help = "functions to be used")
parser$add_argument("--type_data", type = "character", help = "tscore, path_Reactome or path_GO")
parser$add_argument("--type_data_cluster", type = "character", help = "tscore, path_Reactome or path_GO")
parser$add_argument("--type_sim", type = "character", default = 'HK', help = "HK or ED or SNF")
parser$add_argument("--type_input", type = "character", default = 'original', help = "original or zscaled")
parser$add_argument("--pvalresFile", type = "character", default = 'NA', help = "file with pvalue results")
parser$add_argument("--min_genes_path", type = "integer", default = 1, help = "minimum number of genes for a pathway, if > 1 recompute corrected pvalues")
parser$add_argument("--pval_id", type = "integer", default = 0, help = "id to be used on pvalue file")
parser$add_argument("--outFold", type="character", help = "Output file [basename only]")

args <- parser$parse_args()
split_tot <- args$split_tot
sampleAnnFile <- args$sampleAnnFile
inputFile <- args$inputFile
type_cluster <- args$type_cluster
functR <- args$functR
type_data <- args$type_data
type_sim <- args$type_sim
type_input <- args$type_input
clusterFile <- args$clusterFile
pvalresFile <- args$pvalresFile
type_data_cluster <- args$type_data_cluster
min_genes_path <- args$min_genes_path
pval_id <- args$pval_id
outFold <- args$outFold

# ###################################################################################################################
# inputFile <- '/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/OUTPUT_GTEx/predict_CAD/Liver/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/predictedTscores_splitGenes'
# split_tot <- 100
# sampleAnnFile <- '/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/INPUT_DATA_GTEx/CAD/Covariates/UKBB/CAD_HARD_clustering/covariateMatrix_CADHARD_All_phenoAssoc.txt'
# clusterFile <- '/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/OUTPUT_GTEx/predict_CAD/Liver/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/CAD_HARD_clustering/tscore_zscaled_clusterCases_PGmethod_HKmetric.RData'
# type_cluster <- 'Cases'
# type_data <- 'tscore'
# type_data_cluster <- 'tscore'
# type_sim <- 'HK'
# min_genes_path <- 2
# pval_id <- 1
# pvalresFile <- '/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/OUTPUT_GTEx/predict_CAD/Liver/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/pval_CAD_pheno_covCorr.RData'
# outFold <- '/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/OUTPUT_GTEx/predict_CAD/Liver/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/CAD_HARD_clustering/'
# functR <- '/psycl/g/mpsziller/lucia/priler_project/Software/model_clustering/clustering_functions.R'
# type_input <- 'zscaled'
# ####################################################################################################################

source(functR)
cluster_output <- get(load(clusterFile))

sampleAnn <- read.table(sampleAnnFile, h=T, stringsAsFactors = F, check.names = F)
sampleAnn <- sampleAnn[match(cluster_output$samples_id, sampleAnn$Individual_ID), ]

identical(sampleAnn$Individual_ID, cluster_output$samples_id)

# load pval res
res_pval <- get(load(pvalresFile))
if(type_data == 'tscore'){
  res_pval <- res_pval$tscore[[pval_id]]
  id_pval <- 8
  id_info <- 2
}else{
  if(type_data == 'path_Reactome'){
    res_pval <- res_pval$pathScore_reactome[[pval_id]]
    id_pval <- 13
    id_info <- 1
  }else{
    if(type_data == 'path_GO'){
      res_pval <- res_pval$pathScore_GO[[pval_id]]
      id_pval <- 15
      id_info <- 1
    }else{
      stop('unknown pathway called')
    }
  }
}

# recompute pvalue if ngenes_tscore > 1
if(min_genes_path > 1 & grepl('path',type_data)){
  res_pval <- res_pval[res_pval$ngenes_tscore >= min_genes_path, ]
  res_pval[,id_pval+1] <- qvalue(res_pval[,id_pval])$qvalues
  res_pval[,id_pval+2] <- p.adjust(res_pval[,id_pval], method = 'BH')
}

# load input matrix 
if(split_tot == 0){
  
  scoreMat <- get(load(inputFile))
  # filter out based on samples and ids
  id_el <- intersect(scoreMat[,1], res_pval[, id_info])
  scoreMat <- scoreMat[match(id_el,scoreMat[,1]), ]
  
  common_samples <- sampleAnn$Individual_ID
  scoreMat <- t(scoreMat[,match(common_samples,colnames(scoreMat))])
  
  rownames(scoreMat) <- common_samples
  colnames(scoreMat) <- id_el
  res_pval <- res_pval[match(id_el, res_pval[, id_info]),]
  
}else{
  
  ###### load score Mat #######
  scoreMat_list <- vector(mode = 'list', length = split_tot)
  samplesID <- vector(mode = 'list', length = split_tot)
  elementID <- NULL
  
  for(i in 1:split_tot){
    
    print(i)
    if(file.exists(sprintf('%s%i.RData', inputFile, i))){
      tmp <- get(load(sprintf('%s%i.RData', inputFile, i)))
      elementID <- c(elementID,tmp[,1])
      samplesID[[i]] <- intersect(sampleAnn$Individual_ID, colnames(tmp))
      scoreMat_list[[i]] <- t(tmp[,match(samplesID[[i]],colnames(tmp))])
    }else{
      print(sprintf('split %i does not exist', i))
      split_tot <- split_tot - 1
    }
  }
  
  print(split_tot)
  # check samplesID always the same
  if(!all(table(unlist(samplesID)) == split_tot)){print('ERROR: wrong name annotations')}
  
  scoreMat <- do.call(cbind, scoreMat_list)
  colnames(scoreMat) <- elementID
  rm(scoreMat_list)
  
  # filter out elements that are repeated twice:
  id_dup <- names(which(table(colnames(scoreMat)) > 1)) 
  scoreMat <- scoreMat[, !colnames(scoreMat) %in% id_dup]
  
  id_el <- intersect(colnames(scoreMat),  res_pval[, id_info])
  scoreMat <- scoreMat[, match(id_el, colnames(scoreMat))]
  
  rownames(scoreMat) <- samplesID[[1]]
  # remove sample that have NAs
  id_s <- rowSums(is.na(scoreMat)) == 0
  if(!all(id_s)){scoreMat <- scoreMat[id_s, ]}
  
  common_samples <- sampleAnn$Individual_ID
  scoreMat <- scoreMat[match(common_samples,rownames(scoreMat)),]
  res_pval <- res_pval[match(id_el, res_pval[, id_info]),]
  
}

print(identical(colnames(scoreMat), res_pval[, id_info]))
print(identical(rownames(scoreMat), sampleAnn$Individual_ID))

covDat <- sampleAnn[,!colnames(sampleAnn) %in% c('Individual_ID', 'Dx', 'genoSample_ID')]
output <- list(data = scoreMat, cl = cluster_output$cl_best, cov = covDat)

cl <-  cluster_output$cl_best$gr
scale_data <- scale(scoreMat)
attr(scale_data, "scaled:scale") <- NULL
attr(scale_data, "scaled:center") <- NULL

output <- list(inputData = scoreMat, scaleData = scale_data, res_pval = res_pval, cl = cluster_output$cl_best)

##########################
#### check covariates ####
##########################

gr_names <- sort(unique(cl))
P <- length(gr_names)
covDat <- sampleAnn[, !colnames(sampleAnn) %in% c('Individual_ID', 'genoSample_ID', 'Dx')]
output$covDat = covDat

test_cov <- vector(mode = 'list', length = length(gr_names))
for(i in 1:length(gr_names)){
  
  print(paste0('group', gr_names[i], '_vs_all'))
  
  # j vs all
  gr_id <- as.factor(as.numeric(cl == gr_names[i]))
  test_cov[[i]] <- data.frame(cov = colnames(covDat), comp = rep(sprintf('gr%i_vs_all', i), ncol(covDat)))
  test_cov[[i]]$test_type<- NA
  test_cov[[i]]$pval<- NA
  test_cov[[i]]$estimates<- NA
  test_cov[[i]]$CI_low<- NA
  test_cov[[i]]$CI_up <- NA
  
  for(l in 1:ncol(covDat)){
    
    if((is.integer(covDat[, l]) | is.character(covDat[, l])) & colnames(covDat)[l] != 'Age'){
      
      tmp <- as.data.frame(rstatix::chisq_test(table(gr_id, covDat[,l])))
      test_cov[[i]]$pval[l]<-  tmp$p
      test_cov[[i]]$estimates[l] <- tmp$statistic
      test_cov[[i]]$test_type[l] <- tmp$method
      
    }else{
      tmp_data <- data.frame(f = covDat[,l], g = gr_id)
      tmp <- as.data.frame(rstatix::wilcox_test(f ~ g, data = tmp_data, detailed = T))
      test_cov[[i]]$pval[l] <- tmp$p
      test_cov[[i]]$estimates[l] <- tmp$estimate
      test_cov[[i]][l,c('CI_low', 'CI_up')] <- c(tmp$conf.low, tmp$conf.high)
      test_cov[[i]]$test_type[l] <- tmp$method
      
    }
    
  }
  
  test_cov[[i]]$pval_corr <- p.adjust(test_cov[[i]]$pval, method = 'BH')

}


test_cov <- do.call(rbind, test_cov)
test_cov$pval_corr_overall <-  p.adjust(test_cov$pval, method = 'BH')


output$test_cov = test_cov


####################################
#### check features (gene/path) ####
####################################


test_feat <- vector(mode = 'list', length = length(gr_names))
for(i in 1:length(gr_names)){
  
  print(paste0('group', gr_names[i], '_vs_all'))
  
  # j vs all
  gr_id <- as.factor(as.numeric(cl == gr_names[i]))
  test_feat[[i]] <- data.frame(feat = colnames(scale_data), comp = rep(sprintf('gr%i_vs_all', i), ncol(scale_data)))
  test_feat[[i]]$pval<- NA
  test_feat[[i]]$estimates<- NA
  test_feat[[i]]$CI_low<- NA
  test_feat[[i]]$CI_up <- NA
  
  for(l in 1:ncol(scale_data)){

      print(l) 
      tmp_data <- data.frame(f = scale_data[,l], g = gr_id)
      tmp <- as.data.frame(wilcox_test(f~g,data = tmp_data, detailed = T))
      test_feat[[i]]$pval[l] <- tmp$p
      test_feat[[i]]$estimates[l] <- tmp$estimate
      test_feat[[i]][l,c('CI_low', 'CI_up')] <- c(tmp$conf.low, tmp$conf.high)
  }
  
  test_feat[[i]]$pval_corr <- p.adjust(test_feat[[i]]$pval, method = 'BH')
  
}


test_feat <- do.call(rbind, test_feat)
test_feat$pval_corr_overall <-  p.adjust(test_feat$pval, method = 'BH')

output$test_feat = test_feat

# Save
save(output, file = sprintf('%s%sOriginal_%sCluster%s_featAssociation.RData', outFold, type_data_cluster, type_data, type_cluster))


