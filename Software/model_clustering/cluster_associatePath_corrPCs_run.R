#!/usr/bin/env Rscript

# find relevant pathways to cluster structure 
# use wilcoxon test
# possible include multiple tissues

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
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(rlist))
suppressPackageStartupMessages(library(doParallel))
options(bitmapType = 'cairo', device = 'png')

parser <- ArgumentParser(description="Find relevant genes for a cluster comparison")
parser$add_argument("--sampleAnnFile", type = "character", help = "file with samples to be used")
parser$add_argument("--clusterFile", type = "character", help = "file with clustering structure")
parser$add_argument("--inputFold", type = "character", nargs = '*', help = "folder containing pathway scores")
parser$add_argument("--tissues", type = "character", nargs = '*', help = "tissues name")
parser$add_argument("--type_cluster", type = "character", default = 'All', help = "All, Cases, Controls")
parser$add_argument("--functR", type = "character", help = "functions to be used")
parser$add_argument("--type_data", type = "character", help = "tscore, path_Reactome or path_GO")
parser$add_argument("--type_data_cluster", type = "character", help = "tscore, path_Reactome or path_GO")
parser$add_argument("--type_sim", type = "character", default = 'HK', help = "HK or ED or SNF")
parser$add_argument("--type_input", type = "character", default = 'original', help = "original or zscaled")
parser$add_argument("--pvalresFile", type = "character", nargs = '*',  help = "file with pvalue results")
parser$add_argument("--geneInfoFile", type = "character", nargs = '*', default = NULL, help = "file with info model genes")
parser$add_argument("--pval_id", type = "integer", default = 0, help = "id to be used on pvalue file")
parser$add_argument("--pvalcorr_thr", type = "double", default = 0.05, help = "group sign")
parser$add_argument("--ncores", type = "integer", default = 5, help = "n, cores parallelization per tissue")
parser$add_argument("--thr_js", type = "double", default = 0.2, help = "jaccard similarity threshold for filtering pathways")
parser$add_argument("--outFold", type="character", help = "Output file [basename only]")

args <- parser$parse_args()
sampleAnnFile <- args$sampleAnnFile
inputFold <- args$inputFold
type_cluster <- args$type_cluster
functR <- args$functR
type_data <- args$type_data
type_sim <- args$type_sim
type_input <- args$type_input
clusterFile <- args$clusterFile
pvalresFile <- args$pvalresFile
type_data_cluster <- args$type_data_cluster
pval_id <- args$pval_id
tissues <- args$tissues
pvalcorr_thr <- args$pvalcorr_thr
ncores <- args$ncores
thr_js <- args$thr_js
outFold <- args$outFold

###################################################################################################################
# tissues <- c('Liver', 'Adipose_Subcutaneous')
# inputFold <- sprintf('OUTPUT_GTEx/predict_CAD/%s/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/', tissues)
# sampleAnnFile <- 'INPUT_DATA_GTEx/CAD/Covariates/UKBB/CAD_HARD_clustering/covariateMatrix_CADHARD_All_phenoAssoc.txt'
# clusterFile <- 'OUTPUT_GTEx/predict_CAD/Liver/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/CAD_HARD_clustering/update_corrPCs/tscore_corrPCs_zscaled_clusterCases_PGmethod_HKmetric.RData'
# type_cluster <- 'Cases'
# type_data <- 'path'
# type_data_cluster <- 'tscore'
# type_sim <- 'HK'
# pval_id <- 1
# pvalresFile <- sprintf('OUTPUT_GTEx/predict_CAD/%s/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/pval_CAD_pheno_covCorr.RData', tissues)
# outFold <- 'OUTPUT_GTEx/predict_CAD/Liver/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/CAD_HARD_clustering/'
# functR <- '/psycl/g/mpsziller/lucia/castom-igex/Software/model_clustering/clustering_functions.R'
# type_input <- 'zscaled'
# thr_js <- 0.2
####################################################################################################################

source(functR)
# combine function for parallelization
comb <- function(x, ...) {
  lapply(seq_along(x),
         function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
}
###############################

cluster_output <- get(load(clusterFile))

# used in case of prediction
if(!'samples_id' %in% names(cluster_output)){
  cluster_output$samples_id <- cluster_output$sampleAnn$Individual_ID
}
if(!'cl_best' %in% names(cluster_output)){
  cluster_output$cl_best <- cluster_output$cl_new
}

# if(type_data == 'tscore'){
#   id_pval <- 8
#   id_info <- 2
#   id_geno_summ <- 3
# }else{
#   if(type_data == 'path_Reactome'){
#     id_pval <- 13
#     id_info <- 1
#     id_geno_summ <- 4
#   }else{
#     if(type_data == 'path_GO'){
#       id_pval <- 15
#       id_info <- 1
#       id_geno_summ <- 6
#     }else{
#       stop('unknown pathway called')
#     }
#   }
# }

print(inputFold)
print(pvalresFile)

sampleAnn <- read.table(sampleAnnFile, h=T, stringsAsFactors = F, check.names = F)
sampleAnn$Individual_ID <- as.character(sampleAnn$Individual_ID)
sampleAnn <- sampleAnn[match(cluster_output$samples_id, sampleAnn$Individual_ID), ]
name_cov <- setdiff(colnames(sampleAnn),c('Individual_ID', 'genoSample_ID', 'Dx', 'Sex', 'Age', 'Gender', 'Array'))
covDat <- sampleAnn[,!colnames(sampleAnn) %in% c('Individual_ID', 'Dx', 'genoSample_ID')]

identical(as.character(sampleAnn$Individual_ID), as.character(cluster_output$samples_id))
# scale_data_t <- scoreMat_t <- res_pval_t <- list()

# load score matrix, load pval matrix, rescale for p-value
registerDoParallel(cores=min(ncores, length(tissues)))

res <- foreach(id_t=1:length(tissues), .combine='comb', 
               .multicombine=TRUE, 
               .init=list(list(), list(), list()))%dopar%{

#res <- list()                 
#for(id_t in 1:length(tissues)){
                 print(tissues[id_t])
                 # load pval res
                 res_pval <- get(load(pvalresFile[id_t]))
                 res_pval1 <- res_pval$pathScore_reactome[[pval_id]]
                 res_pval2 <- res_pval$pathScore_GO[[pval_id]]
     
                 
                 # filter pathways based on shared structure
                 path_ann <- read.table(sprintf('%sselected_pathways_JSthr%s.txt', inputFold[id_t], as.character(thr_js)), 
                                        h=T, stringsAsFactors = F, sep = '\t')
                 
                 res_pval1 <- res_pval1[res_pval1$path %in% path_ann$path[path_ann$name == 'Reactome'], ]
                 res_pval2 <- res_pval2[res_pval2$path %in% path_ann$path[path_ann$name == 'GO'], ]
                 
                 #### load input matrix ####
                 if(file.exists(sprintf('%sPathway_Reactome_scores.RData', inputFold[id_t]))){
                   pathR_file <- sprintf('%sPathway_Reactome_scores.RData', inputFold[id_t])
                 }else{
                   pathR_file <- sprintf('%sPathway_Reactome_scores.txt', inputFold[id_t])
                 }
                 load_output_reactome <- load_input_matrix(inputFile = pathR_file, 
                                                  sampleAnn = sampleAnn, 
                                                  res_pval = res_pval1, 
                                                  split_tot = 0, 
                                                  id_info = 1)
                 
                 if(file.exists(sprintf('%sPathway_GO_scores.RData', inputFold[id_t]))){
                   pathGO_file <- sprintf('%sPathway_GO_scores.RData', inputFold[id_t])
                 }else{
                   pathGO_file <- sprintf('%sPathway_GO_scores.txt', inputFold[id_t])
                 }
                 load_output_GO <- load_input_matrix(inputFile = pathGO_file, 
                                                           sampleAnn = sampleAnn, 
                                                           res_pval = res_pval2, 
                                                           split_tot = 0, 
                                                           id_info = 1)
                 colnames(load_output_GO$scoreMat) <- load_output_GO$res_pval$path
                 
                 
                 load_output_reactome$res_pval$name <- 'Reactome'
                 load_output_GO$res_pval$name <- 'GO'
                 load_output_GO$res_pval <- load_output_GO$res_pval[, -c(1,3)]
                 res_pval_t <- rbind(load_output_reactome$res_pval, load_output_GO$res_pval)
                 
                 scoreMat_t <- cbind(load_output_reactome$scoreMat, load_output_GO$scoreMat)
                 input_data_notcorr <- scale(scoreMat_t)
                 attr(input_data_notcorr, "scaled:scale") <- NULL
                 attr(input_data_notcorr, "scaled:center") <- NULL
                 
                 # remove PCs1-10 for each genes
                 input_data <- matrix(ncol = ncol(input_data_notcorr), nrow = nrow(input_data_notcorr))
                 rownames(input_data) <- rownames(input_data_notcorr)
                 colnames(input_data) <- colnames(input_data_notcorr)
                 fmla <- as.formula(paste('g ~', paste0(name_cov, collapse = '+')))
                 for(i in 1:ncol(input_data_notcorr)){
                   # print(i)
		   tmp <- data.frame(g = input_data_notcorr[,i], sampleAnn[, name_cov])
		   if(any(is.na(tmp$g))){
                    # only needed when imputed gene expression is computed in non-harmonized data
                    # some genes might have variance zero
                    input_data[,i] <- 0
                   }else{
                    reg <- lm(fmla, data = tmp)
                    input_data[,i] <- reg$residuals
                  }
                 }
                 print("corrected for PCs")
                 
                 print(identical(colnames(input_data), res_pval_t[, 1]))
                 print(identical(rownames(input_data), sampleAnn$Individual_ID))
                 
                 scale_data_t <- input_data
                 res_pval_t$tissue <- tissues[id_t]
                 list(scoreMat_t, scale_data_t, res_pval_t)
                 # res[[id_t]] <- list(scoreMat_t, scale_data_t, res_pval_t)
                 
}

# scoreMat_t <- lapply(res, function(x) x[[1]])
# scale_data_t <- lapply(res, function(x) x[[2]])
# res_pval_t <- lapply(res, function(x) x[[3]])
scoreMat_t <- res[[1]]
scale_data_t <- res[[2]]
res_pval_t <- res[[3]]

output <- list(inputData = scoreMat_t, scaleData = scale_data_t, res_pval = res_pval_t, 
               cl = cluster_output$cl_best, tissues = tissues)

print(str(output))

##########################
#### check covariates ####
##########################

cl <-  cluster_output$cl_best$gr
gr_names <- sort(unique(cl))
P <- length(gr_names)
output$covDat <- covDat

test_cov <- vector(mode = 'list', length = length(gr_names))
for(i in 1:length(gr_names)){
  
  print(paste0('group', gr_names[i], '_vs_all'))
  
  # j vs all
  gr_id <- as.factor(as.numeric(cl == gr_names[i]))
  test_cov[[i]] <- data.frame(cov = colnames(covDat), comp = rep(sprintf('gr%i_vs_all', gr_names[i]), ncol(covDat)))
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

output$test_cov <- test_cov

####################################
#### check features (gene/path) ####
####################################

registerDoParallel(cores=min(ncores, length(tissues)))

test_feat_t <- foreach(id_t=1:length(tissues))%dopar%{
  
  # for(id_t in 1:length(tissues)){
  test_feat <- vector(mode = 'list', length = length(gr_names))
  for(i in 1:length(gr_names)){
    
    print(paste0('group', gr_names[i], '_vs_all, ', tissues[id_t]))
    
    # j vs all
    gr_id <- as.factor(as.numeric(cl != gr_names[i]))
    test_feat[[i]] <- data.frame(feat = colnames(scale_data_t[[id_t]]), 
                                 comp = rep(sprintf('gr%i_vs_all', gr_names[i]), ncol(scale_data_t[[id_t]])))
    test_feat[[i]]$pval<- NA
    test_feat[[i]]$estimates<- NA
    test_feat[[i]]$CI_low<- NA
    test_feat[[i]]$CI_up <- NA
    
    for(l in 1:ncol(scale_data_t[[id_t]])){
      # print(l) 
      tmp_data <- data.frame(f = scale_data_t[[id_t]][,l], g = gr_id)
      tmp <- as.data.frame(wilcox_test(f~g,data = tmp_data, detailed = T))
      test_feat[[i]]$pval[l] <- tmp$p
      test_feat[[i]]$estimates[l] <- tmp$estimate
      test_feat[[i]][l,c('CI_low', 'CI_up')] <- c(tmp$conf.low, tmp$conf.high)
    }
    
    test_feat[[i]]$pval_corr <- p.adjust(test_feat[[i]]$pval, method = 'BH')
    
  }
  
  test_feat <- do.call(rbind, test_feat)
  test_feat$pval_corr_overall <-  p.adjust(test_feat$pval, method = 'BH')
  test_feat$tissue <- tissues[id_t]
  test_feat
  
}

output$test_feat <- test_feat_t

# Save
save(output, file = sprintf('%s%sOriginal_filtJS%s_corrPCs_%sCluster%s_featAssociation.RData', outFold, type_data, 
                            thr_js, type_data_cluster, type_cluster))


####################################################################################
# summarise results for Reactome annotation


