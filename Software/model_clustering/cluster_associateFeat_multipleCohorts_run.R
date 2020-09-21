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
suppressPackageStartupMessages(library(lme4))
suppressPackageStartupMessages(library(pryr))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(rstatix))
suppressPackageStartupMessages(library(umap))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(RColorBrewer))
options(bitmapType = 'cairo', device = 'png')

parser <- ArgumentParser(description="Find relevant genes for a cluster comparison")
parser$add_argument("--sampleAnnFile", type = "character",  nargs = '*', help = "file with samples to be used")
parser$add_argument("--clusterFile", type = "character", help = "file with clustering structure")
parser$add_argument("--name_cohorts", type = "character", nargs = '*', help = "name of the single cohorts")
parser$add_argument("--split_tot", type = "integer", default = 0, help = "if 0 then inpuntFile load alone, otherwise splitted version")
parser$add_argument("--inputFile", type = "character", default = 'NA', nargs = '*', help = "file to be loaded (predicted tscore or pathScore)")
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
name_cohorts <- args$name_cohorts
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

###################################################################################################################
# name_cohorts <- read.table('/home/luciat/eQTL_PROJECT/INPUT_DATA/SCZ_cohort_names_CLUST', header = F, stringsAsFactors = F)$V1[1:2]
# inputFile <- paste0('/home/luciat/eQTL_PROJECT/OUTPUT_CMC/predict_PGC/200kb/',name_cohorts, '/devgeno0.01_testdevgeno0/Pathway_Reactome_scores.txt')
# split_tot <- 0
# sampleAnnFile <- paste0('/home/luciat/eQTL_PROJECT/INPUT_DATA/Covariates/',name_cohorts,'.covariateMatrix_old.txt')
# clusterFile <- '/home/luciat/eQTL_PROJECT/OUTPUT_CMC/predict_PGC/200kb/Meta_Analysis_SCZ/devgeno0.01_testdevgeno0/tscore_zscaled_clusterCases_PGmethod_HKmetric.RData'
# type_cluster <- 'Cases'
# type_data <- 'path_Reactome'
# type_data_cluster <- 'tscore'
# type_sim <- 'HK'
# min_genes_path <- 2
# pval_id <- 1
# pvalresFile <-  '/home/luciat/eQTL_PROJECT/OUTPUT_all/Meta_Analysis_SCZ/DLPC_CMC/pval_Dx_pheno_covCorr.RData'
# outFold <- '/home/luciat/eQTL_PROJECT/OUTPUT_CMC/predict_PGC/200kb/Meta_Analysis_SCZ/devgeno0.01_testdevgeno0/'
# functR <- '/home/luciat/priler_project/Software/model_clustering/clustering_functions.R'
# type_input <- 'zscaled'
# ####################################################################################################################

source(functR)
cluster_output <- get(load(clusterFile))


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

### load score data ###
sampleAnn <- vector(mode = 'list', length = length(name_cohorts))
scoreMat <- vector(mode = 'list', length = length(name_cohorts))

for(c_id in 1:length(name_cohorts)){
  
  print(name_cohorts[c_id])
  
  sampleAnn[[c_id]] <- read.table(sampleAnnFile[[c_id]], h = T, stringsAsFactors = F)
  if(type_cluster == 'Cases'){
    sampleAnn[[c_id]] <- sampleAnn[[c_id]][sampleAnn[[c_id]]$Dx == 1,]
  }else{
    if(type_cluster == 'Controls'){
      sampleAnn[[c_id]] <- sampleAnn[[c_id]][sampleAnn[[c_id]]$Dx == 0,]
    }else{
      if(type_cluster != 'All')
        stop('type_cluster must be either Cases or Controls or All')
    }
  }
  
  sampleAnn[[c_id]]$Temp_ID <- sampleAnn[[c_id]]$Individual_ID
  sampleAnn[[c_id]]$cohort <- rep(name_cohorts[c_id], nrow(sampleAnn[[c_id]]))
  
  ### load score ###
  if(substr(inputFile[c_id], nchar(inputFile[c_id])-3, nchar(inputFile[c_id])) == '.txt'){
    tmp <- read.delim(inputFile[c_id], h=T, stringsAsFactors = F, check.names = F)
    if(type_data == 'tscore'){
      # correct sample names
      sampleID <- unname(sapply(colnames(tmp)[-1], function(x) strsplit(x, split = '.vs')[[1]][1]))
      colnames(tmp)[-1] <- sampleID
    }
    
    # filter out elements that are repeated twice:
    id_dup <- names(which(table(tmp[,1] ) > 1)) 
    tmp <- tmp[!tmp[,1] %in% id_dup, ]
    tmp <- tmp[tmp[,1] %in% res_pval[, id_info],]
    elementID <- tmp[,1]
    tmp <- tmp[, colnames(tmp) %in% sampleAnn[[c_id]]$Temp_ID]
    scoreMat[[c_id]] <- as.matrix(t(tmp))
    colnames(scoreMat[[c_id]]) <- elementID
  }
}

try(if(any(!sapply(scoreMat, function(x) all(colnames(x) %in% res_pval[, id_info])))) stop("ERROR: wrong filtering elements"))
scoreMat <- do.call(rbind, scoreMat)
res_pval <- res_pval[match(colnames(scoreMat), res_pval[, id_info]),]

print(identical(colnames(scoreMat), res_pval[,id_info]))
# remove sample that have NAs
id_s <- rowSums(is.na(scoreMat)) == 0
if(!all(id_s)){scoreMat <- scoreMat[id_s, ]}
sampleAnn_list <- sampleAnn
sampleAnn <- do.call(rbind, sampleAnn_list)

# match to have the same samples and same order with annotation
sampleAnn <- sampleAnn[match(rownames(scoreMat), sampleAnn$Temp_ID), ]
sampleAnn$cohort_id <- as.numeric(as.factor(sampleAnn$cohort))

# consider only samples used for clustering
common_s <- intersect(cluster_output$samples_id, sampleAnn$Individual_ID)
sampleAnn <- sampleAnn[match(common_s, sampleAnn$Individual_ID), ]
cl_res <- cluster_output$cl_best[match(common_s, cluster_output$cl_best$id),]
scoreMat <- scoreMat[match(common_s, rownames(scoreMat)),]

identical(sampleAnn$Individual_ID, cl_res$id)

print(identical(colnames(scoreMat), res_pval[, id_info]))
print(identical(rownames(scoreMat), sampleAnn$Individual_ID))

covDat <- sampleAnn[,!colnames(sampleAnn) %in% c('Individual_ID', 'Dx', 'genoSample_ID', 'Temp_ID', 'cohort_id')]
output <- list(data = scoreMat, cl = cluster_output$cl_best, cov = covDat)

cl <- cl_res$gr
scale_data <- scale(scoreMat)
attr(scale_data, "scaled:scale") <- NULL
attr(scale_data, "scaled:center") <- NULL

output <- list(inputData = scoreMat, scaleData = scale_data, res_pval = res_pval, cl = cl_res)

##########################
#### check covariates ####
##########################

gr_names <- sort(unique(cl))
P <- length(gr_names)

test_cov <- vector(mode = 'list', length = length(gr_names))
for(i in 1:length(gr_names)){
  
  print(paste0('group', gr_names[i], '_vs_all'))
  
  # j vs all
  gr_id <- as.factor(as.numeric(cl == gr_names[i]))
  test_cov[[i]] <- data.frame(cov = colnames(covDat), comp = rep(sprintf('gr%i_vs_all', i), ncol(covDat)))
  test_cov[[i]]$pval<- NA
  test_cov[[i]]$estimates<- NA
  test_cov[[i]]$CI_low<- NA
  test_cov[[i]]$CI_up <- NA
  
  # test for cohort effect:
  tmp <- as.data.frame(rstatix::chisq_test(table(gr_id, covDat[,'cohort'])))
  test_cov[[i]]$pval[test_cov[[i]]$cov == 'cohort']<-  tmp$p
  test_cov[[i]]$estimates[test_cov[[i]]$cov == 'cohort'] <- tmp$statistic
  test_cov[[i]]$test_type[test_cov[[i]]$cov == 'cohort'] <- tmp$method
  
  # use linear mixed models
  df <- cbind(data.frame(gr_id = gr_id), covDat)
  df$cohort <- factor(df$cohort, levels = name_cohorts)
  fmla <- as.formula(paste('gr_id ~ ', paste0(colnames(df[, !colnames(df) %in% c('gr_id', 'cohort')]), collapse = '+'), '+ (1|cohort)'))
  mixed_mod <- glmer(fmla, df, family  = 'binomial')
  tmp <- coef(summary(mixed_mod))
  
  test_cov[[i]]$test_type[match( rownames(tmp)[-1], test_cov[[i]]$cov)] <- 'generalized lmm (1|cohort)'
  test_cov[[i]]$pval[match( rownames(tmp)[-1], test_cov[[i]]$cov)] <- tmp[-1,4]
  test_cov[[i]]$estimates[match(rownames(tmp)[-1], test_cov[[i]]$cov)] <- exp(tmp[-1,1])
  test_cov[[i]]$CI_low[match(rownames(tmp)[-1], test_cov[[i]]$cov)] <- exp(tmp[-1,1] + qnorm(0.025)*tmp[-1,2])
  test_cov[[i]]$CI_up [match( rownames(tmp)[-1], test_cov[[i]]$cov)]<- exp(tmp[-1,1] + qnorm(0.975)*tmp[-1,2])
  
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
save(output, file = sprintf('%s%sOriginal_%sCluster%s_featAssociation.RData', outFold, type_data, type_data_cluster, type_cluster))



