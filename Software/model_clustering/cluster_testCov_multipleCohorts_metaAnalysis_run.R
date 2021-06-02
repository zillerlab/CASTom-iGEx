# associate endophenotypes to cluster structure 
# use GLM, multiplt cohorts

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
suppressPackageStartupMessages(library(coin))
suppressPackageStartupMessages(library(umap))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(lme4))
options(bitmapType = 'cairo', device = 'png')

parser <- ArgumentParser(description="Associate cluster to covariates, multiple cohorts analysed together (meta analysis)")
parser$add_argument("--name_cohorts", type = "character", nargs = '*', help = "name cohorts")
parser$add_argument("--covDatFile", type = "character", nargs = '*', help = "file to be loaded (endophenotypes, must be a unique matrix)")
parser$add_argument("--sampleAnnFile", type = "character", nargs = '*', help = "file with samples to be used")
parser$add_argument("--clusterFile", type = "character", nargs = '*', help = "file with clustering structure")
parser$add_argument("--type_cluster", type = "character", default = 'All', help = "All, Cases, Controls")
parser$add_argument("--functR", type = "character", help = "functions to be used")
parser$add_argument("--type_data", type = "character", help = "tscore, path_Reactome or path_GO")
parser$add_argument("--type_sim", type = "character", default = 'HK', help = "HK or ED or SNF")
parser$add_argument("--type_input", type = "character", default = 'original', help = "original or zscaled")
parser$add_argument("--outFold", type="character", help = "Output file [basename only]")

args <- parser$parse_args()
name_cohorts <- args$name_cohorts
covDatFile <- args$covDatFile
sampleAnnFile <- args$sampleAnnFile
type_cluster <- args$type_cluster
functR <- args$functR
type_data <- args$type_data
type_sim <- args$type_sim
type_input <- args$type_input
clusterFile <- args$clusterFile
outFold <- args$outFold

####################################################################################################################
# name_cohorts <- read.table('INPUT_DATA/SCZ_cohort_names_CLUST')$V1
# covDatFile <- paste0('INPUT_DATA/Covariates/',name_cohorts,'.covariateMatrix_old.txt')
# sampleAnnFile <- paste0('INPUT_DATA/Covariates/',name_cohorts,'.covariateMatrix_old.txt')
# clusterFile <- 'OUTPUT_CMC/predict_PGC/200kb/Meta_Analysis_SCZ/devgeno0.01_testdevgeno0/matchUKBB_tscore_zscaled_clusterCases_PGmethod_HKmetric.RData'
# type_cluster <- 'Cases'
# type_data <- 'tscore'
# type_sim <- 'HK'
# outFold <- 'OUTPUT_CMC/predict_PGC/200kb/Meta_Analysis_SCZ/devgeno0.01_testdevgeno0/matchUKBB_'
# functR <- '/home/luciat/priler_project/Software/model_clustering/clustering_functions.R'
# type_input <- 'zscaled'
####################################################################################################################

source(functR)
covDat <- list()
sampleAnn <- list()
cl_cohort <- list()

if(length(clusterFile) == 1){
  # cluster computed combining multiple cohorts
  cluster_output <- get(load(clusterFile))
}else{
  cluster_output <- list()  
}

for(i in 1:length(name_cohorts)){
  
  print(name_cohorts[i])
  
  covDat[[i]] <- fread(covDatFile[i], h=T, stringsAsFactor = F, data.table = F)
  sampleAnn[[i]] <- read.table(sampleAnnFile[i], h=T, stringsAsFactors = F, check.names = F)
  
  if(type_cluster == 'Cases'){
    sampleAnn[[i]] <- sampleAnn[[i]][sampleAnn[[i]]$Dx == 1,]
  }else{
    if(type_cluster == 'Controls'){
      sampleAnn[[i]] <- sampleAnn[[i]][sampleAnn[[i]]$Dx == 0,]
    }else{
      if(type_cluster != 'All')
        stop('type_cluster must be either Cases or Controls or All')
    }
  }
  
  sampleAnn[[i]]$cohort <- name_cohorts[i]
  
  if(length(clusterFile) == 1){
    
    intersect_samples <- intersect(sampleAnn[[i]]$Individual_ID, cluster_output$sampleInfo$Individual_ID[cluster_output$sampleInfo$cohort == name_cohorts[i]])
    sampleAnn[[i]] <- sampleAnn[[i]][match(intersect_samples, sampleAnn[[i]]$Individual_ID),]
    name_cl <- ifelse('cl_best' %in% names(cluster_output), 'cl_best', 'cl_new')
    tmp <- cluster_output[which(names(cluster_output) == name_cl)][[1]]
    cl_cohort[[i]] <- tmp[match(sampleAnn[[i]]$Individual_ID,tmp$id),]
    
  }else{
    tmp <- get(load(clusterFile[i]))
    name_cl <- ifelse('cl_best' %in% names(tmp), 'cl_best', 'cl_new')
    cluster_output[[i]] <- tmp[which(names(tmp) == name_cl)][[1]]
    sampleAnn[[i]] <- sampleAnn[[i]][match(cluster_output[[i]]$id,sampleAnn[[i]]$Individual_ID),]
  }
  
  covDat[[i]] <- covDat[[i]][match(sampleAnn[[i]]$Individual_ID, covDat[[i]]$Individual_ID),]
  covDat[[i]] <- covDat[[i]][, !colnames(covDat[[i]]) %in% c('Individual_ID', 'genoSample_ID', 'Dx')]
  
  # rescale for each cohort separately (only continuous)
  covDat[[i]] <- scale(covDat[[i]])
  attr(covDat[[i]], "scaled:center") <- NULL
  attr(covDat[[i]], "scaled:scale") <- NULL
  covDat[[i]] <- as.data.frame(covDat[[i]])
  
}

# combine results
common_cov <- names(which(table(unlist(lapply(sampleAnn, colnames))) == length(name_cohorts)))
sampleAnn_tot <- do.call(rbind, lapply(sampleAnn, function(x) x[, match(common_cov,colnames(x))]))

if(length(clusterFile) > 1){
  cl <- do.call(rbind, cluster_output)
  cl_cohort <- cluster_output
}else{
  cl <- cluster_output[which(names(cluster_output) == name_cl)][[1]]
  cl <- cl[match(sampleAnn_tot$Individual_ID,cl$id),]
}
cl$cohort <- sampleAnn_tot$cohort
gr_names <- sort(unique(cl$gr))
P <- length(gr_names)

### for each cohort perform test separately, then perform meta analysis ###
res_cohort <- list()
for(n in 1:length(name_cohorts)){
  
  print(name_cohorts[n])
  # only includes continuous variables (PCs)
  res_test <- apply(covDat[[n]], 2, function(x) kruskal.test(x = x, g = factor(cl_cohort[[n]]$gr))$p.value)
  res_test <- data.frame(cov = names(res_test), pvalue = res_test)
  res_test$cohort <- name_cohorts[n]
  res_cohort[[n]] <- res_test
  
  fmla  <- as.formula(paste('cov~gr_id'))
  gr_names_cohort <- as.numeric(names(which(table(cl_cohort[[n]]$gr) >= 10)))
  bin_reg <- vector(mode = 'list', length = length(gr_names_cohort))
  for(i in 1:length(gr_names_cohort)){
    
    cov_tmp <- list(covDat[[n]][cl_cohort[[n]]$gr == gr_names_cohort[i],], covDat[[n]][cl_cohort[[n]]$gr != gr_names_cohort[i],])
    gr_id <- factor(c(rep(1, nrow(cov_tmp[[1]])), rep(0, nrow(cov_tmp[[2]]))))
    new <- do.call(rbind, cov_tmp)
    res_comp <- matrix(nrow = ncol(new), ncol = 7)
    
    for(l in 1:ncol(new)){
      tmp_dat <- data.frame(cov = new[, l], gr_id = gr_id)
      res <- tryCatch(glm(fmla, data = tmp_dat, family = 'gaussian'),warning=function(...) NA, error=function(...) NA)
      output <- coef(summary(res))[rownames(coef(summary(res))) == 'gr_id1',1:4]
      output[5:7] <- c(output[1], output[1] + qnorm(0.025)*output[2], output[1] + qnorm(0.975)*output[2])
      res_comp[l, ] <- output 
    }
   
    colnames(res_comp) <- c('beta', 'se_beta', 'z', 'pvalue', 'OR_or_Beta', 'CI_low', 'CI_up')
    res_comp <- as.data.frame(res_comp)
    bin_reg[[i]] <- cbind(data.frame(cov_id = colnames(covDat[[n]])), res_comp)
    bin_reg[[i]]$comp <- sprintf('gr%i_vs_all', gr_names_cohort[i])
  }
  
  res_cohort[[n]] <- do.call(rbind, bin_reg)
  res_cohort[[n]]$cohort <- name_cohorts[n]
}

# summary results from metaAnalysis
cov_id <- unique(unlist(lapply(covDat, function(x) colnames(x))))
comp <- unique(unlist(lapply(res_cohort, function(x) unique(x$comp))))
meta_res <- list()
for(j in 1:length(comp)){
  print(comp[j])
  meta_res[[j]] <- data.frame()
  for(i in 1:length(cov_id)){
    
    tmp <- do.call(rbind,lapply(res_cohort, function(x) x[x$comp == comp[j] & x$cov_id == cov_id[i],]))
    tmp$beta[is.na(tmp$beta) | is.nan(tmp$se_beta) | is.na(tmp$se_beta)] <- 0
    tmp$se_beta[is.nan(tmp$se_beta)] <- NA
    
    tmp_meta <- meta_analysis_res(beta = tmp$beta, se_beta = tmp$se_beta, type_pheno = 'CONTINUOUS')
    tmp_meta <- cbind(data.frame(tmp[1,c('cov_id'), drop = F]), tmp_meta)
    tmp_meta$comp <- comp[j]
    meta_res[[j]] <- rbind(meta_res[[j]], tmp_meta)
  }
  meta_res[[j]]$pval_corr <- p.adjust(meta_res[[j]]$pvalue, method = 'BH')
}

tot_meta_res <- do.call(rbind, meta_res)

### test for cohort: ###
test_cohort <- sapply(gr_names, function(x) chisq.test(x = as.numeric(cl$gr == x), y = cl$cohort)$p.value)
test_cohort <- data.frame(comp = comp, pvalue = test_cohort, test = 'chisq')
test_cohort_tot <- chisq.test(x = cl$gr, y = cl$cohort)

output <- list(meta_analysis = tot_meta_res, single_cohorts = res_cohort, 
               cohort_analysis = list(comp = test_cohort, all = test_cohort_tot),
               covDat = covDat, cl = cl, sampleAnn = sampleAnn)
# save
write.table(x = tot_meta_res, sprintf('%s%s_%s_cluster%s_PGmethod_%smetric_covTest_metaAnalysis.txt', outFold, type_data, type_input, type_cluster, type_sim), col.names = T, row.names = F, sep = '\t', quote = F)
write.table(x = test_cohort, sprintf('%s%s_%s_cluster%s_PGmethod_%smetric_cohortTest.txt', outFold, type_data, type_input, type_cluster, type_sim), col.names = T, row.names = F, sep = '\t', quote = F)
save(output, file = sprintf('%s%s_%s_cluster%s_PGmethod_%smetric_covTest_total.RData', outFold, type_data, type_input, type_cluster, type_sim))




