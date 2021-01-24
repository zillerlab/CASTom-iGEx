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

parser <- ArgumentParser(description="Associate cluster to endophenotype, multiple cohorts analysed together (meta analysis)")
parser$add_argument("--name_cohorts", type = "character", nargs = '*', help = "name cohorts")
parser$add_argument("--phenoDatFile", type = "character", nargs = '*', help = "file to be loaded (endophenotypes, must be a unique matrix)")
parser$add_argument("--phenoDescFile", type = "character", help = "description endophenotype")
parser$add_argument("--sampleAnnFile", type = "character", nargs = '*', help = "file with samples to be used")
parser$add_argument("--clusterFile", type = "character", nargs = '*', help = "file with clustering structure")
parser$add_argument("--type_cluster", type = "character", default = 'All', help = "All, Cases, Controls")
parser$add_argument("--functR", type = "character", help = "functions to be used")
parser$add_argument("--type_data", type = "character", help = "tscore, path_Reactome or path_GO")
parser$add_argument("--type_sim", type = "character", default = 'HK', help = "HK or ED or SNF")
parser$add_argument("--type_input", type = "character", default = 'original', help = "original or zscaled")
parser$add_argument("--risk_score", type = "logical", default = F, help = "if true, phenotype is risk score")
parser$add_argument("--outFold", type="character", help = "Output file [basename only]")

args <- parser$parse_args()
name_cohorts <- args$name_cohorts
phenoDatFile <- args$phenoDatFile
phenoDescFile <- args$phenoDescFile
sampleAnnFile <- args$sampleAnnFile
type_cluster <- args$type_cluster
functR <- args$functR
type_data <- args$type_data
type_sim <- args$type_sim
type_input <- args$type_input
clusterFile <- args$clusterFile
risk_score <- args$risk_score
outFold <- args$outFold

####################################################################################################################
# name_cohorts <- c('German1','German2', 'German3', 'German4', 'German5')
# phenoDatFile <- paste0('OUTPUT_GTEx/predict_CAD/Liver/200kb/CAD_GWAS_bin5e-2/',name_cohorts,'/devgeno0.01_testdevgeno0/CAD_HARD_clustering/tscore_corrThr0.5_risk_score_relatedPhenotypes.txt')
# phenoDescFile <- '/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/INPUT_DATA_GTEx/CAD/Covariates/UKBB/CAD_HARD_clustering/phenotypeDescription.txt'
# sampleAnnFile <- paste0('/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/INPUT_DATA_GTEx/CAD/Covariates/',name_cohorts,'/covariateMatrix.txt')
# clusterFile <- paste0('/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/OUTPUT_GTEx/predict_CAD/Liver/200kb/CAD_GWAS_bin5e-2/',name_cohorts,'/devgeno0.01_testdevgeno0/CAD_HARD_clustering/tscore_zscaled_predictClusterCases_PGmethod_HKmetric.RData')
# type_cluster <- 'Cases'
# type_data <- 'tscore'
# type_sim <- 'HK'
# outFold <- '/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/OUTPUT_GTEx/predict_CAD/Liver/200kb/CAD_GWAS_bin5e-2/Meta_Analysis_CAD/CAD_HARD_clustering/'
# functR <- '/psycl/g/mpsziller/lucia/priler_project/Software/model_clustering/clustering_functions.R'
# type_input <- 'zscaled'
####################################################################################################################

####################################################################################################################
# name_cohorts <- read.table('INPUT_DATA/SCZ_cohort_names_CLUST')$V1
# phenoDatFile <- paste0('OUTPUT_CMC/predict_PGC/200kb/',name_cohorts,'/devgeno0.01_testdevgeno0/tscore_corrThr0.5_risk_score_relatedPhenotypes.txt')
# phenoDescFile <- '/home/luciat/UKBB_SCZrelated/phenotypeDescription_rsSCZ.txt'
# sampleAnnFile <- paste0('INPUT_DATA/Covariates/',name_cohorts,'.covariateMatrix_old.txt')
# clusterFile <- 'OUTPUT_CMC/predict_PGC/200kb/Meta_Analysis_SCZ/devgeno0.01_testdevgeno0/tscore_zscaled_clusterCases_PGmethod_HKmetric.RData'
# type_cluster <- 'Cases'
# type_data <- 'tscore'
# type_sim <- 'HK'
# outFold <- './'
# functR <- '/home/luciat/priler_project/Software/model_clustering/clustering_functions.R'
# type_input <- 'zscaled'
####################################################################################################################

source(functR)
phenoDat <- list()
sampleAnn <- list()
cl_cohort <- list()
phenoInfo <- read.delim(phenoDescFile, h=T, stringsAsFactors = F, sep = '\t')

if(length(clusterFile) == 1){
  # cluster computed combining multiple cohorts
  cluster_output <- get(load(clusterFile))
}else{
  cluster_output <- list()  
}

for(i in 1:length(name_cohorts)){
  
  print(name_cohorts[i])
  
  phenoDat[[i]] <- fread(phenoDatFile[i], h=T, stringsAsFactor = F, data.table = F)
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
  
  phenoDat[[i]] <- phenoDat[[i]][match(sampleAnn[[i]]$Individual_ID, phenoDat[[i]]$Individual_ID),]
  phenoDat[[i]] <- phenoDat[[i]][, -1]
  
  # get common pheno id (phenoDescFile can be used to filter)
  common_pheno <- intersect(phenoInfo$pheno_id, colnames(phenoDat[[i]]))
  phenoDat[[i]] <- phenoDat[[i]][,match(common_pheno, colnames(phenoDat[[i]]))]
  
  if(risk_score){
    phenoDat[[i]] <- scale(phenoDat[[i]])
    attr(phenoDat[[i]], "scaled:center") <- NULL
    attr(phenoDat[[i]], "scaled:scale") <- NULL
    phenoDat[[i]] <- as.data.frame(phenoDat[[i]])
  }
  
}

if(risk_score){
  phenoInfo$transformed_type <- 'CONTINUOUS'
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

if(any(phenoInfo$Path %in% 'Online follow-up > Cognitive function online > Fluid intelligence')){
  phenoInfo$Field[phenoInfo$Path == 'Online follow-up > Cognitive function online > Fluid intelligence'] <- paste(phenoInfo$Field[phenoInfo$Path == 'Online follow-up > Cognitive function online > Fluid intelligence'], '(Online)')
}

### for each cohort perform endophenotype association separately, then perform meta analysis ###
reg_cohort <- list()
for(n in 1:length(name_cohorts)){
  
  print(name_cohorts[n])
  phenoInfo_cohort <- phenoInfo[match(colnames(phenoDat[[n]]), phenoInfo$pheno_id),]
  covDat <- sampleAnn[[n]][, !colnames(sampleAnn[[n]]) %in% c('Individual_ID', 'genoSample_ID', 'Dx')]
  fmla  <- as.formula(paste('pheno~gr_id+', paste0(colnames(covDat)[!colnames(covDat) %in% 'cohort'], collapse = '+')))
  
  gr_names_cohort <- as.numeric(names(which(table(cl_cohort[[n]]$gr) >= 10)))
  
  ########################################
  #### binary regression (gi vs all) #####
  ########################################
  
  bin_reg <- vector(mode = 'list', length = length(gr_names_cohort))
  for(i in 1:length(gr_names_cohort)){
    
    print(paste0('group', gr_names_cohort[i], '_vs_all'))
    
    # j vs all
    pheno_case_tmp <- list(phenoDat[[n]][cl_cohort[[n]]$gr == gr_names_cohort[i],], phenoDat[[n]][cl_cohort[[n]]$gr != gr_names_cohort[i],])
    covDat_tmp <- list(covDat[cl_cohort[[n]]$gr == gr_names_cohort[i],], covDat[cl_cohort[[n]]$gr != gr_names_cohort[i],]) 
    
    new <- do.call(rbind, pheno_case_tmp)
    colnames(new) <- paste0('p', colnames(new))
    gr_id <- factor(c(rep(1, nrow(pheno_case_tmp[[1]])), rep(0, nrow(pheno_case_tmp[[2]]))))
    
    new_cov <- do.call(rbind, covDat_tmp)
    
    res_glm <- matrix(nrow = ncol(new), ncol = 7)
    for(l in 1:ncol(new)){
      # print(l)
      type_pheno <- phenoInfo_cohort$transformed_type[paste0('p',phenoInfo_cohort$pheno_id) == colnames(new)[l]]
      tmp_dat <- cbind(data.frame(pheno = new[, l], gr_id = gr_id), new_cov)
      res_glm[l,] <- compute_reg_endopheno(mat = tmp_dat, fmla = fmla, type_pheno = type_pheno)
    }
    
    colnames(res_glm) <- c('beta', 'se_beta', 'z', 'pvalue', 'OR_or_Beta', 'CI_low', 'CI_up')
    res_glm <- as.data.frame(res_glm)
    res_glm$type_pheno <- phenoInfo_cohort$transformed_type[match(colnames(new), paste0('p',phenoInfo_cohort$pheno_id))]
    
    phenoInfo_tmp <- phenoInfo_cohort[match(colnames(new), paste0('p',phenoInfo_cohort$pheno_id)),]
    
    bin_reg[[i]] <- cbind(data.frame(pheno_id = phenoInfo_tmp$pheno_id, Field = phenoInfo_tmp$Field, meaning = phenoInfo_tmp$Coding_meaning), res_glm)
    bin_reg[[i]]$comp <- sprintf('gr%i_vs_all', gr_names_cohort[i])
    
  }
  
  reg_cohort[[n]] <- do.call(rbind, bin_reg)
  reg_cohort[[n]]$cohort <- name_cohorts[n]
  
}

# summary results from metaAnalysis
pheno_id <- unique(unlist(lapply(phenoDat, function(x) colnames(x))))
comp <- unique(unlist(lapply(reg_cohort, function(x) unique(x$comp))))
meta_res <- list()
for(j in 1:length(comp)){
  print(comp[j])
  meta_res[[j]] <- data.frame()
  for(i in 1:length(pheno_id)){
    
    tmp <- do.call(rbind,lapply(reg_cohort, function(x) x[x$comp == comp[j] & x$pheno_id == pheno_id[i],]))
    tmp$beta[is.na(tmp$beta) | is.nan(tmp$se_beta) | is.na(tmp$se_beta)] <- 0
    tmp$se_beta[is.nan(tmp$se_beta)] <- NA
    
    tmp_meta <- meta_analysis_res(beta = tmp$beta, se_beta = tmp$se_beta, type_pheno = unique(tmp$type_pheno))
    tmp_meta <- cbind(data.frame(tmp[1,c('pheno_id','Field','meaning')]), tmp_meta)
    tmp_meta$type_pheno <- unique(tmp$type_pheno)
    tmp_meta$comp <- comp[j]
    meta_res[[j]] <- rbind(meta_res[[j]], tmp_meta)
  }
  meta_res[[j]]$pval_corr <- p.adjust(meta_res[[j]]$pvalue, method = 'BH')
}

tot_meta_res <- do.call(rbind, meta_res)

output <- list(meta_analysis = tot_meta_res, single_cohorts = reg_cohort, phenoDat = phenoDat, phenoInfo = phenoInfo, cl = cl, sampleAnn = sampleAnn)
# save
write.table(x = tot_meta_res, sprintf('%s%s_%s_cluster%s_PGmethod_%smetric_phenoAssociation_GLM_metaAnalysis.txt', outFold, type_data, type_input, type_cluster, type_sim), col.names = T, row.names = F, sep = '\t', quote = F)
save(output, file = sprintf('%s%s_%s_cluster%s_PGmethod_%smetric_phenoAssociation_GLM_metaAnalysis.RData', outFold, type_data, type_input, type_cluster, type_sim))

####
reg_cohort <- list()
for(n in 1:length(name_cohorts)){
  
  print(name_cohorts[n])
  phenoInfo_cohort <- phenoInfo[match(colnames(phenoDat[[n]]), phenoInfo$pheno_id),]
  covDat <- sampleAnn[[n]][, !colnames(sampleAnn[[n]]) %in% c('Individual_ID', 'genoSample_ID', 'Dx')]
  fmla  <- as.formula(paste('pheno~gr_id+', paste0(colnames(covDat)[!colnames(covDat) %in% 'cohort'], collapse = '+')))
  
  gr_names_cohort <- as.numeric(names(which(table(cl_cohort[[n]]$gr) >= 10)))
  
  #######################################
  #### binary regression (gi vs gj) #####
  #######################################
  bin_reg <- vector(mode = 'list', length = length(gr_names_cohort)-1)
  
  for(i in 1:(length(gr_names_cohort)-1)){
    
    print(paste0('group', gr_names_cohort[i], '_vs_groupj'))
    
    # j vs all
    pheno_case_tmp <- lapply(gr_names_cohort[i:length(gr_names_cohort)], function(x) phenoDat[[n]][cl_cohort[[n]]$gr == x,])
    covDat_tmp <- lapply(gr_names_cohort[i:length(gr_names_cohort)], function(x) covDat[cl_cohort[[n]]$gr == x,])
    bin_reg[[i]] <-  vector(mode = 'list', length = length(pheno_case_tmp)-1)
    
    for(j in 2:length(pheno_case_tmp)){
      
      print(j)
      
      new <- rbind(pheno_case_tmp[[1]], pheno_case_tmp[[j]])
      colnames(new) <- paste0('p', colnames(new))
      gr_id <- factor(c(rep(0, nrow(pheno_case_tmp[[1]])), rep(1, nrow(pheno_case_tmp[[j]]))))
      
      new_cov <- rbind(covDat_tmp[[1]], covDat_tmp[[j]])
      res_glm <- matrix(nrow = ncol(new), ncol = 7)
      for(l in 1:ncol(new)){
        # print(l)
        type_pheno <- phenoInfo_cohort$transformed_type[paste0('p',phenoInfo_cohort$pheno_id) == colnames(new)[l]]
        tmp_dat <- cbind(data.frame(pheno = new[, l], gr_id = gr_id), new_cov)
        res_glm[l,] <- compute_reg_endopheno(mat = tmp_dat, fmla = fmla, type_pheno = type_pheno)
      }
      
      colnames(res_glm) <- c('beta', 'se_beta', 'z', 'pvalue', 'OR_or_Beta', 'CI_low', 'CI_up')
      res_glm <- as.data.frame(res_glm)
      res_glm$type_pheno <- phenoInfo_cohort$transformed_type[match(colnames(new), paste0('p',phenoInfo_cohort$pheno_id))]
      phenoInfo_tmp <- phenoInfo_cohort[match(colnames(new), paste0('p',phenoInfo_cohort$pheno_id)),]
      
      bin_reg[[i]][[j-1]] <- cbind(data.frame(pheno_id = phenoInfo_tmp$pheno_id, Field = phenoInfo_tmp$Field, meaning = phenoInfo_tmp$Coding_meaning), res_glm)
      bin_reg[[i]][[j-1]]$comp <- sprintf('gr%i_vs_gr%i', gr_names_cohort[i:length(gr_names_cohort)][j], gr_names_cohort[i])
      
    }
    
  }
  
  reg_cohort[[n]] <- do.call(rbind, do.call(c,bin_reg))
  reg_cohort[[n]]$cohort <- name_cohorts[n]
}  


# summary results from metaAnalysis
pheno_id <- unique(unlist(lapply(phenoDat, function(x) colnames(x))))
comp <- unique(unlist(lapply(reg_cohort, function(x) unique(x$comp))))
meta_res <- list()
for(j in 1:length(comp)){
  print(comp[j])
  meta_res[[j]] <- data.frame()
  for(i in 1:length(pheno_id)){
    
    tmp <- do.call(rbind,lapply(reg_cohort, function(x) x[x$comp == comp[j] & x$pheno_id == pheno_id[i],]))
    tmp$beta[is.na(tmp$beta) | is.nan(tmp$se_beta) | is.na(tmp$se_beta)] <- 0
    tmp$se_beta[is.nan(tmp$se_beta)] <- NA
    
    tmp_meta <- meta_analysis_res(beta = tmp$beta, se_beta = tmp$se_beta, type_pheno = unique(tmp$type_pheno))
    tmp_meta <- cbind(data.frame(tmp[1,c('pheno_id','Field','meaning')]), tmp_meta)
    tmp_meta$type_pheno <- unique(tmp$type_pheno)
    tmp_meta$comp <- comp[j]
    meta_res[[j]] <- rbind(meta_res[[j]], tmp_meta)
  }
  meta_res[[j]]$pval_corr <- p.adjust(meta_res[[j]]$pvalue, method = 'BH')
}

tot_meta_res <- do.call(rbind, meta_res)

output <- list(meta_analysis = tot_meta_res, single_cohorts = reg_cohort, phenoDat = phenoDat, phenoInfo = phenoInfo, cl = cl, sampleAnn = sampleAnn)
# save
write.table(x = tot_meta_res, sprintf('%s%s_%s_cluster%s_PGmethod_%smetric_phenoAssociation_GLMpairwise_metaAnalysis.txt', outFold, type_data, type_input, type_cluster, type_sim), col.names = T, row.names = F, sep = '\t', quote = F)
save(output, file = sprintf('%s%s_%s_cluster%s_PGmethod_%smetric_phenoAssociation_GLMpairwise_metaAnalysis.RData', outFold, type_data, type_input, type_cluster, type_sim))



