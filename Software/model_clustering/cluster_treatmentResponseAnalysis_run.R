# treatment response analysis
# perform groupwise

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
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(RColorBrewer))
options(bitmapType = 'cairo', device = 'png')

parser <- ArgumentParser(description="Treatment response analysis")
parser$add_argument("--phenoDatFile", type = "character", help = "file to be loaded (endophenotypes, must be a unique matrix)")
parser$add_argument("--phenoDescFile", type = "character", help = "description endophenotype")
parser$add_argument("--phenoDescCovFile", type = "character", help = "")
parser$add_argument("--covDatFile", type = "character", help = "")
parser$add_argument("--clusterFile", type = "character", help = "file with clustering structure")
parser$add_argument("--type_cluster", type = "character", default = 'All', help = "All, Cases, Controls")
parser$add_argument("--functR", type = "character", help = "functions to be used")
parser$add_argument("--type_data", type = "character", help = "tscore, path_Reactome or path_GO")
parser$add_argument("--type_sim", type = "character", default = 'HK', help = "HK or ED or SNF")
parser$add_argument("--type_input", type = "character", default = 'original', help = "original or zscaled")
parser$add_argument("--outFold", type="character", help = "Output file [basename only]")

args <- parser$parse_args()
phenoDatFile <- args$phenoDatFile
phenoDescFile <- args$phenoDescFile
covDatFile <- args$covDatFile
type_cluster <- args$type_cluster
functR <- args$functR
type_data <- args$type_data
type_sim <- args$type_sim
type_input <- args$type_input
clusterFile <- args$clusterFile
phenoDescCovFile <- args$phenoDescCovFile
outFold <- args$outFold

###################################################################################################################
# tissue_name <- 'Whole_Blood'
# type_cluster <- 'Cases'
# type_data <- 'tscore'
# type_sim <- 'HK'
# type_input <- 'zscaled'
# clusterFile <- sprintf('OUTPUT_GTEx/predict_UKBB/%s/200kb/noGWAS/devgeno0.01_testdevgeno0/Asthma_pheno/Asthma_clustering/tscore_zscaled_clusterCases_PGmethod_HKmetric.RData', tissue_name)
# functR <- '/psycl/g/mpsziller/lucia/priler_project/Software/model_clustering/clustering_functions.R'
# outFold <- sprintf('OUTPUT_GTEx/predict_CAD/%s/200kb/noGWAS/devgeno0.01_testdevgeno0/Asthma_clustering/', tissue_name)
# covDatFile <- 'INPUT_DATA/Covariates/Asthma_clustering/covariateMatrix_Asthma_All_phenoAssoc_FEV1pred.txt'
# phenoDatFile <- 'INPUT_DATA/Covariates/Asthma_clustering/phenotypeMatrix_Asthma_All_phenoAssoc_FEV1pred.txt'
# phenoDescCovFile <- 'INPUT_DATA/Covariates/Asthma_clustering/phenotypeDescription_covariateMatrix_FEV1pred.txt'
# phenoDescFile <- 'INPUT_DATA/Covariates/Asthma_clustering/phenotypeDescription_FEV1pred.txt'
####################################################################################################################

source(functR)
res_tscore <- get(load(clusterFile))
cl_res <- res_tscore$cl_best
P <- length(unique(cl_res$gr))

covDat <- fread(covDatFile, h=T, stringsAsFactors = F, data.table = F)
covDat <- covDat[match(res_tscore$samples_id, covDat$Individual_ID), ]
treat_pheno <- colnames(covDat)[!colnames(covDat) %in% c('Individual_ID', 'Dx', paste0('PC',1:10), 'Age', 'Gender')]
print(treat_pheno)

phenoDat <- fread(phenoDatFile, h=T, stringsAsFactors = F, data.table = F)
phenoDat <- phenoDat[match(res_tscore$samples_id, phenoDat$Individual_ID), ]
phenoInfo <- fread(phenoDescFile, h=T, stringsAsFactors = F, data.table = F)
# consider only certain class of phenotypes
phenoInfo <- phenoInfo[phenoInfo$pheno_type %in% c('Arterial_stiffness', 'Blood_biochemistry', 'Blood_count', 'Blood_pressure', 
                                                   'Blood_count_ratio','Body_size_measures', 'Hand_grip_strength', 'Impedance_measures', 
                                                   'Residential_air_pollution', 'Spirometry'),]
# exclude Nucleated red blood cell
phenoInfo <- phenoInfo[!grepl('Nucleated red blood cell',phenoInfo$Field), ]
phenoInfo <- phenoInfo[!grepl('Traffic intensity on the nearest road',phenoInfo$Field), ]

if(length(treat_pheno) == 1){
  
  phenoInfo_treat <- fread(phenoDescCovFile, h=T, stringsAsFactors = F, data.table = F)
  phenoInfo_treat <- phenoInfo_treat[paste0("p", phenoInfo_treat$pheno_id) %in% treat_pheno, ]
  
  fmla  <- as.formula(paste('pheno~', treat_pheno, '+', paste0(c(paste0('PC',1:10), 'Age', 'Gender'),  collapse = '+')))
  pheno_id <- as.character(phenoInfo$pheno_id)
    
}else{
  
  phenoInfo_treat <- fread(phenoDescCovFile, h=T, stringsAsFactors = F, data.table = F)
  phenoInfo_treat <- rbind(phenoInfo_treat,
                           data.frame(pheno_id = c('6153_6177_1', '6153_6177_2', '6153_6177_3'), FieldID = c('6153_6177', '6153_6177', '6153_6177'),
                                      Field = rep('Medication for cholesterol, blood pressure or diabetes', 3), Path = rep(NA, 3), Strata  = rep(NA, 3),
                                      Sexed = rep('Unisex', 3), Coding = rep(NA, 3), Coding_meaning = c('Cholesterol lowering medication', 'Blood pressure medication', 'Insulin'),
                                      original_type = rep('CAT_MULTIPLE', 3), transformed_type = rep('CAT_MUL_BINARY_VAR',3), nsamples = rep(NA,3), nsamples_T= rep(NA,3),
                                      nsamples_F= rep(NA,3), pheno_type = rep('Medication', 3)))
  
  fmla  <- as.formula(paste('pheno~', paste0(treat_pheno, collapse = '+'), '+', paste0(c(paste0('PC',1:10), 'Age', 'Gender'),  collapse = '+')))
  pheno_id <- as.character(phenoInfo$pheno_id)
  if('12144der' %in% phenoInfo_treat$pheno_id){
    phenoInfo_treat$Coding_meaning[phenoInfo_treat$pheno_id == '12144der'] <- 'Height derived'
  }
  
}

res_diff <- list()
for(p in 1:P){
  
  print(paste0('#################### gr', p, ' ####################'))
  
  id <- which(cl_res$gr == p)
  id_not <- which(cl_res$gr != p)
  
  df <- data.frame(treat_id = c(), z_diff = c(), pvalue_diff = c(), gr_beta = c(), gr_se_beta = c(), gr_ORorBeta = c(), gr_CI_low = c(), gr_CI_up = c(), 
                   gr_pvalue = c(),  notgr_beta = c(), notgr_se_beta = c(), notgr_ORorBeta = c(), notgr_CI_low = c(), notgr_CI_up = c(), 
                   notgr_pvalue = c(),  comp =c(), pheno_id = c())
  
  for(l in 1:length(pheno_id)){
    
    type_pheno <- phenoInfo$transformed_type[phenoInfo$pheno_id == pheno_id[l]]
    tmp_dat <- cbind(data.frame(pheno = phenoDat[id, pheno_id[l]]), covDat[id,!colnames(covDat) %in% c('Individual_ID', 'Dx')])
    
    if(sum(!is.na(tmp_dat$pheno))> 300){
      print(l)
      # standardize Y if continous and PCs cov
      tmp_dat[, colnames(tmp_dat) %in% paste0('PC', 1:10)] <- scale(tmp_dat[, colnames(tmp_dat) %in% paste0('PC', 1:10)])
      if(type_pheno == 'CONTINUOUS'){
        tmp_dat$pheno <- scale(tmp_dat$pheno)
      }
      
      tmp_gr <- compute_reg_endopheno_multi(mat = tmp_dat, fmla = fmla, cov_int = treat_pheno, type_pheno = type_pheno)
      
      tmp_dat <- cbind(data.frame(pheno = phenoDat[id_not, pheno_id[l]]), covDat[id_not,!colnames(covDat) %in% c('Individual_ID', 'Dx')])
      # standardize Y if continous and PCs cov
      tmp_dat[, colnames(tmp_dat) %in% paste0('PC', 1:10)] <- scale(tmp_dat[, colnames(tmp_dat) %in% paste0('PC', 1:10)])
      if(type_pheno == 'CONTINUOUS'){
        tmp_dat$pheno <- scale(tmp_dat$pheno)
      }
      tmp_notgr <- compute_reg_endopheno_multi(mat = tmp_dat, fmla = fmla, cov_int = treat_pheno, type_pheno = type_pheno)
      
      # keep only treat that are significant in at least 1 (and regression computed)
      if(!(is.null(tmp_gr) | is.null(tmp_notgr))){
        test_diff <- (tmp_gr$beta - tmp_notgr$beta)/sqrt(tmp_gr$se_beta^2 + tmp_notgr$se_beta^2)
        pval_test <- 2*pnorm(-abs(test_diff))
        tmp <- data.frame(treat_id = rownames(tmp_gr), z_diff = test_diff,pvalue_diff = pval_test, 
                          gr_beta = tmp_gr$beta, gr_se_beta = tmp_gr$se_beta, 
                          gr_ORorBeta = tmp_gr$OR_or_beta, gr_CI_low = tmp_gr$CI_low, gr_CI_up = tmp_gr$CI_up, 
                          gr_pvalue = tmp_gr$pvalue, 
                          notgr_beta = tmp_notgr$beta, notgr_se_beta = tmp_notgr$se_beta, 
                          notgr_ORorBeta = tmp_notgr$OR_or_beta, notgr_CI_low = tmp_notgr$CI_low, notgr_CI_up = tmp_notgr$CI_up, 
                          notgr_pvalue = tmp_notgr$pvalue)
        tmp$comp <- sprintf('gr%i_vs_all', p)
        tmp$pheno_id <- pheno_id[l]
        df <- rbind(df, tmp)
      }
    }
  }
  
  # df$pvalue_diff_corr_overall <- p.adjust(df$pvalue_diff, method = 'BH')
  # tmp <- lapply(treat_pheno, function(x) p.adjust(df$pvalue_diff[df$treat_id == x], method = 'BH'))
  # tmp_id <- lapply(treat_pheno, function(x) which(df$treat_id == x))
  # df$pvalue_diff_corr <- unlist(tmp[unlist(tmp_id)])
  df$treat_Field <- phenoInfo_treat$Field[match(sapply(df$treat_id, function(x) paste0(strsplit(x, split = 'p')[[1]][-1], collapse='p')), phenoInfo_treat$pheno_id)]
  df$treat_meaning <- phenoInfo_treat$Coding_meaning[match(sapply(df$treat_id, function(x) paste0(strsplit(x, split = 'p')[[1]][-1], collapse='p')), phenoInfo_treat$pheno_id)]
  df$pheno_Field <- phenoInfo$Field[match(df$pheno_id, phenoInfo$pheno_id)]
  df$pheno_meaning <- phenoInfo$Coding_meaning[match(df$pheno_id, phenoInfo$pheno_id)]
  df$pheno_type <- phenoInfo$transformed_type[match(df$pheno_id, phenoInfo$pheno_id)]
  df$pheno_class <- phenoInfo$pheno_type[match(df$pheno_id, phenoInfo$pheno_id)]
  res_diff[[p]] <- df
  
}

res_tot_diff <- do.call(rbind, res_diff)
# how to correct pvalues?
# separetly for each treatment and comparison
res_tot_diff$pvalue_corr_diff <- NA
for(i in 1:length(unique(res_tot_diff$treat_meaning))){
  for(j in 1:length(unique(res_tot_diff$comp))){
    res_tot_diff$pvalue_corr_diff[res_tot_diff$treat_meaning == unique(res_tot_diff$treat_meaning)[i] 
                                  & res_tot_diff$comp == unique(res_tot_diff$comp)[j]] <- p.adjust(res_tot_diff$pvalue_diff[res_tot_diff$treat_meaning == unique(res_tot_diff$treat_meaning)[i] 
                                                                                                                            & res_tot_diff$comp == unique(res_tot_diff$comp)[j]], method = 'BH')
  }
}

# save results
write.table(x = res_tot_diff, file = sprintf('%s%s_%s_cluster%s_TreatResponse.txt',outFold, type_data, type_input, type_cluster), 
            col.names = T, row.names = F, sep = '\t', quote = F)

# test groupwise differences
res_tot <- list()
for(p in 1:(P-1)){
  
  res_tot[[p]] <- list()
  tmp_2 <- res_tot_diff[res_tot_diff$comp == sprintf('gr%i_vs_all', p),]
  tmp_2$new_id <- paste0(tmp_2$treat_id, '_and_', tmp_2$pheno_id)
  
  for(q in (p+1):P){
    
    tmp_1 <- res_tot_diff[res_tot_diff$comp == sprintf('gr%i_vs_all', q),]
    tmp_1$new_id <- paste0(tmp_1$treat_id, '_and_', tmp_1$pheno_id)
    # common set:
    common_id <- intersect( tmp_1$new_id ,  tmp_2$new_id)
    test_diff <- (tmp_1$gr_beta[match(common_id, tmp_1$new_id)] - tmp_2$gr_beta[match(common_id, tmp_2$new_id)])/sqrt(tmp_1$gr_se_beta[match(common_id, tmp_1$new_id)]^2 + tmp_2$gr_se_beta[match(common_id, tmp_2$new_id)]^2)
    pval_test <- 2*pnorm(-abs(test_diff))
    res_tot[[p]][[q-p]] <- data.frame(treat_id = tmp_1$treat_id[match(common_id, tmp_1$new_id)] , z_diff = test_diff, pvalue_diff = pval_test, 
                                      gr1_beta = tmp_1$gr_beta[match(common_id, tmp_1$new_id)], gr1_se_beta = tmp_1$gr_se_beta[match(common_id, tmp_1$new_id)],
                                      gr1_ORorBeta = tmp_1$gr_ORorBeta[match(common_id, tmp_1$new_id)], 
                                      gr1_CI_low = tmp_1$gr_CI_low[match(common_id, tmp_1$new_id)],  
                                      gr1_CI_up = tmp_1$gr_CI_up[match(common_id, tmp_1$new_id)],  
                                      gr1_pvalue = tmp_1$gr_pvalue[match(common_id, tmp_1$new_id)], 
                                      gr2_beta = tmp_2$gr_beta[match(common_id, tmp_2$new_id)], gr2_se_beta = tmp_2$gr_se_beta[match(common_id, tmp_2$new_id)],
                                      gr2_ORorBeta = tmp_2$gr_ORorBeta[match(common_id, tmp_2$new_id)], 
                                      gr2_CI_low = tmp_2$gr_CI_low[match(common_id, tmp_2$new_id)],  
                                      gr2_CI_up = tmp_2$gr_CI_up[match(common_id, tmp_2$new_id)],  
                                      gr2_pvalue = tmp_2$gr_pvalue[match(common_id, tmp_2$new_id)])
    res_tot[[p]][[q-p]]$comp <- sprintf('gr%i_vs_gr%i', q, p)
    res_tot[[p]][[q-p]]$treat_Field <- tmp_1$treat_Field[match(common_id, tmp_1$new_id)]
    res_tot[[p]][[q-p]]$treat_meaning <- tmp_1$treat_meaning[match(common_id, tmp_1$new_id)]
    res_tot[[p]][[q-p]]$pheno_id <- tmp_1$pheno_id[match(common_id, tmp_1$new_id)]
    res_tot[[p]][[q-p]]$pheno_Field <- tmp_1$pheno_Field[match(common_id, tmp_1$new_id)]
    res_tot[[p]][[q-p]]$pheno_meaning <- tmp_1$pheno_meaning[match(common_id, tmp_1$new_id)]
    res_tot[[p]][[q-p]]$pheno_type <- tmp_1$pheno_type[match(common_id, tmp_1$new_id)]
    res_tot[[p]][[q-p]]$pheno_class <- tmp_1$pheno_class[match(common_id, tmp_1$new_id)]
    
  }
}

res_tot_gr <- do.call(rbind, do.call(c,res_tot))
res_tot_gr$pvalue_corr_diff <- NA
for(i in 1:length(unique(res_tot_gr$treat_meaning))){
  for(j in 1:length(unique(res_tot_gr$comp))){
    res_tot_gr$pvalue_corr_diff[res_tot_gr$treat_meaning == unique(res_tot_gr$treat_meaning)[i] 
                                & res_tot_gr$comp == unique(res_tot_gr$comp)[j]] <- p.adjust(res_tot_gr$pvalue_diff[res_tot_gr$treat_meaning == unique(res_tot_gr$treat_meaning)[i] 
                                                                                                                    & res_tot_gr$comp == unique(res_tot_gr$comp)[j]], method = 'BH')
  }
}
write.table(x = res_tot_gr, file = sprintf('%s%s_%s_cluster%s_TreatResponse_pairwise.txt',outFold, type_data, type_input, type_cluster), 
            col.names = T, row.names = F, sep = '\t', quote = F)

