# Endophenotype when phenoInfo not available

options(stringsAsFactors=F)
options(max.print=1000)
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(qvalue))
suppressPackageStartupMessages(library(pROC))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(SparseM))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(MASS))
suppressPackageStartupMessages(library(lmtest))
suppressPackageStartupMessages(library(ggsci))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(RNOmni))
options(bitmapType = 'cairo', device = 'png')

parser <- ArgumentParser(description="Endophenotype differences on a predicted clustering, phenoInfo not available")
parser$add_argument("--cohort_name", nargs = '*', type = "character", help = "")
parser$add_argument("--phenoNew_file", nargs = '*', type = "character", help = "")
parser$add_argument("--covNew_file", nargs = '*', default = NULL, type = "character", help = "")
parser$add_argument("--type_cluster", type = "character",default = 'All', help = "All, Cases, Controls")
parser$add_argument("--type_data", type = "character", help = "pathway or tscore")
parser$add_argument("--model_name", type = "character", help = "")
parser$add_argument("--clustFile_new", type = "character", nargs = '*', help = "file cluster results")
parser$add_argument("--functR", type = "character", help = "functions to be used")
parser$add_argument("--type_input", type = "character", default = 'original', help = "original or zscaled")
parser$add_argument("--outFold", type="character", help = "Output file [basename only]")

args <- parser$parse_args()
cohort_name <- args$cohort_name
clustFile_new <- args$clustFile_new
functR <- args$functR
type_data <- args$type_data
type_input <- args$type_input
type_cluster <- args$type_cluster
phenoNew_file <- args$phenoNew_file
covNew_file <- args$covNew_file
model_name <- args$model_name
outFold <- args$outFold

###################################################################################################################
# functR <- '/psycl/g/mpsziller/lucia/castom-igex/Software/model_clustering/clustering_functions.R'
# cohort_name = 'SHIP-TREND'
# type_data <- 'tscore_corrPCs'
# type_input <- 'zscaled'
# type_cluster <- 'Cases'
# clustFile_new <- 'Results/PriLer/Liver/tscore_corrPCs_zscaled_predictClusterCases_PGmethod_HKmetric.RData'
# phenoNew_file <-  'Results/PriLer/pheno_test.txt'
# covNew_file <- 'Results/PriLer/SHIP-TREND_gPC_SHIP_2022_27_withSex.txt'
# outFold <- 'Results/PriLer/Liver/'
#################################################################################################################

source(functR)

sampleAnn_new <- list()
phenoDat_new <- list()
covDat_new <- list()
clust_new <- list()
phenoInfo_new <- list()
pairwise_bin_reg <- list()
tot_bin_reg <- list()


for(i in 1:length(cohort_name)){
  
  print(cohort_name[i])
  tmp <- get(load(clustFile_new[i]))
  
  if(any(table(tmp$cl_new$gr) < 3)){
    rm_gr <- names(which(table(tmp$cl_new$gr) < 3))
    tmp$cl_new <- tmp$cl_new[!tmp$cl_new$gr %in% rm_gr, ]
    tmp$samples_id <- tmp$cl_new$id
    tmp$sampleAnn <- tmp$sampleAnn[match(tmp$samples_id, tmp$sampleAnn$Individual_ID),]
    tmp$data_new <- tmp$data_new[match(tmp$samples_id, rownames(tmp$data_new)),]
  }
  
  sampleAnn_new[[i]] <- tmp$sampleAnn
  clust_new[[i]] <- tmp$cl_new
  P <- length(unique(clust_new[[i]]$gr))
  
  
  #### endophenotype association ####
  ################
  ## gri vs grj ##
  ################
  
  phenoDat_new[[i]] <- read.table(phenoNew_file[i], h=T, stringsAsFactors = F)
  common_samples <- intersect(phenoDat_new[[i]]$Individual_ID, sampleAnn_new[[i]]$Individual_ID)
  if(!is.null(covNew_file[[i]])){
    covDat_new[[i]] <- read.table(covNew_file[i], h=T, stringsAsFactors = F)
    common_samples <- intersect(common_samples, covDat_new[[i]]$Individual_ID)
    covDat_new[[i]] <- covDat_new[[i]][match(common_samples, covDat_new[[i]]$Individual_ID),]
    covDat_new[[i]] <- covDat_new[[i]][, !colnames(covDat_new[[i]]) %in% 
                                            c('Individual_ID','RNASample_ID', 'genoSample_ID', 'Dx', 
                                              'Age', 'Gender', 'Sex', 'Array')]
  }
  
  phenoDat_new[[i]] <- phenoDat_new[[i]][match(common_samples, phenoDat_new[[i]]$Individual_ID),]
  sampleAnn_new[[i]] <- sampleAnn_new[[i]][match(common_samples, sampleAnn_new[[i]]$Individual_ID),]
  
  sampleAnn_new[[i]] <- sampleAnn_new[[i]][, colnames(sampleAnn_new[[i]]) %in% 
                                             c('Individual_ID', paste0('C', 1:10), paste0('PC', 1:10), 
                                               'Gender', 'Sex','Age', 'Array')]
  if('Sex' %in% colnames(sampleAnn_new[[i]])){
    sampleAnn_new[[i]]$Sex <- as.character(sampleAnn_new[[i]]$Sex)
    }
  if('Gender' %in% colnames(sampleAnn_new[[i]])){
    sampleAnn_new[[i]]$Gender <- as.character(sampleAnn_new[[i]]$Gender)
    }
  
  
  
  gr_names <- sort(unique(clust_new[[i]]$gr))
  cl <- clust_new[[i]]$gr
  
  covDat <- sampleAnn_new[[i]][, !colnames(sampleAnn_new[[i]]) %in% c('Individual_ID', 'genoSample_ID', 'Dx')]
  if(!is.null(covNew_file[[i]])){
    covDat <- cbind(covDat, covDat_new[[i]])
    covDat <- covDat[, !duplicated(colnames(covDat))]
  }
  
  fmla  <- as.formula(paste('pheno~gr_id+', paste0(colnames(covDat), collapse = '+')))
  phenoDat <- phenoDat_new[[i]][, !colnames(phenoDat_new[[i]]) %in%  
                                  c('Individual_ID', 'RNASample_ID', 'genoSample_ID', 'Dx',
                                    'Age', 'Gender', 'Sex', 'Array')]
  
  if(any(table(cl)<=10)){
    rm_id <- names(which(table(cl)<=10))
    P <- P-length(rm_id)
    gr_names <- gr_names[!gr_names %in% rm_id]
    covDat <- covDat[!cl %in% rm_id,]
    phenoDat <- phenoDat[!cl %in% rm_id,]
    cl <- cl[!cl %in% rm_id]
  }
   
  phenoInfo_new[[i]] <- data.frame(pheno_id = colnames(phenoDat))
  phenoInfo_new[[i]]$type_pheno <- 'CONTINUOUS'
  phenoInfo_new[[i]]$type_pheno[sapply(1:ncol(phenoDat), function(x) is.integer(phenoDat[,x]) & length(unique(na.omit(phenoDat[,x]))) == 2)] <- 'CAT_SINGLE_BINARY'
  phenoInfo_new[[i]]$type_pheno[sapply(1:ncol(phenoDat), function(x) is.integer(phenoDat[,x]) & length(unique(na.omit(phenoDat[,x]))) > 2)] <- 'CAT_ORD'
  for(j in 1:ncol(phenoDat)){
    if(phenoInfo_new[[i]]$type_pheno[j] == 'CONTINUOUS'){
      tmp <- phenoDat[!is.na(phenoDat[,j]),j]
      phenoDat[!is.na(phenoDat[,j]),j] <- RankNorm(tmp)
    }
  }
  
  bin_reg <- vector(mode = 'list', length = length(gr_names)-1)
  for(k in 1:(length(gr_names)-1)){
    
    print(paste0('group', gr_names[k], '_vs_groupj'))
    
    # j vs all
    pheno_case_tmp <- lapply(gr_names[k:length(gr_names)], function(x) phenoDat[cl == x,])
    covDat_tmp <- lapply(gr_names[k:length(gr_names)], function(x) covDat[cl == x,])
    bin_reg[[k]] <-  vector(mode = 'list', length = length(pheno_case_tmp)-1)
    
    for(j in 2:length(pheno_case_tmp)){
      
      print(j)
      
      new <- rbind(pheno_case_tmp[[1]], pheno_case_tmp[[j]])
      gr_id <- factor(c(rep(0, nrow(pheno_case_tmp[[1]])), rep(1, nrow(pheno_case_tmp[[j]]))))
      
      # remove pheno with constant values
      p_rm <- apply(new, 2, function(x) length(unique(x)) == 1)
      new <- new[,!p_rm]
      
      new_cov <- rbind(covDat_tmp[[1]], covDat_tmp[[j]])
      res_glm <- matrix(nrow = ncol(new), ncol = 7)
      for(l in 1:ncol(new)){
        type_pheno <- phenoInfo_new[[i]]$type_pheno[phenoInfo_new[[i]]$pheno_id == colnames(new)[l]]
        tmp_dat <- cbind(data.frame(pheno = new[, l], gr_id = gr_id), new_cov)
        res_glm[l,] <- compute_reg_endopheno(mat = tmp_dat, fmla = fmla, type_pheno = type_pheno)
      }
      colnames(res_glm) <- c('beta', 'se_beta', 'z', 'pvalue', 'OR_or_beta', 'CI_low', 'CI_up')
      res_glm <- as.data.frame(res_glm)
      
      phenoInfo_tmp <- phenoInfo_new[[i]][match(colnames(new), phenoInfo_new[[i]]$pheno_id),]
      
      bin_reg[[k]][[j-1]] <- cbind(data.frame(pheno_id = phenoInfo_tmp$pheno_id, type_pheno = phenoInfo_tmp$type_pheno), res_glm)
      bin_reg[[k]][[j-1]]$pval_corr <- p.adjust(bin_reg[[k]][[j-1]]$pvalue, method = 'BH')
      bin_reg[[k]][[j-1]]$comp <- sprintf('gr%i_vs_gr%i',  gr_names[k:length(gr_names)][j], gr_names[k])
      
    }
    
  }
  pairwise_bin_reg[[i]] <- do.call(rbind, do.call(c,bin_reg))
  pairwise_bin_reg[[i]]$pval_corr_overall <-  p.adjust(pairwise_bin_reg[[i]]$pvalue, method = 'BY')
  pairwise_bin_reg[[i]]$cohort <- cohort_name
  
  ################
  ## gri vs all ##
  ################
  
  bin_reg <- vector(mode = 'list', length = length(gr_names))
  for(k in 1:length(gr_names)){
    
    print(paste0('group', gr_names[k], '_vs_all'))
    
    # j vs all
    pheno_case_tmp <- list(phenoDat[cl == gr_names[k],], phenoDat[cl != gr_names[k],])
    covDat_tmp <- list(covDat[cl == gr_names[k],], covDat[cl != gr_names[k],]) 
    
    new <- do.call(rbind, pheno_case_tmp)
    gr_id <- factor(c(rep(1, nrow(pheno_case_tmp[[1]])), rep(0, nrow(pheno_case_tmp[[2]]))))
    
    # remove pheno with constant values
    p_rm <- apply(new, 2, function(x) length(unique(x)) == 1)
    new <- new[,!p_rm]
    
    new_cov <- do.call(rbind, covDat_tmp)
    res_glm <- matrix(nrow = ncol(new), ncol = 7)
    for(l in 1:ncol(new)){
      # print(l)  
      type_pheno <- phenoInfo_new[[i]]$type_pheno[phenoInfo_new[[i]]$pheno_id == colnames(new)[l]]
      tmp_dat <- cbind(data.frame(pheno = new[, l], gr_id = gr_id), new_cov)
      res_glm[l,] <- compute_reg_endopheno(mat = tmp_dat, fmla = fmla, type_pheno = type_pheno)
    }
    
    colnames(res_glm) <- c('beta', 'se_beta', 'z', 'pvalue', 'OR_or_Beta', 'CI_low', 'CI_up')
    res_glm <- as.data.frame(res_glm)
    
    phenoInfo_tmp <- phenoInfo_new[[i]][match(colnames(new), phenoInfo_new[[i]]$pheno_id),]
    bin_reg[[k]] <- cbind(data.frame(pheno_id = phenoInfo_tmp$pheno_id, type_pheno = phenoInfo_tmp$type_pheno), res_glm)
    bin_reg[[k]]$pval_corr <- p.adjust(bin_reg[[k]]$pvalue, method = 'BH')
    bin_reg[[k]]$comp <- sprintf('gr%i_vs_all', gr_names[k])
    
  }
  
  tot_bin_reg[[i]] <- do.call(rbind, bin_reg)
  tot_bin_reg[[i]]$pval_corr_overall <-  p.adjust(tot_bin_reg[[i]]$pvalue, method = 'BY')
  tot_bin_reg[[i]]$cohort <- cohort_name
  
  covDat_new[[i]] <- covDat
}


# save results
output <- list(pairwise = pairwise_bin_reg, 
               tot = tot_bin_reg, 
               cl = clust_new, 
               phenoDat = phenoDat_new, 
               covDat = covDat_new, 
               phenoInfo = phenoInfo_new)
save(output, 
     file = sprintf('%s%s_%s_cluster%s_phenoAssociationGLMall_prediction_model%s.RData', 
                    outFold, type_data, type_input, type_cluster, model_name))

pairwise_save <- do.call(rbind, pairwise_bin_reg)
tot_save <- do.call(rbind, tot_bin_reg)

write.table(pairwise_save, 
            file = sprintf('%s%s_%s_cluster%s_phenoAssociationGLMpairwise_prediction_model%s.txt', 
                           outFold, type_data, type_input, type_cluster, model_name))
write.table(tot_save, 
            file = sprintf('%s%s_%s_cluster%s_phenoAssociationGLM_prediction_model%s.txt', 
                           outFold, type_data, type_input, type_cluster, model_name))


