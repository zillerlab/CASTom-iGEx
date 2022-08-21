# find relevant genes/pathways to cluster structure 
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
parser$add_argument("--split_tot", type = "integer", default = 0, help = "if 0 then inpuntFile load alone, otherwise splitted version")
parser$add_argument("--inputFile", type = "character", nargs = '*', default = 'NA', help = "file to be loaded (predicted tscore or pathScore)")
parser$add_argument("--tissues", type = "character", nargs = '*', help = "tissues name")
parser$add_argument("--type_cluster", type = "character", default = 'All', help = "All, Cases, Controls")
parser$add_argument("--functR", type = "character", help = "functions to be used")
parser$add_argument("--type_data", type = "character", help = "tscore, path_Reactome or path_GO")
parser$add_argument("--type_data_cluster", type = "character", help = "tscore, path_Reactome or path_GO")
parser$add_argument("--type_sim", type = "character", default = 'HK', help = "HK or ED or SNF")
parser$add_argument("--type_input", type = "character", default = 'original', help = "original or zscaled")
parser$add_argument("--pvalresFile", type = "character", nargs = '*',  help = "file with pvalue results")
parser$add_argument("--geneInfoFile", type = "character", nargs = '*', default = NULL, help = "file with info model genes")
parser$add_argument("--min_genes_path", type = "integer", default = 1, help = "minimum number of genes for a pathway, if > 1 recompute corrected pvalues")
parser$add_argument("--pval_id", type = "integer", default = 0, help = "id to be used on pvalue file")
parser$add_argument("--pvalcorr_thr", type = "double", default = 0.05, help = "group sign")
parser$add_argument("--ncores", type = "integer", default = 5, help = "n, cores parallelization per tissue")
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
tissues <- args$tissues
geneInfoFile <- args$geneInfoFile
pvalcorr_thr <- args$pvalcorr_thr
ncores <- args$ncores
outFold <- args$outFold

###################################################################################################################
# tissues <- read.table('OUTPUT_GTEx/Tissue_PGCgwas_red', h=F, stringsAsFactors = F)$V1
# inputFile <- c('OUTPUT_CMC/predict_PGC/200kb/scz_boco_eur/devgeno0.01_testdevgeno0/predictedTscores.txt',
#                sprintf('OUTPUT_GTEx/predict_PGC/%s/200kb/PGC_GWAS_bin1e-2/scz_boco_eur/devgeno0.01_testdevgeno0/predictedTscores.txt', tissues))
# split_tot <- 0
# sampleAnnFile <- 'INPUT_DATA/Covariates/scz_boco_eur.covariateMatrix_old.txt'
# clusterFile <- 'OUTPUT_CMC/predict_PGC/200kb/scz_boco_eur/devgeno0.01_testdevgeno0/update_corrPCs/matchUKBB_filt0.1_tscore_corrPCs_zscaled_predictClusterCases_PGmethod_HKmetric.RData'
# type_cluster <- 'Cases'
# type_data <- 'tscore'
# type_data_cluster <- 'tscore'
# type_sim <- 'HK'
# min_genes_path <- 2
# pval_id <- 1
# pvalresFile <- c('OUTPUT_all/Meta_Analysis_SCZ/DLPC_CMC/pval_Dx_pheno_covCorr.RData',
#                  sprintf('OUTPUT_all/Meta_Analysis_SCZ/%s/pval_Dx_pheno_covCorr.RData', tissues))
# outFold <- 'OUTPUT_GTEx/predict_CAD/Liver/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/CAD_HARD_clustering/'
# functR <- '/home/luciat/castom-igex/Software/model_clustering/clustering_functions.R'
# type_input <- 'zscaled'
# geneInfoFile <- c('OUTPUT_CMC/train_CMC/200kb/resPrior_regEval_allchr.txt',
#                   sprintf('OUTPUT_GTEx/train_GTEx/%s/200kb/PGC_GWAS_bin1e-2/resPrior_regEval_allchr.txt', tissues))
# tissues <- c('DLPC_CMC', tissues)
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

cl <-  cluster_output$cl_best$gr
gr_names <- sort(unique(cl))

if(any(table(cl) < 3)){
  rm_gr <- names(which(table(cl) < 3))
  cluster_output$cl_best <- cluster_output$cl_best[!cluster_output$cl_best$gr %in% rm_gr, ]
  cluster_output$samples_id <- cluster_output$cl_best$id
  cl <-  cluster_output$cl_best$gr
  gr_names <- sort(unique(cl))
}

if(type_data == 'tscore'){
  id_pval <- 8
  id_info <- 2
  id_geno_summ <- 3
}else{
  if(type_data == 'path_Reactome'){
    id_pval <- 13
    id_info <- 1
    id_geno_summ <- 4
  }else{
    if(type_data == 'path_GO'){
      id_pval <- 15
      id_info <- 1
      id_geno_summ <- 6
    }else{
      stop('unknown pathway called')
    }
  }
}

sampleAnn <- read.table(sampleAnnFile, h=T, stringsAsFactors = F, check.names = F)
sampleAnn$Individual_ID <- as.character(sampleAnn$Individual_ID)
sampleAnn <- sampleAnn[match(cluster_output$samples_id, sampleAnn$Individual_ID), ]
name_cov <- setdiff(colnames(sampleAnn),c('Individual_ID', 'genoSample_ID', 'Dx', 'Sex', 'Age', 'Gender', 'Array'))
covDat <- sampleAnn[,!colnames(sampleAnn) %in% c('Individual_ID', 'Dx', 'genoSample_ID')]

identical(sampleAnn$Individual_ID, cluster_output$samples_id)
# scale_data_t <- scoreMat_t <- res_pval_t <- list()

# load score matrix, load pval matrix, rescale for p-value
registerDoParallel(cores=min(ncores, length(tissues)))
print(getDoParWorkers())
#res <- list()

res <- foreach(id_t=1:length(tissues), .combine='comb', 
               .multicombine=TRUE, 
               .init=list(list(), list(), list()))%dopar%{
                 
#                 for(id_t in 1:length(tissues)){
                 print(tissues[id_t])
                 # load pval res
                 res_pval <- get(load(pvalresFile[id_t]))
                 if(type_data == 'tscore'){
                   res_pval <- res_pval$tscore[[pval_id]]
                 }else{
                   if(type_data == 'path_Reactome'){
                     res_pval <- res_pval$pathScore_reactome[[pval_id]]
                   }else{
                     if(type_data == 'path_GO'){
                       res_pval <- res_pval$pathScore_GO[[pval_id]]
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
                 
                 #### load input matrix ####
                 load_output <- load_input_matrix(inputFile = inputFile[id_t], 
                                                  sampleAnn = sampleAnn, 
                                                  res_pval = res_pval, 
                                                  split_tot = split_tot, 
                                                  id_info = id_info)
                 
                 scoreMat_t <- load_output$scoreMat
                 res_pval_t <- load_output$res_pval
                 
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
                 print(identical(colnames(input_data), res_pval_t[, id_info]))
                 print(identical(rownames(input_data), sampleAnn$Individual_ID))
                 
                 scale_data_t <- input_data
                 res_pval_t$tissue <- tissues[id_t]
                 # res[[id_t]] <- list(scoreMat_t, scale_data_t, res_pval_t)
                 list(scoreMat_t, scale_data_t, res_pval_t)
                 
               }

scoreMat_t <- res[[1]]
scale_data_t <- res[[2]]
res_pval_t <- res[[3]]

output <- list(inputData = scoreMat_t, scaleData = scale_data_t, 
               res_pval = res_pval_t, cl = cluster_output$cl_best, 
               tissues = tissues)

##########################
#### check covariates ####
##########################

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

#test_feat_t <- vector(mode = 'list', length = length(tissues))
test_feat_t <- foreach(id_t=1:length(tissues))%dopar%{
  
#for(id_t in 1:length(tissues)){
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
      if(sd(tmp_data$g) > 0){
        tmp <- as.data.frame(wilcox_test(f~g,data = tmp_data, detailed = T))
        test_feat[[i]]$pval[l] <- tmp$p
        test_feat[[i]]$estimates[l] <- tmp$estimate
        test_feat[[i]][l,c('CI_low', 'CI_up')] <- c(tmp$conf.low, tmp$conf.high)
      }
    }
    
    test_feat[[i]]$pval_corr <- p.adjust(test_feat[[i]]$pval, method = 'BH')
    
  }
  
  test_feat <- do.call(rbind, test_feat)
  test_feat$pval_corr_overall <-  p.adjust(test_feat$pval, method = 'BH')
  test_feat$tissue <- tissues[id_t]
  test_feat
  #test_feat_t[[id_t]] <- test_feat
}

output$test_feat <- test_feat_t

# Save
save(output, file = sprintf('%s%sOriginal_corrPCs_%sCluster%s_featAssociation.RData', outFold, type_data, type_data_cluster, type_cluster))


####################################################################################

if(type_data == 'tscore'){
  
  ### combine in loci ###
  # tissue specific
  geneInfo_t <- tissue_spec_loci <- list()
  
  for(id_t in 1:length(tissues)){
    
    tmp <- read.table(geneInfoFile[id_t], h=T, stringsAsFactors = F, sep = '\t')
    geneInfo_t[[id_t]] <- tmp[match(output$res_pval[[id_t]]$ensembl_gene_id, tmp$ensembl_gene_id),]
    geneInfo_t[[id_t]]$tissue <- tissues[id_t]
    
    tmp_feat <- output$test_feat[[id_t]]
    tmp_feat <- tmp_feat[tmp_feat$pval_corr <= pvalcorr_thr,]
    ## remove Y_RNA if any (repeated in different locations of the genome)
    #tmp_feat <- tmp_feat[!tmp_feat$feat %in% c('Y_RNA'),]
    tmp_info <- geneInfo_t[[id_t]][geneInfo_t[[id_t]]$external_gene_name %in% tmp_feat$feat, ] 
    tmp_info$Zstat <- output$res_pval[[id_t]][match(tmp_info$ensembl_gene_id,output$res_pval[[id_t]]$ensembl_gene_id),id_pval-1]
    # divide per chr
    chr_id <- unique(tmp_info$chrom)
    tmp_info_chr <- lapply(chr_id, function(x) tmp_info[tmp_info$chrom == x,])
    
    tmp_loci <- list()
    for(j in 1:length(chr_id)){
      # print(j)
      tmp_loci[[j]] <- merge_locus_pos(tmp_info_chr[[j]], tissue = tissues[id_t])
    }  
    
    tmp_loci <- do.call(rbind, tmp_loci)
    
    tmp_loci$comp_sign <- NA
    tmp_loci$best_WMW_gene <- NA
    tmp_loci$best_WMW_est <- NA
    tmp_loci$best_WMW_pvalue <- NA
    
    # for each loci, find groups that are significantly different
    for(l in 1:nrow(tmp_loci)){
      genes <- strsplit(tmp_loci$gene[l], split = ',')[[1]]
      tmp_gr <- sapply(unique(tmp_feat$comp[tmp_feat$feat %in% genes]), 
                       function(x) strsplit(x, split = '_vs_all')[[1]][1])
      tmp_gr <- sort(tmp_gr)
      tmp_loci$comp_sign[l] <-  paste0(tmp_gr, collapse = ',') 
      tmp_WMW <- lapply(names(tmp_gr), function(x) 
        tmp_feat[tmp_feat$comp == x & tmp_feat$feat %in% genes,])
      # find common genes in all the groups
      genes_update <- names(which(table(unlist(lapply(tmp_WMW, function(x) x$feat))) == length(tmp_gr)))
      if(length(genes_update)>0){
        tmp_WMW <- lapply(tmp_WMW, function(x) x[x$feat %in% genes_update, ])
        tmp_WMW_all <- do.call(rbind, tmp_WMW)
        tmp_loci$best_WMW_gene[l] <- tmp_WMW_all$feat[which.max(abs(tmp_WMW_all$estimates))]
        tmp_loci$best_WMW_est[l] <- paste0(round(sapply(tmp_WMW, function(x) x$estimates[x$feat == tmp_loci$best_WMW_gene[l]]),digits = 5), collapse = ',') 
        tmp_loci$best_WMW_pvalue[l] <- paste0(sapply(tmp_WMW, function(x) x$pval[x$feat == tmp_loci$best_WMW_gene[l]]), collapse = ',')
      }else{
        tmp_WMW <- lapply(names(tmp_gr), function(x) tmp_feat[tmp_feat$comp == x & tmp_feat$feat %in% genes,])
        tmp_loci$best_WMW_gene[l] <- paste0(sapply(tmp_WMW, function(x) x$feat[which.max(abs(x$estimates))]), collapse = ',')
        tmp_loci$best_WMW_est[l] <-  paste0(round(sapply(tmp_WMW, function(x) x$estimates[which.max(abs(x$estimates))]),digits = 5), collapse = ',') 
        tmp_loci$best_WMW_pvalue[l] <-  paste0(sapply(tmp_WMW, function(x) x$pval[which.min(x$pval)]), collapse = ',')
      }
    }
    
    tissue_spec_loci[[id_t]] <- tmp_loci
    
  }
  
  tissue_spec_loci <- do.call(rbind, tissue_spec_loci)
  # save
  write.table(x = tissue_spec_loci, file = sprintf('%s%s_corrPCs_%s_cluster%s_summary_geneLoci_tissueSpec.txt',outFold, type_data, type_input, type_cluster), 
              col.names = T, row.names = F, sep = '\t', quote = F)
  
  # all tissues
  test_feat_tot <- do.call(rbind, output$test_feat)
  test_feat_tot$new_id <- paste(test_feat_tot$feat, test_feat_tot$tissue, sep = '_')
  
  geneInfo_tot <- do.call(rbind, geneInfo_t)
  geneInfo_tot$new_id <- paste(geneInfo_tot$external_gene_name, geneInfo_tot$tissue, sep = '_')
  
  res_pval_tot <- do.call(rbind, output$res_pval)
  res_pval_tot$new_id <- paste(res_pval_tot$external_gene_name, res_pval_tot$tissue, sep = '_')
  
  # all tissues
  tmp_feat <- test_feat_tot[test_feat_tot$pval_corr <= pvalcorr_thr,]
  # # remove Y_RNA if any (repeated in different locations of the genome)
  # tmp_feat <- tmp_feat[!tmp_feat$feat %in% c('Y_RNA'),]

  tmp_info <- geneInfo_tot[geneInfo_tot$new_id %in% tmp_feat$new_id, ] 
  tmp_info$Zstat <- res_pval_tot[match(tmp_info$new_id,res_pval_tot$new_id),id_pval-1]
  
  # divide per chr
  chr_id <- unique(tmp_info$chrom)
  tmp_info_chr <- lapply(chr_id, function(x) tmp_info[tmp_info$chrom == x,])
  
  tmp_loci <- list()
  for(j in 1:length(chr_id)){
    # print(j)
    tmp_loci[[j]] <- merge_locus_pos(tmp_info_chr[[j]], tissue = 'combined')
  }  
  
  alltissues_loci <- do.call(rbind, tmp_loci)
  alltissues_loci$comp_sign <- NA
  # for each loci, find groups that are significantly different
  for(i in 1:nrow(alltissues_loci)){
    genes <- strsplit(alltissues_loci$gene[i], split = ',')[[1]]
    tmp_gr <- sapply(unique(tmp_feat$comp[tmp_feat$feat %in% genes]), function(x) strsplit(x, split = '_vs_all')[[1]][1])
    alltissues_loci$comp_sign[i] <-  paste0(tmp_gr, collapse = ',') 
  }
  alltissues_loci <- alltissues_loci[order(factor(alltissues_loci$chrom, levels = paste0('chr', 1:22)), alltissues_loci$start), ]
  
  # save
  write.table(x = alltissues_loci, file = sprintf('%s%s_corrPCs_%s_cluster%s_summary_geneLoci_allTissues.txt',outFold, type_data, type_input, type_cluster), 
              col.names = T, row.names = F, sep = '\t', quote = F)
  
}

