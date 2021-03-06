# Mendelian randomization (version 2, cosider all phenotypes regardelss correlation)

options(stringsAsFactors=F)
options(max.print=1000)
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggsci))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(qvalue))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(ggsignif))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(MendelianRandomization))
options(bitmapType = 'cairo', device = 'png')

parser <- ArgumentParser(description="Mendelian Randomization for enriched related phenotypes")
parser$add_argument("--phenoFold", type = "character", help = "Folder with results phenotype annotation UKBB")
parser$add_argument("--inputFold_rel", type = "character", nargs = '*',help = "Folder with results")
parser$add_argument("--pval_FDR_rel", type = "double", default = 0.05, help = "pval threshold to filter the genes and pathways (after BH correction) in related phenotypes")
parser$add_argument("--pval_corr", type = "double", default = 0.05, help = "pval threshold to filter phenotypes based on fisher test enrichment")
parser$add_argument("--tissue_name", type = "character", nargs = '*', help = "tissue considered")
parser$add_argument("--pheno_name_comp", type = "character", help = "name phenotype of interest")
parser$add_argument("--pheno_file", type = "character", nargs = '*', help = "")
parser$add_argument("--type_data", type = "character", help = "")
parser$add_argument("--enrich_file", type = "character", help = "")
parser$add_argument("--MR_pheno_file", type = "character", default = NULL, help = "")
parser$add_argument("--cor_file", type = "character", help = "")
parser$add_argument("--outFold", type="character", help = "Output file [basename only]")

args <- parser$parse_args()
phenoFold <- args$phenoFold
MR_pheno_file <- args$MR_pheno_file
tissue_name <- args$tissue_name
inputFold_rel <- args$inputFold_rel
pval_FDR_rel <- args$pval_FDR_rel
pval_corr <- args$pval_corr
pheno_name_comp <- args$pheno_name_comp
pheno_file <- args$pheno_file
type_data <- args$type_data
enrich_file <- args$enrich_file
cor_file <- args$cor_file
outFold <- args$outFold

########################################################################################################################
# phenoFold <- 'INPUT_DATA_GTEx/CAD/Covariates/UKBB/'
# tissue_name <- 'Whole_Blood'
# pheno_name_comp <- 'CADCardioG'
# MR_pheno_file <- 
# inputFold_rel <- paste0('OUTPUT_GTEx/predict_CAD/Whole_Blood/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/')
# pheno_file <- c('OUTPUT_GTEx/predict_CAD/AllTissues/200kb/CAD_GWAS_bin5e-2/Meta_Analysis_CAD/tscore_pval_Dx_covCorr.txt')
# # pheno_file <- c('Meta_Analysis_SCZ/OUTPUT_all/path_Reactome_pval_SCZ_covCorr_filt.txt',
# #                 'Meta_Analysis_SCZ/OUTPUT_all/path_GO_pval_SCZ_covCorr_filt.txt')
# cor_file <-  paste0('OUTPUT_GTEx/predict_CAD/Whole_Blood/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/correlation_estimate_tscore.RData')
# type_data <- 'tscore'
# enrich_file <- 'OUTPUT_GTEx/predict_CAD/Whole_Blood/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/enrichment_CADHARD_res/CADCardioG_perc0.3_dist200000b_correlation_enrich_CADCardioG_relatedPheno.RData'
# outFold <- 'OUTPUT_GTEx/predict_CAD/Whole_Blood/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/enrichment_CADHARD_res/dist250000b_'
# pval_FDR_rel <- 0.05
# pval_corr <- 0.05
##########################################################################################################################

load_data <- function(data_file, pval_FDR, tissue_name, pval_id){
  
  tscore <- list()
  pathR <- list()
  pathGO <- list()
  
  for(i in 1:length(tissue_name)){
    
    # print(tissue_name[i])
    res_pval <- get(load(data_file[i]))
    
    tscore[[i]] <- list()
    pathR[[i]] <- list()
    pathGO[[i]] <- list()
    
    for(j in 1:length(pval_id)){
      
      tscore[[i]][[j]] <- res_pval$tscore[[pval_id[j]]]
      tscore[[i]][[j]] <- res_pval$tscore[[pval_id[j]]]
      pathR[[i]][[j]] <- res_pval$pathScore_reactome[[pval_id[j]]]
      pathGO[[i]][[j]] <- res_pval$pathScore_GO[[pval_id[j]]]
      tscore[[i]][[j]]$tissue <- pathR[[i]][[j]]$tissue <- pathGO[[i]][[j]]$tissue <- tissue_name[i]
      # filter out pathways with only 1 gene and recompute adjusted pvalue
      pathR[[i]][[j]] <- pathR[[i]][[j]][pathR[[i]][[j]]$ngenes_tscore>1, ]
      pathR[[i]][[j]][, 14] <- qvalue(pathR[[i]][[j]][,13])$qvalue
      pathR[[i]][[j]][, 15] <- p.adjust(p = pathR[[i]][[j]][,13], method = 'BH')
      
      pathGO[[i]][[j]] <- pathGO[[i]][[j]][pathGO[[i]][[j]]$ngenes_tscore>1, ]
      pathGO[[i]][[j]][, 16] <- qvalue(pathGO[[i]][[j]][,15])$qvalue
      pathGO[[i]][[j]][, 17] <- p.adjust(p = pathGO[[i]][[j]][,15], method = 'BH')
      
    }
  }
  
  pheno_tscore <- pheno_pathR <- pheno_pathGO <- list()
  tscore_red <- pathR_red <- pathGO_red <- list()
  for(j in 1:length(pval_id)){
    
    pheno_tscore[[j]] <- do.call(rbind, lapply(tscore, function(x) x[[j]]))
    pheno_pathR[[j]] <- do.call(rbind, lapply(pathR, function(x) x[[j]]))
    pheno_pathGO[[j]] <- do.call(rbind, lapply(pathGO, function(x) x[[j]]))
    pheno_tscore[[j]]$new_id <- paste0(pheno_tscore[[j]][,1],'_tissue_',pheno_tscore[[j]]$tissue)
    pheno_pathR[[j]]$new_id <- paste0(pheno_pathR[[j]][,1],'_tissue_',pheno_pathR[[j]]$tissue)
    pheno_pathGO[[j]]$new_id <- paste0(pheno_pathGO[[j]][,2],'_tissue_',pheno_pathGO[[j]]$tissue)
    
    tscore_red[[j]] <- pheno_tscore[[j]][pheno_tscore[[j]][, 10] <= pval_FDR,]
    pathR_red[[j]] <- pheno_pathR[[j]][pheno_pathR[[j]][, 15] <= pval_FDR,]
    pathGO_red[[j]] <- pheno_pathGO[[j]][pheno_pathGO[[j]][, 17] <= pval_FDR,]
    
  }
  return(list(tscore = pheno_tscore, tscore_red = tscore_red, pathR = pheno_pathR, pathR_red = pathR_red, 
              pathGO = pheno_pathGO, pathGO_red = pathGO_red))
  
}

if(type_data == 'tscore'){
  id_name <- 2
  id_keep <- 2
  pval_id <- 8
  beta_id <- 5
  se_beta_id <- 6
}

if(type_data == 'tot_path'){
  if(grepl('CAD_', pheno_name_comp)){id_name <- 22}
  if(grepl('SCZ', pheno_name_comp) | grepl('CADCardioG', pheno_name_comp)){id_name <- 25}
  pval_id <- 13
  beta_id <- 10
  se_beta_id <- 11
}

feat_res <- lapply(pheno_file, function(x) read.delim(x, h=T, stringsAsFactors = F, sep = '\t'))
feat_res <- lapply(feat_res, function(x) x[x$tissue %in% tissue_name, ])
if(type_data == 'tot_path'){
  feat_res[[1]]$path_id <- feat_res[[1]]$path
  feat_res[[1]]$type <- 'Reactome'
  feat_res[[2]]$type <- 'GO'
  common_h <- intersect(colnames(feat_res[[1]]), colnames(feat_res[[2]]))
  feat_res <- lapply(feat_res, function(x) x[, match(common_h, colnames(x))])
  feat_res <- do.call(rbind, feat_res)
  feat_res$new_id <- paste0(feat_res$path_id, '_type_', feat_res$type)
}
if(type_data == 'tscore'){
  feat_res <- do.call(rbind, feat_res)
  feat_res$new_id <- feat_res[,1]
}

# load OR/correlation results
enrich_res <- get(load(enrich_file))
pheno_info <- enrich_res$pheno
if(type_data == 'tscore'){
  enrich_feat <- enrich_res$tscore
  keep_genes <- enrich_res$gene_ann
  # use same genes used to compute correlation 
  feat_res <- feat_res[feat_res$ensembl_gene_id %in% keep_genes$ensembl_gene_id, ]
}
if(type_data == 'tot_path'){
  enrich_feat <- enrich_res$pathScore
  keep_path <- enrich_res$path_ann
  # use same pathways used to compute correlation 
  feat_res <- feat_res[paste0(feat_res$path, '_type_', feat_res$type) %in% paste0(keep_path$path, '_type_', keep_path$type), ]
}

# load correlation data
cor_feat <- get(load(cor_file))
cor_feat <- cor_feat$cor
common_f <- intersect(feat_res[, id_name], colnames(cor_feat))
feat_res <- feat_res[match(common_f, feat_res[,id_name]), ]
cor_feat <- cor_feat[match(common_f, rownames(cor_feat)), match(common_f, colnames(cor_feat))]

# keep only phenotype that are correlated 
enrich_feat <- enrich_feat[!is.na(enrich_feat$cor_pval_BHcorr),]
enrich_feat <- enrich_feat[enrich_feat$cor_pval <= pval_corr, ]
# enrich_feat <- enrich_feat[!is.na(enrich_feat$fisher_pval_BHcorr),]
# enrich_feat <- enrich_feat[enrich_feat$fisher_pval_BHcorr <= pval_FDR_fisher,]
# enrich_feat <- enrich_feat[enrich_feat$k > 2, ]
pheno_type <- unique(enrich_feat$pheno_type)
if(!is.null(MR_pheno_file)){
  MR_pheno_subset <- read.table(MR_pheno_file, h = F, stringsAsFactors = F, sep = '\t')$V1
  pheno_type <- pheno_type[pheno_type %in% MR_pheno_subset]
}
print(pheno_type)

res_MR_Egger <- list()
res_MR_IVW <- list()

for(i in 1:length(pheno_type)){
  
  print(pheno_type[i])
  
  tmp_pheno_red <- enrich_feat[enrich_feat$pheno_type %in% pheno_type[i], ]
  res_MR_Egger[[i]] <- data.frame(nfeat = rep(NA, nrow(tmp_pheno_red)), MREgg_est = rep(NA, nrow(tmp_pheno_red)), MREgg_est_se = rep(NA, nrow(tmp_pheno_red)), 
                                  MREgg_est_CIl = rep(NA, nrow(tmp_pheno_red)), MREgg_est_CIu = rep(NA, nrow(tmp_pheno_red)), 
                                  MREgg_est_pval = rep(NA, nrow(tmp_pheno_red)), MREgg_int = rep(NA, nrow(tmp_pheno_red)), MREgg_int_se = rep(NA, nrow(tmp_pheno_red)), 
                                  MREgg_int_CIl = rep(NA, nrow(tmp_pheno_red)), MREgg_int_CIu = rep(NA, nrow(tmp_pheno_red)), 
                                  MREgg_int_pval = rep(NA, nrow(tmp_pheno_red)), MREgg_heter_stat  = rep(NA, nrow(tmp_pheno_red)),
                                  MREgg_heter_pval = rep(NA, nrow(tmp_pheno_red)), MREgg_Isq  = rep(NA, nrow(tmp_pheno_red)))
  
  res_MR_IVW[[i]] <-  data.frame(nfeat = rep(NA, nrow(tmp_pheno_red)), MRIVW_est = rep(NA, nrow(tmp_pheno_red)), MRIVW_est_se = rep(NA, nrow(tmp_pheno_red)), 
                                 MRIVW_est_CIl = rep(NA, nrow(tmp_pheno_red)), MRIVW_est_CIu = rep(NA, nrow(tmp_pheno_red)), 
                                 MRIVW_est_pval = rep(NA, nrow(tmp_pheno_red)),
                                 MRIVW_heter_stat  = rep(NA, nrow(tmp_pheno_red)),
                                 MRIVW_heter_pval = rep(NA, nrow(tmp_pheno_red)))
  
  if(grepl('CAD', pheno_name_comp) & pheno_type[i] %in% c('Blood_biochemistry', 'Blood_count')){
    file_toload <- sprintf('%s/pval_%s_withMed_pheno_covCorr.RData', inputFold_rel, pheno_type[i])
  }else{
    file_toload <- sprintf('%s/pval_%s_pheno_covCorr.RData', inputFold_rel, pheno_type[i])
  }
  
  tmp <- get(load(file_toload))
  rm(final)
  tmp_pheno <- tmp$pheno
  
  if(type_data == 'tot_path'){
    new <- list()
    for(l in 1:nrow(tmp_pheno)){
      
      tmp$pathScore_reactome[[l]] <- tmp$pathScore_reactome[[l]][tmp$pathScore_reactome[[l]]$ngenes_tscore > 1, ]
      tmp$pathScore_reactome[[l]][, 15] <-  p.adjust(tmp$pathScore_reactome[[l]][, 13], method = 'BH')
      tmp$pathScore_GO[[l]] <- tmp$pathScore_GO[[l]][tmp$pathScore_GO[[l]]$ngenes_tscore > 1, ]
      tmp$pathScore_GO[[l]][, 17] <-  p.adjust(tmp$pathScore_GO[[l]][, 15], method = 'BH')
      
      new[[l]] <- list(tmp$pathScore_reactome[[l]], tmp$pathScore_GO[[l]])
      new[[l]][[1]]$path_id <- new[[l]][[1]]$path
      new[[l]][[1]]$type <- 'Reactome'
      new[[l]][[2]]$type <- 'GO'
      common_h <- intersect(colnames(new[[l]][[1]]), colnames(new[[l]][[2]]))
      new[[l]] <- lapply(new[[l]], function(x) x[, match(common_h, colnames(x))])
      new[[l]] <- do.call(rbind, new[[l]])
      new[[l]]$new_id <- paste0(new[[l]]$path_id, '_type_', new[[l]]$type)
    }
    tmp <- new
  }
  
  #### for pathways, remove 1 pathways with only 1 gene and update pvalues (otherwise problem in matching!!!) ###
  
  if(type_data == 'tscore'){
    tmp <- tmp[[id_keep]] 
    for(j in 1:length(tmp)){
      # tmp[[j]]$new_id <- paste0(tmp[[j]][,1], '_tissue_', tissue_name)
      tmp[[j]]$new_id <- tmp[[j]][,1]
    }
  }
  
  id_keep_rel <- match(tmp_pheno_red$pheno,tmp_pheno$pheno_id)
  tmp <- tmp[id_keep_rel]
  
  common_f <- lapply(tmp, function(x) intersect(x$new_id, feat_res$new_id))
  tmp <- mapply(function(x, y) x[match(y,x$new_id),], x = tmp, y = common_f,SIMPLIFY = F)
  
  # sign_common_rel <- lapply(tmp, function(x) x[x$new_id %in% feat_res$new_id[feat_res[,pval_id+2] <= pval_FDR_pheno] & x[, pval_id +2] <= pval_FDR_rel,])
  # sign_common_pheno <- lapply(tmp, function(x) feat_res[feat_res$new_id %in% x$new_id[x[,pval_id+2] <= pval_FDR_rel] & feat_res[, pval_id +2] <= pval_FDR_pheno,])
  sign_common_rel <- lapply(tmp, function(x) x[x[, pval_id +2] <= pval_FDR_rel,])
  sign_common_pheno <- lapply(sign_common_rel, function(x) feat_res[feat_res$new_id %in% x$new_id,])
  
  for(j in 1:length(tmp)){
    
    common_f <- intersect(sign_common_rel[[j]]$new_id, sign_common_pheno[[j]]$new_id)
    
    if(length(common_f)>2){
      
      sign_common_rel[[j]] <- sign_common_rel[[j]][match(common_f, sign_common_rel[[j]]$new_id),]
      sign_common_pheno[[j]] <- sign_common_pheno[[j]][match(common_f, sign_common_pheno[[j]]$new_id),]
      tmp_corr <- cor_feat[match(sign_common_pheno[[j]][, id_name], rownames(cor_feat)), 
                           match(sign_common_pheno[[j]][, id_name], colnames(cor_feat))]
      
      MRInputObject <- mr_input(bx = sign_common_rel[[j]][, beta_id], bxse = sign_common_rel[[j]][, se_beta_id],
                                by = sign_common_pheno[[j]][, beta_id], byse = sign_common_pheno[[j]][, se_beta_id], 
                                correlation = tmp_corr, 
                                snps = colnames(tmp_corr))
      
      MREggerObject <- tryCatch(mr_egger(MRInputObject, correl = T), error = function(x) NULL)
      # MREggerObject <- mr_egger(MRInputObject, correl = T)
      if(!is.null(MREggerObject)){
        res_MR_Egger[[i]][j, ] <- c(length(common_f), MREggerObject@Estimate, MREggerObject@StdError.Est, MREggerObject@CILower.Est, 
                                    MREggerObject@CIUpper.Est, MREggerObject@Pvalue.Est, MREggerObject@Intercept, 
                                    MREggerObject@StdError.Int, MREggerObject@CILower.Int, MREggerObject@CIUpper.Int,
                                    MREggerObject@Pvalue.Int, MREggerObject@Heter.Stat, MREggerObject@I.sq)
      }else{
        print(paste('computational error for', tmp_pheno_red$names_field[j]))
      }
      
      MRIVWObject <- tryCatch(mr_ivw(MRInputObject, correl = T, model = 'random'), error = function(x) NULL)
      # MREggerObject <- mr_egger(MRInputObject, correl = T)
      if(!is.null(MRIVWObject)){
        res_MR_IVW[[i]][j, ] <- c(length(common_f), MRIVWObject@Estimate, MRIVWObject@StdError, MRIVWObject@CILower, 
                                  MRIVWObject@CIUpper, MRIVWObject@Pvalue, MRIVWObject@Heter.Stat)
      }else{
        print(paste('computational error for', tmp_pheno_red$names_field[j]))
      }
      
    }
    
  }
  
  res_MR_Egger[[i]] <- cbind(res_MR_Egger[[i]], tmp_pheno_red[, c('names_field', 'pheno', 'pheno_type')])
  res_MR_IVW[[i]] <- cbind(res_MR_IVW[[i]], tmp_pheno_red[, c('names_field', 'pheno', 'pheno_type')])
  
}

res_MR_Egger <- do.call(rbind, res_MR_Egger)
res_MR_IVW <- do.call(rbind, res_MR_IVW)
res_MR_Egger <- res_MR_Egger[!is.na(res_MR_Egger$nfeat), ]
res_MR_IVW <- res_MR_IVW[!is.na(res_MR_IVW$nfeat), ]

res_MR_Egger$MREgg_est_pval_FDR <- p.adjust(res_MR_Egger$MREgg_est_pval, method = 'BH')
res_MR_Egger$MREgg_int_pval_FDR <- p.adjust(res_MR_Egger$MREgg_int_pval, method = 'BH')
res_MR_IVW$MRIVW_est_pval_FDR <- p.adjust(res_MR_IVW$MRIVW_est_pval, method = 'BH')


# save results
write.table(x = res_MR_Egger, file = sprintf('%sMendelian_randomization_Egger_%s_pvalFDRrel%s.txt', outFold, type_data, 
                                             as.character(pval_FDR_rel)), col.names = T, row.names = F, quote = F, sep = '\t')

write.table(x = res_MR_IVW, file = sprintf('%sMendelian_randomization_IVW_%s_pvalFDRrel%s.txt', outFold, type_data, 
                                           as.character(pval_FDR_rel)), col.names = T, row.names = F, quote = F, sep = '\t')


