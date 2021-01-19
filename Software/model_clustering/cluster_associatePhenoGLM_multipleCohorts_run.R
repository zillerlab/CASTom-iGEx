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

parser <- ArgumentParser(description="Associate cluster to endophenotype, multiple cohorts analysed together (lmm)")
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

source(functR)
phenoDat <- list()
sampleAnn <- list()

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
    intersect_samples <- intersect(sampleAnn[[i]]$Individual_ID, cluster_output$sampleAnn$Individual_ID[cluster_output$sampleAnn$cohort == name_cohorts[i]])
    sampleAnn[[i]] <- sampleAnn[[i]][match(intersect_samples, sampleAnn[[i]]$Individual_ID),]
    name_cl <- ifelse('cl_best' %in% names(cluster_output), 'cl_best', 'cl_new')
    cluster_output <- cluster_output[which(names(cluster_output) == name_cl)][[1]]
    
  }else{
    tmp <- get(load(clusterFile[i]))
    name_cl <- ifelse('cl_best' %in% names(tmp), 'cl_best', 'cl_new')
    cluster_output[[i]] <- tmp[which(names(tmp) == name_cl)][[1]]
    sampleAnn[[i]] <- sampleAnn[[i]][match(cluster_output[[i]]$id,sampleAnn[[i]]$Individual_ID),]
  }
  
  phenoDat[[i]] <- phenoDat[[i]][match(sampleAnn[[i]]$Individual_ID, phenoDat[[i]]$Individual_ID),]
  phenoDat[[i]] <- phenoDat[[i]][, -1]
  
}

# combine results
common_cov <- names(which(table(unlist(lapply(sampleAnn, colnames))) == length(name_cohorts)))
sampleAnn_tot <- do.call(rbind, lapply(sampleAnn, function(x) x[, match(common_cov,colnames(x))]))

common_pheno <- names(which(table(unlist(lapply(phenoDat, colnames))) == length(name_cohorts)))
phenoDat_tot <- do.call(rbind, lapply(phenoDat, function(x) x[, match(common_pheno,colnames(x))]))
# remove columns with too many NAs
phenoDat_tot <- phenoDat_tot[,colSums(!is.na(phenoDat_tot))>=100]
# remove binary with too few T
id_bin <- rep(F, ncol(phenoDat_tot))
for(i in 1:ncol(phenoDat_tot)){
  id_bin[i] <- is.integer(phenoDat_tot[, i])
}
rm_col <- colnames(phenoDat_tot[, id_bin & colSums(phenoDat_tot != 0 & !is.na(phenoDat_tot))<50])
phenoDat_tot <- phenoDat_tot[,!colnames(phenoDat_tot) %in% rm_col]

phenoInfo <- read.delim(phenoDescFile, h=T, stringsAsFactors = F, sep = '\t')
# get common pheno id (phenoDescFile can be used to filter)
common_pheno <- intersect(phenoInfo$pheno_id, colnames(phenoDat_tot))
phenoInfo <- phenoInfo[match(common_pheno, phenoInfo$pheno_id),]
phenoDat_tot <- phenoDat_tot[,match(common_pheno,colnames(phenoDat_tot))]

if(risk_score){
  phenoInfo$transformed_type <- 'CONTINUOUS'
  phenoDat_tot <- scale(phenoDat_tot)
  attr(phenoDat_tot, "scaled:center") <- NULL
  attr(phenoDat_tot, "scaled:scale") <- NULL
  phenoDat_tot <- as.data.frame(phenoDat_tot)
}

if(any(phenoInfo$Path %in% 'Online follow-up > Cognitive function online > Fluid intelligence')){
  phenoInfo$Field[phenoInfo$Path == 'Online follow-up > Cognitive function online > Fluid intelligence'] <- paste(phenoInfo$Field[phenoInfo$Path == 'Online follow-up > Cognitive function online > Fluid intelligence'], '(Online)')
}

if(length(clusterFile) > 1){
  cl <- do.call(rbind, cluster_output)
}

cl$cohort <- sampleAnn_tot$cohort

output <- list(phenoDat = phenoDat_tot, phenoInfo = phenoInfo, cl = cl)

#######################################
#### binary regression (gi vs gj) #####
#######################################

gr_names <- sort(unique(cl$gr))
P <- length(gr_names)
covDat <- sampleAnn_tot[, !colnames(sampleAnn_tot) %in% c('Individual_ID', 'genoSample_ID', 'Dx')]
fmla  <- as.formula(paste('pheno~gr_id+', paste0(colnames(covDat)[!colnames(covDat) %in% 'cohort'], collapse = '+'), '+ (1 | cohort)'))
output$covDat = covDat

bin_reg <- vector(mode = 'list', length = length(gr_names)-1)
for(i in 1:(length(gr_names)-1)){
  
  print(paste0('group', gr_names[i], '_vs_groupj'))
  
  # j vs all
  pheno_case_tmp <- lapply(gr_names[i:length(gr_names)], function(x) phenoDat_tot[cl$gr == x,])
  covDat_tmp <- lapply(gr_names[i:length(gr_names)], function(x) covDat[cl$gr == x,])
  bin_reg[[i]] <-  vector(mode = 'list', length = length(pheno_case_tmp)-1)
  
  for(j in 2:length(pheno_case_tmp)){
    
    print(j)
    
    new <- rbind(pheno_case_tmp[[1]], pheno_case_tmp[[j]])
    colnames(new) <- paste0('p', colnames(new))
    gr_id <- factor(c(rep(0, nrow(pheno_case_tmp[[1]])), rep(1, nrow(pheno_case_tmp[[j]]))))
    
    # remove pheno with constant values
    p_rm <- apply(new, 2, function(x) length(unique(x)) == 1)
    new <- new[,!p_rm]
    # remove phenotype with few true overall (find binomial)
    id_b <- intersect(colnames(new), paste0('p',phenoInfo$pheno_id[!phenoInfo$transformed_type %in% c('CONTINUOUS', 'CAT_ORD')]))
    id_rm <- names(which(colSums(new[, id_b, drop = F], na.rm = T) < 50))
    if(length(id_rm)>0){
      new <- new[, !colnames(new) %in% id_rm]
    }
    
    new_cov <- rbind(covDat_tmp[[1]], covDat_tmp[[j]])
    res_glm <- data.frame(beta = c(), se_beta = c(), z = c(), CI_low = c(), CI_up = c())
    for(l in 1:ncol(new)){
      print(l)
      type_pheno <- phenoInfo$transformed_type[paste0('p',phenoInfo$pheno_id) == colnames(new)[l]]
      tmp_dat <- cbind(data.frame(pheno = new[, l], gr_id = gr_id), new_cov)
      tmp_dat$cohort <- factor(tmp_dat$cohort)
      res_glm <- rbind(res_glm, compute_reg_endopheno_lmm(mat = tmp_dat, fmla = fmla, type_pheno = type_pheno))
    }
   
    res_glm$type_pheno <- phenoInfo$transformed_type[match(colnames(new), paste0('p',phenoInfo$pheno_id))]
    
    phenoInfo_tmp <- phenoInfo[match(colnames(new), paste0('p',phenoInfo$pheno_id)),]
    
    bin_reg[[i]][[j-1]] <- cbind(data.frame(pheno_id = phenoInfo_tmp$pheno_id, Field = phenoInfo_tmp$Field, meaning = phenoInfo_tmp$Coding_meaning), res_glm)
    # bin_reg[[i]][[j-1]]$pval_corr <- p.adjust(bin_reg[[i]][[j-1]]$pvalue, method = 'BH')
    bin_reg[[i]][[j-1]]$comp <- sprintf('gr%i_vs_gr%i', gr_names[i:length(gr_names)][j], gr_names[i])
    
  }
  
  
}

tot_bin_reg <- do.call(rbind, do.call(c,bin_reg))
# tot_bin_reg$pval_corr_overall <-  p.adjust(tot_bin_reg$pvalue, method = 'BY')

# save
write.table(x = tot_bin_reg, sprintf('%s%s_%s_cluster%s_PGmethod_%smetric_phenoAssociation_GLMpairwise.txt', outFold, type_data, type_input, type_cluster, type_sim), col.names = T, row.names = F, sep = '\t', quote = F)

output$bin_reg = tot_bin_reg

save(output, file = sprintf('%s%s_%s_cluster%s_PGmethod_%smetric_phenoAssociation_GLMpairwise.RData', outFold, type_data, type_input, type_cluster, type_sim))

###############################################################################
if(!risk_score){
  
  test_pheno <- tot_bin_reg
  test_pheno$sign <- 'no'
  test_pheno$sign[test_pheno$pval_corr_overall <= 0.05] <- 'yes'
  test_pheno$logpval <- -log10(test_pheno$pvalue)
  test_pheno$sign <- factor(test_pheno$sign, levels = c('no', 'yes'))
  test_pheno$new_Field <- test_pheno$Field
  test_pheno$new_Field[test_pheno$Field == 'Blood clot, DVT, bronchitis, emphysema, asthma, rhinitis, eczema, allergy diagnosed by doctor'] <- 'Blood clot, etc. diagnosed by doctor'
  test_pheno$new_id <- test_pheno$new_Field
  test_pheno$new_id[!is.na(test_pheno$meaning)] <- paste(test_pheno$meaning[!is.na(test_pheno$meaning)], test_pheno$new_Field[!is.na(test_pheno$meaning)], sep = '\n')
  
  test_pheno <- test_pheno[order(test_pheno$pvalue),]
  if(length(which(test_pheno$sign == 'yes'))>=20){
    test_pheno <- test_pheno[1:length(which(test_pheno$sign == 'yes'))+1,]
    height_pl <- 5 + (nrow(test_pheno)-20)*0.2
  }else{
    test_pheno <- test_pheno[1:20,] 
    height_pl = 5
  }
  
  test_pheno$comp <- factor(test_pheno$comp)
  test_pheno$new_id <- factor(test_pheno$new_id, levels = unique(test_pheno$new_id))
  
  pl <-  ggplot(test_pheno, aes(x = new_id, y = logpval, fill = sign))+
    geom_bar(stat = 'identity', position = position_dodge(), width = 0.5) + theme_bw()+ 
    ylab('-log10(pvalue)')+xlab('')+geom_hline(yintercept = -log10(0.05))+
    facet_wrap(.~comp, scales = 'free_y', ncol = 1, strip.position="right")+
    scale_fill_manual(values=c("#999999", "#E69F00"))+
    theme(legend.position = 'none', plot.title = element_text(size=9), axis.text.y = element_text(size = 7), strip.text = element_text(size=5))+
    ggtitle(paste(type_data, type_input, 'cluster', type_cluster))+
    coord_flip()
  ggsave(filename = sprintf('%s%s_%s_cluster%s_PGmethod_%smetric_phenoAssociation_GLMpairwise.png', outFold, type_data, type_input, type_cluster, type_sim), width = 4.5, height = height_pl, plot = pl, device = 'png')
  ggsave(filename = sprintf('%s%s_%s_cluster%s_PGmethod_%smetric_phenoAssociation_GLMpairwise.pdf', outFold, type_data, type_input, type_cluster, type_sim), width = 4.5, height = height_pl, plot = pl, device = 'pdf')
}
###############################################################################################
# same but not pairwise (1 gropu against all the others)

output <- list(phenoDat = phenoDat_tot, phenoInfo = phenoInfo, cl = cl)

########################################
#### binary regression (gi vs all) #####
########################################

covDat <- sampleAnn_tot[, !colnames(sampleAnn_tot) %in% c('Individual_ID', 'genoSample_ID', 'Dx')]
fmla  <- as.formula(paste('pheno~gr_id+', paste0(colnames(covDat)[!colnames(covDat) %in% 'cohort'], collapse = '+'), '+ (1 | cohort)'))
output$covDat = covDat

bin_reg <- vector(mode = 'list', length = length(gr_names))
for(i in 1:length(gr_names)){
  
  print(paste0('group', gr_names[i], '_vs_all'))
  
  # j vs all
  pheno_case_tmp <- list(phenoDat_tot[cl$gr == gr_names[i],], phenoDat_tot[cl$gr != gr_names[i],])
  covDat_tmp <- list(covDat[cl$gr == gr_names[i],], covDat[cl$gr != gr_names[i],]) 
  
  new <- do.call(rbind, pheno_case_tmp)
  colnames(new) <- paste0('p', colnames(new))
  gr_id <- factor(c(rep(1, nrow(pheno_case_tmp[[1]])), rep(0, nrow(pheno_case_tmp[[2]]))))
  
  # remove pheno with constant values
  p_rm <- apply(new, 2, function(x) length(unique(x)) == 1)
  new <- new[,!p_rm]
  # remove phenotype with few true overall (find binomial)
  id_b <- intersect(colnames(new), paste0('p',phenoInfo$pheno_id[!phenoInfo$transformed_type %in% c('CONTINUOUS', 'CAT_ORD')]))
  id_rm <- names(which(colSums(new[, id_b, drop = F], na.rm = T) < 50))
  if(length(id_rm)>0){
    new <- new[, !colnames(new) %in% id_rm]
  }
  
  new_cov <- do.call(rbind, covDat_tmp)
  
  res_glm <- data.frame(beta = c(), se_beta = c(), z = c(), CI_low = c(), CI_up = c())
  for(l in 1:ncol(new)){
    print(l)
    type_pheno <- phenoInfo$transformed_type[paste0('p',phenoInfo$pheno_id) == colnames(new)[l]]
    tmp_dat <- cbind(data.frame(pheno = new[, l], gr_id = gr_id), new_cov)
    tmp_dat$cohort <- factor(tmp_dat$cohort)
    res_glm <- rbind(res_glm, compute_reg_endopheno_lmm(mat = tmp_dat, fmla = fmla, type_pheno = type_pheno))
  }
  
  res_glm$type_pheno <- phenoInfo$transformed_type[match(colnames(new), paste0('p',phenoInfo$pheno_id))]
  
  phenoInfo_tmp <- phenoInfo[match(colnames(new), paste0('p',phenoInfo$pheno_id)),]

  bin_reg[[i]] <- cbind(data.frame(pheno_id = phenoInfo_tmp$pheno_id, Field = phenoInfo_tmp$Field, meaning = phenoInfo_tmp$Coding_meaning), res_glm)
  # bin_reg[[i]]$pval_corr <- p.adjust(bin_reg[[i]]$pvalue, method = 'BH')
  bin_reg[[i]]$comp <- sprintf('gr%i_vs_all', gr_names[i])
  
}

tot_bin_reg <- do.call(rbind, bin_reg)
# tot_bin_reg$pval_corr_overall <-  p.adjust(tot_bin_reg$pvalue, method = 'BY')

# save
write.table(x = tot_bin_reg, sprintf('%s%s_%s_cluster%s_PGmethod_%smetric_phenoAssociation_GLM.txt', outFold, type_data, type_input, type_cluster, type_sim), col.names = T, row.names = F, sep = '\t', quote = F)

output$bin_reg = tot_bin_reg

save(output, file = sprintf('%s%s_%s_cluster%s_PGmethod_%smetric_phenoAssociation_GLM.RData', outFold, type_data, type_input, type_cluster, type_sim))

###############################################################################
if(!risk_score){
  test_pheno <- tot_bin_reg
  test_pheno$sign <- 'no'
  test_pheno$sign[test_pheno$pval_corr <= 0.05] <- 'yes'
  test_pheno$logpval <- -log10(test_pheno$pvalue)
  test_pheno$sign <- factor(test_pheno$sign, levels = c('no', 'yes'))
  test_pheno$new_Field <- test_pheno$Field
  test_pheno$new_Field[test_pheno$Field == 'Blood clot, DVT, bronchitis, emphysema, asthma, rhinitis, eczema, allergy diagnosed by doctor'] <- 'Blood clot, etc. diagnosed by doctor'
  test_pheno$new_Field[test_pheno$Field == 'Medication for cholesterol, blood pressure, diabetes, or take exogenous hormones'] <- 'Medication for cholesterol etc. (female)'
  test_pheno$new_Field[test_pheno$Field == 'Medication for cholesterol, blood pressure or diabetes'] <- 'Medication for cholesterol etc. (male)'
  test_pheno$new_id <- test_pheno$new_Field
  test_pheno$new_id[!is.na(test_pheno$meaning)] <- paste(test_pheno$meaning[!is.na(test_pheno$meaning)], test_pheno$new_Field[!is.na(test_pheno$meaning)], sep = '\n')
  
  test_pheno <- test_pheno[order(test_pheno$pvalue),]
  if(length(which(test_pheno$sign == 'yes'))>=20){
    test_pheno <- test_pheno[1:length(which(test_pheno$sign == 'yes'))+1,]
    height_pl <- 5 + (nrow(test_pheno)-20)*0.2
  }else{
    test_pheno <- test_pheno[1:20,] 
    height_pl = 5
  }
  
  test_pheno$comp <- factor(test_pheno$comp)
  test_pheno$new_id <- factor(test_pheno$new_id, levels = unique(test_pheno$new_id))
  
  pl <-  ggplot(test_pheno, aes(x = new_id, y = logpval, fill = sign))+
    geom_bar(stat = 'identity', position = position_dodge(), width = 0.5) + theme_bw()+ 
    ylab('-log10(pvalue)')+xlab('')+geom_hline(yintercept = -log10(0.05))+
    facet_wrap(.~comp, scales = 'free_y', ncol = 1, strip.position="right")+
    scale_fill_manual(values=c("#999999", "#E69F00"))+
    theme(legend.position = 'none', plot.title = element_text(size=9), axis.text.y = element_text(size = 7), strip.text = element_text(size=5))+
    ggtitle(paste(type_data, type_input, 'cluster', type_cluster))+
    coord_flip()
  ggsave(filename = sprintf('%s%s_%s_cluster%s_PGmethod_%smetric_phenoAssociation_GLM.png', outFold, type_data, type_input, type_cluster, type_sim), width = 4.5, height = height_pl, plot = pl, device = 'png')
  ggsave(filename = sprintf('%s%s_%s_cluster%s_PGmethod_%smetric_phenoAssociation_GLM.pdf', outFold, type_data, type_input, type_cluster, type_sim), width = 4.5, height = height_pl, plot = pl, device = 'pdf')
}
