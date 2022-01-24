# associate endophenotypes to cluster structure 
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
suppressPackageStartupMessages(library(coin))
suppressPackageStartupMessages(library(umap))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(RColorBrewer))
options(bitmapType = 'cairo', device = 'png')

parser <- ArgumentParser(description="Associate cluster to endophenotype")
parser$add_argument("--phenoDatFile", type = "character", help = "file to be loaded (endophenotypes, must be a unique matrix)")
parser$add_argument("--phenoDescFile", type = "character", help = "description endophenotype")
parser$add_argument("--sampleAnnFile", type = "character", help = "file with samples to be used")
parser$add_argument("--clusterFile", type = "character", help = "file with clustering structure")
parser$add_argument("--type_cluster", type = "character", default = 'All', help = "All, Cases, Controls")
parser$add_argument("--functR", type = "character", help = "functions to be used")
parser$add_argument("--type_data", type = "character", help = "tscore, path_Reactome, path_GO, PCs")
parser$add_argument("--type_sim", type = "character", default = 'HK', help = "HK or ED or SNF")
parser$add_argument("--type_input", type = "character", default = 'original', help = "original or zscaled")
parser$add_argument("--risk_score", type = "logical", default = F, help = "if true, phenotype is risk score")
parser$add_argument("--rescale_pheno", type = "logical", default = F, help = "if true, continous phenotype rescaled")
parser$add_argument("--outFold", type="character", help = "Output file [basename only]")

args <- parser$parse_args()
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
rescale_pheno <- args$rescale_pheno
outFold <- args$outFold

####################################################################################################################
# sampleAnnFile <- 'INPUT_DATA_GTEx/CAD/Covariates/UKBB/CAD_HARD_clustering/covariateMatrix_CADHARD_All_phenoAssoc_withMedication.txt'
# phenoDatFile <- 'INPUT_DATA_GTEx/CAD/Covariates/UKBB/CAD_HARD_clustering/phenotypeMatrix_CADHARD_All_phenoAssoc_withMedication.txt'
# phenoDescFile <- 'INPUT_DATA_GTEx/CAD/Covariates/UKBB/CAD_HARD_clustering/phenotypeDescription_withMedication.txt'
# clusterFile <- 'INPUT_DATA_GTEx/CAD/Covariates/UKBB/CAD_HARD_clustering/PCs_clusterCases_PGmethod_HKmetric.RData'
# type_cluster <- 'Cases'
# type_data <- 'PCs'
# type_sim <- 'HK'
# outFold <- 'INPUT_DATA_GTEx/CAD/Covariates/UKBB/CAD_HARD_clustering/'
# functR <- '/psycl/g/mpsziller/lucia/priler_project/Software/model_clustering/clustering_functions.R'
# type_input <- 'original'
# rescale_pheno = T
###################################################################################################################

source(functR)

phenoDat <- fread(phenoDatFile, h=T, stringsAsFactor = F, data.table = F)
sampleAnn <- read.table(sampleAnnFile, h=T, stringsAsFactors = F, check.names = F)
cluster_output <- get(load(clusterFile))

sampleAnn <- read.table(sampleAnnFile, h=T, stringsAsFactors = F, check.names = F)

if(type_cluster == 'Cases'){
  sampleAnn <- sampleAnn[sampleAnn$Dx == 1,]
}else{
  if(type_cluster == 'Controls'){
    sampleAnn <- sampleAnn[sampleAnn$Dx == 0,]
  }else{
    if(type_cluster != 'All')
      stop('type_cluster must be either Cases or Controls or All')
  }
}

identical(sampleAnn$Individual_ID, cluster_output$samples_id)

phenoDat <- phenoDat[match(sampleAnn$Individual_ID, phenoDat$Individual_ID),]
phenoDat <- phenoDat[, -1, drop = F]
# remove columns with too many NAs
phenoDat <- phenoDat[,colSums(!is.na(phenoDat))>=100, drop = F]
# remove binary with too few T
id_bin <- rep(F, ncol(phenoDat))
for(i in 1:ncol(phenoDat)){
  id_bin[i] <- is.integer(phenoDat[, i])
}
rm_col <- colnames(phenoDat[, id_bin & colSums(phenoDat != 0 & !is.na(phenoDat))<50, drop = F])
phenoDat <- phenoDat[,!colnames(phenoDat) %in% rm_col, drop = F]

phenoInfo <- read.delim(phenoDescFile, h=T, stringsAsFactors = F, sep = '\t')
# get common pheno id (phenoDescFile can be used to filter)
common_pheno <- intersect(phenoInfo$pheno_id, colnames(phenoDat))
phenoInfo <- phenoInfo[match(common_pheno, phenoInfo$pheno_id),]
phenoDat <- phenoDat[,match(common_pheno,colnames(phenoDat)), drop = F]

if(risk_score){
  phenoInfo$transformed_type <- 'CONTINUOUS'
}

if(rescale_pheno){
  id_c <- which(phenoInfo$transformed_type == 'CONTINUOUS')
  if(length(id_c) > 0){
    for(i in id_c){
      tmp <- scale(phenoDat[, i])
      attr(tmp,  "scaled:center") <- NULL
      attr(tmp,  "scaled:scale") <- NULL
      phenoDat[, i] <- tmp[,1]
    }
  }
}

if(any(phenoInfo$Path %in% 'Online follow-up > Cognitive function online > Fluid intelligence')){
  phenoInfo$Field[phenoInfo$Path == 'Online follow-up > Cognitive function online > Fluid intelligence'] <- paste(phenoInfo$Field[phenoInfo$Path == 'Online follow-up > Cognitive function online > Fluid intelligence'], '(Online)')
}
cl <- cluster_output$cl_best$gr

output <- list(phenoDat = phenoDat, phenoInfo = phenoInfo, cl = cluster_output$cl_best)

#######################################
#### binary regression (gi vs gj) #####
#######################################

# function to check with cad ordinal to remove
remove_pheno_ordinal <- function(pheno_df, group, thr){
  
  name_pheno <- colnames(pheno_df)
  pheno_rm <- c()
  
  for(i in 1:ncol(pheno_df)){
    
    min_p <- min(pheno_df[,i], na.rm = T)
    n_base_gr0 <- sum(pheno_df[group == 0,i] == min_p, na.rm = T)
    n_notbase_gr0 <- sum(pheno_df[group == 0,i] > min_p, na.rm = T)
    n_base_gr1 <- sum(pheno_df[group == 1, i] == min_p, na.rm = T)
    n_notbase_gr1 <- sum(pheno_df[group == 1, i] > min_p, na.rm = T)
    
    if(any(c(n_base_gr0, n_base_gr1, n_notbase_gr0, n_notbase_gr1) < thr)){
       pheno_rm <- c(pheno_rm, name_pheno[i])
    }
  }
  
  return(pheno_rm)
  
}


gr_names <- sort(unique(cl))
P <- length(gr_names)
covDat <- sampleAnn[, !colnames(sampleAnn) %in% c('Individual_ID', 'genoSample_ID', 'Dx')]
if(type_data == 'PCs'){
  covDat <- covDat[, !colnames(covDat) %in% c(paste0('PC', 1:10), paste0('C', 1:10))]
}
fmla  <- as.formula(paste('pheno~gr_id+', paste0(colnames(covDat), collapse = '+')))
output$covDat = covDat

bin_reg <- vector(mode = 'list', length = length(gr_names)-1)
for(i in 1:(length(gr_names)-1)){
  
  print(paste0('group', gr_names[i], '_vs_groupj'))
  
  # j vs all
  pheno_case_tmp <- lapply(gr_names[i:length(gr_names)], function(x) phenoDat[cl == x,,drop = F])
  covDat_tmp <- lapply(gr_names[i:length(gr_names)], function(x) covDat[cl == x,,drop = F])
  bin_reg[[i]] <-  vector(mode = 'list', length = length(pheno_case_tmp)-1)
  
  for(j in 2:length(pheno_case_tmp)){
    
    print(j)
    
    new <- rbind(pheno_case_tmp[[1]], pheno_case_tmp[[j]])
    colnames(new) <- paste0('p', colnames(new))
    gr_id <- factor(c(rep(0, nrow(pheno_case_tmp[[1]])), rep(1, nrow(pheno_case_tmp[[j]]))))
    
    # remove pheno with constant values
    p_rm <- apply(new, 2, function(x) length(unique(x)) == 1)
    new <- new[,!p_rm, drop = F]
    # remove phenotype with few true overall (find binomial)
    id_b <- intersect(colnames(new), paste0('p',phenoInfo$pheno_id[!phenoInfo$transformed_type %in% c('CONTINUOUS', 'CAT_ORD')]))
    id_rm <- names(which(colSums(new[, id_b, drop = F], na.rm = T) < 50))
    if(length(id_rm)>0){
      new <- new[, !colnames(new) %in% id_rm, drop = F]
    }
    # remove phenotype categorical ordinal with less than 10 base 
    # class in each pairwise group (or viceversa)
    id_o <- intersect(colnames(new),
                      paste0('p',phenoInfo$pheno_id[phenoInfo$transformed_type %in% c('CAT_ORD')]))
    id_rm <- remove_pheno_ordinal(pheno_df = new[, id_o], group = gr_id, thr = 10)
    if(length(id_rm)>0){
      new <- new[, !colnames(new) %in% id_rm, drop = F]
    }
    
    new_cov <- rbind(covDat_tmp[[1]], covDat_tmp[[j]])
    res_glm <- matrix(nrow = ncol(new), ncol = 7)
    for(l in 1:ncol(new)){
      type_pheno <- phenoInfo$transformed_type[paste0('p',phenoInfo$pheno_id) == colnames(new)[l]]
      tmp_dat <- cbind(data.frame(pheno = new[, l], gr_id = gr_id), new_cov)
      res_glm[l,] <- compute_reg_endopheno(mat = tmp_dat, fmla = fmla, type_pheno = type_pheno)
    }
    colnames(res_glm) <- c('beta', 'se_beta', 'z', 'pvalue', 'OR_or_Beta', 'CI_low', 'CI_up')
    res_glm <- as.data.frame(res_glm)
    res_glm$type_pheno <- phenoInfo$transformed_type[match(colnames(new), paste0('p',phenoInfo$pheno_id))]
    
    phenoInfo_tmp <- phenoInfo[match(colnames(new), paste0('p',phenoInfo$pheno_id)),]
    
    bin_reg[[i]][[j-1]] <- cbind(data.frame(pheno_id = phenoInfo_tmp$pheno_id, Field = phenoInfo_tmp$Field, meaning = phenoInfo_tmp$Coding_meaning), res_glm)
    bin_reg[[i]][[j-1]]$pval_corr <- p.adjust(bin_reg[[i]][[j-1]]$pvalue, method = 'BH')
    bin_reg[[i]][[j-1]]$comp <- sprintf('gr%i_vs_gr%i', gr_names[i:length(gr_names)][j], gr_names[i])
    
  }
  
  
}

tot_bin_reg <- do.call(rbind, do.call(c,bin_reg))
tot_bin_reg$pval_corr_overall <-  p.adjust(tot_bin_reg$pvalue, method = 'BY')

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

output <- list(phenoDat = phenoDat, phenoInfo = phenoInfo, cl = cluster_output$cl_best)

########################################
#### binary regression (gi vs all) #####
########################################

gr_names <- sort(unique(cl))
P <- length(gr_names)
covDat <- sampleAnn[, !colnames(sampleAnn) %in% c('Individual_ID', 'genoSample_ID', 'Dx')]
fmla  <- as.formula(paste('pheno~gr_id+', paste0(colnames(covDat), collapse = '+')))
output$covDat = covDat

bin_reg <- vector(mode = 'list', length = length(gr_names))
for(i in 1:length(gr_names)){
  
  print(paste0('group', gr_names[i], '_vs_all'))
  
  # j vs all
  pheno_case_tmp <- list(phenoDat[cl == gr_names[i],, drop = F], phenoDat[cl != gr_names[i],, drop = F])
  covDat_tmp <- list(covDat[cl == gr_names[i],,drop = F], covDat[cl != gr_names[i],,drop = F]) 
  
  new <- do.call(rbind, pheno_case_tmp)
  colnames(new) <- paste0('p', colnames(new))
  gr_id <- factor(c(rep(1, nrow(pheno_case_tmp[[1]])), rep(0, nrow(pheno_case_tmp[[2]]))))
  
  # remove pheno with constant values
  p_rm <- apply(new, 2, function(x) length(unique(x)) == 1)
  new <- new[,!p_rm, drop = F]
  # remove phenotype with few true overall (find binomial)
  id_b <- intersect(colnames(new), paste0('p',phenoInfo$pheno_id[!phenoInfo$transformed_type %in% c('CONTINUOUS', 'CAT_ORD')]))
  id_rm <- names(which(colSums(new[, id_b, drop = F], na.rm = T) < 50))
  if(length(id_rm)>0){
    new <- new[, !colnames(new) %in% id_rm, drop = F]
  }
  
  # remove phenotype categorical ordinal with less than 10 base 
  # class in each pairwise group (or viceversa)
  id_o <- intersect(colnames(new),
                    paste0('p',phenoInfo$pheno_id[phenoInfo$transformed_type %in% c('CAT_ORD')]))
  id_rm <- remove_pheno_ordinal(pheno_df = new[, id_o], group = gr_id, thr = 10)
  if(length(id_rm)>0){
    new <- new[, !colnames(new) %in% id_rm, drop = F]
  }
  
  new_cov <- do.call(rbind, covDat_tmp)
  res_glm <- matrix(nrow = ncol(new), ncol = 7)
  for(l in 1:ncol(new)){
    print(l)  
    type_pheno <- phenoInfo$transformed_type[paste0('p',phenoInfo$pheno_id) == colnames(new)[l]]
    tmp_dat <- cbind(data.frame(pheno = new[, l], gr_id = gr_id), new_cov)
    res_glm[l,] <- compute_reg_endopheno(mat = tmp_dat, fmla = fmla, type_pheno = type_pheno)
  }
  
  colnames(res_glm) <- c('beta', 'se_beta', 'z', 'pvalue', 'OR_or_Beta', 'CI_low', 'CI_up')
  res_glm <- as.data.frame(res_glm)
  res_glm$type_pheno <- phenoInfo$transformed_type[match(colnames(new), paste0('p',phenoInfo$pheno_id))]
  
  phenoInfo_tmp <- phenoInfo[match(colnames(new), paste0('p',phenoInfo$pheno_id)),]
  
  bin_reg[[i]] <- cbind(data.frame(pheno_id = phenoInfo_tmp$pheno_id, Field = phenoInfo_tmp$Field, meaning = phenoInfo_tmp$Coding_meaning), res_glm)
  bin_reg[[i]]$pval_corr <- p.adjust(bin_reg[[i]]$pvalue, method = 'BH')
  bin_reg[[i]]$comp <- sprintf('gr%i_vs_all', gr_names[i])
  
}


tot_bin_reg <- do.call(rbind, bin_reg)
tot_bin_reg$pval_corr_overall <-  p.adjust(tot_bin_reg$pvalue, method = 'BY')

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

