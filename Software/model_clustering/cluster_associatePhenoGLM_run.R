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
parser$add_argument("--type_data", type = "character", help = "tscore, path_Reactome or path_GO")
parser$add_argument("--type_sim", type = "character", default = 'HK', help = "HK or ED or SNF")
parser$add_argument("--type_input", type = "character", default = 'original', help = "original or zscaled")
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
outFold <- args$outFold

# ####################################################################################################################
# phenoDatFile <- '/psycl/g/mpsziller/lucia/UKBB/eQTL_PROJECT/INPUT_DATA/Covariates/MDD_pheno_def/phenotypeMatrix_endoPh_LifetimeMDD.txt'
# phenoDescFile <- '/psycl/g/mpsziller/lucia/UKBB/eQTL_PROJECT/INPUT_DATA/Covariates/MDD_pheno_def/phenotypeDescription_endoPh_LifetimeMDD.txt'
# sampleAnnFile <- '/psycl/g/mpsziller/lucia/UKBB/eQTL_PROJECT/INPUT_DATA/Covariates/MDD_pheno_def/covariateMatrix_LifetimeMDD.txt'
# clusterFile <- '/psycl/g/mpsziller/lucia/UKBB/eQTL_PROJECT/OUTPUT_all/LifetimeMDD_pheno/tscore_zscaled_clusterCases_PGmethod_SNFmetric.RData'
# type_cluster <- 'Cases'
# type_data <- 'tscore'
# type_sim <- 'SNF'
# outFold <- '/psycl/g/mpsziller/lucia/UKBB/eQTL_PROJECT/OUTPUT_CMC/predict_UKBB/200kb/devgeno0.01_testdevgeno0/LifetimeMDD_pheno/'
# functR <- '/psycl/g/mpsziller/lucia/priler_project/Software/model_clustering/clustering_functions.R'
# type_input <- 'zscaled'
# ####################################################################################################################

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
phenoDat <- phenoDat[, -1]
# remove columns with too many NAs
phenoDat <- phenoDat[,colSums(!is.na(phenoDat))>=100]
# remove binary with too few T
id_bin <- rep(F, ncol(phenoDat))
for(i in 1:ncol(phenoDat)){
  id_bin[i] <- is.integer(phenoDat[, i])
}
rm_col <- colnames(phenoDat[, id_bin & colSums(phenoDat != 0 & !is.na(phenoDat))<50])
phenoDat <- phenoDat[,!colnames(phenoDat) %in% rm_col]

phenoInfo <- read.delim(phenoDescFile, h=T, stringsAsFactors = F, sep = '\t')
phenoInfo <- phenoInfo[match(colnames(phenoDat), phenoInfo$pheno_id),]
if(any(phenoInfo$Path %in% 'Online follow-up > Cognitive function online > Fluid intelligence')){
  phenoInfo$Field[phenoInfo$Path == 'Online follow-up > Cognitive function online > Fluid intelligence'] <- paste(phenoInfo$Field[phenoInfo$Path == 'Online follow-up > Cognitive function online > Fluid intelligence'], '(Online)')
}
cl <- cluster_output$cl_best$gr

output <- list(phenoDat = phenoDat, phenoInfo = phenoInfo, cl = cluster_output$cl_best)

#######################################
#### binary regression (gi vs gj) #####
#######################################

gr_names <- sort(unique(cl))
P <- length(gr_names)
covDat <- sampleAnn[, !colnames(sampleAnn) %in% c('Individual_ID', 'genoSample_ID', 'Dx')]
fmla  <- as.formula(paste('pheno~gr_id+', paste0(colnames(covDat), collapse = '+')))
output$covDat = covDat

bin_reg <- vector(mode = 'list', length = length(gr_names)-1)
for(i in 1:(length(gr_names)-1)){
  
  print(paste0('group', gr_names[i], '_vs_groupj'))
  
  # j vs all
  pheno_case_tmp <- lapply(gr_names[i:length(gr_names)], function(x) phenoDat[cl == x,])
  covDat_tmp <- lapply(gr_names[i:length(gr_names)], function(x) covDat[cl == x,])
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
    id_rm <- names(which(colSums(new[, id_b], na.rm = T) < 50))
    if(length(id_rm)>0){
      new <- new[, !colnames(new) %in% id_rm]
    }
    
    new_cov <- rbind(covDat_tmp[[1]], covDat_tmp[[j]])
    res_glm <- matrix(nrow = ncol(new), ncol = 4)
    for(l in 1:ncol(new)){
      type_pheno <- phenoInfo$transformed_type[paste0('p',phenoInfo$pheno_id) == colnames(new)[l]]
      tmp_dat <- cbind(data.frame(pheno = new[, l], gr_id = gr_id), new_cov)
      res_glm[l,] <- compute_reg_endopheno(mat = tmp_dat, fmla = fmla, type_pheno = type_pheno)
    }
    colnames(res_glm) <- c('beta', 'se_beta', 'z', 'pvalue')
    res_glm <- as.data.frame(res_glm)
    res_glm$type_pheno <- phenoInfo$transformed_type[match(colnames(new), paste0('p',phenoInfo$pheno_id))]
    
    phenoInfo_tmp <- phenoInfo[match(colnames(new), paste0('p',phenoInfo$pheno_id)),]
    
    bin_reg[[i]][[j-1]] <- cbind(data.frame(pheno_id = phenoInfo_tmp$pheno_id, Field = phenoInfo_tmp$Field, meaning = phenoInfo_tmp$Coding_meaning), res_glm)
    bin_reg[[i]][[j-1]]$pval_corr <- p.adjust(bin_reg[[i]][[j-1]]$pvalue, method = 'BY')
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


