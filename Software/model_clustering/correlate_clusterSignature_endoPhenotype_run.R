# correlate z-stat phenotype with group signature
# use only the tissue that was used for clustering
# remove higly correlated (compute externally)?

options(stringsAsFactors=F)
options(max.print=1000)
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggsci))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(apcluster))
suppressPackageStartupMessages(library(qvalue))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(ggsignif))
suppressPackageStartupMessages(library(ggpubr))
options(bitmapType = 'cairo', device = 'png')


parser <- ArgumentParser(description="correlation cluser signature endophenotypes")
parser$add_argument("--clusterFeat_fold", type = "character", help = "")
parser$add_argument("--inputPheno_fold", type = "character", help = "")
parser$add_argument("--pheno_name_comp", type = "character", help = "")
parser$add_argument("--phenoFold", type = "character", help = "")
parser$add_argument("--type_corr", type = "character", default = 'pearson', help = "")
parser$add_argument("--outFold", type="character", help = "Output file [basename only]")

args <- parser$parse_args()
clusterFeat_fold <- args$clusterFeat_fold
inputPheno_fold <- args$inputPheno_fold
pheno_name_comp <- args$pheno_name_comp
phenoFold <- args$phenoFold
type_corr <- args$type_corr
outFold <- args$outFold

########################################################################################################################
# clusterFeat_fold <- '/psycl/g/mpsziller/lucia/SCZ_PGC/eQTL_PROJECT/clustering_res/Brain_Cortex/'
# inputPheno_fold <- paste0('/psycl/g/mpsziller/lucia/UKBB/eQTL_PROJECT/OUTPUT_GTEx/predict_UKBB/Brain_Cortex/200kb/noGWAS/devgeno0.01_testdevgeno0/')
# outFold <- clusterFeat_fold
# pheno_name_comp <- 'SCZ'
# type_corr = 'pearson'
# phenoFold <- '/psycl/g/mpsziller/lucia/UKBB/eQTL_PROJECT/INPUT_DATA/Covariates/'
##########################################################################################################################

permutation_test_corr <-  function(v1, v2, cor_value, nperm = 10000, method_cor = 'spearman'){
  
  n <- length(v2)
  cor_perm <- c()
  for(i in 1:nperm){
    set.seed(56*i)
    cor_perm <- c(cor_perm , cor(v2[sample(1:n, n, replace = F)], v1, method = method_cor))
  }
  
  prob <- sum(cor_perm > abs(cor_value) | cor_perm < -abs(cor_value))/nperm
  return(prob)
  
}


if(grepl('CAD',pheno_name_comp)){
  pheno_name <- c(read.table(sprintf('%smatch_cov_pheno_CADrel_filter.txt', phenoFold), h=F, stringsAsFactors = F)$V1, 'ICD10_Anaemia', 'ICD10_Circulatory_system', 'ICD10_Endocrine', 'ICD10_Respiratory_system')
  pheno_name <- pheno_name[!pheno_name %in% c('Medication', 'Medical_conditions')]
  pheno_input <- read.delim(sprintf('%sphenotypeDescription_PHESANTproc_CADrelatedpheno_annotated.txt', phenoFold), h=T, stringsAsFactors = F, sep = '\t')
}

if(grepl('SCZ', pheno_name_comp)){
  
  pheno_name <- c(read.table(sprintf('%smatch_cov_pheno.txt', phenoFold), h=F, stringsAsFactors = F)$V1, 'mixedpheno_Psychiatric', 'ICD10_Psychiatric', 'ICD9_Psychiatric', 
                  read.table(sprintf('%smatch_cov_pheno_SchunkertApp.txt', phenoFold), h=F, stringsAsFactors = F)$V1, 'Blood_biochemistry')
  pheno_input <- read.delim(sprintf('%sphenotypeDescription_PHESANTproc.txt', phenoFold), h=T, stringsAsFactors = F, sep = '\t')
  tmp <- read.delim(sprintf('%sphenotypeDescription_PHESANTproc_CADrelatedpheno.txt', phenoFold), h=T, stringsAsFactors = F, sep = '\t')
  pheno_input <- rbind(pheno_input, tmp)
  pheno_input <- pheno_input[!duplicated(pheno_input$pheno_id),]
  pheno_name <- pheno_name[!pheno_name %in% c('Numeric_memory','Diffusion_brain_MRI', 'Medical_conditions', 'ICD10_Psychiatric', 'ICD9_Psychiatric', 'Resting_functional_brain_MRI',
                                              'dMRI_skeleton', 'T2-weighted_brain_MRI', 'Estimated_nutrients_yesterday')]
  
}

# load clustering res
cl_genes <- read.delim(sprintf('%stscoreOriginal_tscoreClusterCases_featAssociation.txt', clusterFeat_fold), header = T, stringsAsFactors = F, sep = '\t')
cl_genes$se <- (cl_genes$CI_up - cl_genes$estimates)/1.96
genes_cl <- unique(cl_genes$feat)
gr_tot <- sapply(unique(cl_genes$comp), function(x) strsplit(x, split = '_vs_all')[[1]][1])
mat_est_genes <- sapply(gr_tot, function(x) sapply(genes_cl, function(y) cl_genes$estimates[cl_genes$comp == paste0(x, '_vs_all') & cl_genes$feat == y]))

tscore_pheno_corr <- list()
tscore_pheno_pval <- list()
pheno_info <- list()

for(j in 1:length(pheno_name)){
  
  print(pheno_name[j])
  
  if(pheno_name[j] %in% c('Blood_biochemistry', 'Blood_count') & grepl('CAD', pheno_name_comp)){
    file_toload <- sprintf('%s/pval_%s_withMed_pheno_covCorr.RData', inputPheno_fold, pheno_name[j])
  }else{
    file_toload <- sprintf('%s/pval_%s_pheno_covCorr.RData', inputPheno_fold, pheno_name[j])
  }
  
  tmp <- get(load(file_toload[1]))
  rm(final)
  
  id_keep <- 1:nrow(tmp$pheno)
  if(grepl('SCZ', pheno_name_comp) & pheno_name[j] == 'Alcohol_use'){
    id_keep <-  which(tmp$pheno$Field %in% c('Amount of alcohol drunk on a typical drinking day', 'Frequency of drinking alcohol', 'Frequency of consuming six or more units of alcohol'))
  }
  if(grepl('SCZ', pheno_name_comp) & pheno_name[j] == 'Cannabis_use'){
    id_keep <-  which(tmp$pheno$Field %in% c('Ever taken cannabis', 'Maximum frequency of taking cannabis'))
  }
  if(grepl('SCZ', pheno_name_comp) & pheno_name[j] == 'Anxiety'){
    id_keep <-  which(!(grepl('Recent', tmp$pheno$Field) | grepl('undertaken', tmp$pheno$Field)))
  }
  if(grepl('SCZ', pheno_name_comp) & pheno_name[j] == 'Depression'){
    id_keep <-  which(!(grepl('Recent', tmp$pheno$Field) | grepl('undertaken', tmp$pheno$Field)))
  }
  if(grepl('SCZ', pheno_name_comp) & pheno_name[j] == 'dMRI_weighted_means'){
    id_keep <-  which(grepl('Weighted-mean MD in', tmp$pheno$Field) | grepl('Weighted-mean ICVF in', tmp$pheno$Field))
  }
  if(grepl('SCZ', pheno_name_comp) & pheno_name[j] == 'Fluid_intelligence'){
    id_keep <-  which(grepl('UK Biobank Assessment Centre',tmp$pheno$Path) & !tmp$pheno$Field %in% c('Attempted fluid intelligence (FI) test.', 'Fluid intelligence completion status'))
  }
  if(grepl('SCZ', pheno_name_comp) & pheno_name[j] == 'Pairs_matching'){
    id_keep <- which(grepl('UK Biobank Assessment Centre',tmp$pheno$Path) & !tmp$pheno$Field %in% c('Number of incorrect matches in round'))
  }
  if(grepl('SCZ', pheno_name_comp) & pheno_name[j] == 'Prospective_memory'){
    id_keep <-  which(tmp$pheno$Field %in% c('Time to answer', 'Prospective memory result'))
  }
  if(grepl('SCZ', pheno_name_comp) & pheno_name[j] == 'Sleep'){
    id_keep <-  which(tmp$pheno$Field %in% c('Sleep duration', 'Sleepness / insomnia', 'Daytime dozing / sleeping (narcolepsy)'))
  }
  if(grepl('SCZ', pheno_name_comp) & pheno_name[j] == 'Smoking'){
    id_keep <-  which(tmp$pheno$Field %in% c('Current tobacco smoking', 'Past tobacco smoking', 'Light smokers, at least 100 smokes in lifetime', 'Number of cigarettes currently smoked daily (current cigarette smokers)', 'Ever smoked'))
  }
  if(grepl('SCZ', pheno_name_comp) & pheno_name[j] == 'Susceptibility_weighted_brain_MRI'){
    id_keep <-  which(!tmp$pheno$Field %in% c('Discrepancy between SWI brain image and T1 brain image'))
  }
  if(grepl('SCZ', pheno_name_comp) & pheno_name[j] == 'Symbol_digit_substitution'){
    id_keep <-  which(!tmp$pheno$Field %in% c('Symbol digit completion status'))
  }
  if(grepl('SCZ', pheno_name_comp) & pheno_name[j] =='T1_structural_brain_MRI'){
    id_keep <- which(grepl('Volume',tmp$pheno$Field))
  }
  if(grepl('SCZ', pheno_name_comp) & pheno_name[j] == 'Trail_making'){
    id_keep <-  which(!tmp$pheno$Field %in% c('Trail making completion status'))
  }
  
  pheno_info[[j]] <- data.frame(pheno = tmp$pheno$pheno_id[id_keep], pheno_type = pheno_name[j], Field = tmp$pheno$Field[id_keep], meaning = tmp$pheno$Coding_meaning[id_keep])
  tscore_pheno_corr[[j]] <- tscore_pheno_pval[[j]] <- matrix(nrow = length(id_keep), ncol = length(gr_tot))
  colnames(tscore_pheno_corr[[j]]) <- colnames(tscore_pheno_pval[[j]]) <- gr_tot
  
  tmp_pheno <- tmp$tscore[id_keep]
  inters_genes <- lapply(tmp_pheno, function(x) intersect(genes_cl, x$external_gene_name))
  tmp_pheno <- mapply(function(x, y) x[match(y,x$external_gene_name),], x = tmp_pheno, y = inters_genes, SIMPLIFY = F)
  mat_est_genes_tmp <- lapply(inters_genes, function(y) mat_est_genes[match(y,genes_cl),])
  
  for(i in 1:length(id_keep)){
    tscore_pheno_corr[[j]][i, ] <- apply(mat_est_genes_tmp[[i]], 2, function(x) cor(tmp_pheno[[i]][, 5], x,  method = type_corr)) # compare with beta and not zeta statistics
    tscore_pheno_pval[[j]][i, ] <- sapply(1:length(gr_tot), function(x) permutation_test_corr(v1 = mat_est_genes_tmp[[i]][,x], v2 = tmp_pheno[[i]][, 5], 
                                                                                              cor_value = tscore_pheno_corr[[j]][i, x], nperm = 10000, method = type_corr))
  }
}

pheno_info <- do.call(rbind, pheno_info)
pheno_info$names_field <- NA
pheno_info$names_field[is.na(pheno_info$meaning)] <- paste(pheno_info$Field, paste0('(',pheno_info$pheno,')'))[is.na(pheno_info$meaning)]
pheno_info$names_field[!is.na(pheno_info$meaning)] <- paste(pheno_info$Field, ':', pheno_info$meaning, paste0('(',pheno_info$pheno,')'))[!is.na(pheno_info$meaning)]

tscore_pheno_corr <- as.data.frame(do.call(rbind, tscore_pheno_corr))
tscore_pheno_pval <- as.data.frame(do.call(rbind, tscore_pheno_pval))

tscore_pheno_corr$names_field <- pheno_info$names_field
tscore_pheno_pval$names_field <- pheno_info$names_field

# save
res_cor_cl <- list(corr = tscore_pheno_corr, pval = tscore_pheno_pval, pheno = pheno_info)
save(res_cor_cl, file = sprintf('%sclusterCases_%s_groupCorrelation_relatedPheno.RData', outFold, pheno_name_comp))


