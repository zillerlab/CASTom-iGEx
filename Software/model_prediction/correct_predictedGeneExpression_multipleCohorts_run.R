#### code written by Lucia Trastulla, e-mail: lucia_trastulla@psych.mpg.de ####
# correct predicted gene espression across multiple cohorts
# save new predicted gene expression for all cohort combined

options(stringsAsFactors=F)
options(max.print=1000)
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(lme4))
suppressPackageStartupMessages(library(umap))


parser <- ArgumentParser(description="Correct predicted gene expression combining multiple cohorts")

parser$add_argument("--input_file", type = "character",  nargs = '*', help = "total path predicted gene expression (gz file)")
parser$add_argument("--name_cohorts", type = "character", nargs = '*', help = "name cohorts")
parser$add_argument("--covDat_file", type = "character", nargs = '*', help = "file with covariate to use for correction")
parser$add_argument("--sampleAnn_file", type = "character", nargs = '*', help = "file with sample info to use as reference (only Dx == 0)")
parser$add_argument("--thr_reliableGenes", type = "double", nargs = '*', default = c(0.01, 0), help = "threshold for reliable genes: dev_geno_tot and test_dev_geno")
parser$add_argument("--n_comp_umap", type = "integer", default = 2, help = "number of component for umap dim reduction")
parser$add_argument("--n_neigh_umap", type = "integer", default = 30, help = "number of neighbours for umap dim reduction")
parser$add_argument("--min_dist_umap", type = "double", default = 0.01, help = "minimum dist for umap dim reduction")
parser$add_argument("--seed_umap", type = "integer", default = 50, help = "seed for umap dim reduction")
parser$add_argument("--outFold", type="character", help = "Output folder")
parser$add_argument("--outFoldAnn", type="character", help = "Output folder for sample annotation")


args <- parser$parse_args()
name_cohorts <- args$name_cohorts
input_file <- args$input_file
covDat_file <- args$covDat_file
sampleAnn_file <- args$sampleAnn_file
thr_reliableGenes <- args$thr_reliableGenes
n_comp_umap <- args$n_comp_umap
n_neigh_umap <- args$n_neigh_umap
min_dist_umap <- args$min_dist_umap
seed_umap <- args$seed_umap
outFold <- args$outFold
outFoldAnn <- args$outFoldAnn

####################################################################
# name_cohorts <- c('CG', 'German1',  'German2',  'German3',  'German4',  'German5', 'LURIC', 'MG', 'WTCCC')
# thr_reliableGenes <- c(0.01, 0)
# covDat_file <- paste0('/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/INPUT_DATA_GTEx/CAD/Covariates/',name_cohorts,'/covariateMatrix.txt')
# sampleAnn_file <- paste0('/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/INPUT_DATA_GTEx/CAD/Covariates/',name_cohorts,'/covariateMatrix.txt')
# input_file <-  paste0('/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/OUTPUT_GTEx/predict_CAD/Liver/200kb/CAD_GWAS_bin5e-2/',name_cohorts,'/predictedExpression.txt.gz')
# outFold = '/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/OUTPUT_GTEx/predict_CAD/Liver/200kb/CAD_GWAS_bin5e-2/Meta_Analysis_CAD/devgeno0.01_testdevgeno0/'
# n_neigh_umap <- 30
# min_dist_umap <- 0.01
# seed_umap <- 50
# n_comp_umap <- 2
# outFoldAnn <- '/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/INPUT_DATA_GTEx/CAD/Covariates/Meta_Analysis_CAD/'
####################################################################

sampleAnn_c <- list()
covDat_c <- list()
expDat_c <- list()
for(i in 1:length(name_cohorts)){
  
  print(name_cohorts[i])
  
  ### load data ###
  sampleAnn_c[[i]] <- read.table(sampleAnn_file[i], header = T, stringsAsFactors = F)
  expDat_c[[i]] <- read.table(gzfile(input_file[i]), sep="\t", header=T,  check.names = F, stringsAsFactors = F)
  covDat_c[[i]] <- read.table(covDat_file[i], header = T, stringsAsFactors = F)
  common_s <- intersect(intersect(sampleAnn_c[[i]]$Individual_ID, colnames(expDat_c[[i]])), covDat_c[[i]]$Individual_ID)
  
  sampleAnn_c[[i]] <-  sampleAnn_c[[i]][match(common_s, sampleAnn_c[[i]]$Individual_ID),]
  sampleAnn_c[[i]]$cohort <- name_cohorts[i]

  covDat_c[[i]] <-  covDat_c[[i]][match(common_s, covDat_c[[i]]$Individual_ID),]
  covDat_c[[i]]$cohort <- name_cohorts[i]
 
  # filter genes
  expDat_c[[i]] <- expDat_c[[i]][!(is.na(expDat_c[[i]]$dev_geno) | is.na(expDat_c[[i]]$test_dev_geno)), ]
  expDat_c[[i]] <- expDat_c[[i]][expDat_c[[i]]$dev_geno >= thr_reliableGenes[1], ]
  expDat_c[[i]] <- expDat_c[[i]][expDat_c[[i]]$test_dev_geno > thr_reliableGenes[2], ]
  geneInfo <- expDat_c[[i]][, -match(common_s, colnames(expDat_c[[i]]))]
  expDat_c[[i]] <- as.matrix(expDat_c[[i]][,match(common_s, colnames(expDat_c[[i]]))])
  rownames(expDat_c[[i]]) <- geneInfo$external_gene_name

  if(!identical(sampleAnn_c[[i]]$Individual_ID, colnames(expDat_c[[i]])) | !identical(covDat_c[[i]]$Individual_ID, colnames(expDat_c[[i]]))){ stop('Sample names do not match') }
  rownames(covDat_c[[i]]) <- covDat_c[[i]]$Individual_ID
  covDat_c[[i]] <-  covDat_c[[i]][, !colnames(covDat_c[[i]]) %in% c('Individual_ID', 'genoSample_ID', 'Dx')]
  
}

# create new names (some may be repeated across cohorts)
common_column <- names(which(table(unlist(sapply(sampleAnn_c, colnames))) == length(name_cohorts)))
sampleAnn_c <- lapply(sampleAnn_c, function(x) x[,colnames(x) %in% common_column])
sampleAnn <- do.call(rbind, sampleAnn_c)
sampleAnn$cohort_Individual_ID <- sampleAnn$Individual_ID
sampleAnn$Individual_ID <- paste(sampleAnn$cohort, sampleAnn$Individual_ID, sep = '_')
expDat <- t(do.call(cbind, expDat_c))
print(identical(rownames(expDat), sampleAnn$cohort_Individual_ID))
rownames(expDat) <- sampleAnn$Individual_ID

# consider common covariates variables 
common_column <- names(which(table(unlist(sapply(covDat_c, colnames))) == length(name_cohorts)))
covDat_c <- lapply(covDat_c, function(x) x[,colnames(x) %in% common_column])
covDat <- do.call(rbind, covDat_c)
rownames(covDat) <- sampleAnn$Individual_ID
covDat$cohort <- factor(covDat$cohort, levels = name_cohorts)
block_chr <- paste0(colnames(covDat)[!colnames(covDat) %in% c('cohort')], collapse = '+')
# fml <- as.formula(paste('y ~ ', block_chr, ' + (',block_chr,'|cohort)'))
fml <- as.formula(paste('y ~ ', block_chr, ' + ( 1 | cohort )')) # random intercept but fixed slope, random slope unfeasable computational time

corr_expData <- matrix(ncol = ncol(expDat), nrow = nrow(expDat))
var_cohort <- c()
for(g in 1:ncol(expDat)){
  print(g)
  tmp <- cbind(data.frame( y = expDat[,g]), covDat)
  fm <- lmer(fml,  tmp)
  var_cohort <- c(var_cohort, as.numeric(summary(fm)$varcor$cohort))
  corr_expData[,g] <- resid(fm)
}

colnames(corr_expData) <- colnames(expDat)
rownames(corr_expData) <- rownames(expDat)

# is the cohort information still present?
ktest <- apply(expDat,2,function(x) kruskal.test(x = x, g = factor(sampleAnn$cohort))$p.value)
ktest_corr <- apply(corr_expData,2,function(x) kruskal.test(x = x, g = factor(sampleAnn$cohort))$p.value)
print(identical(colnames(expDat),  geneInfo$external_gene_name))
df_cohort <- data.frame(ensembl_gene_id = geneInfo$ensembl_gene_id, external_gene_name = geneInfo$external_gene_name, ktest_pval_exp = ktest, ktest_pval_expcorr = ktest_corr)
# correlation original expDat - new version
df_cohort$corr <- sapply(1:ncol(expDat), function(x) cor(expDat[,x], corr_expData[,x]))
df_cohort$var_cohort <- var_cohort

# save new annotation:
write.table(sampleAnn, file = sprintf('%scovariateMatrix.txt', outFoldAnn), quote = F, row.names = F, col.names=T , sep = '\t')
write.table(sampleAnn[, colnames(sampleAnn) %in% c('Individual_ID', 'Dx')], file = sprintf('%sphenoMatrix.txt', outFoldAnn), quote = F, row.names = F, col.names=T , sep = '\t')

# save new expression:
corr_expData_save <- t(corr_expData)
corr_expData_save <- as.data.frame(cbind(geneInfo, corr_expData_save))
write.table(corr_expData_save, file = sprintf('%scorrected_predictedExpression.txt', outFold), quote = F, row.names = F, col.names=T , sep = '\t')
system(paste("gzip",sprintf('%scorrected_predictedExpression.txt', outFold)))

# save info conversion
write.table(df_cohort, file = sprintf('%scorrected_predictedExpression_info.txt', outFold), quote = F, row.names = F, col.names=T , sep = '\t')

print('conversion concluded and saved')

##############
#### umap ####

custom.settings = umap.defaults
custom.settings$min_dist = min_dist_umap
custom.settings$random_state = seed_umap
custom.settings$n_neighbors = n_neigh_umap
custom.settings$n_components <- n_comp_umap

umap_res <- umap(corr_expData, custom.settings)
output <- list(sampleAnn = sampleAnn, umap = umap_res)
save(output, file =  sprintf('%scorrected_predictedExpression_umap.RData', outFold))

