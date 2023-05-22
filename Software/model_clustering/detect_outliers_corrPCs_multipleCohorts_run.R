#!/usr/bin/env Rscript
# cluster (multiple cohort combined)

options(stringsAsFactors=F)
options(max.print=1000)
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(qvalue))
suppressPackageStartupMessages(library(igraph))
suppressPackageStartupMessages(library(pryr))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(coin))
suppressPackageStartupMessages(library(umap))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(RColorBrewer))
options(bitmapType = 'cairo', device = 'png')


parser <- ArgumentParser(description="remove outliers if needed")
parser$add_argument("--inputFile", type = "character", nargs = '*', help = "file to be loaded (predicted tscore or pathScore)")
parser$add_argument("--name_cohorts", type = "character", nargs = '*', help = "name of the single cohorts")
parser$add_argument("--sampleAnnFile", type = "character", nargs = '*', help = "file with samples to be used")
parser$add_argument("--geneRegionFile", type = "character", default=NULL, help = "used if tscore and exclude_MHC")
parser$add_argument("--tissues_name", type = "character", help = "name tissue")
parser$add_argument("--genes_to_filter", type = "character", default = NULL, help = "additional file to filter genes")
parser$add_argument("--exclude_MHC", type = "logical", default = F, help = "if true, MHC region excluded (only ossible for tscore)")
parser$add_argument("--type_cluster", type = "character", help = "All, Cases, Controls")
parser$add_argument("--split_tot", type = "integer", default = 0, help = "if 0 then inpuntFile load alone, otherwise splitted version")
parser$add_argument("--pvalresFile", type = "character", help = "file with pvalue results")
parser$add_argument("--pval_id", type = "integer", default = 1, help = "id to be used on pvalue file")
parser$add_argument("--corr_thr", type = "double", default = -1, help = "correlation among features threshold")
parser$add_argument("--functR", type = "character", help = "functions to be used")
parser$add_argument("--type_data", type = "character", help = "tscore, path_Reactome or path_GO")
parser$add_argument("--type_sim", type = "character", default = 'HK', help = "HK or ED")
parser$add_argument("--type_input", type = "character", default = 'original', help = "original or zscaled")
parser$add_argument("--min_genes_path", type = "integer", default = 1, help = "minimum number of genes for a pathway, if > 1 recompute corrected pvalues")
parser$add_argument("--outFold", type="character", help = "Output file [basename only]")

args <- parser$parse_args()
genes_to_filter <- args$genes_to_filter
name_cohorts <- args$name_cohorts
pvalresFile <- args$pvalresFile
tissues_name <- args$tissues_name
pval_id <- args$pval_id
inputFile <- args$inputFile
exclude_MHC <- args$exclude_MHC
geneRegionFile <- args$geneRegionFile
sampleAnnFile <- args$sampleAnnFile
type_cluster <- args$type_cluster
split_tot <- args$split_tot
functR <- args$functR
type_data <- args$type_data
corr_thr <- args$corr_thr
min_genes_path <- args$min_genes_path
type_sim <- args$type_sim
type_input <- args$type_input
outFold <- args$outFold

####################################################################################################################
# name_cohorts <- c("German1", "German2")
# geneRegionFile <- NULL
# exclude_MHC <- T
# inputFile <- paste0('OUTPUT_GTEx/predict_CAD/Liver/200kb/CAD_GWAS_bin5e-2/',name_cohorts,'/devgeno0.01_testdevgeno0/predictedTscores.txt')
# sampleAnnFile <- paste0('INPUT_DATA_GTEx/CAD/Covariates/',name_cohorts,'/covariateMatrix.txt')
# covDatFile <- sampleAnnFile
# split_tot <- 0
# pvalresFile <- 'OUTPUT_GTEx/predict_CAD/Liver/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/pval_CAD_pheno_covCorr.RData'
# pval_id <- 1
# min_genes_path <- 2
# type_data <- 'tscore'
# type_cluster <- 'Cases'
# type_sim <- 'HK'
# outFold <- 'OUTPUT_GTEx/predict_CAD/Liver/200kb/CAD_GWAS_bin5e-2/'
# functR <- '/psycl/g/mpsziller/lucia/castom-igex/Software/model_clustering/clustering_functions.R'
# corr_thr <- 0.9
# type_input <- 'zscaled'
# kNN_par <- 30
# tissues_name <- 'Liver'
# genes_to_filter <- NULL
####################################################################################################################

source(functR)

# load pval res
res_pval <- get(load(pvalresFile))
if(type_data == 'tscore'){
  res_pval <- res_pval$tscore[[pval_id]]
  id_pval <- 8
  id_info <- 2
}else{
  if(type_data == 'path_Reactome'){
    res_pval <- res_pval$pathScore_reactome[[pval_id]]
    id_pval <- 13
    id_info <- 1
  }else{
    if(type_data == 'path_GO'){
      res_pval <- res_pval$pathScore_GO[[pval_id]]
      id_pval <- 15
      id_info <- 1
    }else{
      stop('unknown pathway called')
    }
  }
}

# remove features that are not concordant across cohorts (meta analysis random model)
if(any(grepl('_model', colnames(res_pval)))){
  res_pval <- res_pval[res_pval[, grepl('_model', colnames(res_pval))] == 'fixed',]
}

# recompute pvalue if ngenes_tscore > 1
if(min_genes_path > 1 & grepl('path',type_data)){
  res_pval <- res_pval[res_pval$ngenes_tscore >= min_genes_path, ]
  res_pval[,id_pval+1] <- qvalue(res_pval[,id_pval])$qvalues
  res_pval[,id_pval+2] <- p.adjust(res_pval[,id_pval], method = 'BH')
}

if(exclude_MHC & type_data == 'tscore'){
  res_pval$start_position <- NA
  res_pval$chrom <- NA
  tmp <- read.table(geneRegionFile, h=T,stringsAsFactors = F)
  tmp <- tmp[match(res_pval$ensembl_gene_id, tmp$ensembl_gene_id),]
  res_pval$start_position <- tmp$start_position
  res_pval$chrom <- tmp$chrom
  HLA_reg <- c(26000000, 34000000)
  res_pval <- res_pval[!(res_pval$chrom %in% 'chr6' & res_pval$start_position <=HLA_reg[2] & res_pval$start_position >= HLA_reg[1]) , ]
}

if(!is.null(genes_to_filter)){
  genes_filt <- read.table(genes_to_filter, h=T, stringsAsFactors = F, sep = '\t')
  genes_filt <- genes_filt[genes_filt$keep & !is.na(genes_filt$keep),]
  res_pval <- res_pval[res_pval$ensembl_gene_id %in% genes_filt$ensembl_gene_id,]
}

### load score data ###
sampleAnn <- vector(mode = 'list', length = length(name_cohorts))
scoreMat <- vector(mode = 'list', length = length(name_cohorts))

for(c_id in 1:length(name_cohorts)){
  
  print(name_cohorts[c_id])
  
  sampleAnn[[c_id]] <- read.table(sampleAnnFile[[c_id]], h = T, stringsAsFactors = F)
  if(type_cluster == 'Cases'){
    sampleAnn[[c_id]] <- sampleAnn[[c_id]][sampleAnn[[c_id]]$Dx == 1,]
  }else{
    if(type_cluster == 'Controls'){
      sampleAnn[[c_id]] <- sampleAnn[[c_id]][sampleAnn[[c_id]]$Dx == 0,]
    }else{
      if(type_cluster != 'All')
        stop('type_cluster must be either Cases or Controls or All')
    }
  }
  
  sampleAnn[[c_id]]$Temp_ID <- sampleAnn[[c_id]]$Individual_ID
  sampleAnn[[c_id]]$cohort <- rep(name_cohorts[c_id], nrow(sampleAnn[[c_id]]))
  
  load_output <- load_input_matrix(inputFile[c_id], sampleAnn[[c_id]], res_pval, split_tot, id_info)
  scoreMat[[c_id]] <- load_output$scoreMat
  sampleAnn[[c_id]] <- load_output$sampleAnn
  
  # remove sample that have NAs
  id_s <- rowSums(is.na(scoreMat[[c_id]])) == 0
  if(!all(id_s)){scoreMat[[c_id]] <- scoreMat[[c_id]][id_s,]}
  sampleAnn[[c_id]] <- sampleAnn[[c_id]][match(rownames(scoreMat[[c_id]]), sampleAnn[[c_id]]$Individual_ID), ]
  
  print(dim(scoreMat[[c_id]]))
  print(dim(sampleAnn[[c_id]]))
  
}

try(if(any(!sapply(scoreMat, function(x) all(colnames(x) %in% res_pval[, id_info])))) 
  stop("ERROR: wrong filtering elements"))
res_pval <- res_pval[match(colnames(scoreMat[[1]]), res_pval[, id_info]),]
print(identical(colnames(scoreMat[[1]]), res_pval[,id_info]))

### clumping: sort according best gene association ###
scoreMat_tot <- do.call(rbind, scoreMat)
cor_score <- cor(scoreMat_tot, use = "pairwise.complete.obs")
element_rm <- clumping_features(res_pval=res_pval, 
                                id_info = id_info, 
                                corr_feat = cor_score, 
                                id_pval = id_pval, 
                                corr_thr = corr_thr)
print(paste(length(element_rm),'features removed due to high correlation'))

sampleAnn_list <- sampleAnn
sampleAnn <- do.call(rbind, sampleAnn)
scoreMat_tot <- scoreMat_tot[, !colnames(scoreMat_tot) %in% element_rm]
# match to have the same samples and same order with annotation
sampleAnn <- sampleAnn[match(rownames(scoreMat_tot), sampleAnn$Temp_ID), ]
sampleAnn$cohort_id <- as.numeric(as.factor(sampleAnn$cohort))

## correct for PCs ##
input_data_notcorr <- scale(scoreMat_tot)
attr(input_data_notcorr, "scaled:scale") <- NULL
attr(input_data_notcorr, "scaled:center") <- NULL

# remove PCs1-10 for each genes
input_data <- matrix(ncol = ncol(input_data_notcorr), nrow = nrow(input_data_notcorr))
rownames(input_data) <- rownames(input_data_notcorr)
colnames(input_data) <- colnames(input_data_notcorr)

name_cov <- setdiff(colnames(sampleAnn),
                    c('Individual_ID', 'genoSample_ID', 'Dx', 'Sex',
                      'Age', 'Temp_ID', 'cohort', 'cohort_id'))
fmla <- as.formula(paste('g ~', paste0(name_cov, collapse = '+')))

for(i in 1:ncol(input_data_notcorr)){
  # print(i)
  tmp <- data.frame(g = input_data_notcorr[,i], sampleAnn[, name_cov])
  reg <- lm(fmla, data = tmp)
  input_data[,i] <- reg$residuals
}
print("corrected for PCs")
res_pval <- res_pval[match(colnames(input_data), res_pval[, id_info]),]

if(type_input == 'zscaled'){
  input_data <- sapply(1:ncol(input_data), function(x) 
    input_data[, x]*res_pval[res_pval[,id_info] == colnames(input_data)[x], id_pval-1])
  colnames(input_data) <- res_pval[, id_info]
}

### plot: UMAP
# remove samples that are outliers (from UMAP)
n_comp_umap <- 2
n_neigh_umap <- 30
min_dist_umap <- 0.01
seed_umap <- 67

custom.settings = umap.defaults
custom.settings$min_dist = min_dist_umap
custom.settings$n_components = n_comp_umap
custom.settings$n_neighbors = n_neigh_umap
custom.settings$random_state <- seed_umap

umap_tot <- umap::umap(input_data, custom.settings)
id_out <- unique(unlist(apply(umap_tot$layout, 2, function(x) which(abs(x - median(x)) > (6 * sd(x))))))
sampleAnn_out <- NULL

df_umap_tot <- data.frame(component_1=umap_tot$layout[,1], component_2=umap_tot$layout[,2], 
                          outlier = rep('no', nrow(umap_tot$layout)), cohort = sampleAnn$cohort)
if(length(id_out)>0){
  
  print(sprintf('remove outliers (%s)', length(id_out)))
  
  df_umap_tot$outlier[id_out] <- 'yes'
  df_umap_tot$outlier <- factor(df_umap_tot$outlier, levels = c('yes', 'no'))
  
  # plot
  tot_pl <- ggplot(df_umap_tot, aes(x = component_1, y = component_2, color = outlier))+
    geom_point(size = 0.05)+
    theme_bw()+theme(legend.position = 'right')
  width_pl <- 4
  ggsave(filename = sprintf('%s%s_corrPCs_%s_cluster%s_PGmethod_%smetric_umap_oultiers.png', outFold, type_data, type_input, type_cluster, type_sim), width = width_pl, height = 4, plot = tot_pl, device = 'png')
  ggsave(filename = sprintf('%s%s_corrPCs_%s_cluster%s_PGmethod_%smetric_umap_outliers.pdf', outFold, type_data, type_input, type_cluster, type_sim), width = width_pl, height = 4, plot = tot_pl, device = 'pdf')
  
  sampleAnn_out <- sampleAnn[id_out, ]
  write.table(x = sampleAnn_out, 
              file = sprintf('%s%s_corrPCs_%s_cluster%s_PGmethod_umap_oultiers.txt', outFold, type_data, type_input, type_cluster), 
              col.names = T, row.names = F, sep = '\t', quote = F)
  
}else{
  print('no samples to be removed')
}

# save umap output
write.table(df_umap_tot, 
            file = sprintf('%s%s_corrPCs_%s_cluster%s_PGmethod_umap.txt', 
                           outFold, type_data, type_input, type_cluster), 
            col.names = T, row.names = F, sep = '\t', quote = F)



