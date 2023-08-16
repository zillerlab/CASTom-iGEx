#!/usr/bin/env Rscript

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
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(bigmemory))
suppressPackageStartupMessages(library(ggsci))
suppressPackageStartupMessages(library(clustAnalytics))
options(bitmapType = 'cairo', device = 'png')


parser <- ArgumentParser(description="clustering multiple cohorts concatenated using PG method")
parser$add_argument("--inputFile", type = "character",  nargs = '*', help = "file to be loaded (predicted tscore or pathScore)")
parser$add_argument("--cluster_method", type = "character", default = "leiden", help = "leidein or louvain, community detection method")
parser$add_argument("--kNN_par", type = "integer", nargs = '*', default = 20, help = "parameter used for PG method")
parser$add_argument("--genes_to_filter", type = "character", default = NULL, help = "additional file to filter genes")
parser$add_argument("--name_cohorts", type = "character", nargs = '*', help = "name of the single cohorts")
parser$add_argument("--sampleAnnFile", type = "character", nargs = '*', help = "file with samples to be used")
parser$add_argument("--sampleOutFile", type = "character", default = NULL, help = "file with samples to be excluded")
parser$add_argument("--geneRegionFile", type = "character", default = NULL, help = "used if tscore and exclude_MHC")
parser$add_argument("--tissues_name", type = "character", help = "name tissue")
parser$add_argument("--exclude_MHC", type = "logical", default = F, help = "if true, MHC region excluded (only ossible for tscore)")
parser$add_argument("--type_cluster", type = "character", help = "All, Cases, Controls")
parser$add_argument("--split_tot", type = "integer", default = 0, help = "if 0 then inpuntFile load alone, otherwise splitted version")
parser$add_argument("--pvalresFile", type = "character",  help = "file with pvalue results")
parser$add_argument("--pval_id", type = "integer", default = 1, help = "id to be used on pvalue file")
parser$add_argument("--corr_thr", type = "double", default = -1, help = "correlation among features threshold")
parser$add_argument("--functR", type = "character", help = "functions to be used")
parser$add_argument("--type_data", type = "character", default = "tscore", help = "tscore, path_Reactome or path_GO")
parser$add_argument("--type_sim", type = "character", default = 'HK', help = "HK or ED")
parser$add_argument("--type_input", type = "character", default = 'original', help = "original or zscaled")
parser$add_argument("--min_genes_path", type = "integer", default = 1, help = "minimum number of genes for a pathway, if > 1 recompute corrected pvalues")
parser$add_argument("--outFold", type="character", help = "Output file [basename only]")

args <- parser$parse_args()
name_cohorts <- args$name_cohorts
cluster_method <- args$cluster_method
genes_to_filter <- args$genes_to_filter
pvalresFile <- args$pvalresFile
sampleOutFile <- args$sampleOutFile
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
kNN_par <- args$kNN_par
outFold <- args$outFold

####################################################################################################################
# name_cohorts <- read.table('/home/luciat/eQTL_PROJECT/INPUT_DATA/SCZ_cohort_names_CLUST', header = F, stringsAsFactors = F)$V1[c(1:2, 16)]
# geneRegionFile <- '/home/luciat/eQTL_PROJECT/OUTPUT_CMC/train_CMC/200kb/resPrior_regEval_allchr.txt'
# exclude_MHC <- T
# inputFile <- paste0('/home/luciat/eQTL_PROJECT/OUTPUT_CMC/predict_PGC/200kb/',name_cohorts, '/devgeno0.01_testdevgeno0/predictedTscores.txt')
# sampleAnnFile <- paste0('/home/luciat/eQTL_PROJECT/INPUT_DATA/Covariates/',name_cohorts,'.covariateMatrix_old.txt')
# covDatFile <- paste0('/home/luciat/eQTL_PROJECT/INPUT_DATA/Covariates/',name_cohorts,'.covariateMatrix_old.txt')
# split_tot <- 0
# pvalresFile <- '/home/luciat/eQTL_PROJECT/OUTPUT_all/Meta_Analysis_SCZ/DLPC_CMC/pval_Dx_pheno_covCorr.RData'
# pval_id <- 1
# min_genes_path <- 2
# type_data <- 'tscore'
# type_cluster <- 'Cases'
# type_sim <- 'HK'
# outFold <- '/home/luciat/eQTL_PROJECT/OUTPUT_CMC/predict_PGC/200kb/Meta_Analysis_SCZ/devgeno0.01_testdevgeno0/update_corrPCs/'
# functR <- '/home/luciat/castom-igex/Software/model_clustering/clustering_functions.R'
# corr_thr <- 0.9
# type_input <- 'zscaled'
# kNN_par <- 20
# color_file <- '/home/luciat/castom-igex/Figures/color_tissues.txt'
# tissues_name <- 'DLPC_CMC'
# genes_to_filter <- '/home/luciat/eQTL_PROJECT/compare_prediction_UKBB_SCZ-PGC/DLPC_CMC_filter_genes_matched_datasets.txt'
#####################################################################################################################

source(functR)
print(paste("kNN parameters:", kNN_par))

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

if(!is.null(sampleOutFile)){
  rm_samples <- read.table(sampleOutFile, header = T, stringsAsFactors = F, sep = '\t')
}else{
  rm_samples <- NULL
}

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
  # remove outliers
  if(!is.null(rm_samples)){
    scoreMat[[c_id]] <- scoreMat[[c_id]][!rownames(scoreMat[[c_id]]) %in%  rm_samples$Temp_ID, ]
  }
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
sampleAnn_list <- sampleAnn
sampleAnn <- do.call(rbind, sampleAnn)

if(corr_thr < 1){
cor_score <- cor(scoreMat_tot, use = "pairwise.complete.obs")
element_rm <- clumping_features(res_pval=res_pval, 
                                id_info = id_info, 
                                corr_feat = cor_score, 
                                id_pval = id_pval, 
                                corr_thr = corr_thr)
print(paste(length(element_rm),'features removed due to high correlation'))
scoreMat_tot <- scoreMat_tot[, !colnames(scoreMat_tot) %in% element_rm]
}else{
print('All features considered')
}
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

## compute umap and plot cohort ##
# plot: UMAP
n_comp_umap <- 2
n_neigh_umap <- 30
min_dist_umap <- 0.01
seed_umap <- 67

custom.settings = umap.defaults
custom.settings$min_dist = min_dist_umap
custom.settings$n_components = n_comp_umap
custom.settings$n_neighbors = n_neigh_umap
custom.settings$random_state <- seed_umap

umap_res <- umap::umap(input_data, custom.settings)
df_umap <- data.frame(component_1=umap_res$layout[,1], 
                      component_2=umap_res$layout[,2], 
                      cohort = factor(sampleAnn$cohort))

pl_umap <- ggplot(df_umap, 
                  aes(x = component_1, y = component_2, color = cohort))+
  geom_point(size = 0.05, alpha = 0.8)+
  xlab('UMAP component 1')+ ylab('UMAP component 2')+
  theme_bw() +
  theme(legend.position = 'right')

ggsave(filename = sprintf('%s%s_corrPCs_%s_cluster%s_PGmethod_%smetric_umap_cohort.png', 
                          outFold, type_data, type_input, type_cluster, type_sim), 
       width = 6, height = 4, plot = pl_umap, device = 'png', dpi = 200)

print('UMAP computed and plotted')

###############################################################
## pheno graph clustering:
# specific for cluster with low computational power (LISA PGC)
# 1) Euclidian distance
if(nrow(sampleAnn)>10000){
  ed_dist <- big.matrix(ncol = nrow(sampleAnn) , nrow = nrow(sampleAnn), 
                        type = "double", init = 0, dimnames = NULL, shared = F)
  # build by rows
  fun_eucl <- function(x,y){
    sum((x-y)^2)
  }
  for(i in 1:(nrow(input_data)-1)){
    print(i)
    tmp <- apply(input_data[(i+1):nrow(input_data),,drop = F], 1, function(x) fun_eucl(input_data[i,], x))
    ed_dist[i,(i+1):nrow(input_data)] <- tmp
    ed_dist[(i+1):nrow(input_data),i] <- tmp
  }
  rm(tmp)
}else{
  ed_dist <- as.matrix(dist(input_data, method = 'euclidian'))^2
}
print(mem_used())
print('euclidian distance computed')

if(type_sim == 'HK'){fun_cl <- clust_PGmethod_HKsim}
if(type_sim == 'ED'){fun_cl <- clust_PGmethod_EDdist}
if(!type_sim %in% c('ED', 'HK')){stop('type similarity must be ED or HK')}

# NOTE: for cycle can be parallelized
PG_cl <- vector(mode = 'list', length = length(kNN_par))
test_cov <- vector(mode = 'list', length = length(kNN_par))
for(i in 1:length(kNN_par)){
  
  PG_cl[[i]] <- fun_cl(kNN = kNN_par[i], score = input_data, 
                       type_Dx = type_cluster, 
                       sample_info=sampleAnn,
                       euclDist=ed_dist[,], 
                       multiple_cohorts = T, 
                       cluster_method = cluster_method)
  
  print(PG_cl[[i]]$info)
  # cluster depend on PC?
  id <- PG_cl[[i]]$cl$membership
  # test only for cohort and PCs (they are computed combining all cohorts)
  df <- cbind(data.frame(cl = id), 
              sampleAnn[,colnames(sampleAnn) %in% c('cohort', name_cov)])
  test_cov[[i]] <- data.frame(cov_id = colnames(df)[-(1)])
  test_cov[[i]]$test_type <- test_cov[[i]]$statistic <- test_cov[[i]]$pval <- NA
  for(j in 1:(ncol(df)-1)){
    
    if(is.integer(df[, j+1]) | is.character(df[, j+1])){
      tmp <- chisq.test(table(df$cl, df[, j+1]))
      test_cov[[i]]$pval[j] <- tmp$p.value
      test_cov[[i]]$statistic[j] <-tmp$statistic
      test_cov[[i]]$test_type[j] <- 'chisq'
    }else{
      tmp <- kruskal_test(df[, j+1] ~ factor(df$cl))
      test_cov[[i]]$pval[j] <- pvalue(tmp)
      test_cov[[i]]$statistic[j] <- tmp@statistic@teststatistic
      test_cov[[i]]$test_type[j] <- 'kruskal'
    }
  }
  # test_cov[[i]]$pval_BHcorr <- p.adjust(test_cov[[i]]$pval, method = 'BH')
  test_cov[[i]]$kNN <- kNN_par[i]
  print(test_cov[[i]])
  
}

test_cov <- do.call(rbind, test_cov)
info_hyperParam <- do.call(rbind, lapply(PG_cl, function(x) x$info))
opt_k <- kNN_par[which.max(info_hyperParam$coverage_and_conductance)]

# if type_clster == 'All' compute percentage for each group
df_perc <- df_perc_test <- list()
if(type_cluster == 'All'){
  
  for(i in 1:length(kNN_par)){
    perc <- table(PG_cl[[i]]$cl$membership, sampleAnn$Dx)/rowSums(table(PG_cl[[i]]$cl$membership, sampleAnn$Dx))
    # test fisher for each group
    cl_id <- sort(unique(PG_cl[[i]]$cl$membership))
    df_perc[[i]] <- data.frame(gr = rep(cl_id,2), Dx = c(rep(0,nrow(perc)), rep(1,nrow(perc))), 
                               perc = as.vector(perc), 
                               count = as.vector(table(PG_cl[[i]]$cl$membership, sampleAnn$Dx)))
    df_perc_test[[i]] <- data.frame(gr = c(cl_id, 'all'), 
                                    fisher_test = c(sapply(cl_id, function(x) fisher.test(table(PG_cl[[i]]$cl$membership == x, sampleAnn$Dx))$p.value), 
                                                    chisq.test(table(PG_cl[[i]]$cl$membership, sampleAnn$Dx))$p.value))
  }
}

output <- list(best_k = opt_k, cl_res = PG_cl, test_cov = test_cov, 
               info_tune = info_hyperParam, feat = colnames(input_data), res_pval = res_pval,
              cl_best = data.frame(id = sampleAnn$Individual_ID, 
               gr = PG_cl[[which.max(info_hyperParam$coverage_and_conductance)]]$cl$membership))
output$Dx_perc <- list(perc = df_perc, test = df_perc_test)
output$samples_id <- rownames(input_data)

# most significant elements
test_diff <- data.frame(id = colnames(input_data), pval = apply(input_data, 2, function(x) kruskal.test(x = x, g = factor(output$cl_best$gr))$p.value))
test_diff$pval_corr <- p.adjust(test_diff$pval, method = 'BH')
output$test_diff_gr <- test_diff
# compute mean for each gr
df_gr_mean <- matrix(ncol = length(unique(output$cl_best$gr)), nrow = nrow(test_diff))
df_gr_sd <- matrix(ncol = length(unique(output$cl_best$gr)), nrow = nrow(test_diff))

df_gr_mean[,] <- t(apply(input_data, 2, function(x) 
  sapply(sort(unique(output$cl_best$gr)), function(y) mean(x[output$cl_best$gr == y]) )))

df_gr_sd[,] <- t(apply(input_data, 2, function(x) 
  sapply(sort(unique(output$cl_best$gr)), function(y) sd(x[output$cl_best$gr == y]) )))

df_gr_cv <- df_gr_mean/df_gr_sd
colnames(df_gr_cv) <- colnames(df_gr_mean) <- colnames(df_gr_sd) <- paste0('gr_', sort(unique(output$cl_best$gr)))
df_gr_cv <- as.data.frame(df_gr_cv)
df_gr_mean <- as.data.frame(df_gr_mean)
df_gr_sd <- as.data.frame(df_gr_sd)
df_gr_mean$id <- df_gr_sd$id <- df_gr_cv$id <- colnames(input_data)

output$gr_input <- list(mean = df_gr_mean, sd = df_gr_sd, cv = df_gr_cv)
output$input_data <- input_data
output$sampleInfo <- sampleAnn

# save also euclidian distance
output$ed_dist <- ed_dist[,]
output$sampleOutliers <- list(sample = rm_samples)

# plot umap clustering
df_umap$gr <- factor(output$cl_best$gr)
P <- length(unique(df_umap$gr))
gr_color <- pal_d3(palette = 'category20')(P)

# save
output$umap <- df_umap
save(output, file = sprintf('%s%s_corrPCs_%s_cluster%s_PGmethod_%smetric.RData', 
                            outFold, type_data, type_input, type_cluster, type_sim))

# save reduced output:
output_red <- list(cl_best = output$cl_best,
                   samples_id = output$samples_id, 
                   feat = output$feat, 
                   test_cov = output$test_cov)

save(output_red, file = sprintf('%s%s_corrPCs_%s_cluster%s_PGmethod_%smetric_minimal.RData', 
                            outFold, type_data, type_input, type_cluster, type_sim))

# group
tot_pl <- ggplot(df_umap, aes(x = component_1, y = component_2, color = gr))+
  geom_point(size = 0.05, alpha = 0.8)+
  xlab('UMAP component 1')+ ylab('UMAP component 2')+
  scale_color_manual(values = gr_color)+
  theme_bw()+theme(legend.position = 'right')
width_pl <- 4.5

ggsave(filename = sprintf('%s%s_corrPCs_%s_cluster%s_PGmethod_%smetric_umap.png', 
                          outFold, type_data, type_input, type_cluster, type_sim), 
       dpi=200, width = width_pl, height = 4, plot = tot_pl, device = 'png')



