# cluster using all Tissues
# merge using SNF

options(stringsAsFactors=F)
options(max.print=1000)
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(qvalue))
suppressPackageStartupMessages(library(igraph))
suppressPackageStartupMessages(library(pryr))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(SNFtool))
suppressPackageStartupMessages(library(coin))
suppressPackageStartupMessages(library(uwot))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(RColorBrewer))
options(bitmapType = 'cairo', device = 'png')


parser <- ArgumentParser(description="clustering using SNF + louvain cluster")
parser$add_argument("--tissues_name", type = "character", nargs = '*', help = "tissues")
parser$add_argument("--color_file", type = "character", help = "file with color based on phenotype")
parser$add_argument("--covDatFile", type = "character", default = 'NA', help = "additional cov to test")
parser$add_argument("--inputFile", type = "character", nargs = '*', help = "file to be loaded (predicted tscore or pathScore)")
parser$add_argument("--sampleAnnFile", type = "character", help = "file with samples to be used")
parser$add_argument("--type_cluster", type = "character", default = 'All', help = "All, Cases, Controls")
parser$add_argument("--split_tot", type = "integer", default = 0, help = "if 0 then inpuntFile load alone, otherwise splitted version")
parser$add_argument("--pvalresFile", type = "character", nargs = '*', help = "file with pvalue results")
parser$add_argument("--pval_id", type = "integer", default = 0, help = "id to be used on pvalue file")
parser$add_argument("--corr_thr", type = "double", default = -1, help = "correlation among features threshold")
parser$add_argument("--functR", type = "character", help = "functions to be used")
parser$add_argument("--type_data", type = "character", help = "tscore, path_Reactome or path_GO")
parser$add_argument("--type_input", type = "character", default = 'original', help = "original or zscaled")
parser$add_argument("--min_genes_path", type = "integer", default = 1, help = "minimum number of genes for a pathway, if > 1 recompute corrected pvalues")
parser$add_argument("--kNN_par", type = "integer", nargs = '*',  default = 30, help = "nearest neighbour used to compute sim")
parser$add_argument("--alphapar_SNF", type = "double", default = 0.5, help = "sd HK scaling parameter")
parser$add_argument("--tpar_SNF", type = "integer", default = 20, help = "iteration parameters SNF")
parser$add_argument("--outFold", type="character", help = "Output file [basename only]")

args <- parser$parse_args()
tissues_name <- args$tissues_name
pvalresFile <- args$pvalresFile
color_file <- args$color_file
covDatFile <- args$covDatFile
pval_id <- args$pval_id
inputFile <- args$inputFile
sampleAnnFile <- args$sampleAnnFile
type_cluster <- args$type_cluster
split_tot <- args$split_tot
functR <- args$functR
type_data <- args$type_data
corr_thr <- args$corr_thr
min_genes_path <- args$min_genes_path
type_input <- args$type_input
kNN_par <- args$kNN_par
alphapar_SNF <- args$alphapar_SNF
tpar_SNF <- args$tpar_SNF
outFold <- args$outFold


# ####################################################################################################################
# tissues_name <- c('DLPC_CMC', 'Whole_Blood')
# inputFile <- c('/psycl/g/mpsziller/lucia/UKBB/eQTL_PROJECT/OUTPUT_CMC/predict_UKBB/200kb/devgeno0.01_testdevgeno0/MDDRecur_pheno/predictedTscores_splitGenes',
# '/psycl/g/mpsziller/lucia/UKBB/eQTL_PROJECT/OUTPUT_GTEx/predict_UKBB/Whole_Blood/200kb/noGWAS/devgeno0.01_testdevgeno0/MDDRecur_pheno/predictedTscores_splitGenes')
# sampleAnnFile <- '/psycl/g/mpsziller/lucia/UKBB/eQTL_PROJECT/INPUT_DATA/Covariates/MDD_pheno_def/covariateMatrix_MDDRecur.txt'
# covDatFile <- '/psycl/g/mpsziller/lucia/UKBB/eQTL_PROJECT/INPUT_DATA/Covariates/covariatesMatrix_batchInfo.txt'
# split_tot <- 0
# pvalresFile <- c('/psycl/g/mpsziller/lucia/UKBB/eQTL_PROJECT/OUTPUT_CMC/predict_UKBB/200kb/devgeno0.01_testdevgeno0/MDDRecur_pheno/pval_MDDRecur_pheno_covCorr.RData',
# '/psycl/g/mpsziller/lucia/UKBB/eQTL_PROJECT/OUTPUT_GTEx/predict_UKBB/Whole_Blood/200kb/noGWAS/devgeno0.01_testdevgeno0/MDDRecur_pheno/pval_MDDRecur_pheno_covCorr.RData')
# pval_id <- 1
# min_genes_path <- 2
# type_data <- 'tscore'
# type_cluster <- 'Cases'
# outFold <- '/psycl/g/mpsziller/lucia/UKBB/eQTL_PROJECT/OUTPUT_all/MDDRecur_pheno/'
# functR <- '/psycl/g/mpsziller/lucia/priler_project/Software/model_clustering/clustering_functions.R'
# corr_thr <- 0.9
# type_input <- 'zscaled'
# color_file <- '/psycl/g/mpsziller/lucia/priler_project/Figures/color_tissues.txt'
# ####################################################################################################################

source(functR)

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

if(covDatFile != 'NA'){
  covDat <- read.table(covDatFile, h=T, stringsAsFactors = F, sep = '\t')
  covDat <- covDat[match(sampleAnn$Individual_ID, covDat$Individual_ID), ]
  sampleAnn$Batch <- covDat$Batch
  sampleAnn$Array <- covDat$Array
  sampleAnn$initial_assessment_centre <- covDat$initial_assessment_centre
}  

res_pval <- list()
scoreMat <- list()
for(i in 1:length(tissues_name)){
  
  print(tissues_name[i])
  
  # load pval res
  res_pval[[i]] <- get(load(pvalresFile[[i]]))
  if(type_data == 'tscore'){
    res_pval[[i]] <- res_pval[[i]]$tscore[[pval_id]]
    id_pval <- 8
    id_info <- 2
  }else{
    if(type_data == 'path_Reactome'){
      res_pval[[i]] <- res_pval[[i]]$pathScore_reactome[[pval_id]]
      id_pval <- 13
      id_info <- 1
    }else{
      if(type_data == 'path_GO'){
        res_pval[[i]] <- res_pval[[i]]$pathScore_GO[[pval_id]]
        id_pval <- 15
        id_info <- 1
      }else{
        stop('unknown pathway called')
      }
    }
  }
  
  # recompute pvalue if ngenes_tscore > 1
  if(min_genes_path > 1 & grepl('path',type_data)){
    res_pval[[i]] <- res_pval[[i]][res_pval[[i]]$ngenes_tscore >= min_genes_path, ]
    res_pval[[i]][,id_pval+1] <- qvalue(res_pval[[i]][,id_pval])$qvalues
    res_pval[[i]][,id_pval+2] <- p.adjust(res_pval[[i]][,id_pval], method = 'BH')
  }

  # load input matrix 
  if(split_tot == 0){
    
    scoreMat[[i]] <- get(load(inputFile[[i]]))
    # filter out based on samples and ids
    scoreMat[[i]] <- scoreMat[[i]][match(res_pval[[i]][, id_info],scoreMat[[i]][,1]), ]
    common_samples <- intersect(sampleAnn$Individual_ID, colnames(scoreMat[[i]]))
    
    scoreMat[[i]] <- t(scoreMat[[i]][,match(common_samples,colnames(scoreMat[[i]]))])
    rownames(scoreMat[[i]]) <- common_samples
    colnames(scoreMat[[i]]) <- res_pval[[i]][, id_info]
    
  }else{
    
    ###### load score Mat #######
    scoreMat_list <- vector(mode = 'list', length = split_tot)
    samplesID <- vector(mode = 'list', length = split_tot)
    elementID <- NULL
    
    for(j in 1:split_tot){
      
      print(j)
      if(file.exists(sprintf('%s%i.RData', inputFile[[i]], j))){
        tmp <- get(load(sprintf('%s%i.RData', inputFile[[i]], j)))
        elementID <- c(elementID,tmp[,1])
        samplesID[[j]] <- intersect(sampleAnn$Individual_ID, colnames(tmp))
        scoreMat_list[[j]] <- t(tmp[,match(samplesID[[j]],colnames(tmp))])
      }else{
        print(sprintf('split %i for tissue %s does not exist', j, tissues_name[i]))
      }  
    }
    
    # check samplesID always the same
    if(!all(table(unlist(samplesID)) == split_tot)){print('ERROR: wrong name annotations')}
    
    scoreMat[[i]] <- do.call(cbind, scoreMat_list)
    colnames(scoreMat[[i]]) <- elementID
    rm(scoreMat_list)
    
    # filter out elements that are repeated twice:
    id_dup <- names(which(table(colnames(scoreMat[[i]])) > 1)) 
    scoreMat[[i]] <- scoreMat[[i]][, !colnames(scoreMat[[i]]) %in% id_dup]
    
    id_el <- intersect(colnames(scoreMat[[i]]), res_pval[[i]][,id_info])
    scoreMat[[i]] <- scoreMat[[i]][, match(id_el, colnames(scoreMat[[i]]))]
    res_pval[[i]] <- res_pval[[i]][match(id_el, res_pval[[i]][,id_info]), ]
    
    rownames(scoreMat[[i]]) <- samplesID[[1]]
    # remove sample that have NAs
    id_s <- rowSums(is.na(scoreMat[[i]])) == 0
    if(!all(id_s)){scoreMat[[i]] <- scoreMat[[i]][id_s, ]}
    
    common_samples <- intersect(sampleAnn$Individual_ID, rownames(scoreMat[[i]]))
    scoreMat[[i]] <- scoreMat[[i]][match(common_samples,rownames(scoreMat[[i]])),]
    
  }
  
  print(identical(colnames(scoreMat[[i]]), res_pval[[i]][, id_info]))
  
  # remove higly correlated features, keep highest association
  cor_score <- cor(scoreMat[[i]])
  element_rm <- c()
  for(j in 1:(nrow(cor_score)-1)){
    id <-  which(abs(cor_score[j:nrow(cor_score),j])>corr_thr)
    if(length(id)>1){
      element_rm <- c(element_rm, names(id)[-which.min(res_pval[[i]][match(names(id), res_pval[[i]][,id_info]),id_pval])])
    }
  }
  element_rm <- unique(element_rm)
  print(paste(length(element_rm),'features removed due to high correlation'))
  scoreMat[[i]] <- scoreMat[[i]][,!colnames(scoreMat[[i]]) %in% element_rm]
  res_pval[[i]] <- res_pval[[i]][match(colnames(scoreMat[[i]]), res_pval[[i]][, id_info]),]
  
}

res_pval_tot <- cbind(do.call(rbind, res_pval), data.frame(tissue = unlist(lapply(1:length(tissues_name), function(x) rep(tissues_name[x], nrow(res_pval[[x]]))))))
res_pval_tot$new_id <- paste(res_pval_tot[, id_info], 'tissue', res_pval_tot$tissue, sep = '_')
  
# filter out genes with the same name, if > corr_thr remove
id_rep <- res_pval_tot[duplicated(res_pval_tot[, id_info]), id_info]
element_rm <- c()
for(i in 1:length(id_rep)){
  # print(i)
  tmp <- lapply(scoreMat, function(x) x[, colnames(x) %in% id_rep[i]])
  tmp <- do.call(cbind, tmp)
  colnames(tmp) <- res_pval_tot$tissue[res_pval_tot[,id_info] %in% id_rep[i]]
  cor_gene <- cor(tmp)
  tmp_pval <- res_pval_tot[res_pval_tot[,id_info] %in% id_rep[i],]
  if(any(cor_gene[lower.tri(cor_gene)]>corr_thr)){
    #if(length(which(cor_gene[lower.tri(cor_gene)] > corr_thr)) >=2){
      #print(cor_gene)
    #}   
    for(j in 1:(nrow(cor_gene)-1)){
      id <-  which(abs(cor_gene[j:nrow(cor_gene),j])>corr_thr)
      new <- tmp_pval[tmp_pval$tissue %in% names(id), ]
      element_rm <- c(element_rm, new$new_id[-which.min(new[,id_pval])])    
    }
  }
}
element_rm <- unique(element_rm)
res_pval_tot <- res_pval_tot[!res_pval_tot$new_id %in% element_rm, ]
for(i in 1:length(tissues_name)){
  tmp <- res_pval_tot[res_pval_tot$tissue %in% tissues_name[i], ]
  res_pval[[i]] <- tmp
  scoreMat[[i]] <- scoreMat[[i]][,match(tmp[,id_info], colnames(scoreMat[[i]]))]
}

# remove features highly correlated across tissues ### computationally expensive ####
# for(i in 1:(length(tissues_name)-1)){
#   print(tissues_name[i])
#   tmp <- list()
#   for(j in (i+1):length(tissues_name)){
#     tmp[[j-i]] <- apply(scoreMat[[i]], 2, function(x) apply(scoreMat[[j]], 2, function(y) cor(x, y)))
#     rownames(tmp[[i-j]])
#   }
#   tmp <- do.call(rbind, tmp)
# }

input_data <- list()
for(i in 1:length(tissues_name)){
  
  print(tissues_name[i])
  input_data[[i]] <- scale(scoreMat[[i]])
  attr(input_data[[i]], "scaled:scale") <- NULL
  attr(input_data[[i]], "scaled:center") <- NULL

  if(type_input == 'zscaled'){
    input_data[[i]] <- sapply(1:ncol(input_data[[i]]), function(x) input_data[[i]][, x]*res_pval[[i]][res_pval[[i]][,id_info] == colnames(input_data[[i]])[x], id_pval-1])
    colnames(input_data[[i]]) <- colnames(scoreMat[[i]])
  }
}

rm(scoreMat)
print(mem_used())

###############################################################
## pheno graph clustering using SNF similarity matrix
# NOTE: for cycle can be parallelized
dist_scoreMat <-  lapply(input_data, function(x) (dist2(as.matrix(x),as.matrix(x)))^(1/2))

PG_cl <- vector(mode = 'list', length = length(kNN_par))
test_cov <- vector(mode = 'list', length = length(kNN_par))
SNF_matrix <- vector(mode = 'list', length = length(kNN_par))

for(i in 1:length(kNN_par)){
  
  print(i)
  
  W_single <- lapply(dist_scoreMat, function(x) affinityMatrix(x, kNN_par[i], alphapar_SNF))
  # overall similarity matrix
  system.time(W <- SNF(W_single, kNN_par[i], tpar_SNF))
  PG_cl[[i]] <- clust_PGmethod_sim(kNN =  kNN_par[i], sim = W, type_Dx = type_cluster, sample_info = sampleAnn)
  print(PG_cl[[i]]$info)
  
  # cluster depend on PC?
  id <- PG_cl[[i]]$cl$membership
  df <- cbind(data.frame(cl = id), sampleAnn[,! colnames(sampleAnn) %in% c('Individual_ID', 'Dx')])
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
  test_cov[[i]]$pval_BHcorr <- p.adjust(test_cov[[i]]$pval, method = 'BH')
  test_cov[[i]]$kNN <- kNN_par[i]
  print(test_cov[[i]])
  
  SNF_matrix[[i]] <- W
  
}

test_cov <- do.call(rbind, test_cov)
info_hyperParam <- do.call(rbind, lapply(PG_cl, function(x) x$info))
opt_k <- kNN_par[which.max(info_hyperParam$DB_mean)]

# if type_clster == 'All' compute percentage for each group
df_perc <- df_perc_test <- list()
if(type_cluster == 'All'){
  
  for(i in 1:length(kNN_par)){
    perc <- table(PG_cl[[i]]$cl$membership, sampleAnn$Dx)/rowSums(table(PG_cl[[i]]$cl$membership, sampleAnn$Dx))
    # test fisher for each group
    cl_id <- sort(unique(PG_cl[[i]]$cl$membership))
    df_perc[[i]] <- data.frame(gr = rep(cl_id,2), Dx = c(rep(0,nrow(perc)), rep(1,nrow(perc))), 
                               perc = as.vector(perc), count = as.vector(table(PG_cl[[i]]$cl$membership, sampleAnn$Dx)))
    df_perc_test[[i]] <- data.frame(gr = c(cl_id, 'all'), 
                                    fisher_test = c(sapply(cl_id, function(x) fisher.test(table(PG_cl[[i]]$cl$membership == x, sampleAnn$Dx))$p.value), 
                                                    chisq.test(table(PG_cl[[i]]$cl$membership, sampleAnn$Dx))$p.value))
  }
}

output <- list(best_k = opt_k, cl_res = PG_cl, test_cov = test_cov, info_tune = info_hyperParam, 
               cl_best = data.frame(id = sampleAnn$Individual_ID, gr = PG_cl[[which.max(info_hyperParam$DB_mean)]]$cl$membership), 
               feat = lapply(input_data, colnames),  res_pval = res_pval_tot, tissues = tissues_name)
output$Dx_perc <- list(perc = df_perc, test = df_perc_test)
output$samples_id <- sampleAnn$Individual_ID
output$SNF_matrix <- SNF_matrix

# most significant element (for each tissue)
test_diff <- data.frame(tissue = unlist(mapply(function(x,y) rep(x, ncol(y)), x = tissues_name, y = input_data, SIMPLIFY = F)), 
                        id = unlist(lapply(input_data, colnames)), pval = NA)
for(i in 1:length(tissues_name)){
  test_diff$pval[test_diff$tissue == tissues_name[i]] <- apply(input_data[[i]], 2, function(x) kruskal.test(x = x, g = factor(output$cl_best$gr))$p.value)
}
test_diff$pval_corr <- p.adjust(test_diff$pval, method = 'BH')
output$test_diff_gr <- test_diff
# compute mean for each gr
df_gr_mean <- matrix(ncol = length(unique(output$cl_best$gr)), nrow = nrow(test_diff))
df_gr_sd <- matrix(ncol = length(unique(output$cl_best$gr)), nrow = nrow(test_diff))

for(i in 1:length(tissues_name)){
  df_gr_mean[test_diff$tissue == tissues_name[i],] <- t(apply(input_data[[i]], 2, function(x) 
    sapply(sort(unique(output$cl_best$gr)), function(y) mean(x[output$cl_best$gr == y]) )))
  
  df_gr_sd[test_diff$tissue == tissues_name[i],] <- t(apply(input_data[[i]], 2, function(x) 
    sapply(sort(unique(output$cl_best$gr)), function(y) sd(x[output$cl_best$gr == y]) )))
}
df_gr_cv <- df_gr_mean/df_gr_sd
colnames(df_gr_cv) <- colnames(df_gr_mean) <- colnames(df_gr_sd) <- paste0('gr_', sort(unique(output$cl_best$gr)))
df_gr_cv <- as.data.frame(df_gr_cv)
df_gr_mean <- as.data.frame(df_gr_mean)
df_gr_sd <- as.data.frame(df_gr_sd)
df_gr_mean$id <- df_gr_sd$id <- df_gr_cv$id <- unlist(lapply(input_data, colnames))
df_gr_mean$tissue <- df_gr_sd$tissue <- df_gr_cv$tissue <- unlist(mapply(function(x, y) rep(x, ncol(y)), x = tissues_name, y = input_data, SIMPLIFY = F))

output$gr_input <- list(mean = df_gr_mean, sd = df_gr_sd, cv = df_gr_cv)

# save results:
save(output, file = sprintf('%s%s_%s_cluster%s_PGmethod_SNFmetric.RData', outFold, type_data, type_input, type_cluster))


### plot: UMAP
n_comp_umap <- 2
n_neigh_umap <- opt_k
min_dist_umap <- 0.01
seed_umap <- 67
dist_mat <- 0.5 - output$SNF_matrix[[which(kNN_par == opt_k)]]
dist_mat <- as.dist(dist_mat)

set.seed(seed_umap)
umap_res <- uwot::umap(dist_mat, n_neighbors = n_neigh_umap, n_components = n_comp_umap, min_dist = min_dist_umap)

df <- data.frame(component_1=umap_res[,1], component_2=umap_res[,2], gr = output$cl_best$gr)
df$gr <- factor(df$gr)

tot_pl <- ggplot(df, aes(x = component_1, y = component_2, color = gr))+
  geom_point(size = 0.05)+
  theme_bw()+theme(legend.position = 'right')
width_pl <- 4
if(type_cluster == 'All'){
  df$Dx <- factor(sampleAnn$Dx)
  pl_extra <- ggplot(df, aes(x = component_1, y=component_2, color = Dx))+
    geom_point(size = 0.05)+
    theme_bw()+theme(legend.position = 'right')
  tot_pl <- ggarrange(plotlist = list(tot_pl, pl_extra), ncol = 2, nrow = 1, align='h')
  width_pl <- 8
}
ggsave(filename = sprintf('%s%s_%s_cluster%s_PGmethod_SNFmetric_umap.png', outFold, type_data, type_input, type_cluster), width = width_pl, height = 4, plot = tot_pl, device = 'png')
ggsave(filename = sprintf('%s%s_%s_cluster%s_PGmethod_SNFmetric_umap.pdf', outFold, type_data, type_input, type_cluster), width = width_pl, height = 4, plot = tot_pl, device = 'pdf')

df$Gender <- factor(sampleAnn$Gender)
df$Batch <- factor(sampleAnn$Batch)
df$Array <- factor(sampleAnn$Array)
df$Centre <- factor(sampleAnn$initial_assessment_centre)
df$Age <- sampleAnn$Age

myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
sc <- scale_colour_gradientn(colours = myPalette(100), limits=c(min(sampleAnn$Age), max(sampleAnn$Age)))

pl1 <- ggplot(df, aes(x = component_1, y=component_2, color = Gender))+
  geom_point(size = 0.03)+ggtitle('Gender')+
  theme_bw()+theme(legend.position = 'none')
pl2 <- ggplot(df, aes(x = component_1, y=component_2, color = Batch))+
  geom_point(size = 0.03)+ggtitle('Batch')+
  theme_bw()+theme(legend.position = 'none')
pl3 <- ggplot(df, aes(x = component_1, y=component_2, color = Array))+
  geom_point(size = 0.03)+ggtitle('Array')+
  theme_bw()+theme(legend.position = 'none')
pl4 <- ggplot(df, aes(x = component_1, y=component_2, color = Centre))+
  geom_point(size = 0.03)+ggtitle('Centre')+
  theme_bw()+theme(legend.position = 'none')
pl5 <- ggplot(df, aes(x = component_1, y=component_2, color = Age))+
  geom_point(size = 0.03)+ggtitle('Age')+sc+
  theme_bw()+theme(legend.position = 'right')
tot_pl <- ggarrange(plotlist = list(pl1, pl2, pl3, pl4, pl5), ncol = 3, nrow = 2)
ggsave(filename = sprintf('%s%s_%s_cluster%s_PGmethod_SNFmetric_umap_cov.png', outFold, type_data, type_input, type_cluster), width = 12, height = 8, plot = tot_pl, device = 'png')
ggsave(filename = sprintf('%s%s_%s_cluster%s_PGmethod_SNFmetric_umap_cov.pdf', outFold, type_data, type_input, type_cluster), width = 12, height = 8, plot = tot_pl, device = 'pdf')


### plot: heatmap
file_name <- sprintf('%s%s_%s_cluster%s_PGmethod_SNFmetric', outFold, type_data, type_input, type_cluster)
tmp <- test_diff[test_diff$pval_corr < 0.05,]
if(nrow(tmp)>50){
  keep_feat <- test_diff[order(test_diff$pval)[1:50],]
}else{
  keep_feat <- tmp
}

for(i in 1:length(tissues_name)){
  
  pheat_pl(mat = input_data[[i]][,keep_feat$id[keep_feat$tissue == tissues_name[i]]], type_mat = type_data, cl = output$cl_best, height_pl = 7, width_pl = 5, 
           outFile = paste0(file_name, '_heatmap', tissues_name[i]))
  if(type_input == 'zscaled'){
    pheat_pl(mat = scale(input_data[[i]][,keep_feat$id[keep_feat$tissue == tissues_name[i]]]), type_mat = type_data, cl = output$cl_best, height_pl = 7, width_pl = 5, 
             outFile = paste0(file_name, '_heatmap_scaled', tissues_name[i]))
  }
}

color_tissues <- read.table(color_file, h=T, stringsAsFactors = F)
color_tissues <- color_tissues[match(tissues_name, color_tissues$tissues),]

mat <- output$gr_input$cv[paste(output$gr_input$cv$tissue, output$gr_input$cv$id, sep = '_') %in% paste(keep_feat$tissue, keep_feat$id, sep = '_'), ]
width_pl <- 7
if(type_data == 'path_GO'){
  for(i in 1:length(tissues_name)){
    mat$id[mat$tissue == tissues_name[i]] <- res_pval[[i]]$path[match(mat$id[mat$tissue == tissues_name[i]], res_pval[[i]][, id_info])]
  }
  width_pl <- 9
  print(str(mat))
}
if(type_data == 'path_Reactome'){
  width_pl <- 9
}

# remove duplicated
mat <- mat[!duplicated(mat$id),]

pheat_pl_gr(mat, type_mat = type_data, height_pl = 7, width_pl = width_pl, color_df = color_tissues, outFile = paste0(file_name, '_heatmap_gr'))




