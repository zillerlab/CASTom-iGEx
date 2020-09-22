# predict cluster probability based on phenograp approach

options(stringsAsFactors=F)
options(max.print=1000)
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(e1071))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(qvalue))
suppressPackageStartupMessages(library(pROC))
suppressPackageStartupMessages(library(pryr))
suppressPackageStartupMessages(library(umap))
suppressPackageStartupMessages(library(igraph))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(SparseM))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(ggplot2))
options(bitmapType = 'cairo', device = 'png')

parser <- ArgumentParser(description="predict cluster probability for new samples")
parser$add_argument("--sampleAnnNew_file", type = "character", help = "")
parser$add_argument("--name_cohort", type = "character", help = "")
parser$add_argument("--type_cluster", type = "character",default = 'All', help = "All, Cases, Controls")
parser$add_argument("--type_data", type = "character", help = "pathway or tscore")
parser$add_argument("--inputFile", type = "character", help = "input files (scores)")
parser$add_argument("--clustFile", type = "character", help = "file cluster results")
parser$add_argument("--split_tot", type = "integer", default = 0, help = "if 0 then inpuntFile load alone, otherwise splitted version")
parser$add_argument("--functR", type = "character", help = "functions to be used")
parser$add_argument("--type_input", type = "character", default = 'original', help = "original or zscaled")
parser$add_argument("--tissues_name", type = "character", help = "name tissue")
parser$add_argument("--color_file", type = "character", help = "file with color based on phenotype")
parser$add_argument("--outFold", type="character", help = "Output file [basename only]")

args <- parser$parse_args()
inputFile <- args$inputFile
sampleAnnNew_file <- args$sampleAnnNew_file
clustFile <- args$clustFile
functR <- args$functR
type_data <- args$type_data
type_input <- args$type_input
split_tot <- args$split_tot
type_cluster <- args$type_cluster
tissues_name <- args$tissues_name
color_file <- args$color_file
name_cohort <- args$name_cohort
outFold <- args$outFold

###################################################################################################################
# sampleAnnNew_file <- '/home/luciat/eQTL_PROJECT/INPUT_DATA/Covariates/scz_boco_eur.covariateMatrix_old.txt'
# name_cohort <- 'scz_boco_eur'
# functR <- '/home/luciat/priler_project/Software/model_clustering/clustering_functions.R'
# type_data <- 'tscore'
# type_input <- 'zscaled'
# type_cluster <- 'Cases'
# clustFile <- '/home/luciat/eQTL_PROJECT/OUTPUT_CMC/predict_PGC/200kb/Meta_Analysis_SCZ/devgeno0.01_testdevgeno0/tscore_zscaled_clusterCases_PGmethod_HKmetric.RData'
# outFold <-  '/home/luciat/eQTL_PROJECT/OUTPUT_CMC/predict_PGC/200kb/scz_boco_eur/devgeno0.01_testdevgeno0/'
# tissues_name <- 'DLPC_CMC'
# color_file <-  '/home/luciat/priler_project/Figures/color_tissues.txt'
# split_tot = 0
# inputFile <- '/home/luciat/eQTL_PROJECT/OUTPUT_CMC/predict_PGC/200kb/scz_boco_eur/devgeno0.01_testdevgeno0/predictedTscores.txt'
###################################################################################################################

source(functR)

clust_res <- get(load(clustFile))
sampleAnn_new <- read.table(sampleAnnNew_file, h=T, stringsAsFactors = F)
sampleAnn_new$cohort <- name_cohort
sampleAnn <- clust_res$sampleInfo

if(type_cluster == 'Cases'){
  sampleAnn <- sampleAnn[sampleAnn$Dx == 1,]
  sampleAnn_new <- sampleAnn_new[sampleAnn_new$Dx == 1,]
}else{
  if(type_cluster == 'Controls'){
    sampleAnn <- sampleAnn[sampleAnn$Dx == 0,]
    sampleAnn_new <- sampleAnn_new[sampleAnn_new$Dx == 0,]
  }else{
    if(type_cluster != 'All')
      stop('type_cluster must be either Cases or Controls or All')
  }
}

# load model
cl <- clust_res$cl_best
sampleAnn_tot <- rbind(sampleAnn[, intersect(colnames(sampleAnn), colnames(sampleAnn_new))], sampleAnn_new[, intersect(colnames(sampleAnn), colnames(sampleAnn_new))])

res_pval <- clust_res$res_pval
if(type_data == 'tscore'){
  id_pval <- 8
  id_info <- 2
}else{
  if(type_data == 'path_Reactome'){
    id_pval <- 13
    id_info <- 1
  }else{
    if(type_data == 'path_GO'){
      id_pval <- 15
      id_info <- 1
    }else{
      stop('unknown pathway called')
    }
  }
}

# load input matrix (new data, old data input saved in clustering result)
if(split_tot == 0){
  
  if(grepl('.txt', inputFile, fixed = TRUE)){
    scoreMat <- read.delim(inputFile, h=T, stringsAsFactors = F, check.names = F)
    scoreMat <- scoreMat[match(res_pval[, id_info],scoreMat[,1]), ]
    rownames(scoreMat) <- scoreMat[,1]
    scoreMat <- scoreMat[,-1]
    if(type_data == 'tscore'){
      new_id <- unname(sapply(colnames(scoreMat), function(x) strsplit(x, split = ' vs reference')[[1]][1]))
      colnames(scoreMat) <- new_id  
    }
  }else{
    scoreMat <- get(load(inputFile))  
    scoreMat <- scoreMat[match(res_pval[, id_info],scoreMat[,1]), ]
    rownames(scoreMat) <- scoreMat[,1]
    scoreMat <- scoreMat[,-1]
  }
  
  common_samples <- intersect(sampleAnn_new$Individual_ID, colnames(scoreMat))
  sampleAnn_new <- sampleAnn_new[match(common_samples, sampleAnn_new$Individual_ID),]
  scoreMat <- t(scoreMat[,match(common_samples,colnames(scoreMat))])
  
}else{
  
  ###### load score Mat #######
  scoreMat_list <- vector(mode = 'list', length = split_tot)
  samplesID <- vector(mode = 'list', length = split_tot)
  elementID <- NULL
  
  for(i in 1:split_tot){
    
    print(i)
    tmp <- get(load(sprintf('%s%i.RData', inputFile, i)))
    elementID <- c(elementID,tmp[,1])
    samplesID[[i]] <- intersect(sampleAnn_new$Individual_ID, colnames(tmp))
    scoreMat_list[[i]] <- t(tmp[,match(samplesID[[i]],colnames(tmp))])
    
  }
  
  # check samplesID always the same
  if(!all(table(unlist(samplesID)) == split_tot)){print('ERROR: wrong name annotations')}
  
  scoreMat <- do.call(cbind, scoreMat_list)
  colnames(scoreMat) <- elementID
  rm(scoreMat_list)
  
  # filter out elements that are repeated twice:
  id_dup <- names(which(table(colnames(scoreMat)) > 1)) 
  scoreMat <- scoreMat[, !colnames(scoreMat) %in% id_dup]
  
  scoreMat <- scoreMat[, match(id_el, colnames(scoreMat))]
  
  rownames(scoreMat) <- samplesID[[1]]
  # remove sample that have NAs
  id_s <- rowSums(is.na(scoreMat)) == 0
  if(!all(id_s)){scoreMat <- scoreMat[id_s, ]}
  
  common_samples <- intersect(sampleAnn_new$Individual_ID, rownames(scoreMat))
  sampleAnn_new <- sampleAnn_new[match(common_samples, sampleAnn_new$Individual_ID),]
  scoreMat <- scoreMat[match(common_samples,rownames(scoreMat)),]
  
}
print(identical(colnames(scoreMat), res_pval[, id_info]))

### scale data ###
input_data <- scale(scoreMat)
attr(input_data, "scaled:scale") <- NULL
attr(input_data, "scaled:center") <- NULL

if(type_input == 'zscaled'){
  input_data <- sapply(1:ncol(input_data), function(x) input_data[, x]*res_pval[res_pval[,id_info] == colnames(input_data)[x], id_pval-1])
  colnames(input_data) <- colnames(scoreMat)
}

data_tot <- rbind(clust_res$input_data, input_data)
print(identical(rownames(data_tot), sampleAnn_tot$Individual_ID))

#### project new data clustering ####
# compute total ED (use already computed parts)
ed_dist_old <- clust_res$ed_dist
ed_dist_new <-  as.matrix(dist(input_data, method = 'euclidian'))^2
ed_dist_new_old <- matrix(NA, nrow = nrow(ed_dist_old), ncol = nrow(ed_dist_new))
# build by rows
fun_eucl <- function(x,y){
  sum((x-y)^2)
}
for(i in 1:nrow(input_data)){
  print(i)
  tmp <- apply(clust_res$input_data, 1, function(x) fun_eucl(input_data[i,], x))
  ed_dist_new_old[,i] <- tmp
}
rm(tmp)

tmp <- cbind(ed_dist_old, ed_dist_new_old)
ed_dist_tot <- rbind(tmp, cbind(t(ed_dist_new_old), ed_dist_new))
print(paste0('Block ED matrix symmetric:', isSymmetric(ed_dist_tot)))

res_pred <- project_clust_PGmethod_HKsim(kNN = clust_res$best_k, score = data_tot, data_mod = clust_res$input_data, 
                                         sample_info = sampleAnn_new, euclDist = ed_dist_tot, cl_mod = clust_res$cl_best)

gr_id <- colnames(res_pred$probability)[-(1:2)]
gr_id <- as.numeric(sapply(gr_id, function(x) strsplit(x, split='gr_')[[1]][2]))
cl_new <- data.frame(id = res_pred$probability$Individiual_ID, gr = gr_id[apply(res_pred$probability[, -(1:2)], 1, function(x) which.max(x))])
output <- list(probability = res_pred$probability, tot_W_sNN = res_pred$tot_W_sNN, sampleAnn = sampleAnn_new, data_new = input_data, cl_new = cl_new, res_pval=res_pval, 
               ed_dist_tot = ed_dist_tot)

# compute mean for each gr
test_diff <- clust_res$test_diff_gr
df_gr_mean <- matrix(ncol = length(unique(output$cl_new$gr)), nrow = nrow(test_diff))
df_gr_sd <- matrix(ncol = length(unique(output$cl_new$gr)), nrow = nrow(test_diff))

df_gr_mean[,] <- t(apply(input_data, 2, function(x) 
  sapply(sort(unique(output$cl_new$gr)), function(y) mean(x[output$cl_new$gr == y]) )))

df_gr_sd[,] <- t(apply(input_data, 2, function(x) 
  sapply(sort(unique(output$cl_new$gr)), function(y) sd(x[output$cl_new$gr == y]) )))

df_gr_cv <- df_gr_mean/df_gr_sd
colnames(df_gr_cv) <- colnames(df_gr_mean) <- colnames(df_gr_sd) <- paste0('gr_', sort(unique(output$cl_new$gr)))
df_gr_cv <- as.data.frame(df_gr_cv)
df_gr_mean <- as.data.frame(df_gr_mean)
df_gr_sd <- as.data.frame(df_gr_sd)
df_gr_mean$id <- df_gr_sd$id <- df_gr_cv$id <- colnames(input_data)

output$gr_input <- list(mean = df_gr_mean, sd = df_gr_sd, cv = df_gr_cv)

### plot UMAP and project new samples ###
### plot: UMAP
n_comp_umap <- 2
n_neigh_umap <- clust_res$best_k
min_dist_umap <- 0.01
seed_umap <- 67

custom.settings = umap.defaults
custom.settings$min_dist = min_dist_umap
custom.settings$n_components = n_comp_umap
custom.settings$n_neighbors = n_neigh_umap
custom.settings$random_state <- seed_umap
custom.settings$transform_state <- seed_umap+10

umap_res <- umap::umap(clust_res$input_data, custom.settings)

df_old <- data.frame(component_1=umap_res$layout[,1], component_2=umap_res$layout[,2], gr = clust_res$cl_best$gr)
df_old$gr <- factor(df_old$gr)
df_old$type <- 'model'
# predict
umap_pred <- predict(umap_res, input_data)
df_new <- data.frame(component_1=umap_pred[,1], component_2=umap_pred[,2], gr = output$cl_new$gr)
df_new$gr <- factor(df_new$gr)
df_new$type <- 'predict'
df <- rbind(df_old, df_new)
df$type <- factor(df$type, levels = c('predict', 'model'))

tot_pl <- ggplot(df, aes(x = component_1, y = component_2, color = gr, alpha = type))+
  geom_point(size = 0.05)+
  scale_alpha_manual(values = c(1, 0.05))+
  theme_bw()+theme(legend.position = 'right')
width_pl <- 4

ggsave(filename = sprintf('%s%s_%s_predictCluster%s_PGmethod_HKmetric_umap.png',  outFold, type_data, type_input, type_cluster), width = width_pl, height = 4, plot = tot_pl, device = 'png')
ggsave(filename = sprintf('%s%s_%s_predictCluster%s_PGmethod_HKmetric_umap.pdf',  outFold, type_data, type_input, type_cluster), width = width_pl, height = 4, plot = tot_pl, device = 'pdf')

# save
output$umap_tot <- df
save(output, file = sprintf('%s%s_%s_predictCluster%s_PGmethod_HKmetric.RData',  outFold, type_data, type_input, type_cluster))

### plot: heatmap
file_name <- sprintf('%s%s_%s_predictCluster%s_PGmethod_HKmetric',  outFold, type_data, type_input, type_cluster)

tmp <- test_diff[test_diff$pval_corr < 0.05,]
if(nrow(tmp)>50){
  keep_feat <- test_diff[order(test_diff$pval)[1:50],]
}else{
  keep_feat <- tmp
}

color_tissues <- read.table(color_file, h=T, stringsAsFactors = F)
color_tissues <- color_tissues[match(tissues_name, color_tissues$tissues),]

mat <- output$gr_input$cv[output$gr_input$cv$id %in%  keep_feat$id, ]
width_pl <- 7
if(type_data == 'path_GO'){
  mat$id <- res_pval$path[match(mat$id, res_pval[, id_info])]
  width_pl <- 9
  print(str(mat))
}
if(type_data == 'path_Reactome'){
  width_pl <- 9
}

# remove duplicated
mat <- mat[!duplicated(mat$id),]
mat$tissue <- tissues_name

pheat_pl_gr(mat, type_mat = type_data, height_pl = 7, width_pl = width_pl, color_df = color_tissues, outFile = paste0(file_name, '_heatmap_gr'))

