# plot for tissue specific clustering (multiple cohorts)
# genome wide genes
# associated pathways (> 5 genes)
# comparison gr vs all: endophenotype association
# treatment response

# script working with R/4.0.2

options(stringsAsFactors=F)
options(max.print=1000)
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(ggsci))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(PGSEA))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(cowplot))
options(bitmapType = 'cairo', device = 'png')


parser <- ArgumentParser(description="clustering using PG method, show genes/pathways distribution and endophenotype division")
parser$add_argument("--clustFile", type = "character", help = "")
parser$add_argument("--tissue_name", type = "character", help = "")
parser$add_argument("--pheno_name", type = "character", help = "")
parser$add_argument("--tscore_featRelFile", type = "character", help = "")
parser$add_argument("--pathR_featRelFile", type = "character", default  = 'NA', help = "")
parser$add_argument("--pathGO_featRelFile", type = "character",  default  = 'NA', help = "")
parser$add_argument("--geneInfoFile", type = "character",  help = "")
parser$add_argument("--pval_feat", default = 1, type = "double",  help = "")
parser$add_argument("--type_cluster", type = "character", default = 'All', help = "All, Cases, Controls")
parser$add_argument("--type_cluster_data", type = "character", help = "")
parser$add_argument("--min_genes_path", type = "integer", default = 1, help = "minimum number of genes for a pathway, if > 1 recompute corrected pvalues")
parser$add_argument("--functR", type = "character", help = "functions to be used")
parser$add_argument("--reactome_file", type = "character", help = "reactome pathway anntation (.gmt)")
parser$add_argument("--GOterms_file", type = "character", help = "GO pathway anntation (.RData)")
parser$add_argument("--outFold", type="character", help = "Output file [basename only]")

args <- parser$parse_args()
clustFile <- args$clustFile
tissue_name <- args$tissue_name
tscore_featRelFile <- args$tscore_featRelFile
pathR_featRelFile <- args$pathR_featRelFile
pathGO_featRelFile <- args$pathGO_featRelFile
functR <- args$functR
reactome_file <- args$reactome_file
GOterms_file <- args$GOterms_file
geneInfoFile <- args$geneInfoFile
min_genes_path <- args$min_genes_path
type_cluster_data <- args$type_cluster_data
type_cluster <- args$type_cluster
pval_feat <- args$pval_feat
pheno_name <- args$pheno_name
outFold <- args$outFold

#####################################################################################################################
# functR <- '/home/luciat/priler_project/Software/model_clustering/clustering_functions.R'
# type_cluster <- 'Cases'
# pval_feat <- 0.05
# type_cluster_data <- 'tscore'
# tissue_name <- 'Brain_Hypothalamus'
# min_genes_path <- 5
# pheno_name <- 'SCZ'
# geneInfoFile <- sprintf('/home/luciat/eQTL_PROJECT/OUTPUT_GTEx/train_GTEx/%s/200kb/PGC_GWAS_bin1e-2/resPrior_regEval_allchr.txt', tissue_name)
# 
# clustFile <- sprintf('/home/luciat/eQTL_PROJECT/OUTPUT_GTEx/predict_PGC/%s/200kb/PGC_GWAS_bin1e-2/Meta_Analysis_SCZ/devgeno0.01_testdevgeno0/excludeMHC_tscore_zscaled_clusterCases_PGmethod_HKmetric.RData', tissue_name)
# 
# tscore_featRelFile <- sprintf('/home/luciat/eQTL_PROJECT/OUTPUT_GTEx/predict_PGC/%s/200kb/PGC_GWAS_bin1e-2/Meta_Analysis_SCZ/devgeno0.01_testdevgeno0/excludeMHC_tscoreOriginal_tscoreClusterCases_featAssociation.RData', tissue_name)
# 
# pathR_featRelFile<- sprintf('/home/luciat/eQTL_PROJECT/OUTPUT_GTEx/predict_PGC/%s/200kb/PGC_GWAS_bin1e-2/Meta_Analysis_SCZ/devgeno0.01_testdevgeno0/excludeMHC_path_ReactomeOriginal_tscoreClusterCases_featAssociation.RData', tissue_name)
# 
# pathGO_featRelFile<- sprintf('/home/luciat/eQTL_PROJECT/OUTPUT_GTEx/predict_PGC/%s/200kb/PGC_GWAS_bin1e-2/Meta_Analysis_SCZ/devgeno0.01_testdevgeno0/excludeMHC_path_GOOriginal_tscoreClusterCases_featAssociation.RData', tissue_name)
# 
# outFold <- sprintf('/home/luciat/eQTL_PROJECT/OUTPUT_GTEx/predict_PGC/%s/200kb/PGC_GWAS_bin1e-2/Meta_Analysis_SCZ/devgeno0.01_testdevgeno0/excludeMHC_', tissue_name)
# 
# GOterms_file <- '/home/luciat/priler_project/refData/GOterm_geneAnnotation_allOntologies.RData'
# reactome_file <- '/home/luciat/priler_project/refData/ReactomePathways.gmt'
#####################################################################################################################

source(functR)

# get clustering structure
res_cl <- get(load(clustFile))
cl <- res_cl$cl_best
print(identical(res_cl$sampleInfo$Individual_ID, cl$id))
cl$cohort <- res_cl$sampleInfo$cohort

##################################################
######## block T-score + path_Reactome ###########
##################################################

tmp <- get(load(tscore_featRelFile))
tscore_input <- tmp$scaleData
print(identical(rownames(tscore_input), cl$id))

keep_gene <- unique(tmp$test_feat$feat[tmp$test_feat$pval_corr <= pval_feat])
geneInfo <- read.table(geneInfoFile, h=T, stringsAsFactors = F, sep ='\t')
geneInfo <- geneInfo[geneInfo$dev_geno > 0.01 & geneInfo$test_dev_geno > 0, ]
tscore_info <- geneInfo
# remove duplicates (already done for input data)
id_dup <- geneInfo$external_gene_name[duplicated(geneInfo$external_gene_name)]
if(length(id_dup)>0){
  geneInfo <- geneInfo[!geneInfo$external_gene_name %in% id_dup, ]
}
geneInfo_keep <- geneInfo[match(keep_gene, geneInfo$external_gene_name),]
geneInfo_keep$Zstat <- tmp$res_pval[match(keep_gene,tmp$res_pval$external_gene_name),7]
test_feat_tscore <- tmp$test_feat

# remove genes that are from the same locus (use 250kb region) and correlated (LD structure)
# matrix of distance based on start position
chr <- sort(as.numeric(unique(sapply(geneInfo_keep$chrom, function(x) strsplit(x, split  = 'chr')[[1]][2]))))
geneInfo_keep_chr <- list()
for(i in 1:length(chr)){
  print(i)
  tmp_gene <- geneInfo_keep[geneInfo_keep$chrom == paste0('chr', chr[i]),]
  cor_mat <- cor(tscore_input[, match(tmp_gene$external_gene_name, colnames(tscore_input)), drop = F])
  
  if(any(abs(cor_mat[upper.tri(cor_mat)]) > 0.6)){
    locus_list <- apply(cor_mat, 1, function(x) abs(x)>0.6)
    len_gene <- c()
    keep_gene <- c()
    for(j in 1:nrow(locus_list)){
      tmp_sel <-  tmp_gene[locus_list[j,],]
      tmp_sel <- tmp_sel[!tmp_sel$external_gene_name %in% len_gene, ]
      len_gene <- unique(c(len_gene, tmp_sel$external_gene_name))
      keep_gene <- unique(c(keep_gene, tmp_sel$external_gene_name[which.max(abs(tmp_sel$Zstat))]))
    }
    #print(keep_gene)
    geneInfo_keep_chr[[i]] <- tmp_gene[tmp_gene$external_gene_name %in% keep_gene, ]
    dist_mat <- sapply( geneInfo_keep_chr[[i]]$start_position, function(x) abs(x -  geneInfo_keep_chr[[i]]$start_position))
    
    if(any(dist_mat[upper.tri(dist_mat)] < 250000)){
      locus_list <- apply(dist_mat, 1, function(x) x<250000)
      len_gene <- c()
      keep_gene <- c()
      for(j in 1:nrow(locus_list)){
        tmp_sel <-  geneInfo_keep_chr[[i]][locus_list[j,],]
        tmp_sel <- tmp_sel[!tmp_sel$external_gene_name %in% len_gene, ]
        len_gene <- unique(c(len_gene, tmp_sel$external_gene_name))
        keep_gene <- unique(c(keep_gene, tmp_sel$external_gene_name[which.max(abs(tmp_sel$Zstat))]))
      }
      #print(keep_gene)
      geneInfo_keep_chr[[i]] <- geneInfo_keep_chr[[i]][geneInfo_keep_chr[[i]]$external_gene_name %in% keep_gene, ]
    }
    
  }else{
    geneInfo_keep_chr[[i]] <- tmp_gene
  }
}

geneInfo_keep <- do.call(rbind, geneInfo_keep_chr)

# save total feat-gene:
write.table(x = tmp$test_feat, file = sprintf('%stscoreOriginal_%sCluster%s_featAssociation.txt', outFold, type_cluster_data,type_cluster), 
            col.names = T, row.names = F, sep = '\t', quote = F)
# save gene info
geneInfo_save <- geneInfo
geneInfo_save$Zstat <- tmp$res_pval[match(geneInfo_save$external_gene_name,tmp$res_pval$external_gene_name),7]
geneInfo_save$relevant_cluster <- F
geneInfo_save$relevant_cluster[geneInfo_save$external_gene_name %in% unique(tmp$test_feat$feat[tmp$test_feat$pval_corr <= pval_feat])] <- T
geneInfo_save$relevant_cluster_plot <- F
geneInfo_save$relevant_cluster_plot[geneInfo_save$external_gene_name %in% geneInfo_keep$external_gene_name] <- T
write.table(x = geneInfo_save, file = sprintf('%stscoreOriginal_%sCluster%s_infoGenes.txt', outFold, type_cluster_data,type_cluster), 
            col.names = T, row.names = F, sep = '\t', quote = F)


#### pathR ####
tmp <- get(load(pathR_featRelFile))
pathR_input <- tmp$scaleData
identical(rownames(pathR_input), cl$id)
feat_rel <- tmp$test_feat
keep_path <- unique(tmp$test_feat$feat[tmp$test_feat$pval_corr <= pval_feat])
pathInfo <- tmp$res_pval
pathInfo_keep <- pathInfo[match(keep_path, pathInfo[,1]),]
pathInfo_keep <- pathInfo_keep[pathInfo_keep$ngenes_tscore >= min_genes_path,]
pathInfo_keep$Zstat <- pathInfo_keep[,12]

if(nrow(pathInfo_keep)>0){
  
  test_feat_pathR <- tmp$test_feat
  # filter out pathway based on shared percentage
  gs <- readGmt(reactome_file)
  
  # annotate gs with the used genes
  for(i in 1:length(gs)){
    gs[[i]]@data <- gs[[i]]@ids[gs[[i]]@ids %in% tscore_info$external_gene_name]
  }
  pathInfo_save <- pathInfo[, colnames(pathInfo) %in% c('path',"ngenes_tscore","ngenes_path", 
                                                        "mean_dev_geno","sd_dev_geno","mean_test_dev_geno",
                                                        "sd_test_dev_geno","mean_gene_corr","sd_gene_corr")]
  pathInfo_save$genes_id <- NA
  for(i in 1:nrow(pathInfo_save)){
    id_tmp <- which(sapply(gs, function(x) x@reference == pathInfo_save$path[i]))
    pathInfo_save$genes_id[i] <- paste0(unique(gs[[id_tmp]]@data), collapse = ',')
  }
  pathInfo_save$Zstat <- pathInfo[match(pathInfo_save$path,pathInfo$path),12]
  
  gs <- gs[which(sapply(gs, function(x) x@reference %in% pathInfo_keep[,1]))]
  cor_mat <- cor(pathR_input[, match(pathInfo_keep[,1], colnames(pathR_input)), drop = F])
  # remove higly correlated (same gene or genes in LD) and pathways that share high percentage of genes
  if(any(abs(cor_mat[upper.tri(cor_mat)]) > 0.6)){
    path_list <- apply(cor_mat, 1, function(x) abs(x)>0.6)
    len_path <- c()
    keep_path <- c()
    for(j in 1:nrow(path_list)){
      tmp_sel <- pathInfo_keep[pathInfo_keep[,1] %in% pathInfo_keep[path_list[j,],1],]
      tmp_sel <-  tmp_sel[!tmp_sel[,1] %in% len_path, ]
      len_path <- unique(c(len_path, tmp_sel[,1]))
      keep_path <- unique(c(keep_path, tmp_sel[which.max(abs(tmp_sel$Zstat)), 1]))
    }
    #print(keep_gene)
    pathInfo_keep <- pathInfo_keep[pathInfo_keep[,1] %in% keep_path, ]
    gs <- gs[which(sapply(gs, function(x) x@reference %in% pathInfo_keep[,1]))]
    perc_mat <- sapply(gs, function(y) sapply(gs, function(x) length(intersect(y@data, x@data))/length(union(y@data, x@data))))
    if(any(perc_mat[upper.tri(perc_mat)] > 0.1)){
      path_list <- apply(perc_mat, 1, function(x) x>0.1)
      len_path <- c()
      keep_path <- c()
      for(j in 1:nrow(path_list)){
        tmp_sel <-  pathInfo_keep[pathInfo_keep[,1] %in% sapply(gs[path_list[j,]], function(x) x@reference),]
        tmp_sel <- tmp_sel[!tmp_sel[,1] %in% len_path, ]
        len_path <- unique(c(len_path, tmp_sel[,1]))
        keep_path <- unique(c(keep_path, tmp_sel[which.max(abs(tmp_sel$Zstat)), 1]))
      }
      #print(keep_path)
      pathInfo_keep <- pathInfo_keep[pathInfo_keep[,1] %in% keep_path, ]
    }
  }
  
  # any pathway have a better performance than the single genes reported?
  P <- length(unique(cl$gr))
  comp <- unique(test_feat_tscore$comp)
  impr_path <- c()
  for(i in 1:nrow(pathInfo_keep)){
    tmp <- gs[[which(sapply(gs, function(x) x@reference == pathInfo_keep[i,1]))]]
    test_feat_tscore_tmp <- test_feat_tscore[test_feat_tscore$feat %in% tmp@data,]
    test_feat_path_tmp <- test_feat_pathR[test_feat_pathR$feat == pathInfo_keep[i,1], ]
    test_feat_path_tmp <- test_feat_path_tmp[match(comp, test_feat_path_tmp$comp),]
    tmp_gr_tscore <- lapply(comp, function(x) test_feat_tscore_tmp$pval[test_feat_tscore_tmp$comp == x])
    tscore_path_comp <- mapply(function(x, y, z) all(x < y & z <= 0.05), y = tmp_gr_tscore, x = test_feat_path_tmp$pval, z = test_feat_path_tmp$pval_corr)
    if(any(tscore_path_comp)){
      print(i)
      impr_path <- c(impr_path, pathInfo_keep[i,1])
    }
  }
  pathInfo_keep$impr <- F
  pathInfo_keep$impr[pathInfo_keep[,1] %in% impr_path] <- T
  
  # save total feat-path:
  write.table(x = feat_rel, file = sprintf('%spath_ReactomeOriginal_%sCluster%s_featAssociation.txt', outFold, type_cluster_data,type_cluster), 
              col.names = T, row.names = F, sep = '\t', quote = F)
  
  # save path info
  pathInfo_save$relevant_cluster <- F
  pathInfo_save$relevant_cluster[pathInfo_save$path %in%  unique(feat_rel$feat[feat_rel$pval_corr <= pval_feat])] <- T
  pathInfo_save$relevant_cluster_plot <- F
  pathInfo_save$relevant_cluster_plot[pathInfo_save$path %in% pathInfo_keep$path] <- T
  write.table(x = pathInfo_save, file = sprintf('%spath_ReactomeOriginal_%sCluster%s_infopath.txt', outFold, type_cluster_data,type_cluster), 
              col.names = T, row.names = F, sep = '\t', quote = F)
  
  
  ### plot ###
  pheat_pl_tot(pheno_name = pheno_name, mat_tscore = tscore_input, mat_path = pathR_input, pval_thr_est = 0.05,
               info_feat_tscore = geneInfo_keep, info_feat_path = pathInfo_keep, cl = cl, 
               test_feat_tscore = test_feat_tscore, test_feat_path = test_feat_pathR, 
               cap = 3, res_pl = 250, width_pl = 14 + round(length(unique(cl$gr))*0.25), 
               height_pl = round(7 + (nrow(geneInfo_keep)+nrow(pathInfo_keep))*0.1), 
               outFile = sprintf('%stscoreOriginal_path_ReactomeOriginal_%sCluster%s',outFold,type_cluster_data,type_cluster))
  
  print('tscore-path_Reactome heatmap finished')
}else{
  print('no significant tscore-path_Reactome')
}
############################################
######## block T-score + path_GO ###########
############################################

tmp <- get(load(pathGO_featRelFile))
pathGO_input <- tmp$scaleData
identical(rownames(pathGO_input), cl$id)

test_feat_pathGO <- tmp$test_feat
print(identical(tmp$res_pval$path_id, colnames(pathGO_input)))
print(identical(tmp$res_pval$path_id, tmp$test_feat$feat))
colnames(pathGO_input) <- tmp$res_pval$path
test_feat_pathGO$feat <- tmp$res_pval$path


keep_path <- unique(tmp$test_feat$feat[tmp$test_feat$pval_corr <= pval_feat])
pathInfo <- tmp$res_pval
pathInfo_keep <- pathInfo[match(keep_path, pathInfo[,1]),]
pathInfo_keep <- pathInfo_keep[pathInfo_keep$ngenes_tscore >= min_genes_path + 5,]
pathInfo_keep$Zstat <- pathInfo_keep[,14]
pathInfo_keep[,1] <- pathInfo_keep$path

if(nrow(pathInfo_keep)>0){
  # filter out pathway based on shared percentage
  go <- get(load(GOterms_file))
  # annotate go with the used genes
  for(i in 1:length(go)){
    go[[i]]$new <- go[[i]]$geneIds[go[[i]]$geneIds %in% tscore_info$external_gene_name]
  }
  pathInfo_save <- pathInfo[, colnames(pathInfo) %in% c('path_id','path',"ngenes_tscore","ngenes_path", 
                                                        "mean_dev_geno","sd_dev_geno","mean_test_dev_geno",
                                                        "sd_test_dev_geno","mean_gene_corr","sd_gene_corr")]
  pathInfo_save$genes_id <- NA
  for(i in 1:nrow(pathInfo_save)){
    id_tmp <- which(sapply(go, function(x) x$Term == pathInfo_save$path[i]))
    pathInfo_save$genes_id[i] <- paste0(unique(go[[id_tmp]]$new), collapse = ',')
  }
  pathInfo_save$Zstat <- pathInfo[match(pathInfo_save$path,pathInfo$path),14]
  
  go <- go[which(sapply(go, function(x) x$Term %in% pathInfo_keep[,1]))]
  cor_mat <- cor(pathGO_input[, match(pathInfo_keep[,1], colnames(pathGO_input)), drop = F])
  # remove higly correlated (same gene or genes in LD) and pathways that share high percentage of genes
  if(any(abs(cor_mat[upper.tri(cor_mat)]) > 0.6)){
    path_list <- apply(cor_mat, 1, function(x) abs(x)>0.6)
    len_path <- c()
    keep_path <- c()
    for(j in 1:nrow(path_list)){
      tmp_sel <- pathInfo_keep[pathInfo_keep[,1] %in% pathInfo_keep[path_list[j,],1],]
      tmp_sel <-  tmp_sel[!tmp_sel[,1] %in% len_path, ]
      len_path <- unique(c(len_path, tmp_sel[,1]))
      keep_path <- unique(c(keep_path, tmp_sel[which.max(abs(tmp_sel$Zstat)), 1]))
    }
    #print(keep_gene)
    pathInfo_keep <- pathInfo_keep[pathInfo_keep[,1] %in% keep_path, ]
    go <- go[which(sapply(go, function(x) x$Term %in% pathInfo_keep[,1]))]
    
    perc_mat <- sapply(go, function(y) sapply(go, function(x) length(intersect(y$new, x$new))/length(union(y$new, x$new))))
    if(any(perc_mat[upper.tri(perc_mat)] > 0.1)){
      path_list <- apply(perc_mat, 1, function(x) x>0.1)
      len_path <- c()
      keep_path <- c()
      for(j in 1:nrow(path_list)){
        tmp_sel <-  pathInfo_keep[pathInfo_keep[,1] %in% sapply(go[path_list[j,]], function(x) x$Term),]
        tmp_sel <- tmp_sel[!tmp_sel[,1] %in% len_path, ]
        len_path <- unique(c(len_path, tmp_sel[,1]))
        keep_path <- unique(c(keep_path, tmp_sel[which.max(abs(tmp_sel$Zstat)), 1]))
      }
      #print(keep_path)
      pathInfo_keep <- pathInfo_keep[pathInfo_keep[,1] %in% keep_path, ]
    }
  }
  
  # any pathway have a better performance than the single genes reported?
  P <- length(unique(cl$gr))
  comp <- unique(test_feat_tscore$comp)
  impr_path <- c()
  for(i in 1:nrow(pathInfo_keep)){
    tmp <- go[[which(sapply(go, function(x) x$Term == pathInfo_keep[i,1]))]]
    test_feat_tscore_tmp <- test_feat_tscore[test_feat_tscore$feat %in% tmp$new,]
    test_feat_path_tmp <- test_feat_pathGO[test_feat_pathGO$feat == pathInfo_keep[i,1], ]
    test_feat_path_tmp <- test_feat_path_tmp[match(comp, test_feat_path_tmp$comp),]
    tmp_gr_tscore <- lapply(comp, function(x) test_feat_tscore_tmp$pval[test_feat_tscore_tmp$comp == x])
    tscore_path_comp <- mapply(function(x, y, z) all(x < y & z <= 0.05), y = tmp_gr_tscore, x = test_feat_path_tmp$pval, z = test_feat_path_tmp$pval_corr)
    if(any(tscore_path_comp)){
      print(i)
      impr_path <- c(impr_path, pathInfo_keep[i,1])
    }
  }
  pathInfo_keep$impr <- F
  pathInfo_keep$impr[pathInfo_keep[,1] %in% impr_path] <- T
  
  # save total feat-path:
  write.table(x = test_feat_pathGO, file = sprintf('%spath_GOOriginal_%sCluster%s_featAssociation.txt', outFold, type_cluster_data, type_cluster), 
              col.names = T, row.names = F, sep = '\t', quote = F)
  
  # save path info
  pathInfo_save$relevant_cluster <- F
  pathInfo_save$relevant_cluster[pathInfo_save$path %in%  unique(test_feat_pathGO$feat[test_feat_pathGO$pval_corr <= pval_feat])] <- T
  pathInfo_save$relevant_cluster_plot <- F
  pathInfo_save$relevant_cluster_plot[pathInfo_save$path %in% pathInfo_keep$path] <- T
  write.table(x = pathInfo_save, file = sprintf('%spath_GOOriginal_%sCluster%s_infopath.txt', outFold, type_cluster_data,type_cluster), 
              col.names = T, row.names = F, sep = '\t', quote = F)
  
  ### plot ###
  pheat_pl_tot(pheno_name = pheno_name, mat_tscore = tscore_input, mat_path = pathGO_input, pval_thr_est = 0.05,
               info_feat_tscore = geneInfo_keep, info_feat_path = pathInfo_keep, cl = cl, 
               test_feat_tscore = test_feat_tscore, test_feat_path = test_feat_pathGO, 
               cap = 3, res_pl = 250,  width_pl = round(13 + length(unique(cl$gr))*0.25), 
               height_pl = round(7 + (nrow(geneInfo_keep)+nrow(pathInfo_keep))*0.1), 
               outFile = sprintf('%stscoreOriginal_path_GOOriginal_%sCluster%s', outFold,type_cluster_data,type_cluster))
  
  print('tscore-path_GO heatmap finished')
  
}else{
  print('no significant tscore-path_GO')
}



