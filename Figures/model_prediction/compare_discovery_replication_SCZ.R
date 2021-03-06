# compare discovery (PGC-SCZ) and replication (CMC) dataset
# consider only DLPC_CMC

options(stringsAsFactors=F)
options(max.print=1000)
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(jaccard))
suppressPackageStartupMessages(library(qvalue))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggrepel))
Sys.setlocale("LC_NUMERIC", "C")
options(bitmapType = 'cairo', device = 'png')


parser <- ArgumentParser(description="compare association results discovery and replication dataset")

parser$add_argument("--discovery_res", type = "character", nargs = '*', help = "file .RData for discovery res")
parser$add_argument("--color_file", type = "character", help = "file with tissues color code")
parser$add_argument("--replication_res", type = "character", nargs = '*', help = "file .RData for replication res")
parser$add_argument("--pheno_name", type = "character", help = "pheno name")
parser$add_argument("--tissues_name", type = "character", nargs = '*', help = "tissues considered")
parser$add_argument("--pval_thr_FDR", type = "double", default = 0.05, help = "thrshold to consider signficant elements (BH correction)")
parser$add_argument("--pval_thr_corr", type = "double", default = 0.01, help = "thrshold to consider signficant elements (no correction) correlation analysis")
parser$add_argument("--outFold", type="character", help = "Output file [basename only]")

args <- parser$parse_args()
discovery_res <- args$discovery_res
color_file <- args$color_file
replication_res <- args$replication_res
pheno_name <- args$pheno_name
tissues_name <- args$tissues_name
pval_thr_corr <- args$pval_thr_corr
pval_thr_FDR <- args$pval_thr_FDR
outFold <- args$outFold

# ####################################################################
# tissues_name <- 'DLPC_CMC'
# discovery_res <- '/ziller/lucia/PGC-SCZ/PriLer_PROJECT/Meta_Analysis_SCZ/DLPC_CMC/pval_Dx_pheno_covCorr.RData'
# replication_res <- '/ziller/lucia/PriLer_PROJECT_CMC/OUTPUT_SCZ-PGC_SCRIPTS_v2/predict_All/DLPC_CMC/devgeno0.01_testdevgeno0/pval_Dx_covCorr.RData'
# outFold <- '/ziller/lucia/PriLer_PROJECT_CMC/OUTPUT_SCZ-PGC_SCRIPTS_v2/predict_All/DLPC_CMC/devgeno0.01_testdevgeno0/'
# pval_thr_FDR <- 0.05
# pval_thr_corr <- 0.01
# color_file <- '/ziller/lucia/priler_project/Figures/color_tissues.txt'
# ####################################################################

#################################
##### version 1: PGC method #####
#################################

# load res
df_tscore_disc <- vector(mode = 'list', length = length(tissues_name))
df_pathR_disc <- vector(mode = 'list', length = length(tissues_name))
df_pathGO_disc <- vector(mode = 'list', length = length(tissues_name))

df_tscore_rep <- vector(mode = 'list', length = length(tissues_name))
df_pathR_rep <- vector(mode = 'list', length = length(tissues_name))
df_pathGO_rep <- vector(mode = 'list', length = length(tissues_name))

dim_mat <- matrix(ncol = 4, nrow = length(tissues_name))
dim_tissues <- matrix(ncol = 4, nrow = length(tissues_name)) 
int_mat <- matrix(ncol = 4, nrow = length(tissues_name))
pval_mat <- matrix(ncol = 4, nrow = length(tissues_name))

for(i in 1:length(tissues_name)){
  
  t <- tissues_name[i]
  print(t)
  
  # discovery
  tmp <- get(load(discovery_res[i]))
  id_pheno <- 1
  
  # remove pathway with only 1 gene, recompute correction
  tmp$pathScore_reactome[[id_pheno]] <- tmp$pathScore_reactome[[id_pheno]][tmp$pathScore_reactome[[id_pheno]]$ngenes_tscore >1,]
  tmp$pathScore_reactome[[id_pheno]][,14] <- qvalue(tmp$pathScore_reactome[[id_pheno]][,13])$qvalues
  tmp$pathScore_reactome[[id_pheno]][,15] <- p.adjust(tmp$pathScore_reactome[[id_pheno]][,13], method = 'BH')
  
  tmp$pathScore_GO[[id_pheno]] <- tmp$pathScore_GO[[id_pheno]][tmp$pathScore_GO[[id_pheno]]$ngenes_tscore >1,]
  tmp$pathScore_GO[[id_pheno]][,16] <- qvalue(tmp$pathScore_GO[[id_pheno]][,15])$qvalues
  tmp$pathScore_GO[[id_pheno]][,17] <- p.adjust(tmp$pathScore_GO[[id_pheno]][,15], method = 'BH')
  
  dim_tissues[i,1] <- nrow(tmp$tscore[[id_pheno]])
  dim_tissues[i,2] <- nrow(tmp$pathScore_reactome[[id_pheno]])
  dim_tissues[i,3] <- nrow(tmp$pathScore_GO[[id_pheno]])
  dim_tissues[i,4] <- nrow(tmp$pathScore_GO[[id_pheno]]) + nrow(tmp$pathScore_reactome[[id_pheno]])
  
  # significance: FDR 0.05
  df_tscore_disc[[i]] <- tmp$tscore[[id_pheno]][tmp$tscore[[id_pheno]][,10] <= pval_thr_FDR, ]
  df_pathR_disc[[i]] <- tmp$pathScore_reactome[[id_pheno]][tmp$pathScore_reactome[[id_pheno]][,15] <= pval_thr_FDR, ]
  df_pathGO_disc[[i]] <- tmp$pathScore_GO[[id_pheno]][tmp$pathScore_GO[[id_pheno]][,17] <= pval_thr_FDR, ]
  
  # replication
  tmp <- get(load(replication_res[i]))
  id_pheno=which(tmp$pheno$pheno_id == pheno_name)
  
  id <- match(df_tscore_disc[[i]][,1],tmp$tscore[[id_pheno]][,1])
  df_tscore_rep[[i]] <- tmp$tscore[[id_pheno]][id, ]
  id <- match(df_pathR_disc[[i]][,1],tmp$pathScore_reactome[[id_pheno]][,1])
  df_pathR_rep[[i]] <- tmp$pathScore_reactome[[id_pheno]][id, ]
  id <- match(df_pathGO_disc[[i]][,1],tmp$pathScore_GO[[id_pheno]][,1])
  df_pathGO_rep[[i]] <- tmp$pathScore_GO[[id_pheno]][id, ]
  
  dim_mat[i,1] <- nrow(df_tscore_disc[[i]])
  dim_mat[i,2] <- nrow(df_pathR_disc[[i]])
  dim_mat[i,3] <- nrow(df_pathGO_disc[[i]])
  dim_mat[i,4] <-  nrow(df_pathR_disc[[i]]) + nrow(df_pathGO_disc[[i]])
  
  int_mat[i,1] <- length(which(sign(df_tscore_disc[[i]][,7])*sign(df_tscore_rep[[i]][,7]) == 1))
  int_mat[i,2] <- length(which(sign(df_pathR_disc[[i]][,12])*sign(df_pathR_rep[[i]][,12]) == 1))
  int_mat[i,3] <- length(which(sign(df_pathGO_disc[[i]][,14])*sign(df_pathGO_rep[[i]][,14]) == 1))
  int_mat[i,4] <- int_mat[i,3] + int_mat[i,2]
  
  pval_mat[i,1] <- ifelse(dim_mat[i,1]>0, binom.test(int_mat[i,1], dim_mat[i,1], p = 0.5, alternative = c("greater"))$p.value, NA)
  pval_mat[i,2] <- ifelse(dim_mat[i,2]>0, binom.test(int_mat[i,2], dim_mat[i,2], p = 0.5, alternative = c("greater"))$p.value, NA)
  pval_mat[i,3] <- ifelse(dim_mat[i,3]>0, binom.test(int_mat[i,3], dim_mat[i,3], p = 0.5, alternative = c("greater"))$p.value, NA)
  pval_mat[i,4] <- ifelse(dim_mat[i,4]>0, binom.test(int_mat[i,4], dim_mat[i,4], p = 0.5, alternative = c("greater"))$p.value, NA)
  
  if(nrow(df_tscore_rep[[i]])>0){
    df_tscore_disc[[i]]$SCZCMC_rep_z <- df_tscore_rep[[i]][,7]
    df_tscore_disc[[i]]$SCZCMC_rep_pval <- df_tscore_rep[[i]][,8]
    df_tscore_disc[[i]]$tissue <- tissues_name[i]
  }
  if(nrow(df_pathR_rep[[i]])>0){
    df_pathR_disc[[i]]$SCZCMC_rep_z <- df_pathR_rep[[i]][,12]
    df_pathR_disc[[i]]$SCZCMC_rep_pval <- df_pathR_rep[[i]][,13]
    df_pathR_disc[[i]]$tissue <- tissues_name[i]
  }
  if(nrow(df_pathGO_rep[[i]])>0){
    df_pathGO_disc[[i]]$SCZCMC_rep_z <- df_pathGO_rep[[i]][,14]
    df_pathGO_disc[[i]]$SCZCMC_rep_pval <- df_pathGO_rep[[i]][,15]
    df_pathGO_disc[[i]]$tissue <- tissues_name[i]
  }
  
}

perc_mat <- int_mat/dim_mat
colnames(perc_mat) = colnames(pval_mat) = colnames(int_mat) = colnames(dim_mat) = colnames(dim_tissues) = c('tscore', 'pathScore_reactome', 'pathScore_GO', 'pathScore_tot')
dim_tissues <- cbind(data.frame(tissue = tissues_name, stringsAsFactors = F), as.data.frame(dim_tissues))
int_mat <- cbind(data.frame(tissue = tissues_name, stringsAsFactors = F), as.data.frame(int_mat))
dim_mat <- cbind(data.frame(tissue = tissues_name, stringsAsFactors = F), as.data.frame(dim_mat))
perc_mat <- cbind(data.frame(tissue = tissues_name, stringsAsFactors = F), as.data.frame(perc_mat))
pval_mat <- cbind(data.frame(tissue = tissues_name, stringsAsFactors = F), as.data.frame(pval_mat))

# add total 
dim_tissues <- rbind(dim_tissues, data.frame(tissue = 'All Tissues', tscore = sum(dim_tissues$tscore), pathScore_reactome = sum(dim_tissues$pathScore_reactome), 
                                             pathScore_GO = sum(dim_tissues$pathScore_GO),  pathScore_tot = sum(dim_tissues$pathScore_tot)))
int_mat <- rbind(int_mat, data.frame(tissue = 'All Tissues', tscore = sum(int_mat$tscore), pathScore_reactome = sum(int_mat$pathScore_reactome),
                                     pathScore_GO = sum(int_mat$pathScore_GO), pathScore_tot = sum(int_mat$pathScore_tot)))
dim_mat <- rbind(dim_mat, data.frame(tissue = 'All Tissues', tscore = sum(dim_mat$tscore), pathScore_reactome = sum(dim_mat$pathScore_reactome), 
                                     pathScore_GO = sum(dim_mat$pathScore_GO), pathScore_tot = sum(dim_mat$pathScore_tot)))
perc_mat <- rbind(perc_mat, data.frame(tissue = 'All Tissues', tscore = int_mat$tscore[int_mat$tissue == 'All Tissues']/dim_mat$tscore[int_mat$tissue == 'All Tissues'], 
                                       pathScore_reactome = int_mat$pathScore_reactome[int_mat$tissue == 'All Tissues']/dim_mat$pathScore_reactome[int_mat$tissue == 'All Tissues'], 
                                       pathScore_GO = int_mat$pathScore_GO[int_mat$tissue == 'All Tissues']/dim_mat$pathScore_GO[int_mat$tissue == 'All Tissues'], 
                                       pathScore_tot = int_mat$pathScore_tot[int_mat$tissue == 'All Tissues']/dim_mat$pathScore_tot[int_mat$tissue == 'All Tissues']))
pval_mat <- rbind(pval_mat, data.frame(tissue = 'All Tissues', tscore = binom.test(int_mat$tscore[int_mat$tissue == 'All Tissues'], dim_mat$tscore[dim_mat$tissue == 'All Tissues'], p = 0.5, alternative = c("greater"))$p.value, 
                                       pathScore_reactome = binom.test(int_mat$pathScore_reactome[int_mat$tissue == 'All Tissues'], dim_mat$pathScore_reactome[dim_mat$tissue == 'All Tissues'], p = 0.5, alternative = c("greater"))$p.value, 
                                       pathScore_GO = binom.test(int_mat$pathScore_GO[int_mat$tissue == 'All Tissues'], dim_mat$pathScore_GO[dim_mat$tissue == 'All Tissues'], p = 0.5, alternative = c("greater"))$p.value, 
                                       pathScore_tot = binom.test(int_mat$pathScore_tot[int_mat$tissue == 'All Tissues'], dim_mat$pathScore_tot[dim_mat$tissue == 'All Tissues'], p = 0.5, alternative = c("greater"))$p.value))

res_tab <- list(dim_tissues = dim_tissues, dim_disc_sign = dim_mat, int_disc_rep = int_mat, perc_disc_rep = perc_mat, pval_disc_rep = pval_mat)

# save
save(res_tab,file = sprintf('%ssignConcordance_discoverySCZPGC_replication%sCMC_signFDR%.2f.RData', outFold, pheno_name, pval_thr_FDR))

df_tscore_disc_allt <- do.call(rbind,df_tscore_disc)
df_pathR_disc_allt <- do.call(rbind,df_pathR_disc)
df_pathGO_disc_allt <- do.call(rbind,df_pathGO_disc)
dim_mat_allt <- c(nrow(df_tscore_disc_allt), nrow(df_pathR_disc_allt),nrow(df_pathGO_disc_allt))
int_mat_allt <- c(length(which(sign(df_tscore_disc_allt[,7])*sign(df_tscore_disc_allt$SCZCMC_rep_z) == 1)), 
                  length(which(sign(df_pathR_disc_allt[,12])*sign(df_pathR_disc_allt$SCZCMC_rep_z) == 1)),
                  length(which(sign(df_pathGO_disc_allt[,14])*sign(df_pathGO_disc_allt$SCZCMC_rep_z) == 1)))
pval_mat_allt <- sapply(1:3, function(x) binom.test(int_mat_allt[x], dim_mat_allt[x], p = 0.5, alternative = c("greater"))$p.value)
df <- data.frame(type = c('tscore', 'pathScore_reactome', 'pathScore_GO'), dim_discSign = dim_mat_allt, int_discSign_rep = int_mat_allt, perc_discSign_rep = int_mat_allt/dim_mat_allt,
                 pval_discSign_rep = pval_mat_allt)

write.table(file = sprintf('%s/tscore_discoverySCZPGC_signFDR%.2f_replication%sCMCInfo_alltissues.txt', outFold, pval_thr_FDR, pheno_name), x = df_tscore_disc_allt, quote = F, sep = '\t', col.names = T, row.names = F)
write.table(file = sprintf('%s/pathR_discoverySCZPGC_signFDR%.2f_replication%sCMCInfo_alltissues.txt', outFold, pval_thr_FDR, pheno_name), x = df_pathR_disc_allt, quote = F, sep = '\t', col.names = T, row.names = F)
write.table(file = sprintf('%s/pathGO_discoverySCZPGC_signFDR%.2f_replication%sCMCInfo_alltissues.txt', outFold,  pval_thr_FDR, pheno_name), x = df_pathGO_disc_allt, quote = F, sep = '\t', col.names = T, row.names = F)
# write.table(file = sprintf('%s/signConcordance_discovery%s_replication_signFDR%.2f_alltissues.txt', outFold, pheno_name, pval_thr_FDR), x = df, quote = F, sep = '\t', col.names = T, row.names = F)


# # plot percentage of concordance with significance
# color_tissues <- read.table(color_file, h=T, stringsAsFactors = F)
# color_tissues <- color_tissues$color[match(tissues_name, color_tissues$tissues)]
# color_tissues <- c(color_tissues, '#666666')
# 
# df <- data.frame(tissue = rep(res_tab$perc_disc_rep$tissue, 3),
#                  type = c(rep('T-score', nrow(res_tab$perc_disc_rep)),rep('Path-score (Reactome)', nrow(res_tab$perc_disc_rep)), rep('Path-score (GO)', nrow(res_tab$perc_disc_rep))), 
#                  perc = c(res_tab$perc_disc_rep$tscore, res_tab$perc_disc_rep$pathScore_reactome, res_tab$perc_disc_rep$pathScore_GO))
# df$pval <- c(res_tab$pval_disc_rep$tscore, res_tab$pval_disc_rep$pathScore_reactome, res_tab$pval_disc_rep$pathScore_GO)
# df$sign_symbol <- ''
# df$sign_symbol[df$pval<=0.05 & df$pval>0.01] <- '*'
# df$sign_symbol[df$pval<=0.01 & df$pval>0.001] <- '**'
# df$sign_symbol[df$pval<=0.001 & df$pval>0.0001] <- '***'
# df$sign_symbol[df$pval<=0.0001] <- '****'
# df$pos <- df$perc+0.03
# df$tissue <- factor(df$tissue, levels = res_tab$perc_disc_rep$tissue)
# df$type <- factor(df$type, levels = c('T-score', 'Path-score (Reactome)', 'Path-score (GO)'))
# 
# pl_bar <- ggplot(data = df, mapping = aes(x = tissue, y = perc, fill = tissue))+
#   geom_bar(stat = 'identity',color = 'black', alpha = 0.7, width = 0.8)+
#   geom_hline(yintercept = 0.5, linetype = 'dashed', size = 1)+
#   geom_text(data = df, aes(x =tissue, y = pos, label = sign_symbol), size = 4, angle = 90)+
#   facet_wrap(.~type, nrow = 1)+
#   theme_bw()+theme(legend.position = 'none', axis.text.y = element_text(colour = color_tissues), axis.title.y = element_blank())+
#   ylab('fraction concordance z')+
#   scale_fill_manual(values = color_tissues)+
#   coord_flip()
# 
# ggsave(plot = pl_bar, filename = sprintf('%s/signConcordance_discovery%s_replication_signFDR%.2f.png', outFold, pheno_name, pval_thr_FDR), width = 9, height = 3.7, dpi = 500)
# ggsave(plot = pl_bar, filename = sprintf('%s/signConcordance_discovery%s_replication_signFDR%.2f.pdf', outFold, pheno_name, pval_thr_FDR), width = 9, height = 3.7,  dpi = 500, compress = F)
# 
# 
# # dot plot
# df$nsign <- c(res_tab$dim_disc_sign$tscore, res_tab$dim_disc_sign$pathScore_reactome, res_tab$dim_disc_sign$pathScore_GO)
# pl_dot <- ggplot(data = df, mapping = aes(x = tissue, y = perc, color = tissue, size = nsign))+
#   geom_point(alpha = 0.7)+
#   geom_hline(yintercept = 0.5, linetype = 'dashed', size = 0.5)+
#   geom_text(data = df, aes(x =tissue, y = pos, label = sign_symbol), size = 4, angle = 90)+
#   facet_wrap(.~type, nrow = 1)+
#   theme_bw()+theme(legend.position = 'bottom', axis.text.y = element_text(colour = color_tissues), axis.title.y = element_blank())+
#   ylab('fraction concordance z')+
#   scale_color_manual(values = color_tissues)+
#   scale_size_continuous(breaks = round(seq(min(df$nsign), max(df$nsign), length.out = 6)))+
#   labs(size = sprintf('n. of genes/pathway significant\n (FDR %.2f)', pval_thr_FDR))+
#   guides(color = FALSE, size=guide_legend(nrow=1,byrow=TRUE)) + coord_flip()
# 
# 
# ggsave(plot = pl_dot, filename = sprintf('%s/signConcordance_discovery%s_replication_signFDR%.2f_dotplot.png', outFold, pheno_name, pval_thr_FDR), width = 9, height = 3.7, dpi = 500)
# ggsave(plot = pl_dot, filename = sprintf('%s/signConcordance_discovery%s_replication_signFDR%.2f_dotplot.pdf', outFold, pheno_name, pval_thr_FDR), width = 9, height = 3.7,  dpi = 500, compress = F)



#################################################
##### version 2: spearman correlation union #####
#################################################

# load res
df_tscore_disc <- vector(mode = 'list', length = length(tissues_name))
df_pathR_disc <- vector(mode = 'list', length = length(tissues_name))
df_pathGO_disc <- vector(mode = 'list', length = length(tissues_name))

df_tscore_rep <- vector(mode = 'list', length = length(tissues_name))
df_pathR_rep <- vector(mode = 'list', length = length(tissues_name))
df_pathGO_rep <- vector(mode = 'list', length = length(tissues_name))

jac_mat <- matrix(ncol = 3, nrow = length(tissues_name))
corP_mat <- matrix(ncol = 3, nrow = length(tissues_name))
pval_corP_mat <- matrix(ncol = 3, nrow = length(tissues_name))
corS_mat <- matrix(ncol = 3, nrow = length(tissues_name))
pval_corS_mat <- matrix(ncol = 3, nrow = length(tissues_name))
dim_tissues <- matrix(ncol = 3, nrow = length(tissues_name)) 
dim_union <- matrix(ncol = 3, nrow = length(tissues_name))

for(i in 1:length(tissues_name)){
  
  t <- tissues_name[i]
  print(t)
  
  # discovery
  tmp <- get(load(discovery_res[i]))
  id_pheno <- 1
  
  dim_tissues[i,1] <- nrow(tmp$tscore[[id_pheno]])
  dim_tissues[i,2] <- nrow(tmp$pathScore_reactome[[id_pheno]])
  dim_tissues[i,3] <- nrow(tmp$pathScore_GO[[id_pheno]])
  
  df_tscore_disc[[i]] <- tmp$tscore[[id_pheno]]
  df_pathR_disc[[i]] <- tmp$pathScore_reactome[[id_pheno]]
  df_pathGO_disc[[i]] <- tmp$pathScore_GO[[id_pheno]]
  
  # replication
  tmp <- get(load(replication_res[i]))
  id_pheno=which(tmp$pheno$pheno_id == pheno_name)
  df_tscore_rep[[i]] <- tmp$tscore[[id_pheno]]
  df_pathR_rep[[i]] <- tmp$pathScore_reactome[[id_pheno]]
  df_pathGO_rep[[i]] <- tmp$pathScore_GO[[id_pheno]]
  
  gene_id <- union(df_tscore_disc[[i]]$ensembl_gene_id[df_tscore_disc[[i]][,8]<=pval_thr_corr], df_tscore_rep[[i]]$ensembl_gene_id[df_tscore_rep[[i]][,8]<=pval_thr_corr])
  pathR_id <- union(df_pathR_disc[[i]]$path[df_pathR_disc[[i]][,13]<=pval_thr_corr], df_pathR_rep[[i]]$path[df_pathR_rep[[i]][,13]<=pval_thr_corr])
  pathGO_id <- union(df_pathGO_disc[[i]]$path[df_pathGO_disc[[i]][,15]<=pval_thr_corr], df_pathGO_rep[[i]]$path[df_pathGO_rep[[i]][,15]<=pval_thr_corr])
  
  tmp_disc <- tmp_rep <- rep(0, dim_tissues[i,1])
  tmp_disc[which(df_tscore_disc[[i]][,8]<=pval_thr_corr)] <- 1
  tmp_rep[which(df_tscore_rep[[i]][,8]<=pval_thr_corr)] <- 1
  jac_mat[i,1] <- jaccard(tmp_disc, tmp_rep)
  
  tmp_disc <- tmp_rep <- rep(0, dim_tissues[i,2])
  tmp_disc[which(df_pathR_disc[[i]][,13]<=pval_thr_corr)] <- 1
  tmp_rep[which(df_pathR_rep[[i]][,13]<=pval_thr_corr)] <- 1
  jac_mat[i,2] <- jaccard(tmp_disc, tmp_rep)
  
  tmp_disc <- tmp_rep <- rep(0, dim_tissues[i,3])
  tmp_disc[which(df_pathGO_disc[[i]][,15]<=pval_thr_corr)] <- 1
  tmp_rep[which(df_pathGO_rep[[i]][,15]<=pval_thr_corr)] <- 1
  jac_mat[i,3] <- jaccard(tmp_disc, tmp_rep)
  
  # filter for union
  df_tscore_disc[[i]] <- df_tscore_disc[[i]][match(gene_id, df_tscore_disc[[i]]$ensembl_gene_id), ]
  df_pathR_disc[[i]] <- df_pathR_disc[[i]][match(pathR_id, df_pathR_disc[[i]]$path), ]
  df_pathGO_disc[[i]] <- df_pathGO_disc[[i]][match(pathGO_id, df_pathGO_disc[[i]]$path), ]
  
  df_tscore_rep[[i]] <- df_tscore_rep[[i]][match(gene_id, df_tscore_rep[[i]]$ensembl_gene_id), ]
  df_pathR_rep[[i]] <- df_pathR_rep[[i]][match(pathR_id, df_pathR_rep[[i]]$path), ]
  df_pathGO_rep[[i]] <- df_pathGO_rep[[i]][match(pathGO_id, df_pathGO_rep[[i]]$path), ]
  
  dim_union[i,1] <- length(gene_id)
  dim_union[i,2] <- length(pathR_id)
  dim_union[i,3] <- length(pathGO_id)
  
  corP_mat[i,1] <- cor.test(df_tscore_disc[[i]][,7], df_tscore_rep[[i]][,7], method = 'pearson')$estimate
  pval_corP_mat[i,1] <- cor.test(df_tscore_disc[[i]][,7], df_tscore_rep[[i]][,7], method = 'pearson')$p.value
  corP_mat[i,2] <- cor.test(df_pathR_disc[[i]][,12], df_pathR_rep[[i]][,12], method = 'pearson')$estimate
  pval_corP_mat[i,2] <- cor.test(df_pathR_disc[[i]][,12], df_pathR_rep[[i]][,12], method = 'pearson')$p.value
  corP_mat[i,3] <- cor.test(df_pathGO_disc[[i]][,14], df_pathGO_rep[[i]][,14], method = 'pearson')$estimate
  pval_corP_mat[i,3] <- cor.test(df_pathGO_disc[[i]][,14], df_pathGO_rep[[i]][,14], method = 'pearson')$p.value
  
  corS_mat[i,1] <- cor.test(df_tscore_disc[[i]][,7], df_tscore_rep[[i]][,7], method = 'spearman')$estimate
  pval_corS_mat[i,1] <- cor.test(df_tscore_disc[[i]][,7], df_tscore_rep[[i]][,7], method = 'spearman')$p.value
  corS_mat[i,2] <- cor.test(df_pathR_disc[[i]][,12], df_pathR_rep[[i]][,12], method = 'spearman')$estimate
  pval_corS_mat[i,2] <- cor.test(df_pathR_disc[[i]][,12], df_pathR_rep[[i]][,12], method = 'spearman')$p.value
  corS_mat[i,3] <- cor.test(df_pathGO_disc[[i]][,14], df_pathGO_rep[[i]][,14], method = 'spearman')$estimate
  pval_corS_mat[i,3] <- cor.test(df_pathGO_disc[[i]][,14], df_pathGO_rep[[i]][,14], method = 'spearman')$p.value
  
}

colnames(jac_mat) = colnames(corS_mat) = colnames(pval_corS_mat) = colnames(corP_mat) = colnames(pval_corP_mat)  = colnames(dim_union) = colnames(dim_tissues) = c('tscore', 'pathScore_reactome', 'pathScore_GO')
dim_tissues <- cbind(data.frame(tissue = tissues_name, stringsAsFactors = F), as.data.frame(dim_tissues))
dim_union <- cbind(data.frame(tissue = tissues_name, stringsAsFactors = F), as.data.frame(dim_union))
corS_mat <- cbind(data.frame(tissue = tissues_name, stringsAsFactors = F), as.data.frame(corS_mat))
pval_corS_mat <- cbind(data.frame(tissue = tissues_name, stringsAsFactors = F), as.data.frame(pval_corS_mat))
corP_mat <- cbind(data.frame(tissue = tissues_name, stringsAsFactors = F), as.data.frame(corP_mat))
pval_corP_mat <- cbind(data.frame(tissue = tissues_name, stringsAsFactors = F), as.data.frame(pval_corP_mat))
jac_mat <- cbind(data.frame(tissue = tissues_name, stringsAsFactors = F), as.data.frame(jac_mat))

res_tab <- list(dim_tissues = dim_tissues, dim_union_sign = dim_union, corS_union = corS_mat, pval_corS_union = pval_corS_mat,
                corP_union = corP_mat, pval_corP_union = pval_corP_mat, jac_sim = jac_mat)
# save
save(res_tab,file = sprintf('%s/corr_discoverySCZPGCSign_replication%sCMCSign_pval%.2f.RData', outFold, pheno_name, pval_thr_corr))


