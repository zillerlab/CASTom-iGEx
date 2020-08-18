# associate endophenotypes to cluster structure 

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

#####################################################################################################################
# phenoDatFile <- '/psycl/g/mpsziller/lucia/UKBB/eQTL_PROJECT/INPUT_DATA/Covariates/MDD_pheno_def/phenotypeMatrix_endoPh_LifetimeMDD.txt'
# phenoDescFile <- '/psycl/g/mpsziller/lucia/UKBB/eQTL_PROJECT/INPUT_DATA/Covariates/MDD_pheno_def/phenotypeDescription_endoPh_LifetimeMDD.txt'
# sampleAnnFile <- '/psycl/g/mpsziller/lucia/UKBB/eQTL_PROJECT/INPUT_DATA/Covariates/MDD_pheno_def/covariateMatrix_LifetimeMDD.txt'
# clusterFile <- '/psycl/g/mpsziller/lucia/UKBB/eQTL_PROJECT/OUTPUT_GTEx/predict_UKBB/Whole_Blood/200kb/noGWAS/devgeno0.01_testdevgeno0/LifetimeMDD_pheno/path_Reactome_zscaled_clusterCases_PGmethod_HKmetric.RData'
# type_cluster <- 'Cases'
# type_data <- 'path_Reactome'
# type_sim <- 'HK'
# outFold <- '/psycl/g/mpsziller/lucia/UKBB/eQTL_PROJECT/OUTPUT_CMC/predict_UKBB/200kb/devgeno0.01_testdevgeno0/LifetimeMDD_pheno/'
# functR <- '/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/RSCRIPTS/SCRIPTS_v2/clustering_functions.R'
# type_input <- 'zscaled'
#####################################################################################################################

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
phenoInfo$Field[phenoInfo$Path == 'Online follow-up > Cognitive function online > Fluid intelligence'] <- paste(phenoInfo$Field[phenoInfo$Path == 'Online follow-up > Cognitive function online > Fluid intelligence'], '(Online)')
cl <- cluster_output$cl_best$gr

output <- list(phenoDat = phenoDat, phenoInfo = phenoInfo, cl = cluster_output$cl_best)

# test phenotype
test_pheno <- data.frame(pheno_id = colnames(phenoDat), Field = phenoInfo$Field,meaning= phenoInfo$Coding_meaning, stringsAsFactors = F)
test_pheno$pval <- NA
test_pheno$pval_BHcorr <- NA
test_pheno$statistic <- NA
test_pheno$test_type <- NA
for(i in 1:nrow(test_pheno)){
  print(i)
  if(is.integer(phenoDat[,i]) | is.character(phenoDat[,i])){
    tmp <- chisq.test(table(phenoDat[,i], cl))
    test_pheno$pval[i] <- tmp$p.value
    test_pheno$statistic[i] <-tmp$statistic
    test_pheno$test_type[i] <- 'chisq'
  }else{
    tmp <- kruskal_test(phenoDat[,i] ~ factor(cl))
    test_pheno$pval[i] <- pvalue(tmp)
    test_pheno$statistic[i] <- tmp@statistic@teststatistic
    test_pheno$test_type[i] <- 'kruskal'
    
  }
}
test_pheno$pval_BHcorr <- p.adjust(test_pheno$pval, method = 'BH')

# save
write.table(x = test_pheno, sprintf('%s%s_%s_cluster%s_PGmethod_%smetric_phenoAssociation.txt', outFold, type_data, type_input, type_cluster, type_sim), col.names = T, row.names = F, sep = '\t', quote = F)

output$test_pheno = test_pheno

###############################################################################
# plot association to covariates, to endophenotypes and Dx_nmi if available

df_cov <- cluster_output$test_cov
df_cov <- df_cov[df_cov$kNN == cluster_output$best_k, ]
df_cov$sign <- 'no'
df_cov$sign[df_cov$pval_BHcorr <= 0.05] <- 'yes'
df_cov$logpval <- -log10(df_cov$pval)
df_cov$sign <- factor(df_cov$sign, levels = c('no', 'yes'))
df_cov$cov_id <- factor(df_cov$cov_id, levels = df_cov$cov_id)

pl <-  ggplot(df_cov, aes(x = cov_id, y = logpval, fill = sign))+
  geom_bar(stat = 'identity', position = position_dodge(), width = 0.5) + theme_bw()+ 
  ylab('-log10(pvalue)')+xlab('')+geom_hline(yintercept = -log10(0.05))+
  scale_fill_manual(values=c("#999999", "#E69F00"))+
  theme(legend.position = 'none', plot.title = element_text(size=9), axis.text.y = element_text(size = 7))+ggtitle(paste(type_data, type_input, 'cluster', type_cluster))+
  coord_flip()
ggsave(filename = sprintf('%s%s_%s_cluster%s_PGmethod_%smetric_covAssociation.png', outFold, type_data, type_input, type_cluster, type_sim), width = 4.5, height = 3, plot = pl, device = 'png')
ggsave(filename = sprintf('%s%s_%s_cluster%s_PGmethod_%smetric_covAssociation.pdf', outFold, type_data, type_input, type_cluster, type_sim), width = 4.5, height = 3, plot = pl, device = 'pdf')


if(type_cluster == 'All'){
  
  perc_info <- cluster_output$Dx_perc$perc[[1]]
  perc_info$gr <- factor(perc_info$gr, levels = sort(unique(perc_info$gr)))
  perc_info$Dx <-  factor(perc_info$Dx)
  df_sign <- data.frame(val = cluster_output$Dx_perc$test[[1]]$fisher_test[-nrow(cluster_output$Dx_perc$test[[1]])])
  df_sign$symb <- ''
  df_sign$symb[df_sign$val<=0.1 & df_sign$val> 0.05] <- '.'
  df_sign$symb[df_sign$val<=0.05 & df_sign$val> 0.01] <- '*'
  df_sign$symb[df_sign$val<=0.01 & df_sign$val> 0.001] <- '**'
  df_sign$symb[df_sign$val<=0.001] <- '***'
  
  pl <-  ggplot(perc_info, aes(x = gr, y = count, fill = Dx))+
    geom_bar(stat = 'identity', width = 0.7) + theme_bw()+ 
    ylab('n. of samples')+xlab('group')+
    annotate("text", x = sort(unique(perc_info$gr)), y = sapply(sort(unique(perc_info$gr)),function(x) sum(perc_info$count[perc_info$gr == x])) + 50, label = df_sign$symb, size = 3) +
    theme(legend.position = 'right', plot.title = element_text(size=9), axis.text.y = element_text(size = 7))+ggtitle(paste(type_data, type_input))
  ggsave(filename = sprintf('%s%s_%s_cluster%s_PGmethod_%smetric_DxPerc.png', outFold, type_data, type_input, type_cluster, type_sim), width = 4.5, height = 3, plot = pl, device = 'png')
  ggsave(filename = sprintf('%s%s_%s_cluster%s_PGmethod_%smetric_DxPerc.pdf', outFold, type_data, type_input, type_cluster, type_sim), width = 4.5, height = 3, plot = pl, device = 'pdf')
  
}

save(output, file = sprintf('%s%s_%s_cluster%s_PGmethod_%smetric_phenoAssociation.RData', outFold, type_data, type_input, type_cluster, type_sim))

test_pheno$sign <- 'no'
test_pheno$sign[test_pheno$pval_BHcorr <= 0.05] <- 'yes'
test_pheno$logpval <- -log10(test_pheno$pval)
test_pheno$sign <- factor(test_pheno$sign, levels = c('no', 'yes'))
test_pheno$new_Field <- test_pheno$Field
test_pheno$new_Field[test_pheno$Field == 'Blood clot, DVT, bronchitis, emphysema, asthma, rhinitis, eczema, allergy diagnosed by doctor'] <- 'Blood clot, etc. diagnosed by doctor'
test_pheno$new_id <- test_pheno$new_Field

test_pheno$new_id[!is.na(test_pheno$meaning)] <- paste(test_pheno$meaning[!is.na(test_pheno$meaning)], test_pheno$new_Field[!is.na(test_pheno$meaning)], sep = '\n')
test_pheno <- test_pheno[order(test_pheno$pval),]

if(length(which(test_pheno$sign == 'yes'))>=20){
  test_pheno <- test_pheno[1:(length(which(test_pheno$sign == 'yes'))+1),]
  height_pl <- 4 + (nrow(test_pheno)-20)*0.2  
}else{
  test_pheno <- test_pheno[1:20,]
  height_pl <- 4
}
test_pheno$new_id <- factor(test_pheno$new_id, levels = test_pheno$new_id)


pl <-  ggplot(test_pheno, aes(x = new_id, y = logpval, fill = sign))+
  geom_bar(stat = 'identity', position = position_dodge(), width = 0.5) + theme_bw()+ 
  ylab('-log10(pvalue)')+xlab('')+geom_hline(yintercept = -log10(0.05))+
  scale_fill_manual(values=c("#999999", "#E69F00"))+
  theme(legend.position = 'none', plot.title = element_text(size=9), axis.text.y = element_text(size = 7))+ggtitle(paste(type_data, type_input, 'cluster', type_cluster))+
  coord_flip()
ggsave(filename = sprintf('%s%s_%s_cluster%s_PGmethod_%smetric_phenoAssociation.png', outFold, type_data, type_input, type_cluster, type_sim), width = 4.5, height = height_pl, plot = pl, device = 'png')
ggsave(filename = sprintf('%s%s_%s_cluster%s_PGmethod_%smetric_phenoAssociation.pdf', outFold, type_data, type_input, type_cluster, type_sim), width = 4.5, height = height_pl, plot = pl, device = 'pdf')


# spider plot for the 10 best phenotypes
library('fmsb')
library('scales')
library('RColorBrewer')

phenoDat_red <- phenoDat[,match(test_pheno$pheno_id[1:10],colnames(phenoDat))]
gr_id <- sort(unique(cl))
perc_gr <- matrix(ncol = length(gr_id),nrow = ncol(phenoDat_red))
for(i in 1:ncol(phenoDat_red)){
  if((is.integer(phenoDat_red[,i]) | is.character(phenoDat_red[,i])) & length(unique(phenoDat_red[,i][!is.na(phenoDat_red[,i])])) == 2){
    print(i)
    perc_gr[i,] <- table(cl, phenoDat_red[,i])[,2]/rowSums(table(cl, phenoDat_red[,i]))
  }else{
    perc_gr[i,] <- sapply(gr_id, function(x) mean(phenoDat_red[cl == x,i], na.rm = T))
  }
}

# transorm to matrix format
df_mat <- t(perc_gr)
colnames(df_mat) <- as.character(test_pheno$new_id[1:10])
rownames(df_mat) <- paste0('group', gr_id)
df_mat <- as.data.frame(df_mat)

# Set graphic colors
coul <- colorRampPalette(brewer.pal(8,"Set2"))(length(gr_id))
colors_border <- alpha(coul, 0.8)
colors_in <- alpha(coul, 0.2)
col_lab <- rep('black', ncol(df_mat))
# col_lab[test_pheno$sign[1:10] == 'yes'] <- 'orange'

# If you remove the 2 first lines, the function compute the max and min of each variable with the available data:
png(sprintf('%s%s_%s_cluster%s_PGmethod_%smetric_phenoAssociation_SP.png', outFold, type_data, type_input, type_cluster, type_sim), width = 10, height = 10, res = 500, units = 'in')
par(xpd = TRUE, mar=c(2,7,2,7))
radarchart(df_mat  , axistype=0, maxmin = F,
           #custom polygon
           pcol=colors_border , pfcol=colors_in , plwd=4 , plty=1,
           #custom the grid
           cglcol="grey", cglty=1, axislabcol="black", cglwd=0.8, 
           #custom labels
           vlcex=1.2)
legend(x=1, y=1.4, legend = rownames(df_mat), bty = "n", pch=20 , col=coul , text.col = "darkgrey", cex=1.2, pt.cex=3)
dev.off()

pdf(sprintf('%s%s_%s_cluster%s_PGmethod_%s_phenoAssociation_SP.pdf', outFold, type_data, type_input, type_cluster, type_sim), width = 10, height = 10, compress = F)
par(xpd = TRUE, mar=c(2,7,2,7))
radarchart(df_mat  , axistype=0, maxmin = F,
           #custom polygon
           pcol=colors_border , pfcol=colors_in , plwd=4 , plty=1,
           #custom the grid
           cglcol="grey", cglty=1, axislabcol="black", cglwd=0.8, 
           #custom labels
           vlcex=1.2)
legend(x=1, y=1.4, legend = rownames(df_mat), bty = "n", pch=20 , col=coul , text.col = "darkgrey", cex=1.2, pt.cex=3)
dev.off()

