# plot prediction cluster

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
suppressPackageStartupMessages(library(SNFtool))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(RNOmni))
options(bitmapType = 'cairo', device = 'png')

parser <- ArgumentParser(description="predict cluster probability for new samples")
parser$add_argument("--cohort_name", nargs = '*', type = "character", help = "")
parser$add_argument("--sampleAnn_file", type = "character", help = "")
parser$add_argument("--phenoNew_file", nargs = '*', type = "character", help = "")
parser$add_argument("--type_cluster", type = "character",default = 'All', help = "All, Cases, Controls")
parser$add_argument("--type_data", type = "character", help = "pathway or tscore")
parser$add_argument("--clustFile", type = "character", help = "file cluster results")
parser$add_argument("--clustFile_new", type = "character", nargs = '*', help = "file cluster results")
parser$add_argument("--functR", type = "character", help = "functions to be used")
parser$add_argument("--type_input", type = "character", default = 'original', help = "original or zscaled")
parser$add_argument("--tissues_name", type = "character", help = "name tissue")
parser$add_argument("--outFold", type="character", help = "Output file [basename only]")

args <- parser$parse_args()
cohort_name <- args$cohort_name
sampleAnn_file <- args$sampleAnn_file
clustFile <- args$clustFile
clustFile_new <- args$clustFile_new
functR <- args$functR
type_data <- args$type_data
type_input <- args$type_input
type_cluster <- args$type_cluster
tissues_name <- args$tissues_name
phenoNew_file <- args$phenoNew_file
outFold <- args$outFold

###################################################################################################################
# functR <- '/psycl/g/mpsziller/lucia/priler_project/Software/model_clustering/clustering_functions.R'
# cohort_name <- c('German1', 'German2', 'German3', 'German4','German5')
# type_data <- 'tscore'
# type_input <- 'zscaled'
# type_cluster <- 'Cases'
# clustFile_new <- paste0('/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/OUTPUT_GTEx/predict_CAD/Adipose_Visceral_Omentum/200kb/CAD_GWAS_bin5e-2/',cohort_name,'/devgeno0.01_testdevgeno0/CAD_HARD_clustering/tscore_zscaled_predictClusterCases_PGmethod_HKmetric.RData')
# clustFile <- '/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/OUTPUT_GTEx/predict_CAD/Adipose_Visceral_Omentum/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/CAD_HARD_clustering/tscore_zscaled_clusterCases_PGmethod_HKmetric.RData'
# tissues_name <- 'Adipose_Visceral_Omentum'
# phenoNew_file <-  paste0('/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/INPUT_DATA_GTEx/CAD/Covariates/',cohort_name,'/phenotypeMatrix_CADrel_Cases.txt')
# outFold <- '/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/OUTPUT_GTEx/predict_CAD/Adipose_Visceral_Omentum/200kb/CAD_GWAS_bin5e-2/Meta_Analysis_CAD/CAD_HARD_clustering'
# ##################################################################################################################

source(functR)

sampleAnn <- read.table(sampleAnn_file, h=T, stringsAsFactors = F)
tmp <- get(load(clustFile)) 
clust <- tmp$cl_best
P <- length(unique(clust$gr))
df <- data.frame(dataset = rep('UKBB', P), type = rep('model', P), gr = paste0('gr_', sort(unique(clust$gr))))
df$n <- sapply(sort(unique(clust$gr)), function(x) length(which(clust$gr == x)))
df$percentage <- sapply(sort(unique(clust$gr)), function(x) length(which(clust$gr == x))/nrow(clust))
data <- tmp$input_data
mean_gr <- tmp$gr_input$mean

sampleAnn_new <- list()
phenoDat_new <- list()
clust_new <- list()
data_new <- list()
df_new <- list()
mean_gr_new <- list()
df_corr <- list()

for(i in 1:length(cohort_name)){
  
  tmp <- get(load(clustFile_new[i])) 
  
  sampleAnn_new[[i]] <- tmp$sampleAnn
  #phenoDat_new[[i]] <- read.table(phenoNew_file[i], h=T, stringsAsFactors = F)
 
  #print(all(phenoDat_new[[i]]$Individual_ID %in% sampleAnn_new[[i]]$Individual_ID))
  #phenoDat_new[[i]] <- phenoDat_new[[i]][match(sampleAnn_new[[i]]$Individual_ID, phenoDat_new[[i]]$Individual_ID),]
  
  #sampleAnn_new[[i]] <- sampleAnn_new[[i]][, colnames(sampleAnn_new[[i]]) %in% c('Individual_ID', paste0('C', 1:10))]
  #sampleAnn_new[[i]]$Age <- phenoDat_new[[i]]$Age
  #sampleAnn_new[[i]]$Gender <- phenoDat_new[[i]]$Gender
  
  clust_new[[i]] <- tmp$cl_new
  data_new[[i]] <- tmp$data_new
  mean_gr_new[[i]] <- tmp$gr_input$mean
    
  df_new[[i]] <- data.frame(dataset = rep(cohort_name[i], P), type = rep('predict', P), gr = df$gr)
  df_new[[i]]$n <- sapply(sort(unique(clust$gr)), function(x) length(which(clust_new[[i]]$gr == x)))
  df_new[[i]]$percentage <- sapply(sort(unique(clust$gr)), function(x) length(which(clust_new[[i]]$gr == x))/nrow(clust_new[[i]]))
  
  df_corr[[i]] <- data.frame(dataset = rep(cohort_name[i], P), gr = df$gr, corr = rep(NA, P), pvalue = rep(NA, P), CI_low= rep(NA, P),  CI_up= rep(NA, P))
  df_corr[[i]] <- df_corr[[i]][df_corr[[i]]$gr %in% colnames(mean_gr_new[[i]]), ]
  df_corr[[i]]$corr <- sapply(df$gr[df$gr %in% colnames(mean_gr_new[[i]])], function(x) cor.test(mean_gr_new[[i]][,x], mean_gr[,x])$estimate)
  df_corr[[i]]$pvalue <- sapply(df$gr[df$gr %in% colnames(mean_gr_new[[i]])], function(x) cor.test(mean_gr_new[[i]][,x], mean_gr[,x])$p.value)
  df_corr[[i]]$CI_low <- sapply(df$gr[df$gr %in% colnames(mean_gr_new[[i]])], function(x) cor.test(mean_gr_new[[i]][,x], mean_gr[,x])$conf.int[1])
  df_corr[[i]]$CI_up <- sapply(df$gr[df$gr %in% colnames(mean_gr_new[[i]])], function(x) cor.test(mean_gr_new[[i]][,x], mean_gr[,x])$conf.int[2])
  
}

df_tot <- rbind(df, do.call(rbind, df_new))
df_corr_tot <- do.call(rbind, df_corr)

# save and plot
write.table(df_tot, file = sprintf('%s%s_%s_cluster%s_percentageGropus_prediction_modelUKBB.txt', outFold, type_data, type_input, type_cluster),quote = F, 
            col.names = T, row.names = T, sep = '\t')

write.table(df_corr_tot, file = sprintf('%s%s_%s_cluster%s_correlationMeanGroups_prediction_modelUKBB.txt', outFold, type_data, type_input, type_cluster),quote = F, 
            col.names = T, row.names = T, sep = '\t')

###
df_tot$new_id <- paste0(df_tot$dataset, '\n(', df_tot$type, ')')
df_tot$new_id <- factor(df_tot$new_id, levels = c("UKBB\n(model)",  paste0(cohort_name, '\n(predict)')))
df_tot$gr <- factor(df_tot$gr, levels = paste0('gr_', sort(unique(clust$gr))))
df$gr <- factor(df$gr, levels = paste0('gr_', sort(unique(clust$gr))))

pl <- ggplot(df_tot, aes(x = new_id, y = percentage, color = gr, group = gr))+
  geom_point(position = position_dodge(width = 0.3))+
  theme_bw()+ 
  geom_segment(data = df, aes(x = 2, y = percentage, xend = length(cohort_name)+1, yend = percentage, group = gr), linetype = 2, alpha = 0.6)+
  ylab('Fraction of Cases')+ 
  theme(legend.position = 'right', axis.title.x = element_blank())
  # scale_shape_manual(values=c(1, 19))+
ggsave(filename = sprintf('%s%s_%s_cluster%s_percentageGropus_prediction_modelUKBB.png', outFold, type_data, type_input, type_cluster), width = 5, height = 3.5, plot = pl, device = 'png')
ggsave(filename = sprintf('%s%s_%s_cluster%s_percentageGropus_prediction_modelUKBB.pdf', outFold, type_data, type_input, type_cluster), width = 5, height = 3.5, plot = pl, device = 'pdf')

###
df_corr_tot$dataset <-factor(df_corr_tot$dataset, levels = cohort_name)
df_corr_tot$gr <- factor(df_corr_tot$gr, levels = paste0('gr_', sort(unique(clust$gr))))

pl <- ggplot(df_corr_tot, aes(x = dataset, y = corr, fill = gr, group = gr))+
  geom_bar(stat = 'identity',width = 0.7, color = 'black', alpha = 0.7, position = position_dodge())+
  geom_errorbar(aes(ymin=CI_low, ymax=CI_up), width=.2, position=position_dodge(.75))+
  theme_bw()+ 
  coord_cartesian(ylim = c(ifelse(min(df_corr_tot$corr)<0.8, 0, 0.8),1)) +
  ylab('correlation mean scores\nwith UKBB (model)')+ 
  theme(legend.position = 'right', axis.title.x = element_blank())
# scale_shape_manual(values=c(1, 19))+
ggsave(filename = sprintf('%s%s_%s_cluster%s_correlationMeanGroups_prediction_modelUKBB.png', outFold, type_data, type_input, type_cluster), width = 6, height = 3.5, plot = pl, device = 'png')
ggsave(filename = sprintf('%s%s_%s_cluster%s_correlationMeanGroups_prediction_modelUKBB.pdf', outFold, type_data, type_input, type_cluster), width = 6, height = 3.5, plot = pl, device = 'pdf')


#### endophenotype association ####
phenoInfo_new <- list()
tot_bin_reg <- list()

for(i in 1:length(cohort_name)){
  
  if(file.exists(phenoNew_file[i])){        
    
    phenoDat_new[[i]] <- read.table(phenoNew_file[i], h=T, stringsAsFactors = F)
  
    print(all(phenoDat_new[[i]]$Individual_ID %in% sampleAnn_new[[i]]$Individual_ID))
    phenoDat_new[[i]] <- phenoDat_new[[i]][match(sampleAnn_new[[i]]$Individual_ID, phenoDat_new[[i]]$Individual_ID),]
  
    sampleAnn_new[[i]] <- sampleAnn_new[[i]][, colnames(sampleAnn_new[[i]]) %in% c('Individual_ID', paste0('C', 1:10))]
    sampleAnn_new[[i]]$Age <- phenoDat_new[[i]]$Age
    sampleAnn_new[[i]]$Gender <- phenoDat_new[[i]]$Gender
    
    P <- length(unique(clust_new[[i]]$gr))
    gr_names <- sort(unique(clust_new[[i]]$gr))
    cl <- clust_new[[i]]$gr
    covDat <- sampleAnn_new[[i]][, !colnames(sampleAnn_new[[i]]) %in% c('Individual_ID', 'genoSample_ID', 'Dx')]
    fmla  <- as.formula(paste('pheno~gr_id+', paste0(colnames(covDat), collapse = '+')))
    phenoDat <- phenoDat_new[[i]][, !colnames(phenoDat_new[[i]]) %in%  c('Individual_ID', 'Dx', 'Age', 'Gender')]
    
    phenoInfo_new[[i]] <- data.frame(pheno_id = colnames(phenoDat))
    phenoInfo_new[[i]]$type_pheno <- 'CONTINUOUS'
    phenoInfo_new[[i]]$type_pheno[sapply(1:ncol(phenoDat), function(x) is.integer(phenoDat[,x]) & length(unique(na.omit(phenoDat[,x]))) == 2)] <- 'CAT_SINGLE_BINARY'
    for(j in 1:ncol(phenoDat)){
      if(phenoInfo_new[[i]]$type_pheno[j] == 'CONTINUOUS'){
        tmp <- phenoDat[!is.na(phenoDat[,j]),j]
        phenoDat[!is.na(phenoDat[,j]),j] <- rankNorm(tmp)
      }
    }
    
    
    bin_reg <- vector(mode = 'list', length = length(gr_names)-1)
    for(k in 1:(length(gr_names)-1)){
      
      print(paste0('group', gr_names[k], '_vs_groupj'))
      
      # j vs all
      pheno_case_tmp <- lapply(gr_names[k:length(gr_names)], function(x) phenoDat[cl == x,])
      covDat_tmp <- lapply(gr_names[k:length(gr_names)], function(x) covDat[cl == x,])
      bin_reg[[k]] <-  vector(mode = 'list', length = length(pheno_case_tmp)-1)
      
      for(j in 2:length(pheno_case_tmp)){
        
        print(j)
        
        new <- rbind(pheno_case_tmp[[1]], pheno_case_tmp[[j]])
        gr_id <- factor(c(rep(0, nrow(pheno_case_tmp[[1]])), rep(1, nrow(pheno_case_tmp[[j]]))))
        
        # remove pheno with constant values
        p_rm <- apply(new, 2, function(x) length(unique(x)) == 1)
        new <- new[,!p_rm]
        
        new_cov <- rbind(covDat_tmp[[1]], covDat_tmp[[j]])
        res_glm <- matrix(nrow = ncol(new), ncol = 4)
        for(l in 1:ncol(new)){
          type_pheno <- phenoInfo_new[[i]]$type_pheno[phenoInfo_new[[i]]$pheno_id == colnames(new)[l]]
          tmp_dat <- cbind(data.frame(pheno = new[, l], gr_id = gr_id), new_cov)
          res_glm[l,] <- compute_reg_endopheno(mat = tmp_dat, fmla = fmla, type_pheno = type_pheno)
        }
        colnames(res_glm) <- c('beta', 'se_beta', 'z', 'pvalue')
        res_glm <- as.data.frame(res_glm)

        phenoInfo_tmp <- phenoInfo_new[[i]][match(colnames(new), phenoInfo_new[[i]]$pheno_id),]
        
        bin_reg[[k]][[j-1]] <- cbind(data.frame(pheno_id = phenoInfo_tmp$pheno_id, type_pheno = phenoInfo_tmp$type_pheno), res_glm)
        bin_reg[[k]][[j-1]]$pval_corr <- p.adjust(bin_reg[[k]][[j-1]]$pvalue, method = 'BH')
        bin_reg[[k]][[j-1]]$comp <- sprintf('gr%i_vs_gr%i',  gr_names[k:length(gr_names)][j], gr_names[k])
        
      }
      
    }
    tot_bin_reg[[i]] <- do.call(rbind, do.call(c,bin_reg))
    tot_bin_reg[[i]]$pval_corr_overall <-  p.adjust(tot_bin_reg[[i]]$pvalue, method = 'BH')
  }  
}

# save results
output <- list(bin_reg = tot_bin_reg, cl = clust_new, phenoDat = phenoDat_new, phenoInfo = phenoInfo_new)
save(output, file = sprintf('%s%s_%s_cluster%s_phenoAssociationGLM_prediction_modelUKBB.RData', outFold, type_data, type_input, type_cluster))

