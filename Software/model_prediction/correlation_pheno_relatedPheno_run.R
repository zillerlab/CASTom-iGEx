# find genes and pathway that are significant both for CAD-UKBB and UKBB CAD related phenotypes
# combine all tissues
# compute correlation 

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

parser <- ArgumentParser(description="correlation tscore and pathway phenotypes")
parser$add_argument("--phenoFold", type = "character", help = "Folder with results phenotype annotation UKBB")
parser$add_argument("--inputFold_rel", type = "character", nargs = '*',help = "Folder with results")
parser$add_argument("--pval_FDR_rel", type = "double", default = 0.05, help = "pval threshold to filter the genes and pathways (after BH correction) in related phenotypes")
parser$add_argument("--pval_FDR_pheno", type = "double", default = 0.05,  help = "pval threshold to filter the genes and pathways (after BH correction) in phenotype")
parser$add_argument("--tissue_name", type = "character", nargs = '*', help = "tissue considered")
parser$add_argument("--pheno_name_comp", type = "character", help = "name phenotype of interest")
parser$add_argument("--perc_par", type = "double", default = 0.1, help = "")
parser$add_argument("--pathR_pheno_file", type = "character", help = "")
parser$add_argument("--pathGO_pheno_file", type = "character", help = "")
parser$add_argument("--tscore_pheno_file", type = "character", help = "")
parser$add_argument("--refFold", type = "character", help = "")
parser$add_argument("--feat_filt", type = "character", default = NULL, nargs = 3, help = "first gene, second Reactome, third GO")
parser$add_argument("--thr_dist_par", type = "integer", default = 250000, help = "")
parser$add_argument("--outFold", type="character", help = "Output file [basename only]")

args <- parser$parse_args()
phenoFold <- args$phenoFold
tissue_name <- args$tissue_name
inputFold_rel <- args$inputFold_rel
pval_FDR_rel <- args$pval_FDR_rel
pval_FDR_pheno <- args$pval_FDR_pheno
pval_id <- args$pval_id
pheno_name_comp <- args$pheno_name_comp
perc_par <- args$perc_par
pathR_pheno_file <- args$pathR_pheno_file
pathGO_pheno_file <- args$pathGO_pheno_file
tscore_pheno_file <- args$tscore_pheno_file
refFold <- args$refFold
thr_dist_par <- args$thr_dist_par
feat_filt <- args$feat_filt
outFold <- args$outFold

########################################################################################################################
# phenoFold <- '/psycl/g/mpsziller/lucia/UKBB/eQTL_PROJECT/INPUT_DATA/Covariates/'
# tissue_name <- c('Brain_Frontal_Cortex_BA9')
# pheno_name_comp <- 'SCZ'
# inputFold_rel <- paste0('/psycl/g/mpsziller/lucia/UKBB/eQTL_PROJECT/OUTPUT_GTEx/predict_UKBB/', tissue_name, '/200kb/noGWAS/devgeno0.01_testdevgeno0/')
# tscore_pheno_file <- 'Meta_Analysis_SCZ/OUTPUT_all/tscore_pval_SCZ_covCorr.txt'
# pathR_pheno_file <- 'Meta_Analysis_SCZ/OUTPUT_all/path_Reactome_pval_SCZ_covCorr_filt.txt'
# pathGO_pheno_file <- 'Meta_Analysis_SCZ/OUTPUT_all/path_GO_pval_SCZ_covCorr_filt.txt'
# outFold <- 'Meta_Analysis_SCZ/Brain_Frontal_Cortex_BA9/enrichment_SCZ-UKBB_res/'
# pval_FDR_pheno <- 0.05
# pval_FDR_rel <- 0.05
# perc_par <- 0.3
# thr_dist_par <- 250000
# corrFile <- '/psycl/g/mpsziller/lucia/UKBB/eQTL_PROJECT/OUTPUT_GTEx/predict_UKBB/Brain_Frontal_Cortex_BA9/200kb/noGWAS/devgeno0.01_testdevgeno0/correlation_estimate_tscore.RData'
# geneAnn_file <- '/psycl/g/mpsziller/lucia/UKBB/eQTL_PROJECT/OUTPUT_GTEx/train_GTEx/Brain_Frontal_Cortex_BA9/200kb/noGWAS/resPrior_regEval_allchr.txt'
# refFold <- '/psycl/g/mpsziller/lucia/priler_project/refData/'
# feat_filt <- paste0('/psycl/g/mpsziller/lucia/compare_prediction_UKBB_SCZ-PGC/' ,tissue_name, c('_filter_genes_matched_datasets.txt', '_filter_path_Reactome_matched_datasets.txt', '_filter_path_GO_matched_datasets.txt'))
# #########################################################################################################################

## create function to load data across all tissues ##
load_data <- function(data_file, pval_FDR, tissue_name, pval_id){
  
  tscore <- list()
  pathR <- list()
  pathGO <- list()
  
  for(i in 1:length(tissue_name)){
    
    # print(tissue_name[i])
    res_pval <- get(load(data_file[i]))
    
    tscore[[i]] <- list()
    pathR[[i]] <- list()
    pathGO[[i]] <- list()
    
    for(j in 1:length(pval_id)){
      
      tscore[[i]][[j]] <- res_pval$tscore[[pval_id[j]]]
      tscore[[i]][[j]] <- res_pval$tscore[[pval_id[j]]]
      pathR[[i]][[j]] <- res_pval$pathScore_reactome[[pval_id[j]]]
      pathGO[[i]][[j]] <- res_pval$pathScore_GO[[pval_id[j]]]
      tscore[[i]][[j]]$tissue <- pathR[[i]][[j]]$tissue <- pathGO[[i]][[j]]$tissue <- tissue_name[i]
      # filter out pathways with only 1 gene and recompute adjusted pvalue
      pathR[[i]][[j]] <- pathR[[i]][[j]][pathR[[i]][[j]]$ngenes_tscore>1, ]
      pathR[[i]][[j]][, 14] <- qvalue(pathR[[i]][[j]][,13])$qvalue
      pathR[[i]][[j]][, 15] <- p.adjust(p = pathR[[i]][[j]][,13], method = 'BH')
      
      pathGO[[i]][[j]] <- pathGO[[i]][[j]][pathGO[[i]][[j]]$ngenes_tscore>1, ]
      pathGO[[i]][[j]][, 16] <- qvalue(pathGO[[i]][[j]][,15])$qvalue
      pathGO[[i]][[j]][, 17] <- p.adjust(p = pathGO[[i]][[j]][,15], method = 'BH')
      
    }
  }
  
  pheno_tscore <- pheno_pathR <- pheno_pathGO <- list()
  tscore_red <- pathR_red <- pathGO_red <- list()
  for(j in 1:length(pval_id)){
    
    pheno_tscore[[j]] <- do.call(rbind, lapply(tscore, function(x) x[[j]]))
    pheno_pathR[[j]] <- do.call(rbind, lapply(pathR, function(x) x[[j]]))
    pheno_pathGO[[j]] <- do.call(rbind, lapply(pathGO, function(x) x[[j]]))
    pheno_tscore[[j]]$new_id <- paste0(pheno_tscore[[j]][,1],'_tissue_',pheno_tscore[[j]]$tissue)
    pheno_pathR[[j]]$new_id <- paste0(pheno_pathR[[j]][,1],'_tissue_',pheno_pathR[[j]]$tissue)
    pheno_pathGO[[j]]$new_id <- paste0(pheno_pathGO[[j]][,2],'_tissue_',pheno_pathGO[[j]]$tissue)
    
    tscore_red[[j]] <- pheno_tscore[[j]][pheno_tscore[[j]][, 10] <= pval_FDR,]
    pathR_red[[j]] <- pheno_pathR[[j]][pheno_pathR[[j]][, 15] <= pval_FDR,]
    pathGO_red[[j]] <- pheno_pathGO[[j]][pheno_pathGO[[j]][, 17] <= pval_FDR,]
    
  }
  return(list(tscore = pheno_tscore, tscore_red = tscore_red, pathR = pheno_pathR, pathR_red = pathR_red, 
              pathGO = pheno_pathGO, pathGO_red = pathGO_red))
  
}

permutation_test_corr <-  function(zstat_pheno, zstat_rel, cor_value, nperm = 10000){
  
  n <- length(zstat_rel)
  cor_perm <- c()
  for(i in 1:nperm){
    set.seed(56*i)
    cor_perm <- c(cor_perm , cor(zstat_rel[sample(1:n, n, replace = F)], zstat_pheno, method = 'spearman'))
  }
  
  prob <- sum(cor_perm > abs(cor_value) | cor_perm < -abs(cor_value))/nperm
  return(prob)
  
}

###################################################
## load pheno results
tscore <- read.table(tscore_pheno_file, h=T, stringsAsFactors = F, sep = '\t')
pathR <- read.delim(pathR_pheno_file, h=T, stringsAsFactors = F, sep = '\t')
pathGO <- read.delim(pathGO_pheno_file, h=T, stringsAsFactors = F, sep = '\t')
tscore <- tscore[tscore$tissue %in% tissue_name, ]
pathR <- pathR[pathR$tissue %in% tissue_name, ]
pathGO <- pathGO[pathGO$tissue %in% tissue_name, ]

if(!is.null(feat_filt)){
  gene_filt <- read.table(feat_filt[1], h=T, stringsAsFactors = F, sep = '\t')
  gene_filt <- gene_filt[gene_filt$keep & !is.na(gene_filt$keep), ]
  tscore <- tscore[tscore$ensembl_gene_id %in% gene_filt$ensembl_gene_id, ]
  
  pathR_filt <- read.delim(feat_filt[2], h=T, stringsAsFactors = F, sep = '\t')
  pathR_filt <- pathR_filt[pathR_filt$keep & !is.na(pathR_filt$keep), ]
  pathR <- pathR[pathR$path %in% pathR_filt$path_id, ]
  
  pathGO_filt <- read.delim(feat_filt[3], h=T, stringsAsFactors = F, sep = '\t')
  pathGO_filt <- pathGO_filt[pathGO_filt$keep & !is.na(pathGO_filt$keep), ]
  pathGO <- pathGO[pathGO$path_id %in% pathGO_filt$path_id, ]
  
}

tscore$new_id <- paste0(tscore$ensembl_gene_id, '_tissue_', tscore$tissue)

tot_path <- rbind(cbind(pathR, data.frame(type = rep('Reactome', nrow(pathR)))), 
                  cbind(pathGO[, !colnames(pathGO) %in% c('path_id', 'path_ont')], data.frame(type = rep('GO', nrow(pathGO)))))
tot_path$new_id <- paste0(tot_path$path, '_tissue_', tot_path$tissue, '_type_', tot_path$type)



# create list with pathway and gene annotation
tot_path_ann <- vector(mode = 'list', length = nrow(tot_path))
for(i in 1:nrow(tot_path)){
  tot_path_ann[[i]] <- list(tissue = tot_path$tissue[i], name = tot_path$path[i], id = paste0(tot_path$path[i], '_tissue_', tot_path$tissue[i], '_type_', tot_path$type[i]), 
                            genes = strsplit(tot_path$genes_path[i], split = '[,]')[[1]])  
}

##
# seperately for each tissue, find shared amount of genes and remove
filt_path <- list()
for(i in 1:length(tissue_name)){
  
  print(tissue_name[i])
  id <- sapply(tot_path_ann, function(x) x$tissue == tissue_name[i])
  tmp <- tot_path_ann[id]
  perc_mat <- sapply(tmp, function(y) sapply(tmp, function(x) length(intersect(y$genes, x$genes))/length(union(y$genes, x$genes))))
  keep_path <- tot_path$new_id
  
  # recursevly until no intersection
  while(any(perc_mat[upper.tri(perc_mat)] > perc_par)){
    
    path_list <- apply(perc_mat, 1, function(x) x>perc_par)
    len_path <- c()
    keep_path <- c()
    for(j in 1:nrow(path_list)){
      tmp_sel <-  tot_path[tot_path$new_id %in% sapply(tmp[path_list[j,]], function(x) x$id),]
      tmp_sel <- tmp_sel[!tmp_sel$new_id %in% len_path, ]
      len_path <- unique(c(len_path, tmp_sel$new_id))
      set.seed(j)
      keep_path <- unique(c(keep_path, tmp_sel$new_id[sample(1:nrow(tmp_sel), size = 1)]))
      # keep_path <- unique(c(keep_path, tmp_sel$new_id[which.max(abs(tmp_sel[,12]))]))
    }
    
    new_filt_path <- tot_path[tot_path$new_id %in% keep_path, ]
    id <- sapply(tot_path_ann, function(x) x$id %in% new_filt_path$new_id)
    tmp <- tot_path_ann[id]
    perc_mat <- sapply(tmp, function(y) sapply(tmp, function(x) length(intersect(y$genes, x$genes))/length(union(y$genes, x$genes))))
    
  }
  
  filt_path[[i]] <- tot_path[tot_path$new_id %in% keep_path, ]
  
}

filt_path <- do.call(rbind, filt_path)

##
# seperately for each tissue, remove genes correlated
filt_genes <- list()
biomart_annTSS <- read.table(sprintf("%s/hg19.ENSEMBL_geneTSS_biomart_correct.txt", refFold), h=T, stringsAsFactors = F, sep = '\t')
for(i in 1:length(tissue_name)){
  
  print(tissue_name[i])
  ## load correlation matrix
  #corrMat <- get(load(corrFile[i]))
  #corrMat <- corrMat$cor
  common_g <- intersect(tscore$ensembl_gene_id[tscore$tissue == tissue_name[i]], biomart_annTSS$ensembl_gene_id)
  geneAnn <- biomart_annTSS[match(common_g, biomart_annTSS$ensembl_gene_id),]
  
  dist_mat <- mapply(function(x, y) abs(x - geneAnn$chromstart) < 250000 & y == geneAnn$chrom, x = geneAnn$chromstart, y = geneAnn$chrom)
  tmp <- tscore[tscore$tissue == tissue_name[i], ]
  
  # recursevly until no intersection
  while(any(dist_mat[upper.tri(dist_mat)])){
    
    len_gene <- c()
    keep_gene <- c()
    for(j in 1:nrow(dist_mat)){
      tmp_sel <-  tscore[tscore$new_id %in% tmp$new_id[dist_mat[j,]],]
      tmp_sel <- tmp_sel[!tmp_sel$new_id %in% len_gene, ]
      len_gene <- unique(c(len_gene, tmp_sel$new_id))
      set.seed(j)
      keep_gene <- unique(c(keep_gene, tmp_sel$new_id[sample(1:nrow(tmp_sel), size = 1)]))
      # keep_path <- unique(c(keep_path, tmp_sel$new_id[which.max(abs(tmp_sel[,12]))]))
    }
    keep_gene <- keep_gene[!is.na(keep_gene)]
    
    new_filt_gene <- tscore[tscore$new_id %in% keep_gene, ]
    common_g <- intersect(new_filt_gene$ensembl_gene_id[new_filt_gene$tissue == tissue_name[i]], biomart_annTSS$ensembl_gene_id)
    geneAnn <- biomart_annTSS[match(common_g, biomart_annTSS$ensembl_gene_id),]
    dist_mat <- mapply(function(x, y) abs(x - geneAnn$chromstart) < 250000 & y == geneAnn$chrom, x = geneAnn$chromstart, y = geneAnn$chrom)
    tmp <- new_filt_gene
  }
  filt_genes[[i]] <- tscore[tscore$new_id %in% keep_gene, ]
}

filt_genes <- do.call(rbind, filt_genes)

# # filter based on BHcorr
# tscore_red <- tscore[tscore[, 10] <= pval_FDR_pheno,]
# path_filt_red <- filt_path[filt_path[, 15] <= pval_FDR_pheno,]

#######################################################
### load related phenotype and intersect with pheno ###
if(grepl('CAD',pheno_name_comp)){
  pheno_name <- c(read.table(sprintf('%smatch_cov_pheno_CADrel_filter.txt', phenoFold), h=F, stringsAsFactors = F)$V1, 'ICD10_Anaemia', 'ICD10_Circulatory_system', 'ICD10_Endocrine', 'ICD10_Respiratory_system')
  pheno_name <- pheno_name[!pheno_name %in% c('Medication', 'Medical_conditions')]
  pheno_input <- read.delim(sprintf('%sphenotypeDescription_PHESANTproc_CADrelatedpheno_annotated.txt', phenoFold), h=T, stringsAsFactors = F, sep = '\t')
}

if(grepl('SCZ', pheno_name_comp)){
  
  pheno_name <- c(read.table(sprintf('%smatch_cov_pheno.txt', phenoFold), h=F, stringsAsFactors = F)$V1, 'mixedpheno_Psychiatric', 'ICD10_Psychiatric', 'ICD9_Psychiatric', 
                  read.table(sprintf('%smatch_cov_pheno_SchunkertApp.txt', phenoFold), h=F, stringsAsFactors = F)$V1, 'Blood_biochemistry', 'Blood_count_ratio')
  
  pheno_input <- read.delim(sprintf('%sphenotypeDescription_PHESANTproc.txt', phenoFold), h=T, stringsAsFactors = F, sep = '\t')
  tmp <- read.delim(sprintf('%sphenotypeDescription_PHESANTproc_CADrelatedpheno.txt', phenoFold), h=T, stringsAsFactors = F, sep = '\t')
  pheno_input <- rbind(pheno_input, tmp)
  pheno_input <- pheno_input[!duplicated(pheno_input$pheno_id),]
  tmp <- read.delim(sprintf('%sphenotypeDescription_ratioBC_PHESANTproc.txt', phenoFold), h=T, stringsAsFactors = F, sep = '\t')
  pheno_input <- rbind(pheno_input, tmp)
  pheno_input <- pheno_input[!duplicated(pheno_input$pheno_id),]
  pheno_name <- pheno_name[!pheno_name %in% c('Numeric_memory','Diffusion_brain_MRI', 'Medical_conditions', 'ICD10_Psychiatric', 'ICD9_Psychiatric', 'Resting_functional_brain_MRI',
                                              'dMRI_skeleton', 'T2-weighted_brain_MRI', 'Estimated_nutrients_yesterday')]
  print(pheno_name)
}

pheno_info <- list()
test_tscore <- test_path <- list()

for(j in 1:length(pheno_name)){
  
  print(pheno_name[j])
  
  if(pheno_name[j] %in% c('Blood_biochemistry', 'Blood_count') & grepl('CAD', pheno_name_comp)){
    file_toload <- sprintf('%s/pval_%s_withMed_pheno_covCorr.RData', inputFold_rel, pheno_name[j])
  }else{
    file_toload <- sprintf('%s/pval_%s_pheno_covCorr.RData', inputFold_rel, pheno_name[j])
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
  test_tscore[[j]] <- test_path[[j]] <- data.frame(N = rep(0, length(id_keep)),
                                                   K =rep(0, length(id_keep)), n =rep(0, length(id_keep)), k= rep(0, length(id_keep)), 
                                                   fisher_pval = rep(NA, length(id_keep)), fisher_OR = rep(NA, length(id_keep)), 
                                                   cor_spearman = rep(NA, length(id_keep)), cor_pval = rep(NA, length(id_keep)))
  
  new_pheno <- load_data(data_file = file_toload, pval_FDR = pval_FDR_rel, tissue_name = tissue_name, pval_id = id_keep)
  tscore_pheno <- new_pheno$tscore
  pathR_pheno <- new_pheno$pathR
  pathGO_pheno <- new_pheno$pathGO
  tot_path_pheno <- lapply(1:length(id_keep), function(x) rbind(cbind(pathR_pheno[[x]], data.frame(type = rep('Reactome', nrow(pathR_pheno[[x]])))), 
                                                                cbind(pathGO_pheno[[x]][, !colnames(pathGO_pheno[[x]]) %in% c('path_id', 'path_ont')], data.frame(type = rep('GO', nrow(pathGO_pheno[[x]]))))))
  
  for(i in 1:length(tot_path_pheno)){
    tot_path_pheno[[i]]$new_id <- paste0(tot_path_pheno[[i]]$new_id, '_type_', tot_path_pheno[[i]]$type)
  }
  
  # tscore_pheno <- lapply(tscore_pheno, function(x) x[x$new_id %in% tscore$new_id,])
  # tot_path_pheno <- lapply(tot_path_pheno, function(x) x[x$new_id %in% filt_path$new_id,])
  
  common_g <- lapply(tscore_pheno, function(x) intersect(x$new_id, filt_genes$new_id))
  test_tscore[[j]]$N <- sapply(common_g, length)
  test_tscore[[j]]$K <- sapply(common_g, function(x) nrow(tscore[tscore$new_id %in% x & tscore[,10] <= pval_FDR_pheno,]))
  test_tscore[[j]]$n <- mapply(function(x, y) sum(x[,10]<=pval_FDR_rel & x$new_id %in% y), x = tscore_pheno, y = common_g)
  test_tscore[[j]]$k <- mapply(function(x, y) nrow(x[x$new_id %in% tscore$new_id[tscore[,10] <= pval_FDR_pheno] & x[, 10] <= pval_FDR_rel & x$new_id %in% y,]), x = tscore_pheno, y=common_g)
  vect_pheno <- rep(0,length(common_g[[1]]))
  vect_pheno[tscore[match(common_g[[1]], tscore$new_id), 10]<= pval_FDR_pheno] <- 1
  vect_rel <- lapply(common_g, function(x) rep(0,length(x)))
  
  for(k in 1:length(common_g)){
    vect_rel[[k]][tscore_pheno[[k]][match(common_g[[k]], tscore_pheno[[k]]$new_id), 10]<= pval_FDR_rel] <-1 
  }
  
  # compute fisher test
  id_notnull <- which(sapply(vect_rel, function(x) any(x == 1)))
  if(length(id_notnull)>0 & any(vect_pheno == 1)){
    test_tscore[[j]]$fisher_pval[id_notnull] <- sapply(vect_rel[id_notnull], function(x) fisher.test(x = vect_pheno, y =x, alternative = 'greater')$p.value)
    test_tscore[[j]]$fisher_OR[id_notnull] <- sapply(vect_rel[id_notnull], function(x) fisher.test(x = vect_pheno, y =x, alternative = 'greater')$estimate)
  }
  
  test_tscore[[j]]$cor_spearman <- mapply(function(x, y) cor(x[match(y, x$new_id), 7], tscore[match(y, tscore$new_id), 7], method = 'spearman', use='pairwise.complete.obs'), x =tscore_pheno, y = common_g)
  test_tscore[[j]]$cor_pval <- mapply(function(x, y, z) permutation_test_corr(zstat_pheno = tscore[match(y, tscore$new_id), 7], zstat_rel =  x[match(y, x$new_id), 7], cor_value = z, nperm = 10000),
                                      x = tscore_pheno, y = common_g, z = test_tscore[[j]]$cor_spearman)
  
  ### path ###
  common_p <- lapply(tot_path_pheno, function(x) intersect(x$new_id, filt_path$new_id))
  test_path[[j]]$N <- sapply(common_p, length)
  test_path[[j]]$K <- sapply(common_p, function(x) nrow(filt_path[filt_path$new_id %in% x & filt_path[,15] <= pval_FDR_pheno,]))
  test_path[[j]]$n <- mapply(function(x, y) sum(x[,15]<=pval_FDR_rel & x$new_id %in% y), x = tot_path_pheno, y = common_p)
  test_path[[j]]$k <- mapply(function(x, y) nrow(x[x$new_id %in% filt_path$new_id[filt_path[,15] <= pval_FDR_pheno] & x[, 15] <= pval_FDR_rel & x$new_id %in% y,]), x = tot_path_pheno, y=common_p)
  vect_pheno <- rep(0,length(common_p[[1]]))
  vect_pheno[filt_path[match(common_p[[1]], filt_path$new_id), 15]<= pval_FDR_pheno] <- 1
  vect_rel <- lapply(common_p, function(x) rep(0,length(x)))
  
  for(k in 1:length(common_p)){
    vect_rel[[k]][tot_path_pheno[[k]][match(common_p[[k]], tot_path_pheno[[k]]$new_id), 15]<= pval_FDR_rel] <-1 
  }
  
  # compute fisher test
  id_notnull <- which(sapply(vect_rel, function(x) any(x == 1)))
  if(length(id_notnull)>0 & any(vect_pheno == 1)){
    test_path[[j]]$fisher_pval[id_notnull] <- sapply(vect_rel[id_notnull], function(x) fisher.test(x = vect_pheno, y =x, alternative = 'greater')$p.value)
    test_path[[j]]$fisher_OR[id_notnull] <- sapply(vect_rel[id_notnull], function(x) fisher.test(x = vect_pheno, y =x, alternative = 'greater')$estimate)
  }
  test_path[[j]]$cor_spearman <- mapply(function(x, y) cor(x[match(y, x$new_id), 12], filt_path[match(y, filt_path$new_id), 12], method = 'spearman', use='pairwise.complete.obs'), x = tot_path_pheno, y = common_p)
  test_path[[j]]$cor_pval <- mapply(function(x, y, z) permutation_test_corr(zstat_pheno = filt_path[match(y, filt_path$new_id), 12], zstat_rel =  x[match(y, x$new_id), 12], cor_value = z, nperm = 10000),
                                    x = tot_path_pheno, y = common_p, z = test_path[[j]]$cor_spearman)
  
  # test_path[[j]]$cor_pval <-mapply(function(x, y) cor.test(x[match(y, x$new_id), 12], filt_path[match(y, filt_path$new_id), 12], method = 'spearman', use='pairwise.complete.obs')$p.value, x = tot_path_pheno, y = common_p)
  
}

pheno_info <- do.call(rbind, pheno_info)
pheno_info$names_field <- NA
pheno_info$names_field[is.na(pheno_info$meaning)] <- paste(pheno_info$Field, paste0('(',pheno_info$pheno,')'))[is.na(pheno_info$meaning)]
pheno_info$names_field[!is.na(pheno_info$meaning)] <- paste(pheno_info$Field, ':', pheno_info$meaning, paste0('(',pheno_info$pheno,')'))[!is.na(pheno_info$meaning)]

test_tscore <- do.call(rbind, test_tscore)
test_path <- do.call(rbind, test_path)

# correct pvalues
# correct pvalues
test_tscore$fisher_pval_BHcorr[!is.na(test_tscore$fisher_pval)] <- p.adjust(test_tscore$fisher_pval[!is.na(test_tscore$fisher_pval)], method = 'BH')
test_path$fisher_pval_BHcorr[!is.na(test_path$fisher_pval)] <- p.adjust(test_path$fisher_pval[!is.na(test_path$fisher_pval)], method = 'BH')
test_tscore$cor_pval_BHcorr[!is.na(test_tscore$cor_pval)] <- p.adjust(test_tscore$cor_pval[!is.na(test_tscore$cor_pval)], method = 'BH')
test_path$cor_pval_BHcorr[!is.na(test_path$cor_pval)] <- p.adjust(test_path$cor_pval[!is.na(test_path$cor_pval)], method = 'BH')

# add pheno name
test_tscore$names_field <- pheno_info$names_field
test_path$names_field <- pheno_info$names_field
test_tscore$pheno <- pheno_info$pheno
test_path$pheno <- pheno_info$pheno
test_tscore$pheno_type <- pheno_info$pheno_type
test_path$pheno_type <- pheno_info$pheno_type

# save test results (to be used together for all the tissues)
res_test_enrichment <- list(tscore = test_tscore, pathScore = test_path, pheno = pheno_info, path_ann = filt_path, gene_ann = filt_genes)
save(res_test_enrichment, file = sprintf('%scorrelation_enrich_%s_relatedPheno.RData', outFold, pheno_name_comp))


