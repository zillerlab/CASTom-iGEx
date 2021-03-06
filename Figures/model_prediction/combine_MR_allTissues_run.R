# combine results overall tissue
# correlation, MR IVW and MREgg

options(stringsAsFactors=F)
options(max.print=1000)
suppressPackageStartupMessages(library(argparse))
options(bitmapType = 'cairo', device = 'png')

parser <- ArgumentParser(description="combine results MR abd correlations across all tissues")
parser$add_argument("--MR_pheno_file", type = "character", default = NULL, help = "")
parser$add_argument("--mrRes_tissue_file", type = "character", nargs = '*', help = "results for MR analysis tissue spec")
parser$add_argument("--comb_corr_file", type = "character", help = "results for corr analysis all tissues")
parser$add_argument("--thr_pval", type = "double", default = 0.05, help = "")
parser$add_argument("--tissue_name", type = "character",nargs = '*', help = "tissues")
parser$add_argument("--mr_type", type = "character", help = "")
parser$add_argument("--type_data", type = "character", help = "")
parser$add_argument("--type_analysis", type = "character", help = "")
parser$add_argument("--outFold", type="character", help = "Output file [basename only]")

args <- parser$parse_args()
mrRes_tissue_file <- args$mrRes_tissue_file
MR_pheno_file <- args$MR_pheno_file
tissue_name <- args$tissue_name
type_data <- args$type_data
thr_pval <- args$thr_pval
mr_type <- args$mr_type
type_analysis <- args$type_analysis
comb_corr_file <- args$comb_corr_file
outFold <- args$outFold


########################################################################################################################
# tissue_name <- c('Adipose_Subcutaneous','Adipose_Visceral_Omentum', 'Adrenal_Gland','Artery_Coronary',
#                  'Artery_Aorta', 'Colon_Sigmoid', 'Colon_Transverse', 'Heart_Atrial_Appendage',
#                  'Heart_Left_Ventricle', 'Liver', 'Whole_Blood')
# mrRes_tissue_file <- paste0('/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/OUTPUT_GTEx/predict_CAD/', tissue_name,'/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/enrichment_CADHARD_res/CADCardioG_dist200000b_Mendelian_randomization_Egger_tscore_pvalFDRrel0.05.txt')
# mr_type <- 'Egger'
# type_data <- 'tscore'
# type_analysis <- 'not reverse'
# comb_corr_file <- 'OUTPUT_GTEx/predict_CAD/AllTissues/200kb/CAD_GWAS_bin5e-2/UKBB/enrichment_CADHARD_res/CADCardioG_dist200000b_tscore_correlation_sign_pval0.05_allTissues.txt'
######################################################################################################################

corr_sign <- read.delim(comb_corr_file, h=T, stringsAsFactors = F, sep = '\t')

if(!is.null(MR_pheno_file)){
    MR_pheno_subset <- read.table(MR_pheno_file, h = F, stringsAsFactors = F, sep = '\t')$V1
    corr_sign <- corr_sign[corr_sign$pheno_type %in% MR_pheno_subset,]
}

feat_mr <- feat_mr_pval <- feat_mr_FDRpval <- corr_sign[, c('names_field', 'pheno', 'pheno_type')]
if(mr_type == 'Egger'){
  feat_mr_int <- feat_mr_int_pval <- feat_mr_int_FDRpval <- corr_sign[, c('names_field', 'pheno', 'pheno_type')]
}

for(i in 1:length(tissue_name)){
  
  t <- tissue_name[i]
  print(t)
  tmp <- read.delim(mrRes_tissue_file[i], h=T, stringsAsFactors = F, sep = '\t')
  id <- match(corr_sign$pheno, tmp$pheno)
  
  if(mr_type == 'Egger'){
    feat_mr <- cbind(feat_mr,tmp$MREgg_est[id])
    feat_mr_pval <- cbind(feat_mr_pval,tmp$MREgg_est_pval[id])
    feat_mr_FDRpval <- cbind(feat_mr_FDRpval,tmp$MREgg_est_pval_FDR[id])
    feat_mr_int <- cbind(feat_mr_int,tmp$MREgg_int[id])
    feat_mr_int_pval <- cbind(feat_mr_int_pval,tmp$MREgg_int_pval[id])
    feat_mr_int_FDRpval <- cbind(feat_mr_int_FDRpval, tmp$MREgg_int_pval_FDR[id])
    colnames(feat_mr)[ncol(feat_mr)] <- tissue_name[i]
    colnames(feat_mr_pval)[ncol(feat_mr_pval)] <- tissue_name[i]
    colnames(feat_mr_FDRpval)[ncol(feat_mr_FDRpval)] <- tissue_name[i]
    colnames(feat_mr_int)[ncol(feat_mr_int)] <- tissue_name[i]
    colnames(feat_mr_int_pval)[ncol(feat_mr_int_pval)] <- tissue_name[i]
    colnames(feat_mr_int_FDRpval)[ncol(feat_mr_int_FDRpval)] <- tissue_name[i]
  }
  if(mr_type == 'IVW'){
    feat_mr <- cbind(feat_mr, tmp$MRIVW_est[id])
    feat_mr_pval <- cbind(feat_mr_pval,tmp$MRIVW_est_pval[id])
    feat_mr_FDRpval <- cbind(feat_mr_FDRpval,tmp$MRIVW_est_pval_FDR[id])
    colnames(feat_mr)[ncol(feat_mr)] <- tissue_name[i]
    colnames(feat_mr_pval)[ncol(feat_mr_pval)] <- tissue_name[i]
    colnames(feat_mr_FDRpval)[ncol(feat_mr_FDRpval)] <- tissue_name[i]
  }
  
}

# convert extimate to signed results
feat_mr_signedPval <- feat_mr[, 1:3]
for(i in 1:length(tissue_name)){
  if(any(!is.na(feat_mr[, i+3]))){
    feat_mr_signedPval <- cbind(feat_mr_signedPval, sign(feat_mr[, i+3])*-log10(feat_mr_pval[, i+3]))  
  }else{
    feat_mr_signedPval <- cbind(feat_mr_signedPval, data.frame(rep(NA, nrow(feat_mr_signedPval))))
  }
}
colnames(feat_mr_signedPval)[-c(1:3)] <- tissue_name
# matrix with only significant res
feat_mr_signedPval_sign <- feat_mr_signedPval
for(i in 1:length(tissue_name)){
  feat_mr_signedPval_sign[feat_mr_pval[,i+3]> 0.05 & !is.na(feat_mr_pval[,i+3]), i+3] <- NA
}


feat_mr_signedPval_sign <- feat_mr_signedPval_sign[rowSums(!is.na(feat_mr_signedPval_sign[, tissue_name])) > 0,]

name_save <- ifelse(type_analysis == 'reverse', paste0('reverse_', mr_type), mr_type)
print(name_save)

# save
write.table(x = feat_mr, file = sprintf('%s%s_Mendelian_randomization_%s_estimates_allTissues.txt', outFold, type_data, name_save), 
            quote = F, sep = '\t', col.names = T, row.names = F)

write.table(x = feat_mr_pval, file = sprintf('%s%s_Mendelian_randomization_%s_estimates_pvalue_allTissues.txt', outFold, type_data, name_save), 
            quote = F, sep = '\t', col.names = T, row.names = F)

write.table(x = feat_mr_FDRpval, file = sprintf('%s%s_Mendelian_randomization_%s_estimates_FDRpvalue_allTissues.txt', outFold, type_data, name_save),
            quote = F, sep = '\t', col.names = T, row.names = F)

write.table(x = feat_mr_signedPval_sign, file = sprintf('%s%s_Mendelian_randomization_%s_sign_pval%s_allTissues.txt', outFold, type_data, name_save, as.character(thr_pval)), 
            quote = F, sep = '\t', col.names = T, row.names = F)


if(mr_type == 'Egger'){
  
  write.table(x = feat_mr_int, file = sprintf('%s%s_Mendelian_randomization_%s_inter_estimates_allTissues.txt', outFold, type_data, name_save), 
              quote = F, sep = '\t', col.names = T, row.names = F)
  
  write.table(x = feat_mr_int_pval, file = sprintf('%s%s_Mendelian_randomization_%s_inter_estimates_pvalue_allTissues.txt', outFold, type_data, name_save), 
              quote = F, sep = '\t', col.names = T, row.names = F)

  write.table(x = feat_mr_int_FDRpval, file = sprintf('%s%s_Mendelian_randomization_%s_inter_estimates_FDRpvalue_allTissues.txt', outFold, type_data, name_save),
              quote = F, sep = '\t', col.names = T, row.names = F)

  
  
}
