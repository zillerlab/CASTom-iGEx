options(stringsAsFactors=F)
options(max.print=1000)
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(data.table))
options(bitmapType = 'cairo', device = 'png')

parser <- ArgumentParser(description="predicted gene expression correlation multiple datasets matching")
parser$add_argument("--geneExpPred_file", nargs = '*', type = "character", help = "2 files to be compared")
parser$add_argument("--tissue_name", type = "character", help = "tissue name")
parser$add_argument("--corr_thr", default = 0.8, type = "double", help = "threshold correlation to filter genes")
parser$add_argument("--outFold", type="character", help = "Output file [basename only]")

args <- parser$parse_args()
geneExpPred_file <- args$geneExpPred_file
corr_thr <- args$corr_thr
tissue_name <- args$tissue_name
outFold <- args$outFold

###########
# corr_thr <- 0.8
# tissue_name <- 'Brain_Cortex'
# outFold <- '/mnt/lucia/PriLer_PROJECT_GTEx/compare_prediction_UKBB_SCZ-PGC/'
# geneExpPred_file <- c('OUTPUT_SCRIPTS_v2_UKBB/Brain_Cortex/200kb/noGWAS/predict/predictedExpression.txt.gz', 
#                       'OUTPUT_SCRIPTS_v2_SCZ-PGC/Brain_Cortex/200kb/PGC_GWAS_bin1e-2/predict/predictedExpression.txt.gz')
# ###########

geneExpPred <- lapply(geneExpPred_file, function(x) read.table(gzfile(x), h=T, stringsAsFactors = F, check.names = F, sep = '\t'))

# check sample and gene order is the same
gene_names <- lapply(geneExpPred, function(x) x$ensembl_gene_id)
sample_names <- lapply(geneExpPred, function(x) colnames(x)[-(1:29)])

try(if(!do.call(identical, gene_names)) stop('different list of genes'))
try(if(!do.call(identical, sample_names)) stop('different list of samples'))

mat_genes <- lapply(geneExpPred, function(x) t(x[, -(1:29)]))

df_cor <- data.frame(ensembl_gene_id = c(), external_gene_name = c() , cor=c(), cor_pval=c())
for(i in 1:ncol(mat_genes[[1]])){
  print(i)
  tmp <- cor.test(mat_genes[[1]][,i], mat_genes[[2]][,i], alternative='greater')
  df_cor <- rbind(df_cor, data.frame(ensembl_gene_id = geneExpPred[[1]]$ensembl_gene_id[i], 
                                     external_gene_name = geneExpPred[[1]]$external_gene_name[i],
                                     cor = tmp$estimate, cor_pval = tmp$p.value))
}

df_cor$keep <- NA
id_rel <- (geneExpPred[[2]]$dev_geno >= 0.01 & geneExpPred[[2]]$test_dev_geno > 0) & 
  (geneExpPred[[1]]$dev_geno >= 0.01 & geneExpPred[[1]]$test_dev_geno > 0)
df_cor$keep[id_rel] <- F
df_cor$keep[df_cor$cor >= corr_thr & id_rel] <- T

write.table(df_cor, file = sprintf('%s%s_filter_genes_matched_datasets.txt', outFold, tissue_name), 
            col.names = T, row.names = F, sep = '\t', quote = F)


