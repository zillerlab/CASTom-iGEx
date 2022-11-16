options(stringsAsFactors=F)
options(max.print=1000)
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(data.table))
options(bitmapType = 'cairo', device = 'png')

parser <- ArgumentParser(description="pathway correlation multiple datasets matching")
parser$add_argument("--pathScore_file", nargs = '*', type = "character", help = "2 files to be compared")
parser$add_argument("--tissue_name", type = "character", help = "tissue name")
parser$add_argument("--type_path", type = "character", help = "path name")
parser$add_argument("--corr_thr", default = 0.8, type = "double", help = "threshold correlation to filter genes")
parser$add_argument("--outFold", type="character", help = "Output file [basename only]")

args <- parser$parse_args()
pathScore_file <- args$pathScore_file
corr_thr <- args$corr_thr
tissue_name <- args$tissue_name
type_path <- args$type_path
outFold <- args$outFold

###########
# corr_thr <- 0.8
# tissue_name <- 'Brain_Cortex'
# outFold <- '/mnt/lucia/PriLer_PROJECT_GTEx/compare_prediction_UKBB_SCZ-PGC/'
# pathScore_file <- c('../../OUTPUT_SCRIPTS_v2_UKBB/Brain_Cortex/200kb/noGWAS/predict/Pathway_Reactome_scores.txt',
#                       '../../OUTPUT_SCRIPTS_v2_SCZ-PGC/Brain_Cortex/200kb/PGC_GWAS_bin1e-2/predict/Pathway_Reactome_scores.txt')
# ###########

pathScore <- lapply(pathScore_file, function(x) read.delim(x, h=T, stringsAsFactors = F, check.names = F, sep = '\t'))

# check sample order is the same
sample_names <- lapply(pathScore, function(x) colnames(x)[-(1)])
try(if(!do.call(identical, sample_names)) stop('different list of samples'))

# filter to get the same pathways id
common_p <- do.call(intersect, lapply(pathScore, function(x) x[,1]))
mat_path <- lapply(pathScore, function(x) t(x[match(common_p, x[, 1]), -1]))

df_cor <- data.frame(path_id = c() , cor=c(), cor_pval=c())
for(i in 1:length(common_p)){
  print(i)
  tmp <- cor.test(mat_path[[1]][,i], mat_path[[2]][,i], alternative='greater')
  df_cor <- rbind(df_cor, data.frame(path_id = common_p[i], 
                                     cor = tmp$estimate, cor_pval = tmp$p.value))
}

df_cor$keep <- F
df_cor$keep[df_cor$cor >= corr_thr] <- T

write.table(df_cor, file = sprintf('%s%s_filter_path_%s_matched_datasets.txt', outFold, tissue_name, type_path), 
            col.names = T, row.names = F, sep = '\t', quote = F)


