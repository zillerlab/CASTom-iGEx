#### code written by Lucia Trastulla, e-mail: lucia_trastulla@psych.mpg.de ####

options(stringsAsFactors=F)
options(max.print=1000)
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(data.table))

parser <- ArgumentParser(description="split data in base and target")

parser$add_argument("--covDat_file", type = "character", help = "file with sample info for new genotype data, must contain Dx column (0 controls 1 cases)")
parser$add_argument("--perc_base", type="double", help = "percentage of samples to be used in base")
parser$add_argument("--perc_ext", type="double", help = "percentage of samples to be used in base")
parser$add_argument("--seed_fixed", type="integer", help = "seed for sample extraction")
parser$add_argument("--outFold", type="character", help = "Output file [basename only]")

args <- parser$parse_args()
covDat_file <- args$covDat_file
perc_base <- args$perc_base
perc_ext <- args$perc_ext
seed_fixed <- args$seed_fixed
outFold <- args$outFold

# ################################################################################################################################
# covDat_file <- 'INPUT_DATA/Covariates/scz_aber_eur.covariateMatrix_old.txt'
# perc_base <- 0.7
# seed_fixed <- 1235
# outFold <- 'INPUT_DATA/Covariates/'
# ################################################################################################################################

sampleAnn <- read.table(covDat_file, header = T, stringsAsFactors = F)
n_cases <- length(which(sampleAnn$Dx == 1))
n_controls <- length(which(sampleAnn$Dx == 0))

# keep percentage case and control
set.seed(seed_fixed)
id_base <-c(sample(sampleAnn$Individual_ID[sampleAnn$Dx == 1], size = round(n_cases*perc_base), replace = F), 
             sample(sampleAnn$Individual_ID[sampleAnn$Dx == 0], size = round(n_controls*perc_base), replace = F))
id_target <- setdiff(sampleAnn$Individual_ID, id_base)

cov_base <- sampleAnn[sampleAnn$Individual_ID %in% id_base, ]
cov_target <- sampleAnn[sampleAnn$Individual_ID %in% id_target, ]

# save
write.table(sprintf('%scovariateMatrix_base_seed%i.txt', outFold, seed_fixed), x = cov_base, col.names = T, row.names = F, sep = '\t', quote = F)
write.table(sprintf('%scovariateMatrix_target_seed%i.txt', outFold, seed_fixed), x = cov_target, col.names = T, row.names = F, sep = '\t', quote = F)

# split into external data
set.seed(9)
id_ext <- sample(1:nrow(cov_target), size = round(nrow(cov_target)*perc_ext),replace = F)
write.table(sprintf('%scovariateMatrix_target_seed%i_externalData.txt', outFold, seed_fixed), x = cov_target[id_ext,], col.names = T, row.names = F, sep = '\t', quote = F)
write.table(sprintf('%scovariateMatrix_target_seed%i_trainData.txt', outFold, seed_fixed), x = cov_target[-id_ext,],col.names = T, row.names = F, sep = '\t', quote = F)


