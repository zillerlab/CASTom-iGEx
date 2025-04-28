#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(argparse)
  library(Matrix)
})


parser <- argparse::ArgumentParser()

parser$add_argument(
  "--VarInfo_file",
  type = "character",
  help = "Full path to the variant information file until 'chr'"
)

parser$add_argument(
  "--outTrain_fold",
  type = "character",
  help = "Path to the folder with the final PriLer model"
)

parser$add_argument(
  "--LDRef_file",
  type = "character",
  help = "Full path to the LD reference file until the chromosome number"
)

parser$add_argument(
  "--modelName",
  type = "character",
  help = "Name for the FUSION-compatible PriLer model"
)

parser$add_argument(
  "--outFold",
  type = "character",
  default = ".",
  help = "Folder to save the results"
)

parser$add_argument(
  "--rsidCol",
  type = "character",
  help = "Name of the column in the variant information file that contains rsIds"
)

parser$add_argument(
  "--N_tot",
  type = "integer",
  help = "Total number of samples used for the training of the model"
)

parser$add_argument(
  "--thr_reliableGenes",
  type = "double",
  nargs = "*",
  default = c(0.01, 0),
  help = "Threshold for reliable genes: dev_geno_tot and test_dev_geno [default %(default)s]"
)


args <- parser$parse_args()

# Create output directories if they don't exist already
if (!dir.exists(sprintf("%s/%s", args$outFold, args$modelName))) {
  dir.create(sprintf("%s/%s", args$outFold, args$modelName), recursive = TRUE)
}

priler_model <- read.table(sprintf("%s/resPrior_regEval_allchr.txt", args$outTrain_fold), header = TRUE)
weights <- get(load(sprintf("%s/resPrior_regCoeffSnps_allchr.RData", args$outTrain_fold)))

N.tot <- args$N_tot

fusion_pos <- list()

for (chr in 1:22) {
  chr_ld_ref <- read.table(sprintf("%s%s.bim", args$LDRef_file, chr))
  chr_variants <- read.table(sprintf("%schr%s.txt", args$VarInfo_file, chr), header = TRUE)

  chr_genes <- subset(priler_model, chrom == paste0("chr", chr))

  reliable_id <- (seq_len(nrow(chr_genes)))[(chr_genes$dev_geno >= args$thr_reliableGenes[1]) & (chr_genes$test_dev_geno > args$thr_reliableGenes[2])]
  chr_genes_reliable <- chr_genes[reliable_id, ]

  chr_weights <- as.matrix(weights[[chr]])[, reliable_id]
  chr_var_sig <- rowSums(chr_weights != 0) != 0

  chr_weights <- chr_weights[chr_var_sig, ]
  chr_variants <- chr_variants[chr_var_sig, ]

  chr_var_ld <- chr_variants[, args$rsidCol, drop = TRUE] %in% chr_ld_ref$V2

  chr_weights <- chr_weights[chr_var_ld, ]
  chr_variants <- chr_variants[chr_var_ld, ]

  for (id in seq_len(nrow(chr_genes_reliable))) {
    chr_weights_id_sig <- chr_weights[, id] != 0

    snps <- as.matrix(data.frame(
      chr = chr,
      rsid = chr_variants[, args$rsidCol, drop = TRUE],
      X3 = 0,
      pos = chr_variants[, "POS", drop = TRUE],
      alt = chr_variants[, "ALT", drop = TRUE],
      ref = chr_variants[, "REF", drop = TRUE]
    ))

    wgt.matrix <- as.matrix(data.frame(
      priler = chr_weights[, id]
    ))

    rownames(wgt.matrix) <- chr_variants[, args$rsidCol, drop = TRUE]

    snps <- snps[chr_weights_id_sig, ]
    wgt.matrix <- wgt.matrix[chr_weights_id_sig, , drop = FALSE]

    # Pval is zero since we define reliable genes differently in PriLer
    cv.performance <- as.matrix(data.frame(
      priler = c(chr_genes_reliable[id, "dev_geno", drop = TRUE], 0)
    ))

    rownames(cv.performance) <- c("rsq", "pval")

    gene_id <- chr_genes_reliable[id, "ensembl_gene_id", drop = TRUE]

    fusion_pos_gene <- list(data.frame(
      PANEL = args$modelName,
      WGT = sprintf("%s/%s.wgt.RDat", args$modelName, gene_id),
      ID = gene_id,
      CHR = chr,
      P0 = chr_genes_reliable[id, "start_position", drop = TRUE],
      P1 = chr_genes_reliable[id, "end_position", drop = TRUE],
      N = args$N_tot,
      priler_r2 = chr_genes_reliable[id, "dev_geno", drop = TRUE],  # Save overall R2 (genetic component only)
      priler_r2cv = chr_genes_reliable[id, "test_dev_geno", drop = TRUE]  # Save average cross-validation R2
    ))

    fusion_pos <- append(fusion_pos, fusion_pos_gene)

    save(
      list = c("cv.performance", "snps", "wgt.matrix", "N.tot"),
      file = sprintf(sprintf("%s/%s/%s.wgt.RDat", args$outFold, args$modelName, gene_id))
    )
  }
}

fusion_pos <- do.call(rbind, fusion_pos)

write.table(fusion_pos, sprintf("%s/%s.pos", args$outFold, args$modelName), sep = "\t", row.names = FALSE)