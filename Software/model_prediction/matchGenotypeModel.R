#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(argparse))


parser <- argparse::ArgumentParser()

parser$add_argument(
  "--varInfoFile",
  type = "character",
  help = "Full path to reference variant annotation files until \"chr\""
)

parser$add_argument(
  "--snpStatsFile",
  type = "character",
  help = "Full path to the output of qctool until \"chr\""
)

parser$add_argument(
  "--cohortName",
  type = "character",
  help = "Short codename for genotype data"
)

parser$add_argument(
  "--outInfoFold",
  type = "character",
  help = "Folder for matched variant annotation files"
)

parser$add_argument(
  "--altFrqColumn",
  type = "character",
  default = NULL,
  help = "Column name in reference model with alternative allele frequence estimates. Default is NULL (skip this part)"
)

parser$add_argument(
  "--altFrqDiff",
  type = "double",
  default = 0.15,
  help = "Maximum allowed difference in the frequency of alternative allele between reference model and new data"
)


args <- parser$parse_args()


info_name <- basename(args$varInfoFile)


match_stats <- data.frame(
  chr = c(as.character(1:22), "total"),
  snps_model = rep(NA, 23),
  snps_cohort = rep(NA, 23),
  common_snps = rep(NA, 23),
  ref_alt_rev = rep(NA, 23),
  alt_frq_fail = rep(NA, 23),
  snps_final = rep(NA, 23)
)


for (chr in 1:22) {
  priler_ref <- read.table(
    sprintf("%schr%s.txt", args$varInfoFile, chr),
    header = TRUE,
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  ext_gen <- read.table(
    sprintf("%schr%s.snps_stats", args$snpStatsFile, chr),
    header = TRUE,
    stringsAsFactors = FALSE,
    row.names = NULL
  )


  # Although matching is done by position, later filtering will be done via rsId, so need to remove duplicates
  ext_gen <- ext_gen[!duplicated(ext_gen$rsid), ]


  # 'alternate_ids' is actually a chromosome number
  ext_gen$temp_id <- with(ext_gen, paste(alternate_ids, position, alleleA, alleleB, sep = "_"))
  ext_gen$temp_id_rev <- with(ext_gen, paste(alternate_ids, position, alleleB, alleleA, sep = "_"))

  priler_ref$temp_id <- with(priler_ref, paste(CHR, POS, REF, ALT, sep = "_"))


  common_var <- intersect(ext_gen$temp_id, priler_ref$temp_id)
  common_var_rev <- intersect(ext_gen$temp_id_rev, priler_ref$temp_id)

  ext_gen_common <- ext_gen[ext_gen$temp_id %in% common_var | ext_gen$temp_id_rev %in% common_var_rev, ]
  priler_ref_common <- priler_ref[priler_ref$temp_id %in% common_var | priler_ref$temp_id %in% common_var_rev, ]


  rev_id <- which(ext_gen_common$temp_id_rev %in% common_var_rev)

  ext_gen_common[rev_id, c("alleleA", "alleleB")] <- rev(ext_gen_common[rev_id, c("alleleA", "alleleB")])
  ext_gen_common$temp_id <- with(ext_gen_common, paste(alternate_ids, position, alleleA, alleleB, sep = "_"))


  ext_gen_common$ALTfrq <- with(
    ext_gen_common,
    # When minor allele is NA, this is a case where both alleles have a frequency of 0.5
    ifelse(alleleB == minor_allele | is.na(minor_allele), minor_allele_frequency, 1 - minor_allele_frequency)
  )

  if (is.null(args$altFrqColumn)) {
    ext_gen_pass <- ext_gen_common
  } else {
    altfrq_pass <- abs(priler_ref_common[, args$altFrqColumn] - ext_gen_common$ALTfrq) <= args$altFrqDiff

    ext_gen_pass <- ext_gen_common[altfrq_pass, ]
  }


  ext_gen_pass <- subset(ext_gen_pass, select = c(alternate_ids, temp_id, rsid, position, alleleA, alleleB, ALTfrq))
  colnames(ext_gen_pass) <- c("CHR", "ID", "rsID", "POS", "REF", "ALT", "ALTfrq")


  match_stats[chr, ] <- c(
    chr,
    nrow(priler_ref),
    nrow(ext_gen),
    nrow(ext_gen_common),
    length(common_var_rev),
    nrow(ext_gen_common) - nrow(ext_gen_pass),
    nrow(ext_gen_pass)
  )

  write.table(
    ext_gen_pass,
    sprintf("%s/%s%s_chr%s.txt", args$outInfoFold, info_name, args$cohortName, chr),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}


match_stats[23, 2:ncol(match_stats)] <- colSums(match_stats[1:22, 2:ncol(match_stats)])
match_stats$model_perc <- with(match_stats, round(snps_final / snps_model * 100, 2))


write.table(
  match_stats,
  sprintf("%s/%s_match_stats.txt", args$outInfoFold, args$cohortName),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


tot_overlap <- match_stats$model_perc[23]


message(sprintf("\nMatching finished - %s%% reference variants were found in the new data", tot_overlap))

if (tot_overlap < 60) message("\nThis is a low number - the selected model might perform poorly")

message(sprintf("\nDetailed stats can be found in %s/%s_match_stats.txt\n", args$outInfoFold, args$cohortName))
