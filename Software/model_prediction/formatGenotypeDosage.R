#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(argparse)
  library(data.table)
})


parser <- argparse::ArgumentParser()

parser$add_argument(
  "--genFile",
  type = "character",
  help = "Full path to the genotype files until \"chr\""
)

parser$add_argument(
  "--dosageThresh",
  type = "double",
  default = 0.1,
  help = "Dosage below this value is set to 0"
)

parser$add_argument(
  "--outDosageFold",
  type = "character",
  help = "Folder to save dosage files"
)


args <- parser$parse_args()


gen_name <- basename(args$genFile)


for (chr in 1:22) {
  gen_file <- data.table::fread(sprintf("%schr%s.gen.gz", args$genFile, chr), header = FALSE)

  sample_file <- read.table(
    sprintf("%schr%s.sample", args$genFile, chr),
    stringsAsFactors = FALSE,
    row.names = NULL
  )


  gen_file <- gen_file[, 6:ncol(gen_file)]

  gen_file <- gen_file[, seq(3, ncol(gen_file), 3), with = FALSE] * 2 +
              gen_file[, seq(2, ncol(gen_file), 3), with = FALSE]

  gen_file <- round(gen_file, 2)

  gen_file[gen_file < args$dosageThresh] <- 0

  colnames(gen_file) <- sample_file$V2[3:nrow(sample_file)]


  out_name <- sprintf("%s/%sdosage_chr%s_matrix.txt.gz", args$outDosageFold, gen_name, chr)

  data.table::fwrite(gen_file, out_name, sep = "\t", compress = "gzip")
}
