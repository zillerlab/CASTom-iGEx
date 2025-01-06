#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(argparse)
  library(data.table)
})


parser <- argparse::ArgumentParser()

parser$add_argument(
  "--trawFile",
  type = "character",
  help = "Full path to the traw files until \"chr\""
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
  default = ".",
  help = "Folder to save dosage files"
)

args <- parser$parse_args()

gen_name <- basename(args$trawFile)

for (chr in 1:22) {
  gen_file <- data.table::fread(sprintf("%schr%s.traw", args$trawFile, chr))

  gen_file <- gen_file[, 7:ncol(gen_file)]

  gen_file <- round(gen_file, 2)

  gen_file[gen_file < args$dosageThresh] <- 0

  out_name <- sprintf("%s/%schr%s_matrix.txt.gz", args$outDosageFold, gen_name, chr)

  data.table::fwrite(gen_file, out_name, sep = "\t", compress = "gzip")
}
