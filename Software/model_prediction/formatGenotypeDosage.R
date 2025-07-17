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
  "--sampleFile",
  type = "character",
  help = "Full path to the sample file, including the suffix"
)

parser$add_argument(
  "--sampleNameColumn",
  type = "integer",
  default = NULL,
  help = "Which column in the sample file to use for column names. Default is a concatenation of FID (1) and IID (2)"
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

sample_file <- read.delim(
  args$sampleFile,
  header = FALSE,
  comment.char = "#",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

gen_name <- basename(args$trawFile)

for (chr in 1:22) {
  gen_file <- data.table::fread(sprintf("%schr%s.traw", args$trawFile, chr))

  gen_file[, (1:6) := NULL]

  gen_file[gen_file < args$dosageThresh] <- 0

  gen_file <- round(gen_file, 2)

  if (is.null(args$sampleNameColumn)) {
    colnames(gen_file) <- paste(sample_file[[1]], sample_file[[2]], sep = "_")
  } else {
    colnames(gen_file) <- sample_file[[args$sampleNameColumn]]
  }

  out_name <- sprintf("%s/%schr%s_matrix.txt.gz", args$outDosageFold, gen_name, chr)

  data.table::fwrite(gen_file, out_name, sep = "\t", compress = "gzip")
}
