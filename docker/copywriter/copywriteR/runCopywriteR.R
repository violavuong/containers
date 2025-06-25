#!/usr/bin/r

# file: runCopywriteR.R
# aim: Detect and analyze CNAs
# version: 0.1
# last update: 31-10-2024

options(warn = -1)

suppressWarnings(suppressMessages(library(optparse)))

# ---- Args definition from the command line ----
options <- list(
  make_option(c("--libDir", "-l"), type = "character", default = "/cwR_lib/", help = "Path to copywriteR library directory."),
  make_option(c("--bamDir", "-b"), type = "character", help = "Path to bam files directory."),
  make_option(c("--normDir", "-n"), type = "character", help = "Path to normal file directory."),
  make_option(c("--refDir", "-r"), type = "character", help = "Path to reference file"), 
  make_option(c("--outDir", "-o"), type = "character", help = "Path to output directory.")
)

parseobj <- OptionParser(option_list=options)
opt <- parse_args(parseobj)
print(opt)
options(scipen=0, stringsAsFactors=F)

# ---- Args definition in the script ----
libDir <- opt$libDir
bamDir <- opt$bamDir
normDir <- opt$normDir
refDir <- opt$refDir
outDir <- opt$outDir


# ---- Main ----
suppressWarnings(suppressMessages(library(CopywriteR)))

## creating output directory
dir.create(paste0(outDir), recursive = TRUE)


## input files
bams <- list.files(bamDir, pattern = ".bam$", full.names = TRUE)
normal <- file.path(paste0(normDir))
samples <- c(bams, normal)

ctrl <- rep(normal, length(samples))
sample_ctrl <- data.frame(samples, ctrl)


## run copywriteR
CopywriteR(sample.control = sample_ctrl,
           destination.folder = file.path(outDir),
           reference.folder = file.path(paste0(libDir, refDir)),
           bp.param = SnowParam(workers = 8, type = "SOCK"))

## plotting
plotCNA(destination.folder = file.path(outDir))
