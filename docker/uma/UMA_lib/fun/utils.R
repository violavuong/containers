#!/usr/bin/r

# file: utils.R
# description: data manipulation functions
# last update: 23-09-24


createGR <- function(df, sample, tool_name){
  filtered_df <- df %>% 
    filter(ID==sample) %>% 
    mutate(tool = paste0(tool_name))
  return(makeGRangesFromDataFrame(filtered_df, keep.extra.columns = T))
}


createVariantID <- function(df){
  df$CHROM <- gsub("chr", "", df$CHROM)
  return(df %>% mutate(variant_name = paste(DNA, CHROM, POS, sep = "_")))
}


createVcf <- function(files, filters = NULL){
  VCFs <- lapply(files, function(x) read.vcfR(x, verbose = F)) 
  tidyVCFs <- lapply(VCFs, function(x) vcfR2tidy(x, single_frame = T, dot_is_NA = T, verbose = F))
  VCF_dfs <- lapply(seq_along(tidyVCFs), function(x) tidyVCFs[[x]]$dat) %>%
    Map(cbind, ., DNA
        = str_extract(files, "[\\d]+_S[\\d]+") %>% str_remove("_S[\\d]+"))
  
  # filtering
  if (is.null(filters)) {
    filtered_VCF <- rbindlist(VCF_dfs)
  } else if (!is.null(filters)) {
    VCF_df <- rbindlist(VCF_dfs)
    filtered_VCF <- rbindlist(lapply(filters, function(x) VCF_df %>% filter(FILTER==x)))
  } else {
    message("Incorrect option.\n")
  }
  
  return(filtered_VCF)
}


detectIgH <- function(df, col_name){
  return(sapply(strsplit(df[, col_name], "/"), function(x) paste(mixedsort(unique(x)), collapse = ",")))
}


searchFiles <- function(wd, regex){
  return(list.files(path = wd, pattern = regex, recursive = T, full.names = T))
}

