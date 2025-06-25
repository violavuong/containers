
#' @title Popeye2
#'
#' @description
#' This function annotates an input SEG file, adding the p or q chromosome to every segment in the seg file. It support both hg19 and hg38 genome version.
#' Differently from Popeye (version 1 from BOBaFIT) it can be used also if some chromosomes are missing.
#'
#' @param segments A data.frame representing a SEG file, mandatory columns: "chr", "start", "end". Chromosomes must be numeric (e.g. 1,2,3...) NOT charachters (e.g. chr1, chr2, chr3..)
#' @param refGenome reference genome. Currently supported "hg19" (default) and "hg38"
#' @param removeXY Logical. Option to remove chromosome X and Y (chromosome 23 and 24)
#'
#' @return
#' A data.frame with the same columns as the input data.frame, plus two new columns: "chrarm" and "arm".
#'
#' @export
#'
#' @examples
#'

Popeye2 <- function(segments, refGen = "hg19", removeXY = TRUE) {

  options(scipen = 999)

  # load chr arm table with arms start/end info - hg19 or hg38
  if(refGen=="hg19") {
    # REF file hg19: http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/cytoBandIdeo.txt.gz
    bands <- read.delim(paste0(libDir, "inst/extdata/cytoBandIdeo_hg19.txt"))
    } else if (refGen=="hg38") {
    # REF file hg38: http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cytoBandIdeo.txt.gz
    bands <- read.delim(paste0(libDir, "inst/extdata/cytoBandIdeo_hg38.txt"))
    } else { stop("Invalid reference genome. Please choose hg19 or hg38")}

  names(bands) <- c("chr", "start", "end", "band", "giemsa")
  bands$arm <- str_extract(bands$band, "[pq]")
  bands$chrarm <- paste0(str_remove(bands$chr, "chr"), bands$arm)
  bands2 <- bands %>% filter(!is.na(arm))
  
  chrtab <- bands2 %>% group_by(chrarm) %>%
    summarise(chr = unique(chr) %>% str_remove("chr"),
              start = min(start),
              end = max(end),
              arm = str_extract(chrarm, "[pq]") %>% unique,
              ref = refGen) %>%
    select(chr, start, end, chrarm, arm, ref)

  chrtabGR <- makeGRangesFromDataFrame(chrtab, keep.extra.columns = TRUE)

  # optional: remove X Y chr
  if(removeXY==TRUE){
    segments <- segments %>% filter(!chr %in% c("X", "Y", "23", "24"))
  }

  segmentsGR <- makeGRangesFromDataFrame(segments, keep.extra.columns = TRUE)
  
  overlaps <- join_overlap_intersect(segmentsGR, chrtabGR)
  
  overlaps <- as.data.frame(overlaps) %>%
    select(-c(arm, ref, strand), chr = seqnames) %>%
    relocate(ID, .before = chr) %>%
    relocate(chrarm, .after = chr)
  
  return(as.data.frame(overlaps))
}
