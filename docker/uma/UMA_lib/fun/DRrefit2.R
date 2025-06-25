#!/usr/bin/r

DRrefit2 <- function(segs_df, IDs, normChr, maxCN = 6, cluster = "ward.D2"){
  BOB_out <- list()
  segs <- reports <- data.frame()

  segs_df$CN[segs_df$CN == Inf] <- maxCN
  
  for(i in seq_along(IDs)){
    #compute CN weighted mean 
    seg_df <- segs_df %>% filter(ID==IDs[i])
    
    weighted_CN_df <- as.data.frame(seg_df %>%
                                      group_by(chrarm) %>%
                                      summarise(weighted_mean_CN = weighted.mean(CN, width), .groups = "drop"))
  
    #diploid region correction
    diploid_reg <- median(weighted_CN_df$weighted_mean_CN)
    CF <- 2 - diploid_reg
  
    seg_df$CN_corrected <- seg_df$CN + CF
    seg_df$CN_corrected <- ifelse(seg_df$CN_corrected < 0, 0.001, seg_df$CN_corrected)
    
    #clusterization
    test <- try({
      invisible(capture.output(cluster_res <- NbClust(weighted_CN_df$weighted_mean_CN, distance = "euclidean", method = cluster, index = "all", min.nc = 2, max.nc = 6), 
                               type = c("output", "message")))
      
      cluster_tb <- data.frame(chr = weighted_CN_df$chrarm, 
                               cluster = cluster_res$Best.partition,
                               stringsAsFactors = FALSE)
      
      cluster_df <- cluster_tb %>%
        filter(chr %in% normChr) %>%
        group_by(cluster) %>%
        summarize(num = dplyr::n()) %>%
        arrange(desc(num))
      
      ref_cluster <- cluster_df$cluster[1]
      ref_df <- cluster_tb %>% filter(cluster %in% ref_cluster)
      cluster_normChr <- ref_df$chr
      
    }, silent = TRUE)
  
    #report
    if (is(test, "try-error")){
      report <- data.frame(ID = IDs[i],
                           clustering = "FAIL",
                           ref_chrs = paste0(normChr, collapse = ","), 
                           n_cluster = NA, 
                           correction_factor = CF, 
                           correction_class = ifelse(abs(CF) > 0.5, "REFITTED", ifelse(abs(CF) <= 0.1, "NO CHANGES", "RECALIBRATED")))
    } else {
      report <- data.frame(ID = IDs[i],
                           clustering = "SUCCESS",
                           ref_chrs = paste0(cluster_normChr, collapse = ","), 
                           n_cluster = max(cluster_res$Best.partition), 
                           correction_factor = CF, 
                           correction_class = ifelse(abs(CF) > 0.5, "REFITTED", ifelse(abs(CF) <= 0.1, "NO CHANGES", "RECALIBRATED")))
    }
    
    #binding
    segs <- bind_rows(segs, seg_df)
    reports <- bind_rows(reports, report)
  }
  
  #output
  BOB_out[["corrected_segments"]] <- segs
  BOB_out[["report"]] <- reports
  
  return(BOB_out)
}
