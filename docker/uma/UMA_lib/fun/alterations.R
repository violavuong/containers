#!/usr/bin/r

# file: alterations.R
# description: MM alterations-specific custom function
# last update: 25-03-25


bcnCNA <- function(broad_calls, IDs, arm = c("1p", "1q", "17p"), alt = c("amp", "del")) {
  alt_df <- broad_calls %>% filter(chrarm==arm & call==alt)
  alt_dt <- rbindlist(sapply(IDs, function(x) alt_df %>% filter(ID==x) %>% mutate(alt = ((cwR_CN + CNV_CN) / 2) %>% round(3)), simplify = FALSE)) %>%
    complete(ID = IDs) %>%
    select(ID, alt)
  
  return(alt_dt)
}


broadConcordance <- function(GR_df, ID, broad_quality, broad_th, th){
  dj <- disjoin(GR_df, with.revmap = T)
  
  # extracting DJ info
  tools_stats <- data.frame()
  for (i in seq_along(dj$revmap)){
    stat <- extractDJInfo(GR_df, dj$revmap[[i]])
    tools_stats <- rbind(tools_stats, stat)
  }
  
  dj$revmap <- sapply(dj$revmap, function(x) paste0(x, collapse = ","))
  
  # define concordance between tools
  conc_df <- as.data.frame(dj) %>%
    mutate(ID, 
           cwR_CN = tools_stats$cwR_CN, cwR_probes = tools_stats$cwR_probes, 
           CNV_CN = tools_stats$CNV_CN, CNV_probes = tools_stats$CNV_probes, 
           broad_quality, broad_th) %>%
    filter(width > 10) %>%
    mutate(CN_diff = abs(cwR_CN - CNV_CN), 
           conc = case_when(cwR_CN > 2 + th & CNV_CN > 2 + th ~ "conc_AMP",
                            cwR_CN < 2 - th & CNV_CN < 2 - th ~ "conc_DEL", 
                            cwR_CN < 2 + th & cwR_CN > 2 - th & CNV_CN < 2 + th & CNV_CN > 2 - th ~ "conc_NORMAL",
                            is.na(CNV_CN) | is.na(cwR_CN) ~ NA, 
                            TRUE ~ "discordant")) %>%
    rename_with(~"chr", seqnames)
  
  conc_df <- Popeye2(conc_df, refGen, removeXY = T) %>%
    relocate(CN_diff, conc, .before = broad_quality)
  
  return(conc_df)
}


broadStats <- function(df, sample){
  df <- df %>% filter(ID==sample)
  
  discordant_mb <- df %>% 
    filter(conc=="discordant") %>% 
    pull(width) %>% sum() %>% `/`(1000000) %>% round(3)
  
  weighted_mean <- weighted.mean(df$CN_diff, df$width, na.rm = T)
  
  return(data.frame(ID = sample, 
                    discordant_mb = discordant_mb,
                    weighted_mean = weighted_mean))
}


correctSegments <- function(df, type = c("broad", "focal"), quality){
  output <- list()
  
  # diploid region correction
  if(type=="broad"){
    df_bob <- DRrefit2(Popeye2(df, refGen, removeXY = F), IDs, normChr)
    df_bob_segs <- df_bob$corrected_segments
    output[["bobafit"]] <- df_bob
  } else if (type=="focal"){
    df_tmp_bob <- Popeye2(df, refGen, removeXY = F)
    
    norm_report <- data.frame()
    for (i in seq_along(IDs)){
      norm_exons <- df_tmp_bob %>% filter(ID==IDs[i], chrarm %in% normChr)
      norm_exons_CN <- median(norm_exons$CN)
      CF <- 2 - norm_exons_CN %>% round(3)
      class <- ifelse(abs(CF) > 0.5, "refitted", ifelse(abs(CF) > 0.1, "recalibrated", "no_changes"))
      
      report <- data.frame(ID = IDs[i],
                           chr_list = paste0(normChr, collapse = ","), 
                           correction_factor = CF, 
                           correction_class = class, 
                           row.names = NULL)
      norm_report <- bind_rows(norm_report, report)
    }
    focal_CFs <- norm_report %>% select(ID, correction_factor)
    df_bob_segs <- merge(df_tmp_bob, focal_CFs, by = "ID") %>% mutate(CN_corrected = CN + correction_factor)
    output[["bobafit"]] <- list(corrected_segments = df_bob_segs, report = norm_report)
  } else {message("Incorrect option.\n")}
  
  # sample purity correction
  df_purity <- merge(df_bob_segs, quality, by = "ID") %>%
    mutate(CN_purity_corrected = ifelse(purity > purityTh, ((CN_corrected - 2)/ purity + 2), CN_corrected), .after = purity)
  output[["purity"]] <- df_purity
  
  return(output)
}


computeFocalCN <- function(df, ID_name){
  return(df %>%
           filter(ID==ID_name) %>% 
           mutate(adj_mean_cov  = as.numeric(mean_coverage) / NORM_FACT, 
                  adj_ratio_cov = adj_mean_cov / median(adj_mean_cov), 
                  logR = log2(adj_ratio_cov), 
                  CN = 2^(logR + 1)))
}


createGenes <- function(df, col_name){
  return(df %>%
           group_by(ID, gene, chr) %>%
           summarise(start = min(start), end = max(end),
                     mean_CN = mean(!!sym(col_name)), median_CN = median(!!sym(col_name)),
                     Q1_CN = quantile(!!sym(col_name), probs = 0.25, na.rm = T), Q3_CN = quantile(!!sym(col_name), probs = 0.75, na.rm = T),
                     CI95_low_CN = smean.cl.boot(!!sym(col_name), conf.int = 0.95, B = 1000) %>% .["Lower"], 
                     CI95_up_CN = smean.cl.boot(!!sym(col_name), conf.int = 0.95, B = 1000) %>% .["Upper"]) %>%
           as.data.frame())
}


extractDJInfo <- function(GR_df, k){
  if(length(k)==2){
    cwR_CN <- GR_df$CN_purity_corrected[k][1]
    CNV_CN <- GR_df$CN_purity_corrected[k][2]
    cwR_probes <- GR_df$probes[k][1]
    CNV_probes <- GR_df$probes[k][2]
  } else if(length(k)==1){
    tool <- GR_df$tool[k]
    try(
      if(tool=="copywriteR"){
        cwR_CN <- GR_df$CN_purity_corrected[k]
        cwR_probes <- GR_df$probes[k]
        CNV_CN <- CNV_probes <- NA
      } else if(tool=="CNVkit"){
        CNV_CN <- GR_df$CN_purity_corrected[k]
        CNV_probes <- GR_df$probes[k]
        cwR_CN <- cwR_probes <- NA
      }
    )
  } else {cwR_CN <- CNV_CN <- cwR_probes <- CNV_probes <- NA}
  
  return(data.frame(cwR_CN, cwR_probes, CNV_CN, CNV_probes))
}


filterForChr14 <- function(df){
  df$CHROM <- gsub("chr", "", df$CHROM)
  df$CHR2 <- gsub("chr", "", df$CHR2)
  
  IgH_df <- df %>% filter((CHROM==14 & !is.na(CHR2)) | (CHR2==14 & !is.na(CHROM)))
  
  IgH_df$POS_partner <- ifelse(IgH_df$CHROM==14, IgH_df$POS2, IgH_df$POS)
  IgH_df$CHROM_partner <- ifelse(IgH_df$CHROM==14, IgH_df$CHR2, IgH_df$CHROM)
  
  IgH_df$POS_IgH <- ifelse(IgH_df$CHROM==14, IgH_df$POS, IgH_df$POS2)
  IgH_df$CHROM_IgH <- ifelse(IgH_df$CHROM==14, IgH_df$CHROM, IgH_df$CHR2)
  
  return(IgH_df)
}


filterForTargets <- function(IgH_df, refGen, partners, IgHLocus, targetGene){
  ## gene loci taken from NCBI
  if (refGen=="hg19") {
    IgH_df <- IgH_df %>% filter(POS_IgH > 106052774 - IgHLocus & POS_IgH < 107288051 + IgHLocus)
    
    mean_chr4 <- mean(mean(1795020, 1810594), mean(1873120, 1983919)) #FGFR3, NSD2, t(4;14)
    mean_chr6 <- mean(41902671, 42016632) #CCND3, t(6;14)
    mean_chr11 <- mean(69455924, 69469242) #CCND1, t(11;14)
    mean_chr16 <- mean(79627735, 79634634) #MAF, t(14;16)
    mean_chr20 <- mean(39314488, 39317876) #MAFB, t(14;20)
    
    means <- c(mean_chr4, mean_chr6, mean_chr11, mean_chr16, mean_chr20)
  } else if (refGen=="hg38") {
    IgH_df <- IgH_df %>% filter(POS_IgH > 105586437 - IgHLocus & POS_IgH < 106879844 + IgHLocus)
    
    mean_chr4 <- mean(mean(1793293, 1808872), mean(1871393, 1982207)) #FGFR3, NSD2, t(4;14)
    mean_chr6 <- mean(41934933, 42050357) #CCND3, t(6;14)
    mean_chr11 <- mean(69641156, 69654474) #CCND1,t(11;14)
    mean_chr16 <- mean(79202622, 79600737) #MAF, t(14;16)
    mean_chr20 <- mean(40685848, 40689236) #MAFB, t(14;20)
    
    means <- c(mean_chr4, mean_chr6, mean_chr11, mean_chr16, mean_chr20)
  } else {
    message("Incorrect reference genome.")
  }
  
  # filter translocations based on target genes
  complete_df <- data.frame()
  
  for (i in seq_along(partners)){
    trans_df <- IgH_df %>% filter(CHROM_partner==partners[i] & (abs(POS_partner - means[i]) < targetGene))
    complete_df <- rbind(trans_df, complete_df)
  }
  
  complete_df <- complete_df %>%
    mutate(t_IgH = ifelse(as.numeric(CHROM_IgH) < as.numeric(CHROM_partner), 
                          paste0("t(",CHROM_IgH,";",CHROM_partner,")"), 
                          paste0("t(",CHROM_partner,";",CHROM_IgH,")")))
  
  return(complete_df)
}

filterSNV <- function(snv_df, tool = c("Mutect2", "Strelka")){
  snv_df <- data.frame(lapply(snv_df, function(x) gsub("\\\\x3b", ";", x)))
  snv_df <- data.frame(lapply(snv_df, function(x) gsub("\\\\x3d", "=", x)))
  
  if (tool=="Mutect2"){
    snv_tmp_df <- snv_df %>% mutate(VAF = as.numeric(gt_AF), DP = as.numeric(DP))
  }
  else if (tool=="Strelka"){
    snv_df$gt_AF <- sub('.*,\\s*','', snv_df$gt_AD)
    snv_tmp_df <- snv_df %>% 
      filter(gt_FT=="PASS") %>%
      mutate(VAF = as.integer(gt_AF)/as.integer(gt_DP), DP = as.integer(gt_DP)) %>% 
      rename_with(~"fathmm.MKL_coding_pred", 67)
  }
    
  snv_filtered_df <- snv_tmp_df %>%
    filter(!Func.refGene %in% c("intergenic", "intronic", "upstream", "downstream","upstream;downstream", "UTR3", "UTR5", "ncRNA_intronic"), 
           !ExonicFunc.refGene %in% c("synonymous_SNV"), #"." 
           VAF >= 0.05 & VAF < 1) %>%
      select(Indiv, CHROM, POS, REF, ALT, DP, FILTER, Func.refGene, Gene.refGene, ExonicFunc.refGene, AAChange.refGene, cytoBand, gnomAD_genome_ALL, 
             avsnp150, SIFT_pred, SIFT4G_pred, Polyphen2_HDIV_pred, Polyphen2_HVAR_pred, LRT_pred, MutationTaster_pred, MutationAssessor_pred, fathmm.MKL_coding_pred, 
             PROVEAN_pred, VEST4_score, MetaSVM_pred, MetaLR_pred, CLNSIG, cosmic70, InterVar_automated, variant_name, VAF)
  
  return(snv_filtered_df)
}


reportIgH <- function(IgH_df, IDs){
  IgHs <- list("4;14", "6;14", "11;14", "14;16", "14;20")
  
  trans <- sapply(IDs, function(x) IgH_df %>% filter(Indiv==x) %>% pull(t_IgH) %>% paste0(collapse="/"))
  call_IDs <- sapply(IDs, function(x) IgH_df %>% filter(Indiv==x) %>% pull(call_ID) %>% paste0(collapse="/"))
  type_IDs <- sapply(IDs, function(x) IgH_df %>% filter(Indiv==x) %>% pull(ID) %>% paste0(collapse="/"))
  
  calls <- sapply(IgHs, function(x) ifelse(trans %>% str_detect(x), 1, 0))
  
  report_df <- data.frame(IDs, trans, call_IDs, type_IDs, calls, row.names = NULL)
  colnames(report_df) <- c("ID", "t_IgH", "call", "type", "t_IgH_4_14", "t_IgH_6_14", "t_IgH_11_14", "t_IgH_14_16", "t_IgH_14_20")
  
  return(report_df)
}

validateSNV <- function(filtered_df, fb_df, var_df, tool = c("Mutect2", "Strelka")){
  df_fb_tmp <- semi_join(filtered_df, fb_df, by = "variant_name")
  df_var_tmp <- semi_join(filtered_df, var_df, by = "variant_name")
  
  merged_df <- rbind(df_fb_tmp, df_var_tmp) %>%
    distinct(.keep_all = T) %>%
    mutate(callers = ifelse(variant_name %in% df_fb_tmp$variant_name & variant_name %in% df_var_tmp$variant_name, paste0(tool, ", FreeBayes, VarScan2"),
                            ifelse(variant_name %in% df_fb_tmp$variant_name, paste0(tool, ", FreeBayes"), paste0(tool, ", VarScan2"))))
  
  return(merged_df)
}
