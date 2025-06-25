#!/usr/bin/r

# file: runUMAPost.R
# aim: Running the entire post-processing analyses of UMA panel experiments
# version: 0.1
# last update: 25-03-25

library(optparse)


# ---- Args definition from the command line ----
#to describe in detail in a separated doc/READ.ME file each option specific format
options <- list(
  make_option(c("--wd"), type = "character", help = "Path to working directory"),
  make_option(c("--libDir", "-l"), type = "character", default = "/UMA_lib/", help = "Path to UMA library directory"),
  #make_option(c("--runID"), type = "character", help = "Run ID"),
  make_option(c("--refGen", "-r"), type = "character", default = "hg19", help = "Reference genome. Default: [%default]"),
  make_option(c("--broadLow"), type = "numeric", default = 2e6, help = "Broad CNAs reads lower limit. Default: [%default]"),
  make_option(c("--broadHigh"), type = "numeric", default = 3e6, help = "Broad CNAs reads upper limit. Default: [%default]"),
  make_option(c("--focalLow"), type = "numeric", default = 1e6, help = "Focal CNAs reads lower limit. Default: [%default]"),
  make_option(c("--focalHigh"), type = "numeric", default = 2.25e6, help = "Focal CNAs reads upper limit. Default: [%default]"),
  make_option(c("--percOff"), type = "numeric", default = 0.3, help = "Percentage of off-target reads cut-off. Default: [%default]"), 
  make_option(c("--MAD"), type = "numeric", default = 0.15, help = "MAD value cut-off. Default: [%default]"),
  make_option(c("--totReads"), type = "numeric", default = 4e6, help = "Total reads cut-off. Default: [%default]"),
  make_option(c("--normChr"), type = "character", default = "c(\"1p\",\"2p\",\"2q\",\"4p\",\"4q\",\"6p\",\"8p\",\"10p\",\"12p\",\"12q\",\"16p\",\"17p\",\"17q\",\"18p\",\"18q\",\"20p\",\"20q\",\"22q\")", help = "Normal chromosomes. Default: [%default]"),
  make_option(c("--clnTh"), type = "numeric", default = 0.2, help = "Broad CNA calling threshold. Default: [%default]"),
  make_option(c("--purity"), type = "numeric", default = 0.3, help = "Sample purity CN value correction threshold. Default: [%default]"),
  make_option(c("--IgH"), type = "numeric", default = 1e6, help = "IgH locus length. Default: [%default]"),
  make_option(c("--target"), type = "numeric", default = 5e6, help = "Maximum distance from the target gene supported. Default: [%default]"),
  make_option(c("--partners"), type = "character", default = "c(4,6,11,16,20)", help = "Partner chromosomes of chromosome 14. Default: [%default]"),
  make_option(c("--outDir", "-o"), type = "character", default = "/post_processing/", help = "Path to output directory. Default: [%default]")
)

parseobj <- OptionParser(option_list=options)
opt <- parse_args(parseobj)
print(opt)
options(scipen = 0, stringsAsFactors = F)

# ---- Args definition in the script ----
wd <- opt$wd 
libDir <- opt$libDir
#runID <- opt$runID
refGen <- opt$refGen
broadLow <- opt$broadLow
broadHigh <- opt$broadHigh
focalLow <- opt$focalLow
focalHigh <- opt$focalHigh
percOff <- opt$percOff
MAD <- opt$MAD
totReads <- opt$totReads
normChr <- as.character(eval(parse(text = opt$normChr)))
clnTh <- opt$clnTh
purityTh <- opt$purity
IgHLocus <- opt$IgH
targetGene <- opt$target
partners <- eval(parse(text = opt$partners))
outDir <- opt$outDir


# ---- Main ----
options(warn = -1)
options(dplyr.summarise.inform = FALSE)

suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(GenomicRanges)))
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(gtools)))
suppressWarnings(suppressMessages(library(Hmisc)))
suppressWarnings(suppressMessages(library(mgsub)))
suppressWarnings(suppressMessages(library(NbClust)))
suppressWarnings(suppressMessages(library(paletteer)))
suppressWarnings(suppressMessages(library(plyranges)))
suppressWarnings(suppressMessages(library(readr)))
suppressWarnings(suppressMessages(library(readxl)))
suppressWarnings(suppressMessages(library(reshape2)))
suppressWarnings(suppressMessages(library(stringr)))
suppressWarnings(suppressMessages(library(tidyr)))
suppressWarnings(suppressMessages(library(vcfR)))


## loading custom function
source(paste0(libDir, "fun/alterations.R"))
source(paste0(libDir, "fun/DRrefit2.R"))
source(paste0(libDir, "fun/plot.R"))
source(paste0(libDir, "fun/Popeye2.R"))
source(paste0(libDir, "fun/utils.R"))
source(paste0(libDir, "fun/VarianThinker.R"))


## creating output directory
dir.create(paste0(outDir), recursive = TRUE)

## defining global variables
regex <- list(hs = "HsMetrics_\\d+\\.txt", cnv = ".call.cns", exons = "targetCoverage", 
              delly = "delly.*vcf", manta = "manta.*vcf",
              mut = "mutect2.*multianno.vcf", str = "strelka.*multianno.vcf", var2 = "varscan.*vcf", fb = "^freebayes.*vcf")

filters <- list(delly = list("PASS"), manta = list("PASS", "MaxDepth", "MaxDepth;MaxMQ0Frac"),
                mut = list("PASS", "germline", "clustered_events", "clustered_events;germline"))

# ---- Sequencing metrics ----
message("Sequencing metrics extraction.\n")
dir.create(paste0(outDir, "metrics/"), recursive = TRUE)


## cwR on- and off-target metrics
cwR_off_on <- fread(file = paste0(wd, "CopywriteR/CNAprofiles/CopywriteR.log"), skip = "off.target") %>% 
  select(-starts_with("unmappable")) %>%
  filter(V1 != "remDup_recal_DONOR.bam") %>%
  rename_with(~"ID", 1) %>%
  mutate(ID = mgsub(ID, c("remDup_recal_", "\\.bam"), c("", "")))
colnames(cwR_off_on) <- gsub("\\.", "_", colnames(cwR_off_on))

IDs <- cwR_off_on$ID #sample IDs

## MAD
cwR_MAD <- fread(file = paste0(wd, "CopywriteR/CNAprofiles/CopywriteR.log"), skip = "MAD", select = c(8, 10)) %>%
  rename_with(~c("ID", "MAD"), c(1, 2)) %>% 
  filter(!(grepl("none|DONOR.bam.vs", ID))) %>%
  mutate(ID = mgsub(ID, c("remDup_recal_", "\\.vs.*", "log2.", "\\.bam"), c("", "", "", "")))

## merging cwR metrics 
cwR_metrics <- left_join(cwR_MAD, cwR_off_on, by = "ID")
cwR_metrics$tot_reads <- rowSums(cwR_metrics[, c("off_target", "on_target")]) #add tot_reads
cwR_metrics$perc_off_target <- (cwR_metrics$off_target / cwR_metrics$tot_reads) %>% round(3) #add percentage of off-target


## HSmetrics coverage metrics
HS_files <- searchFiles(wd, regex$hs)
coverage <- rbindlist(lapply(seq_along(HS_files), function(x) fread(HS_files[[x]], skip = "METRICS CLASS", select = c("MEAN_BAIT_COVERAGE", "MEAN_TARGET_COVERAGE"), dec = "."))) %>%
  mutate(ID = IDs, .before = MEAN_BAIT_COVERAGE)
colnames(coverage)[2:3] <- tolower(colnames(coverage)[2:3])

## merging cwR and HSmetrics
metrics_df <- left_join(cwR_metrics, coverage, by = "ID") %>%
  mutate(DNA = ID %>% str_remove("_.*"), .after = ID)

## purity
purity_df <- suppressMessages(read_xlsx(paste0(wd, "PERCENTUALI POST ARRICCHIMENTO.xlsx")))[, 2:3] %>%
  rename_with(~c("DNA", "purity"), c(n_DNA, `%POST`)) %>%
  mutate(DNA = as.character(DNA), 
         purity = purity / 100)

quality_df <- merge(metrics_df, purity_df, by = "DNA") %>%
  relocate(DNA, purity, .after = ID)

## sample quality and dynamic threshold
quality_df <- quality_df %>%
  mutate(broad_quality = case_when(off_target < broadLow ~ "Low",
                                   off_target < broadHigh ~ "Medium", 
                                   off_target >= broadHigh ~ "High", 
                                   TRUE ~ NA), 
         focal_quality = case_when(on_target < focalLow ~ "Low",
                                   on_target < focalHigh ~ "Medium", 
                                   on_target >= focalHigh ~ "High", 
                                   TRUE ~ NA), 
         broad_th = case_when(broad_quality == "Low" ~ 0.4,
                              broad_quality == "Medium" ~ 0.3,
                              broad_quality == "High" ~ 0.25, 
                              TRUE ~ NA), 
         focal_th = case_when(focal_quality == "Low" ~ 0.4,
                              focal_quality == "Medium" ~ 0.3,
                              focal_quality == "High" ~ 0.2, 
                              TRUE ~ NA))

write_tsv(quality_df, paste0(outDir, "metrics/quality_stats.txt")) #saving


## plotting
cut_offs <- data.frame(variable = c("MAD", "off_target", "off_target", "on_target", "on_target", "tot_reads", "perc_off_target"), 
                       value = c(MAD, broadLow, broadHigh, focalLow, focalHigh, totReads, percOff))

metrics_plot <- plotMetrics(metrics_df %>% select(-ends_with("coverage"), -DNA), cut_offs)
ggsave(plot = metrics_plot, file = paste0(outDir, "metrics/seq_metrics.svg"), width = 8, height = 15, dpi = 300) #saving

message("Sequencing metrics extraction: done!\n\n")


# ---- CNA calling algorithm ----
message("CNAs calling algorithm.\n")
dir.create(paste0(outDir, "CNAs/raw_calls/"), recursive = TRUE)
dir.create(paste0(outDir, "CNAs/refitted_calls/"), recursive = TRUE)


### ANALYZING copywriteR ###
message("Running tool: copywriteR.\n")

## cwR raw segments
load(paste0(wd, "CopywriteR/CNAprofiles/segment.Rdata"))
cwR_raw_segs <- segment.CNA.object$output %>%
  filter(!(grepl("none|DONOR.bam.vs.log2.remDup_recal_DONOR", ID))) %>%
  mutate(ID = str_remove_all(ID, paste(c("log2.","\\.bam.*", "remDup_recal_"), collapse = "|")), 
         width = loc.end - loc.start, .after = loc.end) %>%
  mutate(CN = 2^(seg.mean + 1), .after = seg.mean, 
         chrom = recode(as.character(chrom), "23"="X", "24"="Y")) %>%
  rename_with(~c("chr", "start", "end", "probes"), c(chrom, loc.start, loc.end, num.mark))

write_tsv(cwR_raw_segs, paste0(outDir, "CNAs/raw_calls/cwR_raw_segs.txt")) #saving

## diploid region and purity correction
broad_df <- quality_df %>% select(ID, purity, broad_quality, broad_th)
cwR_corrected_list <- correctSegments(cwR_raw_segs, "broad", broad_df)

write_tsv(cwR_corrected_list$bobafit$corrected_segments, paste0(outDir, "CNAs/refitted_calls/cwR_bobafit_segs.txt")) #saving bobafit corrected segments
write_tsv(cwR_corrected_list$bobafit$report, paste0(outDir, "CNAs/refitted_calls/cwR_bobafit_report.txt")) #saving bobafit report
write_tsv(cwR_corrected_list$purity, paste0(outDir, "CNAs/refitted_calls/cwR_corrected_segs.txt")) #saving purity-corrected segments


### ANALYZING CNVkit ###
message("Running tool: CNVkit.\n")

## CNVkit raw segments
CNV_files <- searchFiles(wd, regex$cnv)
CNV_raw_segs <- rbindlist(lapply(CNV_files, fread) %>%
  Map(cbind, ., ID = str_extract(CNV_files, "remDup.*") %>% str_remove_all("\\..*|remDup_recal_"))) %>%
  mutate(width = end - start, .after = end) %>%
  mutate(CN = 2^(log2 + 1)) %>%
  select(-c(gene, cn, depth, p_ttest, weight)) %>%
  rename_with(~"chr", chromosome) %>%
  relocate(ID, .before = chr) %>%
  relocate(log2, .after = probes)
CNV_raw_segs$chr <- gsub("chr", "", CNV_raw_segs$chr)
  
write_tsv(CNV_raw_segs, paste0(outDir, "CNAs/raw_calls/CNV_raw_segs.txt")) #saving

## diploid region and purity correction
CNV_corrected_list <- correctSegments(CNV_raw_segs, "broad", broad_df)

write_tsv(CNV_corrected_list$bobafit$corrected_segments, paste0(outDir, "CNAs/refitted_calls/CNV_bobafit_segs.txt")) #saving bobafit corrected segments
write_tsv(CNV_corrected_list$bobafit$report, paste0(outDir, "CNAs/refitted_calls/CNV_bobafit_report.txt")) #saving bobafit report
write_tsv(CNV_corrected_list$purity, paste0(outDir, "CNAs/refitted_calls/CNV_corrected_segs.txt")) #saving purity-corrected segments


### CHECKING FOR TOOLS CONCORDANCE - BROAD CNAs ###
cwR_GR_dfs <- sapply(IDs, function(x) createGR(cwR_corrected_list$purity, x, "copywriteR"))
CNV_GR_dfs <- sapply(IDs, function(x) createGR(CNV_corrected_list$purity, x, "CNVkit"))
merged_GR_dfs <- mapply(c, cwR_GR_dfs, CNV_GR_dfs)

## tools concordance - clinical threshold
broad_conc <- rbindlist(mapply(broadConcordance, merged_GR_dfs, IDs, broad_df$broad_quality, broad_df$broad_th, clnTh, SIMPLIFY = FALSE)) #tools concordance
broad_stats <- rbindlist(lapply(IDs, function(x) broadStats(broad_conc, x))) #broad stats per sample

write_tsv(broad_conc, paste0(outDir, "CNAs/refitted_calls/broad_DJ_concordance.txt")) #saving
write_tsv(broad_stats, paste0(outDir, "CNAs/refitted_calls/broad_stats.txt"))

## concordant broad CNA classification
broad_calls <- broad_conc %>%
  mutate(cwR_call = ifelse(cwR_CN >= 2 + clnTh, "amp", ifelse(cwR_CN <= 2 - clnTh, "del", NA)), 
         CNV_call = ifelse(CNV_CN >= 2 + clnTh, "amp", ifelse(CNV_CN <= 2 - clnTh, "del", NA)), 
         call = ifelse(cwR_call==CNV_call & CNV_call=="amp", "amp", ifelse(cwR_call==CNV_call & CNV_call=="del", "del", NA))) %>%
  filter(!is.na(call)) %>%
  select(ID, chrarm, cwR_CN, CNV_CN, cwR_call, CNV_call, call)

write_tsv(broad_calls, paste0(outDir, "CNAs/refitted_calls/broad_calls.txt")) #saving


### ANALYZING FOCAL CNAs ###
message("Running tool: HsMetrics.\n")

## loading target genes normalization factors  
#focal_NF <- fread(paste0(libDir, "inst/extdata/targets_Normalization_Factors.txt"))
focal_NF <- fread(paste0(libDir, "inst/extdata/RTM-NovaSeq_targets_Normalization_Factors.txt"))

## exons data
exon_files <- searchFiles(wd, regex$exons)
exons_tmp_df <- rbindlist(lapply(seq_along(exon_files), function(x) fread(exon_files[[x]], dec = ",")) %>%
  Map(cbind, ., ID = str_extract(exon_files, "(?<=HsMetrics_).*(?=_targetCoverage.txt)"))) %>%
  filter(!name %>% str_detect("IgH"))

## merging exons data with focal normalization factors
exons_tmp_df <- left_join(exons_tmp_df, focal_NF, by = "name") %>%
  mutate(gene = str_extract(name, ".*(?=_exon)")) %>%
  rename_with(~"chr", chrom) %>%
  relocate(ID, .before = chr) 
exons_raw_df <- rbindlist(sapply(IDs, function(x) computeFocalCN(exons_tmp_df, x), simplify = FALSE))

## adding CDKN2C and FAF1 
CDKN2C_FAF1 <- exons_raw_df %>% filter(str_detect(gene, "CDKN2C|FAF1")) 
CDKN2C_FAF1$gene <- unique(CDKN2C_FAF1$gene) %>% paste0(collapse = "_")
exons_raw_df <- bind_rows(exons_raw_df, CDKN2C_FAF1) %>% #binding
  select(-ends_with(c("probe", "gc")), -pct_0x, -read_count, -n, -starts_with("N")) %>%
  mutate(chr = gsub("chr", "", chr))

write_tsv(exons_raw_df, paste0(outDir, "CNAs/raw_calls/exons_raw_data.txt")) #saving

## exons correction
# chr_list <- c("10p", "10q", "12p", "12q", "16p", "17p", "17q", "18p", "18q", "1p", "20p", "20q", "22q", "2p", "2q", "4p", "4q", "6p", "8q")
focal_df <- quality_df %>% select(ID, purity, focal_quality, focal_th)
exons_corrected_list <- correctSegments(exons_raw_df, "focal", focal_df) 

write_tsv(exons_corrected_list$bobafit$report, paste0(outDir, "CNAs/refitted_calls/exons_bobafit_report.txt")) #saving
write_tsv(exons_corrected_list$bobafit$corrected_segments, paste0(outDir, "CNAs/refitted_calls/exons_bobafit_data.txt"))
write_tsv(exons_corrected_list$purity, paste0(outDir, "CNAs/refitted_calls/exons_corrected_data.txt"))


## genes data
genes_raw_df <- createGenes(exons_raw_df, "CN")
write_tsv(genes_raw_df, paste0(outDir, "CNAs/raw_calls/genes_raw_data.txt")) #saving

## genes correction
genes_corrected_df <- createGenes(exons_corrected_list$purity, "CN_purity_corrected")
genes_corrected_df <- Popeye2(genes_corrected_df, refGen, removeXY = F) %>%
  rename_with(~str_replace(., "_CN", "_CN_corrected"))
genes_corrected_df <- merge(genes_corrected_df, focal_df, by = "ID")

write_tsv(genes_corrected_df, paste0(outDir, "CNAs/refitted_calls/genes_corrected_data.txt")) #saving

## focal CNA classification 
focal_calls <- genes_corrected_df %>% 
  mutate(alt = ifelse(Q1_CN_corrected >= 2 + clnTh, "amp", ifelse(Q3_CN_corrected <= 2 - clnTh, "del", NA))) %>%
  select(-c(chr, start, end, width, purity), -starts_with(c("CI95", "focal"))) %>%
  filter(!is.na(alt))

write_tsv(focal_calls, paste0(outDir, "CNAs/refitted_calls/focal_calls.txt"))

message("CNAs calling algorithm: done!\n\n")


# ---- t-IgH calling algorithm ----
message("IgH translocation calling algorithm.\n")
dir.create(paste0(outDir, "t_IgH/"), recursive = TRUE)


### ANALYZING DELLY ###
message("Running tool: DELLY.\n")

## reading vcf files
delly_df <- createVcf(searchFiles(wd, regex$delly), filters$delly)

## filtering
delly_IgH_tmp <- filterForChr14(delly_df)
delly_IgH_df <- filterForTargets(delly_IgH_tmp, refGen, partners, IgHLocus, targetGene) %>%
  mutate(call_ID = paste0("SR=", SR,"_", "PE=", PE,"_", "Q=", QUAL), 
         VAF_split = (gt_DV / (gt_DR + gt_DV)) %>% round(3), 
         VAF_junct = (gt_RV / (gt_RR + gt_RV)) %>% round(3)) %>%
  arrange(CHROM)

write_tsv(delly_IgH_df, paste0(outDir, "t_IgH/delly_IgH.txt")) #saving

## DELLY's t_IgH calls 
delly_calls <- reportIgH(delly_IgH_df, IDs)
write_tsv(delly_calls, paste0(outDir, "t_IgH/delly_calls.txt")) #saving


### ANALYZING Manta ###
message("Running tool: Manta.\n")

## reading vcf files
manta_df <- createVcf(searchFiles(wd, regex$manta), filters$manta) %>%
  mutate(CHR2 = ALT %>% str_extract("[\\[\\]].*:") %>% str_remove_all(":|\\[|\\]"), 
         POS2 = as.numeric(ALT %>% str_extract(":[0-9]+[\\[\\]]") %>% str_remove_all(":|\\[|\\]")), 
         TOT_DEPTH = BND_DEPTH + MATE_BND_DEPTH)

## filtering
manta_IgH_tmp <- filterForChr14(manta_df) %>% filter(TOT_DEPTH > 5)
manta_IgH_df <- filterForTargets(manta_IgH_tmp, refGen, partners, IgHLocus, targetGene) %>%
  filter(CHROM==14) %>%
  mutate(call_ID = paste0("DP=", BND_DEPTH, ";",MATE_BND_DEPTH))

write_tsv(manta_IgH_df, paste0(outDir, "t_IgH/manta_IgH.txt")) #saving

## Manta's t_IgH calls
manta_calls <- reportIgH(manta_IgH_df, IDs)
write_tsv(manta_calls, paste0(outDir, "t_IgH/manta_calls.txt")) #saving


### CHECKING FOR TOOLS CONCORDANCE ###
IgH_calls <- full_join(delly_calls, manta_calls, by = "ID", suffix = c("_delly", "_manta")) %>%
  mutate(conc_4_14 = ifelse((t_IgH_4_14_delly + t_IgH_4_14_manta)>=1, 1, 0), 
         conc_6_14 = ifelse((t_IgH_6_14_delly + t_IgH_6_14_manta)>=1, 1, 0), 
         conc_11_14 = ifelse((t_IgH_11_14_delly + t_IgH_11_14_manta)>=1, 1, 0), 
         conc_14_16 = ifelse((t_IgH_14_16_delly + t_IgH_14_16_manta)>=1, 1, 0), 
         conc_14_20 = ifelse((t_IgH_14_20_delly + t_IgH_14_20_manta)>=1, 1, 0)) %>%
  select(ID, t_IgH_delly, call_delly, type_delly, t_IgH_manta, call_manta, type_manta,
         t_IgH_4_14_delly, t_IgH_4_14_manta, conc_4_14, 
         t_IgH_6_14_delly, t_IgH_6_14_manta, conc_6_14, 
         t_IgH_11_14_delly, t_IgH_11_14_manta, conc_11_14, 
         t_IgH_14_16_delly, t_IgH_14_16_manta, conc_14_16, 
         t_IgH_14_20_delly, t_IgH_14_20_manta, conc_14_20)

write_tsv(IgH_calls, paste0(outDir, "t_IgH/t_IgH_calls.txt")) #saving

message("IgH translocation calling algorithm: done!\n\n")


# ---- Mutations calling algorithm ----
message("Mutations calling algorithm.\n")
dir.create(paste0(outDir, "mutations/"), recursive = TRUE)


### ANALYZING Mutect2 ###
message("Running tool: Mutect2.\n")

mut_files <- searchFiles(wd, regex$mut)

## creating vcf df
mut_tmp_df <- createVariantID(createVcf(mut_files, filters$mut))
mut_df <- filterSNV(mut_tmp_df, "Mutect2")

write_tsv(mut_df, paste0(outDir, "mutations/mutect2_variants.txt"))


### ANALYZING Strelka ###
message("Running tool: Strelka.\n")

str_files <- searchFiles(wd, regex$str)

## creating vcf df
str_tmp_df <- createVariantID(createVcf(str_files))
str_df <- filterSNV(str_tmp_df, "Strelka")

write_tsv(str_df, paste0(outDir, "mutations/strelka_variants.txt"))


### ANALYZING Varscan2 ###
message("Running tool: Varscan2.\n")

var_files <- searchFiles(wd, regex$var2)

## creating vcf df
var_df <- createVariantID(createVcf(var_files))
write_tsv(var_df, paste0(outDir, "mutations/varscan2_variants.txt"))


### ANALYZING Freebayes ###
message("Running tool: FreeBayes.\n")

fb_files <- searchFiles(wd, regex$fb)
message("Checking if Freebayes files are the correct ones... \n")
print(fb_files)
message("\n")

## creating vcf df
fb_df <- createVariantID(createVcf(fb_files))
fb_df$gt_AO <- fb_df$gt_AO %>% as.numeric() #manipulating
fb_df$gt_RO <- fb_df$gt_RO %>% as.numeric()
fb_df <- fb_df %>%  mutate(VAF = gt_AO / (gt_AO + gt_RO))

write_tsv(fb_df, paste0(outDir, "mutations/freebayes_variants.txt"))


### CHECKING FOR TOOLS CONCORDANCE AND ANNOTATING MUTATIONS ###
join_cols <- colnames(mut_df)[-c(6, 7, 31)]

## first validation step
mut_str_tmp <- inner_join(mut_df, str_df, by = join_cols, suffix = c("_mutect2", "_strelka"), relationship = "many-to-many") %>%
  mutate(callers = paste0("Mutect2, Strelka"))

mut_str_tmp$VAF <- (mut_str_tmp$VAF_mutect2 + mut_str_tmp$VAF_strelka) / 2
mut_str_tmp$DP <- (mut_str_tmp$DP_mutect2 + mut_str_tmp$DP_strelka) / 2
mut_str_tmp$FILTER <- paste0(mut_str_tmp$FILTER_mutect2, ", ", mut_str_tmp$FILTER_strelka)

mut_str_df <- mut_str_tmp %>%
  select(-ends_with(c("_mutect2", "_strelka"))) %>%
  relocate(DP, FILTER, .before = Func.refGene) %>%
  relocate(VAF, .before = callers)

## second validation step
str_filtered_df <- str_df[!(str_df$variant_name %in% mut_str_tmp$variant_name),]
strelka_variants_df <- validateSNV(str_filtered_df, fb_df, var_df, "Strelka")

mut_filtered_df <- mut_df[!(mut_df$variant_name %in% mut_str_tmp$variant_name),]
mutect2_variants_df <- validateSNV(mut_filtered_df, fb_df, var_df, "Mutect2")

mut_filtered_df <- mut_df[!(mut_df$variant_name %in% mut_str_tmp$variant_name),]
mutect2_variants_df <- validateSNV(mut_filtered_df, fb_df, var_df, "Mutect2")

variants_df <- rbind(mut_str_df, strelka_variants_df, mutect2_variants_df)

## annotating
annotated_variants <- VarianThinker(variants_df, refGen)
write_tsv(annotated_variants, paste0(outDir, "mutations/annotated_variants.txt"))

message("Mutations calling algorithm: done!\n\n")

# ---- Barcellona criteria ---- 
message("Barcellona criteria.\n")


## del1p
del1p_dt <- bcnCNA(broad_calls, IDs, arm = "1p", alt = "del") %>%
  rename_with(~"del1p", alt)

## amp1q
amp1q_dt <- bcnCNA(broad_calls, IDs, arm = "1q", alt = "amp") %>%
  rename_with(~"amp1q", alt)

## del17p
del17p_dt <- bcnCNA(broad_calls, IDs, arm = "17p", alt = "del") %>%
  rename_with(~"del17p", alt)


## t_IgH
IgH_dt <- IgH_calls %>%
  mutate(alt_delly = detectIgH(IgH_calls, "t_IgH_delly"),
         alt_manta = detectIgH(IgH_calls, "t_IgH_manta"), 
         t_IgH = ifelse(!is.na(alt_delly), paste0(alt_delly), ifelse(!is.na(alt_manta), paste0(alt_manta), NA))) %>%
  select(ID, t_IgH)
IgH_dt[IgH_dt==""] <- NA


## TP53
tp53_dt <- annotated_variants %>% filter(Gene.refGene=="TP53" & InterVar_automated %in% c("Pathogenic", "Likely_pathogenic")) %>%
  select(ID = Indiv, TP53 = POS) %>%
  complete(ID = IDs) 

#binding
bcn_dt <- Reduce(function(x, y) merge(x, y, all = TRUE), list(del1p_dt, amp1q_dt, del17p_dt, IgH_dt, tp53_dt))

write_tsv(bcn_dt, paste0(outDir, "bcn_scriteria.txt"))

message("Barcellona criteria: done!\n")

