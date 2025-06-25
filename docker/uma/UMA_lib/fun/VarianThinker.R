#!/usr/bin/r

VarianThinker <- function(variants, refGen){
  if (refGen=="hg19"){
    annovar_feat <- c("SIFT_pred", "MutationTaster_pred", "fathmm.MKL_coding_pred", "Polyphen2_HDIV_pred", "LRT_pred", "MutationAssessor_pred", "Polyphen2_HVAR_pred",
                      "ExonicFunc.refGene", "Func.refGene", "CLNSIG",
                      "MetaLR_pred", "MetaSVM_pred",
                      "MutationTaster_pred", "VEST4_score", "PROVEAN_pred",
                      "gnomAD_genome_ALL")
    } else if (refGen=="hg38"){
      annovar_feat <- c("SIFT_pred", "MutationTaster_pred", "fathmm.MKL_coding_pred", "Polyphen2_HDIV_pred", "LRT_pred", "MutationAssessor_pred", "Polyphen2_HVAR_pred",
                      "ExonicFunc.refGene", "Func.refGene", "CLNSIG",
                      "MetaLR_pred", "MetaSVM_pred",
                      "MutationTaster_pred", "VEST4_score", "PROVEAN_pred",
                      "gnomad41_genome_AF")
    }

  if (all(annovar_feat %in% names(variants))){
    pseudo_sophia <- notes <- confidence <- preds <- tot_preds <- pred_ratio <- vector() 

    for (i in 1:nrow(variants)){
      variant <- variants[i, ]
      
      variant_pred <- variant %>% select(SIFT_pred, MutationTaster_pred, fathmm.MKL_coding_pred,
                                         Polyphen2_HDIV_pred, LRT_pred, MutationAssessor_pred, Polyphen2_HVAR_pred)
      pred_cols <- c("SIFT_pred", "MutationTaster_pred", "fathmm.MKL_coding_pred", "Polyphen2_HDIV_pred", "LRT_pred", "MutationAssessor_pred", "Polyphen2_HVAR_pred")
      pred_score_count <- c("D", "A", "H", "H/M")

      bio_pred_count <- sum(variant[, pred_cols] %in% pred_score_count)
      bio_pred_tot <- sum(variant[, pred_cols] != "")

      preds <- append(preds, bio_pred_count)
      tot_preds <- append(tot_preds, bio_pred_tot )
      pred_ratio <- append(pred_ratio, bio_pred_count/bio_pred_tot)
      
      # CODING CONSEQUENCE [1st level of confidence]
      if(variant[, "ExonicFunc.refGene"] %in% c("frameshift insertion", "stopgain", "frameshift deletion", "stoploss")){
        pseudo_sophia <- append(pseudo_sophia, "A")
        notes <- append(notes, "1:loss of function")
        confidence <- append(confidence, 1)
      } else if(variant[, "Func.refGene"]=="splicing"){
        pseudo_sophia <- append(pseudo_sophia, "A")
        notes <- append(notes, "1:splicing")
        confidence <- append(confidence, 1)
      #} else if(grepl("^synonymous",variant[, "ExonicFunc.refGene"])){
        #pseudo_sophia <- append(pseudo_sophia, "C")
        #notes <- append(notes, "1:synonymous")
        #confidence <- append(confidence, 1)
      #} else if(grepl("intronic|UTR", variant[, "Func.refGene"])){
        #pseudo_sophia <- append(pseudo_sophia, "C")
        #notes <- append(notes, "1:non coding")
        #confidence <- append(confidence, 1)
        
        # CLINVAR [2nd level of confidence]
      } else if(grepl("^Pathogenic|Likely pathogenic", variant[, "CLNSIG"])==T){
        pseudo_sophia <- append(pseudo_sophia, "A")
        notes <- append(notes, "2:CLINVAR pathogenic")
        confidence <- append(confidence, 2)
      } else if(grepl("^drug response$", variant[, "CLNSIG"])){
        pseudo_sophia <- append(pseudo_sophia, "D")
        notes <- append(notes, "2:CLINVAR drug response")
        confidence <- append(confidence, 2)
      } else if(grepl("^risk factor$", variant[, "CLNSIG"])){
        pseudo_sophia <- append(pseudo_sophia, "D")
        notes <- append(notes, "2:CLINVAR risk factor")
        confidence <- append(confidence, 2)
      } else if(grepl("Benign", variant[, "CLNSIG"])){
        pseudo_sophia <- append(pseudo_sophia, "D")
        notes <- append(notes, "2:CLINVAR Benign")
        confidence <- append(confidence, 2)

        # CLINVAR + META [3rd level of confidence]
      } else if(grepl("Conflicting|Uncertain", variant[, "CLNSIG"]) & (variant[, "MetaLR_pred"]=="D" & variant[, "MetaSVM_pred"]=="D")){
        if((bio_pred_count/bio_pred_tot) > 0.5 | ((variant[, "VEST4_score"] > 0.5) & (variant[, "PROVEAN_pred"]=="D"))){ # confirmation center
          pseudo_sophia <- append(pseudo_sophia, "A")
          notes <- append(notes, "3:CLINVAR uncertain - META pathogenic - Confirmed")
          confidence <- append(confidence, 3)
        } else {
          pseudo_sophia <- append(pseudo_sophia, "B")
          notes <- append(notes, "3:CLINVAR uncertain - META pathogenic- Not Confirmed")
          confidence <- append(confidence, 3)
        }
      } else if(grepl("Conflicting|Uncertain", variant[, "CLNSIG"]) & (variant[, "MetaLR_pred"]=="T" | variant[, "MetaSVM_pred"]=="T")){
        pseudo_sophia <- append(pseudo_sophia, "B")
        notes <- append(notes, "3:CLINVAR uncertain - META tolerated")
        confidence <- append(confidence, 3)
      } else if(grepl("Likely benign",variant[, "CLNSIG"]) & (variant[, "MetaLR_pred"]=="T" & variant[, "MetaSVM_pred"]=="T")){
        pseudo_sophia <- append(pseudo_sophia, "D")
        notes <- append(notes, "3:CLINVAR likely-benign - META benign")
        confidence <- append(confidence, 3)

        # SINGLE BioPred Algorhitm (MutationTaster) [4th level of confidence]
      } else if(variant[, "MutationTaster_pred"] %in% c("N","D")){
        if((bio_pred_count/bio_pred_tot) > 0.5) { # Majority vote
          pseudo_sophia <- append(pseudo_sophia, "B")
          notes <- append(notes, "4:MutTaster Poly-DiseaseCause - Confirmed in center (Majority Vote)")
          confidence <- append(confidence, 4)
        } else if(variant[, "MetaLR_pred"]=="D"){
          pseudo_sophia <- append(pseudo_sophia, "B")
          notes <- append(notes, "4:MutTaster Poly-DiseaseCause - Confirmed in center (Meta)")
          confidence <- append(confidence, 4)
        } else if((variant[, "VEST4_score"] > 0.5) & (variant[, "PROVEAN_pred"]=="D")){
          pseudo_sophia <- append(pseudo_sophia, "B")
          notes <- append(notes, "4:MutTaster Poly-DiseaseCause - Confirmed in center (VEST3 & PROVEAN)")
          confidence <- append(confidence, 4)
        } else {
          pseudo_sophia <- append(pseudo_sophia, "C")
          notes <- append(notes, "4:MutTaster Poly-DiseaseCause - Not Confirmed from center")
          confidence <- append(confidence, 4)
        }
      } else if(variant[, "MutationTaster_pred"]=="P"){
        pseudo_sophia <- append(pseudo_sophia, "D")
        notes <- append(notes, "4:MutTaster automatic polymorphism")
        confidence <- append(confidence, 4)
      } else if(variant[, "MutationTaster_pred"]=="A"){
        pseudo_sophia <- append(pseudo_sophia, "A")
        notes <- append(notes, "4:MutTaster automatic DiseaseCause")
        confidence <- append(confidence, 4)

        # GNOMAD [5th level of confidence]
      } else if(!is.na(variant[, annovar_feat[16]]) < 0.01){
        pseudo_sophia <- append(pseudo_sophia, "B")
        notes <- append(notes, "5: GnomAD not a SNP")
        confidence <- append(confidence, 5)
      } else if(!is.na(variant[, annovar_feat[16]]) >= 0.01){
        pseudo_sophia <- append(pseudo_sophia, "C")
        notes <- append(notes, "5: GnomAD SNP")
        confidence <- append(confidence, 5)
      } else {
        pseudo_sophia <- append(pseudo_sophia, "else")
        notes <- append(notes, "else")
        confidence <- append(confidence, 6)
      }
    }
    
    #vt_cols <- c("VarianThinker_Category", "VT_notes", "VT_confidence", "pred_pathogenic", "pred_tot_calls", "pred_ratio")
    #vt_res <- c(pseudo)
    variants <- variants %>% mutate(VarianThinker_Category = pseudo_sophia, 
                                    VT_notes = notes, 
                                    VT_confidence = confidence, 
                                    pred_pathogenic = preds, 
                                    pred_tot_calls = tot_preds, 
                                    pred_ratio = pred_ratio)
  } else {
    message("Error: wrong variables names.")
  }
  
  return(variants)
}
