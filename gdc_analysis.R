setwd("~/EstagioFCTNOVA2019")

# Calculate MATH
MATH_calculator <- function(interest_data){
  # Calculating MAF for each mutation
  interest_data$MAF<-with(interest_data,t_alt_count/t_depth)
  
  # Calculating MATH for each patient
  MATH_data <-{aggregate(interest_data$MAF,
                  list(interest_data$Tumor_Sample_Barcode),
                  function(MAF){
                    MAD<-median(abs(MAF-median(MAF)))*1.4826
                    MATH<-100*MAD/median(MAF)
                    MATH})}
  
  names(MATH_data)<-c("Tumor_Sample_Barcode","OUTPUT")
  
  return(MATH_data)
}

# Pathways list of lists (Daniel Sobral)
{
# pathways_list<-{list(
#   TF=list("VHL","GATA3","TSHZ3","EP300","CTCF","TAF1","TSHZ2","RUNX1","MECOM","TBX3","SIN3A","WT1","EIF4A2","FOXA1","PHF6","CBFB","SOX9","ELF3","VEZF1","CEBPA","FOXA2"),
#   HM=list("MLL3","MLL2","ARID1A","PBRM1","SETD2","NSD1","SETBP1","KDM5C","KDM6A","MLL4","ARID5B","ASXL1","EZH2"),
#   GI=list("TP53","ATM","ATRX","BRCA2","ATR","STAG2","BAP1","BRCA1","SMC1A","SMC3","CHEK2","RAD21","ERCC2"),
#   RTK=list("EGFR","FLT3","EPHA3","ERBB4","PDGFRA","EPHB6","FGFR2","KIT","FGFR3"),
#   CC=list("CDKN2A","RB1","CDK12","CDKN1B","CCND1","CDKN1A","CDKN2C"),
#   MAPK=list("KRAS","NF1","MAP3K1","BRAF","NRAS","MAP2K4","MAPK8IP1"),
#   PI3K=list("PIK3CA","PTEN","PIK3R1","TLR4","PIK3CG","AKT1"),
#   TGF_B=list("SMAD4","TGFBR2","ACVR1B","SMAD2","ACVR2A"),
#   Wnt_B=list("APC","CTNNB1","AXIN2","TBL1XR1","SOX17"),
#   H=list("HIST1H1C","H3F3C","HIST1H2BD"),
#   P=list("FBXW7","KEAP1","SPOP"),
#   S=list("SF3B1","U2AF1","PCBP1"),
#   HIPPO=list("CDH1","AJUBA"),
#   DNA_M=list("DNMT3A","TET2"),
#   M=list("IDH1","IDH2"),
#   NFE2L=list("NFE2L2","NFE2L3"),
#   PP=list("PPP2R1A","PTPN11"),
#   R=list("RPL22","RPL5"),
#   TOR=list("MTOR","STK11"),
#   Other=list("NAV3","NOTCH1","LRRK2","MALAT1","ARHGAP35","POLQ","NCOR1","USP9X","NPM1","HGF","EPPK1","AR","LIFR","PRX","CRIPAK","EGR3","B4GALT3","MIR142")
# )}
}

# Gene List
pathways_list <- {list(
  Transcription_factor=list("FOXA2","CEBPA","VEZF1","SOX9","PHF6","EIF4A2","WT1","SIN3A","EP300","TBX3","MECOM","RUNX1","TSHZ2","TAF1","CTCF","TSHZ3","GATA3","VHL"),
  EpigeneticMod=list("EZH2","ASXL1","ARID5B","MLL4","KDM6A","KDM5C","SETBP1","NSD1","SETD2","PBRM1","ARID1A","MLL2","MLL3","TET1","TET2","DNMT3A","DNMT3B","DNMT1","HIST1H1C","HIST1H2BD","H3F3C"), #remove "TET1","DNMT3B","DNMT1"??
  Genome_integrity=list("ERCC2","CHEK2","SMC3","SMC1A","BRCA1","BAP1","STAG2","ATR","BRCA2","ATRX","ATM","TP53"),
  RTK_signalling=list("FGFR3","KIT","FGFR2","EPHB6","PDGFRA","ERBB4","EPHA3","FLT3","EGFR"),
  Cell_cycle=list("CDKN2C","CDKN1A","CDK12","RB1","CDKN2A"),
  MAPK_signalling=list("MAPK8IP1","BRAF","MAP3K1","NF1","KRAS"),
  PIK_signalling=list("AKT1","PIK3CG","TLR4","PIK3R1","PTEN","PIK3CA"),
  TGFB_signalling=list("ACVR2A","SMAD2","ACVR1B","TGFBR2","SMAD4"),
  Wnt_BCatenin_signalling=list("TBL1XR1","AXIN2","CTNNB1","APC"),
  Proteolysis=list("SPOP","KEAP1","FBXW7"),
  Splicing=list("PCBP1","SF3B1"),
  HIPPO_signalling=list("CDH1"),
  Metabolism=list("IDH2","IDH1"),
  NFE2L=list("NFE2L3","NFE2L2"),
  Protein_phosphatase=list("PTPN11","PPP2R1A"),
  Ribosome=list("RPL5","RPL22"),
  TOR_signalling=list("STK11","MTOR")
)}

# Creating pathway variables
pathway_variables_func <- function(interest_data){
  aggregated <- aggregate(interest_data$Hugo_Symbol,
                          list(interest_data$Tumor_Sample_Barcode),
                          function(muts){as.character(muts)})
  
  variables <- {t(sapply(apply(aggregated[2],1,function(x){
    sapply(unlist(x),function(y){
      sapply(pathways_list,function(z){y%in%z})})}),function(x){
        apply(x,1,function(x){sum(x)>0})}))}
  
  MATH_data <- MATH_calculator(interest_data)
  
  final_data<-cbind(MATH_data,variables)[,-c(1)]
  return(final_data)
}

# Exclude silent mutations
delete_silent <- function(interest_data){
  silent_consequences<-c("3_prime_UTR_variant","5_prime_UTR_variant","downstream_gene_variant","intergenic_variant","intron_variant","non_coding_transcript_exon_variant","synonymous_variant","upstream_gene_variant")
  not_silent<-interest_data[!interest_data$One_Consequence%in%silent_consequences,]
  not_silent_final <- pathway_variables_func(not_silent)
  return(not_silent_final)
}

# Start to process data
paths <- sapply(read.table("SNV/MANIFEST.txt", header = T)[2],as.character)[1:33]
r_list <- c()
for (path in paths){
  final_data <- delete_silent(read.table(paste("SNV/",path,sep=""), sep ="\t", header = T))
  fit <- lm(OUTPUT ~ ., data=final_data)
  adj_r_sqrd <- (summary(fit)$adj.r.squared)
  cancer_name <- strsplit(path,"[.]")[[1]][2]
  r_list[cancer_name] <- adj_r_sqrd
  coef <- as.data.frame(summary(fit)$coef)
  coef <- coef[-c(1),]
  sig <- coef[coef$`Pr(>|t|)`<.05,]
  sig <- sig[order(sig$`Pr(>|t|)`),]
  }

# Plots
barplot(r_list, las=2, ylab = "Adjusted R squared")
barplot(r_list[order(r_list,decreasing = TRUE)], las=2, ylab = "Adjusted R squared")