setwd("~/EstagioFCTNOVA2019")
#Raw data
coad_patient_data<-read.table("gdc/4060482f-eedf-4959-97f1-f8b6c529c368/nationwidechildrens.org_clinical_patient_coad.txt", sep="\t", header=T)[-c(1,2),]
coad_TCGA_data<-read.table("gdc/03652df4-6090-4f5a-a2ff-ee28a37f9301/TCGA.COAD.mutect.03652df4-6090-4f5a-a2ff-ee28a37f9301.DR-10.0.somatic.maf", sep="\t", header=T)

#Data of interest
coad_interest_data<-coad_TCGA_data[,c("Hugo_Symbol","Tumor_Sample_Barcode","One_Consequence")]

#Calculating MAF for each mutation
coad_interest_data$MAF<-with(coad_TCGA_data,t_alt_count/t_depth)

#Calculating MATH for each patient
coad_MATH_data={aggregate(coad_interest_data$MAF,
                                 list(coad_interest_data$Tumor_Sample_Barcode),
                                 function(MAF){
                                   MAD<-median(abs(MAF-median(MAF)))*1.4826
                                   MATH<-100*MAD/median(MAF)
                                   MATH})}
names(coad_MATH_data)<-c("Tumor_Sample_Barcode","MATH")

#Pathways list of lists (Daniel Sobral)

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


#Genes Lists
pathways_list<-{list(
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

#Creating pathway variables
aggregated<-aggregate(coad_interest_data$Hugo_Symbol,
                      list(coad_interest_data$Tumor_Sample_Barcode),
                      function(muts){as.character(muts)})
variables<-{t(sapply(apply(aggregated[2],1,function(x){
  sapply(unlist(x),function(y){
    sapply(pathways_list,function(z){y%in%z})})}),function(x){
      apply(x,1,function(x){sum(x)>0})}))}
coad_final_data<-cbind(coad_MATH_data,variables)

#Linear Modeling
fit<-lm(MATH~EpigeneticMod+Transcription_factor+Genome_integrity+RTK_signalling+Cell_cycle+MAPK_signalling+PIK_signalling+TGFB_signalling+Wnt_BCatenin_signalling+Proteolysis+Splicing+HIPPO_signalling+Metabolism+NFE2L+Protein_phosphatase+Ribosome+TOR_signalling,
   data=coad_final_data)
coef<-as.data.frame(summary(fit)$coef)
sig<-coef[coef$`Pr(>|t|)`<.05,]
sig<-sig[order(sig$`Pr(>|t|)`),]

#Excluir Silent das pathways
silent_consequences<-c("3_prime_UTR_variant","5_prime_UTR_variant","downstream_gene_variant","intergenic_variant","intron_variant","non_coding_transcript_exon_variant","synonymous_variant","upstream_gene_variant")
not_silent<-!coad_interest_data$One_Consequence%in%silent_consequences
not_silent_aggregated<-aggregate(coad_interest_data$Hugo_Symbol[not_silent],
                      list(coad_interest_data$Tumor_Sample_Barcode[not_silent]),
                      function(muts){as.character(muts)})
not_silent_variables<-{t(sapply(apply(not_silent_aggregated[2],1,function(x){
  sapply(unlist(x),function(y){
    sapply(pathways_list,function(z){y%in%z})})}),function(x){
      apply(x,1,function(x){sum(x)>0})}))}
not_silent_final<-cbind(coad_MATH_data,not_silent_variables)

#Linear Modeling (not silent)
not_silent_fit<-lm(MATH~EpigeneticMod+Transcription_factor+Genome_integrity+RTK_signalling+Cell_cycle+MAPK_signalling+PIK_signalling+TGFB_signalling+Wnt_BCatenin_signalling+Proteolysis+Splicing+HIPPO_signalling+Metabolism+NFE2L+Protein_phosphatase+Ribosome+TOR_signalling,
        data=not_silent_final)
not_silent_coef<-as.data.frame(summary(not_silent_fit)$coef)
not_silent_sig<-not_silent_coef[not_silent_coef$`Pr(>|t|)`<.05,]
not_silent_sig<-not_silent_sig[order(not_silent_sig$`Pr(>|t|)`),]

#Excluir variantes com copy numbers!=0
coad_copy_numbers<-read.table("gdc/7f01e47a-2f4c-4b91-8db6-e1b6b5595390/COAD.focal_score_by_genes.txt", sep="\t")
