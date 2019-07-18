setwd("~/Estágio2019/gdc")
#Raw data
coad_patient_data<-read.table("4060482f-eedf-4959-97f1-f8b6c529c368/nationwidechildrens.org_clinical_patient_coad.txt", sep="\t", header=T)
coad_TCGA_data<-read.table("03652df4-6090-4f5a-a2ff-ee28a37f9301/TCGA.COAD.mutect.03652df4-6090-4f5a-a2ff-ee28a37f9301.DR-10.0.somatic.maf", sep="\t", header=T)
#Treated data
coad_patient_data_filtered<-Filter(function(x)(length(unique(x))>1),coad_patient_data[-c(1,2),])
coad_TCGA_data_filtered<-Filter(function(x)(length(unique(x))>1),coad_TCGA_data)

#Calculating MAF for each mutation
coad_TCGA_data_filtered_interest<-coad_TCGA_data_filtered[,c(1,14,29,27,36)]
coad_TCGA_data_filtered_interest$Tumor_Sample_Barcode<-substr(coad_TCGA_data_filtered_interest$Tumor_Sample_Barcode,9,12)
names(coad_TCGA_data_filtered_interest)[2]<-"patient_id"
coad_TCGA_data_filtered_interest$MAF<-with(coad_TCGA_data,t_alt_count/t_depth)

#Calculating MAD for each patient
coad_final_data=aggregate(coad_TCGA_data_filtered_interest$MAF,
                                 list(coad_TCGA_data_filtered_interest$patient_id),
                                 function(MAF){MAD<-median(abs(MAF-median(MAF)))*1.4826
                                 100*MAD/median(MAF)})
names(coad_final_data)<-c("patient_id","MATH")

#Creating pathways data frame
{pathways_df<-data.frame(
  TF=c("VHL","GATA3","TSHZ3","EP300","CTCF","TAF1","TSHZ2","RUNX1","MECOM","TBX3","SIN3A","WT1","EIF4A2","FOXA1","PHF6","CBFB","SOX9","ELF3","VEZF1","CEBPA","FOXA2"),
  HM=c("MLL3","MLL2","ARID1A","PBRM1","SETD2","NSD1","SETBP1","KDM5C","KDM6A","MLL4","ARID5B","ASXL1","EZH2",rep(NA,8)),
  GI=c("TP53","ATM","ATRX","BRCA2","ATR","STAG2","BAP1","BRCA1","SMC1A","SMC3","CHEK2","RAD21","ERCC2",rep(NA,8)),
  RTK=c("EGFR","FLT3","EPHA3","ERBB4","PDGFRA","EPHB6","FGFR2","KIT","FGFR3",rep(NA,12)),
  CC=c("CDKN2A","RB1","CDK12","CDKN1B","CCND1","CDKN1A","CDKN2C",rep(NA,14)),
  MAPK=c("KRAS","NF1","MAP3K1","BRAF","NRAS","MAP2K4","MAPK8IP1",rep(NA,14)),
  PI3K=c("PIK3CA","PTEN","PIK3R1","TLR4","PIK3CG","AKT1",rep(NA,15)),
  TGF_B=c("SMAD4","TGFBR2","ACVR1B","SMAD2","ACVR2A",rep(NA,16)),
  Wnt_B=c("APC","CTNNB1","AXIN2","TBL1XR1","SOX17",rep(NA,16)),
  H=c("HIST1H1C","H3F3C","HIST1H2BD",rep(NA,18)),
  P=c("FBXW7","KEAP1","SPOP",rep(NA,18)),
  S=c("SF3B1","U2AF1","PCBP1",rep(NA,18)),
  HIPPO=c("CDH1","AJUBA",rep(NA,19)),
  DNA_M=c("DNMT3A","TET2",rep(NA,19)),
  M=c("IDH1","IDH2",rep(NA,19)),
  NFE2L=c("NFE2L2","NFE2L3",rep(NA,19)),
  PP=c("PPP2R1A","PTPN11",rep(NA,19)),
  R=c("RPL22","RPL5",rep(NA,19)),
  TOR=c("MTOR","STK11",rep(NA,19)),
  Other=c("NAV3","NOTCH1","LRRK2","MALAT1","ARHGAP35","POLQ","NCOR1","USP9X","NPM1","HGF","EPPK1","AR","LIFR","PRX","CRIPAK","EGR3","B4GALT3","MIR142",rep(NA,3)))}
#Creating pathways list of lists
{pathways_list<-list(
  TF=list("VHL","GATA3","TSHZ3","EP300","CTCF","TAF1","TSHZ2","RUNX1","MECOM","TBX3","SIN3A","WT1","EIF4A2","FOXA1","PHF6","CBFB","SOX9","ELF3","VEZF1","CEBPA","FOXA2"),
  HM=list("MLL3","MLL2","ARID1A","PBRM1","SETD2","NSD1","SETBP1","KDM5C","KDM6A","MLL4","ARID5B","ASXL1","EZH2",rep(NA,8)),
  GI=list("TP53","ATM","ATRX","BRCA2","ATR","STAG2","BAP1","BRCA1","SMC1A","SMC3","CHEK2","RAD21","ERCC2",rep(NA,8)),
  RTK=list("EGFR","FLT3","EPHA3","ERBB4","PDGFRA","EPHB6","FGFR2","KIT","FGFR3",rep(NA,12)),
  CC=list("CDKN2A","RB1","CDK12","CDKN1B","CCND1","CDKN1A","CDKN2C",rep(NA,14)),
  MAPK=list("KRAS","NF1","MAP3K1","BRAF","NRAS","MAP2K4","MAPK8IP1",rep(NA,14)),
  PI3K=list("PIK3CA","PTEN","PIK3R1","TLR4","PIK3CG","AKT1",rep(NA,15)),
  TGF_B=list("SMAD4","TGFBR2","ACVR1B","SMAD2","ACVR2A",rep(NA,16)),
  Wnt_B=list("APC","CTNNB1","AXIN2","TBL1XR1","SOX17",rep(NA,16)),
  H=list("HIST1H1C","H3F3C","HIST1H2BD",rep(NA,18)),
  P=list("FBXW7","KEAP1","SPOP",rep(NA,18)),
  S=list("SF3B1","U2AF1","PCBP1",rep(NA,18)),
  HIPPO=list("CDH1","AJUBA",rep(NA,19)),
  DNA_M=list("DNMT3A","TET2",rep(NA,19)),
  M=list("IDH1","IDH2",rep(NA,19)),
  NFE2L=list("NFE2L2","NFE2L3",rep(NA,19)),
  PP=list("PPP2R1A","PTPN11",rep(NA,19)),
  R=list("RPL22","RPL5",rep(NA,19)),
  TOR=list("MTOR","STK11",rep(NA,19)),
  Other=list("NAV3","NOTCH1","LRRK2","MALAT1","ARHGAP35","POLQ","NCOR1","USP9X","NPM1","HGF","EPPK1","AR","LIFR","PRX","CRIPAK","EGR3","B4GALT3","MIR142",rep(NA,3)))}
