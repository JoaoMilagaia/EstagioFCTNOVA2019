setwd("~/EstagioFCTNOVA2019")

# Calculate MATH
MATH_calculator <- function(interest_data){#interest_data is the data.frame read from the raw text file
  # Calculating MAF for each mutation
  interest_data$MAF<-with(interest_data,t_alt_count/t_depth)
  
  # Calculating MATH for each patient
  MATH_data <-{aggregate(interest_data$MAF, # Grouping the MAF values
                  list(interest_data$Tumor_Sample_Barcode), # Per patient
                  function(MAF){ # And calculating each patient's MATH from it's MAFs
                    MAD<-median(abs(MAF-median(MAF)))*1.4826
                    MATH<-100*MAD/median(MAF)
                    return(MATH)})}
  
  names(MATH_data)<-c("Tumor_Sample_Barcode","OUTPUT") # The output here is MATH values
  
  return(MATH_data)
}

# Pathways list of lists (Daniel Sobral)
{
# pathways_list<-list(
#   TF=list("VHL","GATA3","TSHZ3","EP300","CTCF","TAF1","TSHZ2","RUNX1","MECOM","TBX3","SIN3A","WT1","EIF4A2","FOXA1","PHF6","CBFB",
#           "SOX9","ELF3","VEZF1","CEBPA","FOXA2"),
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
#   Other=list("NAV3","NOTCH1","LRRK2","MALAT1","ARHGAP35","POLQ","NCOR1","USP9X","NPM1","HGF","EPPK1","AR","LIFR","PRX","CRIPAK",
#              "EGR3","B4GALT3","MIR142")
#   )
}

# Gene List
pathways_list <- {list(
  Transcription_factor=list("FOXA2","CEBPA","VEZF1","SOX9","PHF6","EIF4A2","WT1","SIN3A","EP300","TBX3","MECOM","RUNX1","TSHZ2",
                            "TAF1","CTCF","TSHZ3","GATA3","VHL"),
  EpigeneticMod=list("EZH2","ASXL1","ARID5B","MLL4","KDM6A","KDM5C","SETBP1","NSD1","SETD2","PBRM1","ARID1A","MLL2","MLL3","TET1",
                     "TET2","DNMT3A","DNMT3B","DNMT1","HIST1H1C","HIST1H2BD","H3F3C"), #remove "TET1","DNMT3B","DNMT1"??
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
pathway_variables_func <- function(interest_data, gene_group){
  aggregated <- aggregate(interest_data$Hugo_Symbol, # Aggregate the mutated genes
                          list(interest_data$Tumor_Sample_Barcode),
                          function(muts){as.character(muts)})
  
  variables <- {t(sapply(apply(aggregated[2],1,function(x){
    sapply(unlist(x),function(y){
      sapply(gene_group,function(z){y%in%z})})}),function(x){
        apply(x,1,function(x){sum(x)>0})}))}
  
  MATH_data <- MATH_calculator(interest_data)
  
  final_data<-cbind(MATH_data,variables)[,-c(1)]
  return(final_data)
}

# Exclude silent mutations
delete_silent <- function(interest_data, gene_group){
  silent_consequences<-c("3_prime_UTR_variant","5_prime_UTR_variant","downstream_gene_variant","intergenic_variant",
                         "intron_variant","non_coding_transcript_exon_variant","synonymous_variant","upstream_gene_variant")
  not_silent<-interest_data[!interest_data$One_Consequence%in%silent_consequences,]
  not_silent_final <- pathway_variables_func(not_silent, gene_group)
  return(not_silent_final)
}

# SNV
paths <- sapply(read.table("SNV/MANIFEST.txt", header = T)[2],as.character)[1:33]
# Smaller sample to test
# paths <- paths[1:2]

cancer_names <- character(length(paths))
for (i in 1:length(paths)){
  cancer_names[i] <- strsplit(paths[i],"[.]")[[1]][2]
}

# SNV_data<- list()
# for (i in 1:length(paths)){
#   SNV_data[[i]] <- read.table(paste("SNV/",paths[i],sep=""),sep="\t",header=T,quote="")
# }
# names(SNV_data) <- cancer_names
# save(x=SNV_data, file = "SNV_data.RData")
load("SNV_data.RData")

# Process data
process <- function(data, gene_group){
  processed_data <- list()
  for (i in 1:length(data)){
    processed_data[[i]] <- delete_silent(data[[i]], gene_group)
  }
  names(processed_data) <- cancer_names
  return(processed_data)
}

# R_list
adj_r_sqrd <- function(processed_data){
  r_list <- numeric(length(processed_data)); names(r_list) <- cancer_names
  for (i in 1:length(processed_data)){
    fit <- lm(OUTPUT ~ ., data=processed_data[[i]])
    adj_r_sqrd <- (summary(fit)$adj.r.squared)
    r_list[i] <- adj_r_sqrd
  }
  return(r_list)
  }

# Coef_Matrix
coef <- function(processed_data){
  coef_matrix <- matrix(NA, nrow = length(paths), ncol = length(processed_data[[1]])-1, dimnames = list(cancer_names,
                                                                                                  names(processed_data[[1]])[-1]))
  for (i in 1:length(processed_data)){
    fit <- lm(OUTPUT ~ ., data=processed_data[[i]])
    coef <- as.data.frame(summary(fit)$coef)
    coef <- coef[-c(1),]
    sig <- coef[coef$`Pr(>|t|)`<.05,]
    if (nrow(sig) > 0){
      for (j in 1:nrow(sig)){
        rname <- row.names(sig)[j]
        rname <- substr(rname,1,nchar(rname)-4)
        coef_matrix[i,rname]<-sig[j,1]
      }
    }
  }
  return(coef_matrix)
}

# Plots
library(gplots)
library(pheatmap)

plotting <- function(processed_data){
  Adjusted_R_Squared<-adj_r_sqrd(processed_data)
  CoefMatrix<-coef(processed_data)
  pheatmap(CoefMatrix,
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           color = bluered(9),
           annotation_row = as.data.frame(Adjusted_R_Squared),
           angle_col = 315,
           breaks = c(seq(min(CoefMatrix,na.rm = T),0,length=5),seq(0.000001,max(CoefMatrix,na.rm = T),length=5)),
           na_col = 0
           )
}

# Results
# SNV_pathways <- process(SNV_data, pathways_list)
# save(x=SNV_pathways, file = "SNV_pathways.RData")
load("SNV_pathways.RData")

plotting(SNV_pathways)

# Different grouping of genes
load("geneSets.Rdata")
for (i in 1:length(geneSets)){
  geneSets[[i]]<-as.list(geneSets[[i]][1][,1])
}

# SNV_geneSets <- process(SNV_data, geneSets)
# save(x=SNV_geneSets, file = "SNV_geneSets.RData")
load("SNV_geneSets.RData")

plotting(SNV_geneSets)

# Forward Stepwise Selection
test <- SNV_pathways$COAD
regfit.fwd = regsubsets(OUTPUT~., data=test, nvmax=17, method = "forward")
plot(regfit.fwd,scale="Cp")

# Model selection using a validation set
dim(test)
set.seed(1)
train = sample(seq(399),270,replace=F) # 399 rows total, 270 is about 2/3 of the data
regfit.fwd=regsubsets(OUTPUT~.,data=test[train,],nvmax = 17,method = "forward")

# Predictions on observations not used for training
val.errors=rep(NA,17) # 17 variables so, 17 subset models
x.test=model.matrix(OUTPUT~.,data=test[-train,]) # Validation data set
for(i in 1:17){
  coefi=coefficients(regfit.fwd,id=i) # i is model size, coefi vector has the subset of variables used in that model
  pred=x.test[,names(coefi)]%*%coefi # Predicted response multiplying each used variable by it's coef (matrix multiplication)
  val.errors[i]=mean((test$OUTPUT[-train]-pred)^2) # mean squared error
}
plot(sqrt(val.errors),ylab="Root MSE",pch=19,type="b") # Validation error
points(sqrt(regfit.fwd$rss[-1]/270),col="blue",pch=19,type = "b")
legend("topright",legend=c("Training","Validation"),col=c("blue","black"),pch=19)

predict.regsubsets=fuction(object,newdata,id,...){
  form=as.formula(object$call[[2]])
  mat=model.matrix(form,newdata)
  coefi=coefficients(object,id=id)
  mat[,names(coefi)]%*%coefi
}
# More in gdc_analysis.Rmd
