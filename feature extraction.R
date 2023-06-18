#feature extraction
startTime <- Sys.time()
#######################################################################################################################
#library---------------------------------------------------------------------------------------------------------------
library("readxl")
library("xlsx")
library(foreign)
library(tidyr)
library(ftrCOOL)
#######################################################################################################################
#reading files---------------------------------------------------------------------------------------------------------
#PseudomonassyringaeTable contains Pseudomonas syringae sequences
#ArabidopsisThalianaTable contains Arabidopsis Thaliana sequences
#AT* contains all sequences of proteins that have more than one sequence
PseudomonassyringaeTable <- read_excel("G:/Thesis/RawData/Data-Sequence/Pseudomonassyringa_sequence.xlsx")
ArabidopsisThalianaTable <- read_excel("G:/Thesis/RawData/Data-Sequence/ArabidopsisThaliana_uniprot_sequence.xlsx")
AT4G17680 <- read_excel("G:/Thesis/RawData/Data-Sequence/AT4G17680.xlsx")
AT1G24050 <- read_excel("G:/Thesis/RawData/Data-Sequence/AT1G24050.xlsx")
AT3G11590 <- read_excel("G:/Thesis/RawData/Data-Sequence/AT3G11590.xlsx")
AT3G11720 <- read_excel("G:/Thesis/RawData/Data-Sequence/AT3G11720.xlsx")
AT3G48550 <- read_excel("G:/Thesis/RawData/Data-Sequence/AT3G48550.xlsx")
AT4G02550 <- read_excel("G:/Thesis/RawData/Data-Sequence/AT4G02550.xlsx")
IndependentPseu <- read_excel("G:/Thesis/RawData/Data-Sequence/IndependentTest-Sequence.xlsx", sheet = "Pseudomonas syringae")
IndependentAra <- read_excel("G:/Thesis/RawData/Data-Sequence/IndependentTest-Sequence.xlsx", sheet = "Arabidopsis thaliana")

#these four files are trains and tests file that is constructed by @constructing test & train script
randomTest <- read.csv("G:/Thesis/CodeForThesis-R/Data/randomTest.csv")
randomTrain <- read.csv("G:/Thesis/CodeForThesis-R/Data/randomTrain.csv")
SDDTest <- read.csv("G:/Thesis/CodeForThesis-R/Data/SDDTest.csv")
SDDTrain <- read.csv("G:/Thesis/CodeForThesis-R/Data/SDDTrain.csv")

Independent <- read.csv("G:/Thesis/CodeForThesis-R/Data/IndependentTest.csv")
Independent <- cbind(Independent, rep("interaction",nrow(Independent)))
names(Independent) <- names(randomTrain)

#######################################################################################################################
#feature extraction----------------------------------------------------------------------------------------------------
#ArabidopsisThaliana and PseudomonasSyringae the IDs that participate in interactions
ArabidopsisThalianaID <- c(ArabidopsisThalianaTable$`Gene names` ,"AT4G17680","AT1G24050","AT3G11590","AT3G11720","AT3G48550","AT4G02550")
PseudomonasSyringaeID <- PseudomonassyringaeTable$`Gene names`

#meanFeature calculate the mean of features for proteins that have more than one sequence
meanFeature <- function(MyTable)
{
  rowMyTable <- nrow(MyTable)
  temp <- c(1,2,colSums(MyTable[,-c(1,2)])/rowMyTable)
  return(temp)
}
#----------------------------------------------------------------------------------------------------------------------
#all functions which start with "Get" calculate feature for a table 
#features are extracted by ftrCOOL package
#number 1--------------------------------------------------------------------------------------------------------------
GetKAAComposition <- function(MyTable)
{
  rowMyTable <- nrow(MyTable)
  
  Temp <- kAAComposition(as.character(MyTable[1,2]), rng = 2, normalized = FALSE)
  for(i in 2:rowMyTable)
  {
    temp <- kAAComposition(as.character(MyTable[i,2]), rng = 2, normalized = FALSE)
    Temp <- rbind(Temp,temp)
  }
  MyTable <- cbind(MyTable,Temp)
  return(MyTable)
}
featurePseu1 <- GetKAAComposition(PseudomonassyringaeTable)
featureArab1 <- GetKAAComposition(ArabidopsisThalianaTable)
AT4G17680_1 <- meanFeature(GetKAAComposition(AT4G17680))
AT1G24050_1 <- meanFeature(GetKAAComposition(AT1G24050))
AT3G11590_1 <- meanFeature(GetKAAComposition(AT3G11590))
AT3G11720_1 <- meanFeature(GetKAAComposition(AT3G11720))
AT3G48550_1 <- meanFeature(GetKAAComposition(AT3G48550))
AT4G02550_1 <- meanFeature(GetKAAComposition(AT4G02550))
featureArab1 <- rbind(featureArab1 ,AT4G17680_1 ,AT1G24050_1 ,AT3G11590_1 ,AT3G11720_1 ,AT3G48550_1 ,AT4G02550_1)

IndependentAra1 <- GetKAAComposition(IndependentAra)
IndependentPseu1 <- GetKAAComposition(IndependentPseu)

colnames(IndependentAra1) <- paste("arabidopsis_KAAC", colnames(IndependentAra1), sep = "_")
colnames(IndependentPseu1) <- paste("Sringae_KAAC", colnames(IndependentPseu1), sep = "_")

colnames(featureArab1) <- paste("arabidopsis_KAAC", colnames(featureArab1), sep = "_")
colnames(featurePseu1) <- paste("Sringae_KAAC", colnames(featurePseu1), sep = "_")
#----------------------------------------------------------------------------------------------------------------------
#number 2--------------------------------------------------------------------------------------------------------------
GetPSEAAC <- function(MyTable)
{
  rowMyTable <- nrow(MyTable)
  
  Temp <- PSEAAC(as.character(MyTable[1,2]))
  for(i in 2:rowMyTable)
  {
    temp <- PSEAAC(as.character(MyTable[i,2]))
    Temp <- rbind(Temp,temp)
  }
  MyTable <- cbind(MyTable,Temp)
  return(MyTable)
}
featurePseu2 <- GetPSEAAC(PseudomonassyringaeTable)
featureArab2 <- GetPSEAAC(ArabidopsisThalianaTable)
AT4G17680_2 <- meanFeature(GetPSEAAC(AT4G17680))
AT1G24050_2 <- meanFeature(GetPSEAAC(AT1G24050))
AT3G11590_2 <- meanFeature(GetPSEAAC(AT3G11590))
AT3G11720_2 <- meanFeature(GetPSEAAC(AT3G11720))
AT3G48550_2 <- meanFeature(GetPSEAAC(AT3G48550))
AT4G02550_2 <- meanFeature(GetPSEAAC(AT4G02550))
featureArab2 <- rbind(featureArab2 ,AT4G17680_2 ,AT1G24050_2 ,AT3G11590_2 ,AT3G11720_2 ,AT3G48550_2 ,AT4G02550_2)

IndependentAra2 <- GetPSEAAC(IndependentAra)
IndependentPseu2 <- GetPSEAAC(IndependentPseu)

colnames(IndependentAra2) <- paste("arabidopsis_PAAC", colnames(IndependentAra2), sep = "_")
colnames(IndependentPseu2) <- paste("Sringae_PAAC", colnames(IndependentPseu2), sep = "_")

colnames(featureArab2) <- paste("arabidopsis_PAAC", colnames(featureArab2), sep = "_")
colnames(featurePseu2) <- paste("Sringae_PAAC", colnames(featurePseu2), sep = "_")
#----------------------------------------------------------------------------------------------------------------------
#number 3--------------------------------------------------------------------------------------------------------------
GetAAutoCor <- function(MyTable)
{
  rowMyTable <- nrow(MyTable)
  
  Temp <- AAutoCor(as.character(MyTable[1,2]), type = c("AC"), maxlag = 30,
                   selectedAAidx = list(c("CIDH920105", "BHAR880101", "CHAM820101", "CHAM820102",
                                          "CHOC760101", "BIGC670101", "KLEP840101", "KUHL950101","GRAR740102")))
  for(i in 2:rowMyTable)
  {
    temp <- AAutoCor(as.character(MyTable[i,2]), type = c("AC"), maxlag = 30,
                     selectedAAidx = list(c("CIDH920105", "BHAR880101", "CHAM820101", "CHAM820102",
                                            "CHOC760101", "BIGC670101", "KLEP840101", "KUHL950101","GRAR740102")))
    Temp <- rbind(Temp,temp)
  }
  MyTable <- cbind(MyTable,Temp)
  return(MyTable)
}
featurePseu3 <- GetAAutoCor(PseudomonassyringaeTable)
featureArab3 <- GetAAutoCor(ArabidopsisThalianaTable)
AT4G17680_3 <- meanFeature(GetAAutoCor(AT4G17680))
AT1G24050_3 <- meanFeature(GetAAutoCor(AT1G24050))
AT3G11590_3 <- meanFeature(GetAAutoCor(AT3G11590))
AT3G11720_3 <- meanFeature(GetAAutoCor(AT3G11720))
AT3G48550_3 <- meanFeature(GetAAutoCor(AT3G48550))
AT4G02550_3 <- meanFeature(GetAAutoCor(AT4G02550))
featureArab3 <- rbind(featureArab3 ,AT4G17680_3 ,AT1G24050_3 ,AT3G11590_3 ,AT3G11720_3 ,AT3G48550_3 ,AT4G02550_3)

IndependentAra3 <- GetAAutoCor(IndependentAra)
IndependentPseu3 <- GetAAutoCor(IndependentPseu)

colnames(IndependentAra3) <- paste("arabidopsis_AC", colnames(IndependentAra3), sep = "_")
colnames(IndependentPseu3) <- paste("Sringae_AC", colnames(IndependentPseu3), sep = "_")

colnames(featureArab3) <- paste("arabidopsis_AC", colnames(featureArab3), sep = "_")
colnames(featurePseu3) <- paste("Sringae_AC", colnames(featurePseu3), sep = "_")
#----------------------------------------------------------------------------------------------------------------------
#number 4--------------------------------------------------------------------------------------------------------------
GetSAAC <- function(MyTable)
{
  rowMyTable <- nrow(MyTable)
  
  Temp <- SAAC(as.character(MyTable[1,2]), normalized = FALSE)
  for(i in 2:rowMyTable)
  {
    temp <- SAAC(as.character(MyTable[i,2]), normalized = FALSE)
    Temp <- rbind(Temp,temp)
  }
  MyTable <- cbind(MyTable,Temp)
  return(MyTable)
}
featurePseu4 <- GetSAAC(PseudomonassyringaeTable)
featureArab4 <- GetSAAC(ArabidopsisThalianaTable)
AT4G17680_4 <- meanFeature(GetSAAC(AT4G17680))
AT1G24050_4 <- meanFeature(GetSAAC(AT1G24050))
AT3G11590_4 <- meanFeature(GetSAAC(AT3G11590))
AT3G11720_4 <- meanFeature(GetSAAC(AT3G11720))
AT3G48550_4 <- meanFeature(GetSAAC(AT3G48550))
AT4G02550_4 <- meanFeature(GetSAAC(AT4G02550))
featureArab4 <- rbind(featureArab4 ,AT4G17680_4 ,AT1G24050_4 ,AT3G11590_4 ,AT3G11720_4 ,AT3G48550_4 ,AT4G02550_4)

IndependentAra4 <- GetSAAC(IndependentAra)
IndependentPseu4 <- GetSAAC(IndependentPseu)

colnames(IndependentAra4) <- paste("arabidopsis_SAAC", colnames(IndependentAra4), sep = "_")
colnames(IndependentPseu4) <- paste("Sringae_SAAC", colnames(IndependentPseu4), sep = "_")

colnames(featureArab4) <- paste("arabidopsis_SAAC", colnames(featureArab4), sep = "_")
colnames(featurePseu4) <- paste("Sringae_SAAC", colnames(featurePseu4), sep = "_")
#----------------------------------------------------------------------------------------------------------------------
#number 5--------------------------------------------------------------------------------------------------------------
GetconjointTriad <- function(MyTable)
{
  rowMyTable <- nrow(MyTable)
  
  Temp <- conjointTriad(as.character(MyTable[1,2]), normalized = FALSE)
  for(i in 2:rowMyTable)
  {
    temp <- conjointTriad(as.character(MyTable[i,2]), normalized = FALSE)
    Temp <- rbind(Temp,temp)
  }
  MyTable <- cbind(MyTable,Temp)
  return(MyTable)
}
featurePseu5 <- GetconjointTriad(PseudomonassyringaeTable)
featureArab5 <- GetconjointTriad(ArabidopsisThalianaTable)
AT4G17680_5 <- meanFeature(GetconjointTriad(AT4G17680))
AT1G24050_5 <- meanFeature(GetconjointTriad(AT1G24050))
AT3G11590_5 <- meanFeature(GetconjointTriad(AT3G11590))
AT3G11720_5 <- meanFeature(GetconjointTriad(AT3G11720))
AT3G48550_5 <- meanFeature(GetconjointTriad(AT3G48550))
AT4G02550_5 <- meanFeature(GetconjointTriad(AT4G02550))
featureArab5 <- rbind(featureArab5 ,AT4G17680_5 ,AT1G24050_5 ,AT3G11590_5 ,AT3G11720_5 ,AT3G48550_5 ,AT4G02550_5)

IndependentAra5 <- GetconjointTriad(IndependentAra)
IndependentPseu5 <- GetconjointTriad(IndependentPseu)

colnames(IndependentAra5) <- paste("arabidopsis_CTRI", colnames(IndependentAra5), sep = "_")
colnames(IndependentPseu5) <- paste("Sringae_CTRI", colnames(IndependentPseu5), sep = "_")

colnames(featureArab5 ) <- paste("arabidopsis_CTRI", colnames(featureArab5), sep = "_")
colnames(featurePseu5 ) <- paste("Sringae_CTRI", colnames(featurePseu5), sep = "_")
#----------------------------------------------------------------------------------------------------------------------
#number 6--------------------------------------------------------------------------------------------------------------
GetCTD <- function(MyTable)
{
  rowMyTable <- nrow(MyTable)
  
  Temp <- CTD(as.character(MyTable[1,2]))
  for(i in 2:rowMyTable)
  {
    temp <- CTD(as.character(MyTable[i,2]))
    Temp <- rbind(Temp,temp)
  }
  MyTable <- cbind(MyTable,Temp)
  return(MyTable)
}
featurePseu6 <- GetCTD(PseudomonassyringaeTable)
featureArab6 <- GetCTD(ArabidopsisThalianaTable)
AT4G17680_6 <- meanFeature(GetCTD(AT4G17680))
AT1G24050_6 <- meanFeature(GetCTD(AT1G24050))
AT3G11590_6 <- meanFeature(GetCTD(AT3G11590))
AT3G11720_6 <- meanFeature(GetCTD(AT3G11720))
AT3G48550_6 <- meanFeature(GetCTD(AT3G48550))
AT4G02550_6 <- meanFeature(GetCTD(AT4G02550))
featureArab6 <- rbind(featureArab6 ,AT4G17680_6 ,AT1G24050_6 ,AT3G11590_6 ,AT3G11720_6 ,AT3G48550_6 ,AT4G02550_6)

IndependentAra6 <- GetCTD(IndependentAra)
IndependentPseu6 <- GetCTD(IndependentPseu)

colnames(IndependentAra6) <- paste("arabidopsis_CTD", colnames(IndependentAra6), sep = "_")
colnames(IndependentPseu6) <- paste("Sringae_CTD", colnames(IndependentPseu6), sep = "_")


colnames(featureArab6) <- paste("arabidopsis_CTD", colnames(featureArab6), sep = "_")
colnames(featurePseu6) <- paste("Sringae_CTD", colnames(featurePseu6), sep = "_")
#----------------------------------------------------------------------------------------------------------------------
#writing feature-------------------------------------------------------------------------------------------------------
ArabidopsisThalianaFeature <- cbind(ArabidopsisThalianaID,
                       featureArab1[,-c(1,2)],featureArab2[,-c(1,2)],featureArab3[,-c(1,2)],
                       featureArab4[,-c(1,2)],featureArab5[,-c(1,2)],featureArab6[,-c(1,2)])
write.csv(file="G:/Thesis/CodeForThesis-R/Data-FeatureExtraction/featureArabidopsisThaliana.csv",x=ArabidopsisThalianaFeature, row.names = FALSE)

PseudomonasSyringaeFeature <- cbind(PseudomonasSyringaeID,
                      featurePseu1[,-c(1,2)],featurePseu2[,-c(1,2)],featurePseu3[,-c(1,2)],
                      featurePseu4[,-c(1,2)],featurePseu5[,-c(1,2)],featurePseu6[,-c(1,2)])
write.csv(file="G:/Thesis/CodeForThesis-R/Data-FeatureExtraction/featurePseudomonasSyringae.csv",x=PseudomonasSyringaeFeature, row.names = FALSE)

IndependentAraFeature <- cbind(IndependentAra$`Gene names`,
                                    IndependentAra1[,-c(1,2)],IndependentAra2[,-c(1,2)],IndependentAra3[,-c(1,2)],
                                    IndependentAra4[,-c(1,2)],IndependentAra5[,-c(1,2)],IndependentAra6[,-c(1,2)])
write.csv(file="G:/Thesis/CodeForThesis-R/Data-FeatureExtraction/featureIndependentAra.csv",x=IndependentAraFeature, row.names = FALSE)

IndependentPseuFeature <- cbind(IndependentPseu$`Gene names`,
                                    IndependentPseu1[,-c(1,2)],IndependentPseu2[,-c(1,2)],IndependentPseu3[,-c(1,2)],
                                    IndependentPseu4[,-c(1,2)],IndependentPseu5[,-c(1,2)],IndependentPseu6[,-c(1,2)])
write.csv(file="G:/Thesis/CodeForThesis-R/Data-FeatureExtraction/featureIndependentPseu.csv",x=IndependentPseuFeature, row.names = FALSE)
#######################################################################################################################
#constructing data set-------------------------------------------------------------------------------------------------
#findFeature finds the feature vector corresponding to a protein that is in  myTable
findFeature<- function(myTable, ArabidopsisFeature, PseudomonasFeature)
{
  Nrow <- nrow(myTable)
  tempAra <- ArabidopsisFeature[which(ArabidopsisFeature[,1]==as.character(myTable[1,1])),-1]
  tempPse <- PseudomonasFeature[which(PseudomonasFeature[,1]==as.character(myTable[1,2])),-1]
  myTableset <- cbind(myTable[1,-3],tempAra,tempPse,myTable[1,3])
  for (i in 2:Nrow) 
  {
    tempAra <- ArabidopsisFeature[which(ArabidopsisFeature[,1]==as.character(myTable[i,1])),-1]
    tempPse <- PseudomonasFeature[which(PseudomonasFeature[,1]==as.character(myTable[i,2])),-1]
    temp <- cbind(myTable[i,-3],tempAra,tempPse,myTable[i,3])
    names(temp) <- names(myTableset)
    myTableset <- rbind(myTableset,temp)
  }
  myTableset = unite(myTableset, Interaction, c(Arabidopsis.Thaliana, Pseudomonas.Syringae), remove=TRUE)
  myTableset <- data.frame(myTableset[,-1], row.names=myTableset[,1])
  colnames(myTableset)[ncol(myTableset)] <- "class"
  return(myTableset)
}
#----------------------------------------------------------------------------------------------------------------------
#Normalization data set------------------------------------------------------------------------------------------------
#bring all values into the range [0, 1]--------------------------------------------------------------------------------
normalization <- function(Train, Test,indep)
{
  Ncol <- ncol(Train)-1
  columns <- c()
  for (i in 1:Ncol)
  {
    Max <- max(Train[,i])
    Min <- min(Train[,i])
    if(Max == Min)
    {
      columns <- append(columns, i)
    }
    else
    {
      Train[,i] <- (Train[,i]-Min)/(Max-Min)
      Test[,i] <- (Test[,i]-Min)/(Max-Min)
      indep[,i] <- (indep[,i]-Min)/(Max-Min)
    }
  }
  if(length(columns) > 0)
  {
    Train <- Train[, -columns]
    Test <- Test[, -columns]
    indep <- indep[, -columns]
  }
  
  corilation <- cor(Train[,-ncol(Train)])
  Ncol <- ncol(corilation)
  ccolumns <- c()
  for (i in 1:(Ncol-1))
  {
    k = i+1 
    if(length(which((corilation[i,k:Ncol]) > 0.9)) >= 1)
    {
      ccolumns <- append(ccolumns, rownames(corilation)[i])
    }
  }
  if(length(ccolumns) > 0)
  {
    Train <- Train[,!(names(Train) %in% ccolumns)]
    Test <- Test[, !(names(Test) %in% ccolumns)]
    indep <- indep[, !(names(indep) %in% ccolumns)]
  }
  return(list(Train,Test,indep))
}
#----------------------------------------------------------------------------------------------------------------------
Independent <- findFeature(Independent, IndependentAraFeature, IndependentPseuFeature)

randomTest <- findFeature(randomTest, ArabidopsisThalianaFeature, PseudomonasSyringaeFeature)
randomTrain <- findFeature(randomTrain, ArabidopsisThalianaFeature, PseudomonasSyringaeFeature)

randomResult <- normalization(randomTrain, randomTest, Independent)

write.arff(file="G:/Thesis/CodeForThesis-R/Data-FeatureExtraction/Independent-randomTest.arff", x=randomResult[3], relation = "Independent-Random Test")
write.arff(file="G:/Thesis/CodeForThesis-R/Data-FeatureExtraction/randomTest.arff", x=randomResult[2], relation = "Random Test")
write.arff(file="G:/Thesis/CodeForThesis-R/Data-FeatureExtraction/randomTrain.arff", x=randomResult[1], relation = "Random Train")

write.csv(file="G:/Thesis/CodeForThesis-R/Data-FeatureExtraction/Independent-randomTest.csv", x=randomResult[3])
write.csv(file="G:/Thesis/CodeForThesis-R/Data-FeatureExtraction/randomTest.csv", x=randomResult[2])
write.csv(file="G:/Thesis/CodeForThesis-R/Data-FeatureExtraction/randomTrain.csv", x=randomResult[1])


SDDTest <- findFeature(SDDTest, ArabidopsisThalianaFeature, PseudomonasSyringaeFeature)
SDDTrain <- findFeature(SDDTrain, ArabidopsisThalianaFeature, PseudomonasSyringaeFeature)

SDDResult <- normalization(SDDTrain, SDDTest,Independent)

write.arff(file="G:/Thesis/CodeForThesis-R/Data-FeatureExtraction/Independent-SDDTest.arff", x=SDDResult[3], relation = "Independent-SDD Test")
write.arff(file="G:/Thesis/CodeForThesis-R/Data-FeatureExtraction/SDDTest.arff", x=SDDResult[2], relation = "SDD Test")
write.arff(file="G:/Thesis/CodeForThesis-R/Data-FeatureExtraction/SDDTrain.arff", x=SDDResult[1], relation = "SDD Train")

write.csv(file="G:/Thesis/CodeForThesis-R/Data-FeatureExtraction/Independent-SDDTest.csv", x=SDDResult[3])
write.csv(file="G:/Thesis/CodeForThesis-R/Data-FeatureExtraction/SDDTest.csv", x=SDDResult[2])
write.csv(file="G:/Thesis/CodeForThesis-R/Data-FeatureExtraction/SDDTrain.csv", x=SDDResult[1])
#######################################################################################################################
endTime <- Sys.time()
#end-------------------------------------------------------------------------------------------------------------------
cat("Total Run Time is:" ,endTime-startTime)




