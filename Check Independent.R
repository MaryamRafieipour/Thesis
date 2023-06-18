library("readxl")
library("xlsx")

independentPse<- read_xlsx("D:/Thesis/RawData/Data-Sequence/IndependentTest-Sequence.xlsx", sheet = "Pseudomonas syringae")
compeletpse <- read_xlsx("D:/Thesis/RawData/Data-Sequence/Pseudomonassyringa_sequence.xlsx", sheet = "Sheet1")

for(i in 1:nrow(independentPse))
{
  if(independentPse$Sequence[i] %in% compeletpse$Sequence)
  {
    t <- which(compeletpse$Sequence == independentPse$Sequence[i])
    print(t)
    print(compeletpse$`Gene names`[t])
    print(i)
  }
}

ArabidopsisThalianaTable <- read_excel("D:/Thesis/RawData/Data-Sequence/ArabidopsisThaliana_uniprot_sequence.xlsx")
AT4G17680 <- read_excel("D:/Thesis/RawData/Data-Sequence/AT4G17680.xlsx")
AT1G24050 <- read_excel("D:/Thesis/RawData/Data-Sequence/AT1G24050.xlsx")
AT3G11590 <- read_excel("D:/Thesis/RawData/Data-Sequence/AT3G11590.xlsx")
AT3G11720 <- read_excel("D:/Thesis/RawData/Data-Sequence/AT3G11720.xlsx")
AT3G48550 <- read_excel("D:/Thesis/RawData/Data-Sequence/AT3G48550.xlsx")
AT4G02550 <- read_excel("D:/Thesis/RawData/Data-Sequence/AT4G02550.xlsx")

independentAra<- read_xlsx("D:/Thesis/RawData/Data-Sequence/IndependentTest-Sequence.xlsx", sheet = "Arabidopsis thaliana")
completara <- rbind(ArabidopsisThalianaTable,AT4G17680,AT1G24050,AT3G11590,AT3G11720,AT3G48550,AT4G02550)

for(i in 1:nrow(independentAra))
{
  if(independentAra$Sequence[i] %in% completara$Sequence)
  {
    t <- which(completara$Sequence == independentAra$Sequence[i])
    print(t)
    print(completara$`Gene names`[t])
    print(i)
  }
}

#result pseu: HopAO1_Pto DC3000 
#result Arab: AT3G25070