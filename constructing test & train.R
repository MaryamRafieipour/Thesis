#constructing test & train
startTime <- Sys.time()
#######################################################################################################################
#library---------------------------------------------------------------------------------------------------------------
library("readxl")
library("xlsx")

#######################################################################################################################
#reading files---------------------------------------------------------------------------------------------------------
#first sheet of tables1 contains gene information(name,sequence,....)
#third sheet of tables2 contains interaction data
tableS1 <- read_excel("G:/Thesis/RawData/tableS1.xlsx", sheet = "1")
tableS2 <- read_excel("G:/Thesis/RawData/tableS2.xlsx", sheet = "3")

S1 <- as.data.frame(tableS1)
S2 <- as.data.frame(tableS2)
#omitting * from end of sequence in S1
Row <- nrow(S1)
for(i in 1:Row)
{
  S1$`Protein_ sequence`[i]<-substr(S1$`Protein_ sequence`[i],1,(nchar(S1$`Protein_ sequence`[i])-1))
}


#######################################################################################################################
#finding interactions--------------------------------------------------------------------------------------------------
#all Arabidopsis Thaliana ID's is not in s1$Species--------------------------------------------------------------------
IDs <- S1$Gene_name[which(S1$Species == "Pseudomonas syringae")]
ROW <- nrow(S2)
interactions <- data.frame("Number"=integer(), "Arabidopsis thaliana ID"=character(),"Pseudomonas syringae ID"=character())
#Number in data frame shows the row in S2
for (i in 1:ROW)
{ 
  if(S2[i,1] %in% IDs | S2[i,2] %in% IDs) #finding rows that contains Pseudomonas syringae ID
  { 
    if(!(S2[i,1] %in% IDs & S2[i,2] %in% IDs)) #choosing rows which just one of them is Pseudomonas syringae ID
    {
      interactions <- rbind(interactions, data.frame("Number"=i,"Arabidopsis thaliana ID"=S2[i,1],"Pseudomonas syringae ID"=S2[i,2]))
    }
  }
}

interactions <- cbind(interactions, rep("interaction",nrow(interactions)))
names(interactions) <- c("number" ,"Arabidopsis Thaliana","Pseudomonas Syringae","state")
write.xlsx(file="G:/Thesis/CodeForThesis-R/Data/interactions.xlsx",x=interactions, row.names = FALSE)

#######################################################################################################################
#finding sequences-----------------------------------------------------------------------------------------------------
#finding Pseudomonas Syringae sequences
PseudomonasSyringae <- unique(interactions$`Pseudomonas Syringae`)

Pseudomonassyringae_sequence <- data.frame("Pseudomonas syringae ID"=character(), "sequence"=character())
len <- length(PseudomonasSyringae)
for (i in 1:len)
{ 
  if(PseudomonasSyringae[i] %in% IDs){temp <-  S1$`Protein_ sequence`[which(S1$Gene_name == PseudomonasSyringae[i])]}
  else{temp <- "not found"}
  Pseudomonassyringae_sequence <- rbind( Pseudomonassyringae_sequence, data.frame(PseudomonasSyringae[i], temp ))
}
names(Pseudomonassyringae_sequence) <- c("Gene names", "Sequence")
write.xlsx(file="G:/Thesis/RawData/Data-Sequence/Pseudomonassyringa_sequence.xlsx",x=Pseudomonassyringae_sequence, row.names = FALSE)

#######################################################################################################################
#constructing nagative data--------------------------------------------------------------------------------------------
#CalculateAdjMatrit return adjacency matrix----------------------------------------------------------------------------
CalculateAdjMatrit <- function(myMatrix){
  
  rowNames <- sort(unique(as.character(myMatrix[,2])))
  colNames <- sort(unique(as.character(myMatrix[,3])))
  adjMatrix <-matrix(0,nrow=length(rowNames),ncol = length(colNames))
  rownames(adjMatrix) <- rowNames
  colnames(adjMatrix) <- colNames
  
  rowMyMatrix <- nrow(myMatrix)
  for( i in 1:rowMyMatrix)
  {
    adjMatrix[which(rowNames==myMatrix[i,2]),which(colNames==myMatrix[i,3])]=1
  }
  return(adjMatrix)
}
adjacentMatrix <- CalculateAdjMatrit(interactions)

#A) constructing nagative data-----------------------------------------------------------------------------------------
#constructing negative data that the degree distribution of negative network is as same as degree distribution of positive network.
#adding a column at end of the matrix. This show the degree of each Arabidopsis Thaliana protein
sumArabidopsisThaliana <- rep(0, nrow(adjacentMatrix))
myMatrix <- cbind(adjacentMatrix, sumArabidopsisThaliana)
myMatrix[,"sumArabidopsisThaliana"] <- rowSums(myMatrix)
#adding a row at end of the matrix. 
sumPseudomonasSyringae <- rep(0,ncol(myMatrix))
myMatrix <- rbind(myMatrix, sumPseudomonasSyringae)
rownames(myMatrix)[nrow(myMatrix)] <- "sumPseudomonasSyringae"
#sort descending matrix by the degree of Arabidopsis Thaliana protein 
myMatrix <- myMatrix[order(-myMatrix[,"sumArabidopsisThaliana"]),]
#calculate the degree of each Pseudomonas Syringae protein and number of interactions
myMatrix["sumPseudomonasSyringae",] <- colSums(myMatrix)
#----------------------------------------------------------------------------------------------------------------------
NrowMyMatrix <- nrow(myMatrix)
NcolMyMatrix <- ncol(myMatrix)
NegativeData <- data.frame("Arabidopsis thaliana ID"=character(),"Pseudomonas syringae ID"=character())
for(i in 1:(NrowMyMatrix-1))
{
  probability <- rep(0,NcolMyMatrix-1)
  for(j in 1:(NcolMyMatrix-1))
  {
    if(myMatrix[i,j]==0 && myMatrix[NrowMyMatrix,j]> 0)#there is not a interaction and there is a chance to choose this protein
    {
      probability[j]<- myMatrix[NrowMyMatrix,j]
    }
    
  }
  #sampling according the probability is calculated 
  mySample <- sample(myMatrix[i,-NcolMyMatrix],size=myMatrix[i,"sumArabidopsisThaliana"],prob = probability/myMatrix[i,"sumArabidopsisThaliana"])
  #set -1 to element that is chosen for negative data 
  for(j in 1:length(mySample))
  {
    myMatrix[i,names(mySample)[j]] <- -1
    NegativeData <-  rbind(NegativeData , data.frame(rownames(myMatrix)[i],names(mySample)[j]))
  }
  #updating the amount of last row and column 
  myMatrix[i,"sumArabidopsisThaliana"] <- sum(myMatrix[i,1:NcolMyMatrix-1])
  myMatrix["sumPseudomonasSyringae",] <- colSums(myMatrix[-NrowMyMatrix,])
}
NegativeData <- cbind(NegativeData, rep("noninteraction",nrow(NegativeData)))
names(NegativeData) <- c("Arabidopsis Thaliana","Pseudomonas Syringae","state")
write.xlsx(file="G:/Thesis/CodeForThesis-R/Data/SDDNegativeData.xlsx",x=NegativeData, row.names = FALSE)

#B) constructing negative data-----------------------------------------------------------------------------------------
#constructing randomly negative data
negativetable <- which(adjacentMatrix==0, arr.ind = T)
numberOfInteractions <- sum(adjacentMatrix)
randomNegative <- negativetable[sample(nrow(negativetable),numberOfInteractions),]
randomNegativeData <- data.frame("Arabidopsis thaliana ID"=character(),"Pseudomonas syringae ID"=character())
for(i in 1:numberOfInteractions)
{
  randomNegativeData <- rbind(randomNegativeData , data.frame(rownames(adjacentMatrix)[randomNegative[i,1]],colnames(adjacentMatrix)[randomNegative[i,2]]))
}
randomNegativeData <- cbind(randomNegativeData, rep("noninteraction",nrow(randomNegativeData)))
names(randomNegativeData) <- c("Arabidopsis Thaliana","Pseudomonas Syringae", "state")
write.xlsx(file="G:/Thesis/CodeForThesis-R/Data/randomNegativeData.xlsx",x=randomNegativeData, row.names = FALSE)

#######################################################################################################################
#test and train--------------------------------------------------------------------------------------------------------
# choosing 80 percent of each data set for training and 20 percent for test 
Nrow = nrow(interactions)
pTest <- sample(Nrow, floor(0.2*Nrow))
positiveTest <- interactions[pTest,-1]
positiveTrain <- interactions[-pTest,-1]

nTest <- sample(Nrow, floor(0.2*Nrow))
negativeTest <- NegativeData[nTest,]
negativeTrain <- NegativeData[-nTest,]

nrTest <- sample(Nrow, floor(0.2*Nrow))
randomNegativeTest <- randomNegativeData[nrTest,]
randomNegativeTrain <- randomNegativeData[-nrTest,]
#----------------------------------------------------------------------------------------------------------------------
#randomly Data set-----------------------------------------------------------------------------------------------------
randomTest <- rbind(positiveTest,randomNegativeTest)
randomTrain <- rbind(positiveTrain,randomNegativeTrain)

write.csv(file="G:/Thesis/CodeForThesis-R/Data/randomTest.csv",x=randomTest, row.names = FALSE)
write.csv(file="G:/Thesis/CodeForThesis-R/Data/randomTrain.csv",x=randomTrain, row.names = FALSE)

#same degree distribution Data set-------------------------------------------------------------------------------------
SDDTest <- rbind(positiveTest,negativeTest)
SDDTrain <- rbind(positiveTrain,negativeTrain)

write.csv(file="G:/Thesis/CodeForThesis-R/Data/SDDTest.csv",x=SDDTest, row.names = FALSE)
write.csv(file="G:/Thesis/CodeForThesis-R/Data/SDDTrain.csv",x=SDDTrain, row.names = FALSE)

#######################################################################################################################
endTime <- Sys.time()
#end-------------------------------------------------------------------------------------------------------------------
cat("Total Run Time is:" ,endTime-startTime)

