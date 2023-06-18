#library---------------------------------------------------------------------------------------------------------------
library("readxl")
library("xlsx")
#----------------------------------------------------------------------------------------------------------------------
ID <- read.xlsx("G:/Thesis/RawData/maping Id.xlsx", sheetName = "Sheet1")
intact <- read.xlsx("G:/Thesis/RawData/intact.xlsx", sheetName = "Sheet1")

add_organism <- data.frame("A"=character(),"B"=character())
Nrow <- nrow(intact)
for (i in 1:Nrow)
{
  if(intact$interactor.A[i] %in% ID$From & intact$interactor.B[i] %in% ID$From)
  {
    tempa <- ID$name[which(ID$From == intact$interactor.A[i])]
    tempb <- ID$name[which(ID$From == intact$interactor.B[i])]
    Temp <- cbind(tempa, tempb)
  }
  else {Temp <- c("not", "not")}
  add_organism <- rbind(add_organism, Temp)
}
newIntact <- cbind(intact, add_organism)
write.csv("G:/Thesis/CodeForThesis-R/Data/newIntact.csv", x=newIntact, row.names = FALSE)

IndependentTest <- data.frame("Arabidopsis thaliana ID"=character(), "Pseudomonas syringae ID"=character(), 
                              "Arabidopsis thaliana"=character(), "Pseudomonas syringae"=character())
temp <- data.frame("Arabidopsis thaliana ID"=character(), "Pseudomonas syringae ID"=character(), 
                   "Arabidopsis thaliana"=character(), "Pseudomonas syringae"=character())
name <- c("Arabidopsis thaliana ID", "Pseudomonas syringae ID", "Arabidopsis thaliana", "Pseudomonas syringae")
Nrow <- nrow(newIntact)
for (i in 1:Nrow)
{
  if(newIntact$tempa[i] == "Arabidopsis thaliana" & newIntact$tempb[i] == "Pseudomonas syringae") 
  {
    temp <- data.frame(newIntact$interactor.A[i],newIntact$interactor.B[i],newIntact$tempa[i],newIntact$tempb[i])
    names(temp) <- name
    IndependentTest <- rbind(IndependentTest,data.frame(temp))
  }
  else if(newIntact$tempb[i]== "Arabidopsis thaliana" & newIntact$tempa[i] == "Pseudomonas syringae") 
  {
    temp <- data.frame(newIntact$interactor.B[i],newIntact$interactor.A[i],newIntact$tempb[i],newIntact$tempa[i])
    names(temp) <- name
    IndependentTest <- rbind(IndependentTest,data.frame(temp))
  }
}

IndependentTest <- unique(IndependentTest[,c(1,2)])
write.csv("G:/Thesis/CodeForThesis-R/Data/IndependentTest.csv", x=IndependentTest, row.names = FALSE)