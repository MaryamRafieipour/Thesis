library(foreign)
library(caret)
s <- Sys.time()

allData <- read.arff("D:/Thesis/CodeForThesis-R/randomTrain1.arff")

x <- allData[,-ncol(allData)]
y <- as.factor(allData[,ncol(allData)])


subsets <- c(1:100)

ctrl <- rfeControl(functions = rfFuncs,
                   method = "repeatedcv",
                   repeats = 5,
                   verbose = FALSE)

lmProfile <- rfe(x, y,
                 sizes = subsets,
                 rfeControl = ctrl)

e <- Sys.time()
cat("total run time is:", e-s)

