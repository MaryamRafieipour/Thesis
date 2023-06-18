library(caret)
library(Hmisc)
library(randomForest)
library(e1071)
library(ipred)
library(foreign)
library(adabag)
library(plyr)
library(Rcpp)
library(MASS)
library(klaR)

#SDDTrain---------------------------------------------------------------------------------------------------------------------------------
s <- Sys.time()

allData <- read.arff("D:/Thesis/CodeForThesis-R/Data-FeatureExtraction/SDDTrain.arff")

X <- allData[,-ncol(allData)]
Y <- as.factor(allData$class)

featureImportance <-data.frame('name'= names(X),'score'=rep(0,ncol(X)))
control <- trainControl(method="repeatedcv", number=10, repeats=3)


#adabag
search_grid_adabag <- expand.grid(mfinal=c(5,10,15,20,25), maxdepth = c(3,5,7,10))
model_adabag <- train(x = X, y = Y, method="AdaBag", trControl=control, tuneGrid = search_grid_adabag)
importance_adabag <- varImp(model_adabag, scale=FALSE)
write.csv("D:/Thesis/CodeForThesis-R/feature selection/SDD_importance_adabag.csv", x= importance_adabag$importance)

#knn
search_grid_knn <- expand.grid(k=c(1,3,5,7,11))
model_knn <- train(x = X, y = Y, method="knn", trControl=control, tuneGrid = search_grid_knn)
importance_knn <- varImp(model_knn, scale=FALSE)
write.csv("D:/Thesis/CodeForThesis-R/feature selection/SDD_importance_knn.csv", x= importance_knn$importance)

#randomforest
search_grid_rf <- expand.grid(mtry=c(10,15,20))
model_rf <- train(x = X, y = Y, method="rf", trControl=control,tuneGrid = search_grid_rf)
importance_rf <- varImp(model_rf, scale=FALSE)
write.csv("D:/Thesis/CodeForThesis-R/feature selection/SDD_importance_rf.csv", x= importance_rf$importance)

#naiveBayes
search_grid_nb <- expand.grid(usekernel=c(TRUE),fL=0:5,adjust=seq(0,5,by=1) )
model_nb <- train(x = X, y = Y, method="nb", trControl=control,tuneGrid = search_grid_nb)
importance_nb <- varImp(model_nb, scale=FALSE)
write.csv("D:/Thesis/CodeForThesis-R/feature selection/SDD_importance_nb.csv", x= importance_nb$importance)

#treebag
model_treebag <- train(x = X, y = Y, method="treebag", trControl=control)
importance_treebag <- varImp(model_treebag, scale=FALSE)
write.csv("D:/Thesis/CodeForThesis-R/feature selection/SDD_importance_treebag.csv", x= importance_treebag$importance)


for (i in featureImportance$name)
{

featureImportance$score[which(featureImportance$name==i)] <- (importance_adabag$importance$Overall[which(rownames(importance_adabag$importance)==i)] +
                                                              importance_knn$importance$interaction [which(rownames(importance_knn$importance)==i)] +
                                                              importance_rf$importance$Overall [which(rownames(importance_rf$importance)==i)]+
                                                              importance_nb$importance$interaction [which(rownames(importance_nb$importance)==i)]+
                                                              importance_treebag$importance$Overall [which(rownames(importance_treebag$importance)==i)])
}

featureImportanceSort <- featureImportance[order(-featureImportance$score),]

write.csv("D:/Thesis/CodeForThesis-R/feature selection/importance_feature_SDD.csv", x= featureImportanceSort)

e<-Sys.time()

#randomTrain----------------------------------------------------------------------------------------------------------------------------------------

s1 <- Sys.time()

allData <- read.arff("D:/Thesis/CodeForThesis-R/Data-FeatureExtraction/randomTrain.arff")

X <- allData[,-ncol(allData)]
Y <- as.factor(allData$class)

featureImportance <-data.frame('name'= names(X),'score'=rep(0,ncol(X)))
control <- trainControl(method="repeatedcv", number=10, repeats=3)

#adabag
search_grid_adabag <- expand.grid(mfinal=c(5,10,15,20,25), maxdepth = c(3,5,7,10))
model_adabag <- train(x = X, y = Y, method="AdaBag", trControl=control, tuneGrid = search_grid_adabag)
importance_adabag <- varImp(model_adabag, scale=FALSE)
write.csv("D:/Thesis/CodeForThesis-R/feature selection/random_importance_adabag.csv", x= importance_adabag$importance)

#knn
search_grid_knn <- expand.grid(k=c(1,3,5,7,11))
model_knn <- train(x = X, y = Y, method="knn", trControl=control, tuneGrid = search_grid_knn)
importance_knn <- varImp(model_knn, scale=FALSE)
write.csv("D:/Thesis/CodeForThesis-R/feature selection/random_importance_knn.csv", x= importance_knn$importance)

#randomforest
search_grid_rf <- expand.grid(mtry=c(10,15,20))
model_rf <- train(x = X, y = Y, method="rf", trControl=control,tuneGrid = search_grid_rf)
importance_rf <- varImp(model_rf, scale=FALSE)
write.csv("D:/Thesis/CodeForThesis-R/feature selection/random_importance_rf.csv", x= importance_rf$importance)

#naiveBayes
search_grid_nb <- expand.grid(usekernel=c(TRUE), fL=0:5,adjust=seq(0,5,by=1) )
model_nb <- train(x = X, y = Y, method="nb", trControl=control,tuneGrid = search_grid_nb)
importance_nb <- varImp(model_nb, scale=FALSE)
write.csv("D:/Thesis/CodeForThesis-R/feature selection/random_importance_nb.csv", x= importance_nb$importance)

#treebag
model_treebag <- train(x = X, y = Y, method="treebag", trControl=control)
importance_treebag <- varImp(model_treebag, scale=FALSE)
write.csv("D:/Thesis/CodeForThesis-R/feature selection/random_importance_treebag.csv", x= importance_treebag$importance)




for (i in featureImportance$name)
{
  
  featureImportance$score[which(featureImportance$name==i)] <- (importance_adabag$importance$Overall[which(rownames(importance_adabag$importance)==i)] +
                                                                  importance_knn$importance$interaction [which(rownames(importance_knn$importance)==i)] +
                                                                  importance_rf$importance$Overall [which(rownames(importance_rf$importance)==i)]+
                                                                  importance_nb$importance$interaction [which(rownames(importance_nb$importance)==i)]+
                                                                  importance_treebag$importance$Overall [which(rownames(importance_treebag$importance)==i)])
}

featureImportanceSort <- featureImportance[order(-featureImportance$score),]

write.csv("D:/Thesis/CodeForThesis-R/feature selection/importance_feature_random.csv", x= featureImportanceSort)

e1 <- Sys.time()

cat("Total Run Time is SDD", (e-s))
cat("Total Run Time is Random", (e1-s1))