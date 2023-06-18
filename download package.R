#' Get package dependencies
#'
#' @param packs A string vector of package names
#'
#' @return A string vector with packs plus the names of any dependencies

S <- Sys.time()
getDependencies <- function(packs){
  dependencyNames <- unlist(
    tools::package_dependencies(packages = packs, db = available.packages(), 
                                which = c("Depends", "Imports"),
                                recursive = TRUE))
  packageNames <- union(packs, dependencyNames)
  packageNames
}
# Calculate dependencies
packages <- getDependencies(c("readxl", "xlsx","foreign","tidyr","ftrCOOL","caret","AppliedPredictiveModeling",
                              "ggplot2","tidyverse","e1071","networkD3","doParallel","pROC","ROSE","recipes",
                              "dplyr","kernlab","caTools","mboost","pls","randomForest","reshape2","C50",
                              "Hmisc","desirability","rsample", "fastAdaboost","adabag","plyr","adaptDA",
                              "ipred","earth","mda","logicFS","bartMachine","arm","binda","ada","import",
                              "bst","party","partykit","RWeka","rpart","rpartScore","CHAID","VGAM",
                              "deepboost","sparsediscrim","kerndwd","randomGLM","xgboost","elmNN","HiDimDA",
                              "frbs","gam","mgcv","MASS","gpls","glmnet","Matrix","h2o","proxy","protoclass","hda",
                              "HDclassif","kknn","LiblineaR","class","klaR","LogicReg","bnclassify","RRF",
                              "nnet","monmlp","RSNNS","msaenet","FCNN4R","keras","naivebayes","pamr","mxnet",
                              "obliqueRF","snn","foreach","partDSA","plsRglm","supervisedPRIM",
                              "penalizedLDA","stepPlr","ordinalNet","rFerns","ranger","ordinalForest","Rborist",
                              "extraTrees","inTrees","rrcov","robustDA","rrlda", "rrcovHD","rocc","rotationForest",
                              "kohonen","sda","sdwd","sparseLDA","spls","spikeslab","deepnet","gbm","evtree",
                              "nodeHarvest","vbmp","wsrf"))

# Download the packages to the working directory.
# Package names and filenames are returned in a matrix.
setwd("E:/packages/")
pkgInfo <- download.packages(pkgs = packages, destdir = getwd(), type = "win.binary")
# Save just the package file names (basename() strips off the full paths leaving just the filename)
write.csv(file = "pkgFilenames.csv", basename(pkgInfo[, 2]), row.names = FALSE)
E <- Sys.time()
cat("total Run time is", E-S)

# Set working directory to the location of the package files
setwd("E:/packages/")
# Read the package filenames and install
pkgFilenames <- read.csv("pkgFilenames.csv", stringsAsFactors = FALSE)[, 1]
install.packages(pkgFilenames, repos = NULL, type = "win.binary")