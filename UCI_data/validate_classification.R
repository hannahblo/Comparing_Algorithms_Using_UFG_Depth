###
# input:  train, test
# output: metrics to be evaluated
# better inside the CV-loop to guarantee that same splits
library(dplyr)
library(pROC)
library(stringr)
library(glmnet)
library(randomForest)
library(ggplot2)
library(tidyr)
library(gbm)
library(rpart)
# validation function for xgboost


validate_rf = function(xtrain, xtest, ytrain, ytest){
  rfmod            = randomForest(x = xtrain, 
                                  y=factor(ytrain, levels = c("0", "1")))
  ypred            = predict(rfmod,
                             xtest,
                             type = "prob")[,2]
  auc_rf           = pROC::auc(ytest, ypred)
  acc_rf           = mean(ifelse(ypred > 0.5, 1, 0) ==  ytest)
  brier_rf         = brier_score(ytest, ypred)
  
  list(data.frame(method = c("rf", 
                             "rf", 
                             "rf"),
                  metric = c("auc", "acc", "brier_score"),
                  val = c(auc_rf,
                          acc_rf,
                          brier_rf)
                  ),
       model = NULL)
}

validate_gbm = function(xtrain, xtest, ytrain, ytest){
  gbmod            = gbm(y~., data = data.frame(xtrain,y=as.numeric(factor(ytrain, levels = c("0", "1")))-1), distribution = "bernoulli",
                         n.trees = 300,
                         interaction.depth = 3,
                         shrinkage = 0.02)
  ypred            = predict(gbmod, newdata = data.frame(xtest), type = "response", n.trees = 300)
  auc_rf           = pROC::auc(ytest, ypred)
  acc_rf           = mean(ifelse(ypred > 0.5, 1, 0) ==  ytest)
  brier_rf         = brier_score(ytest, ypred)
  
  list(data.frame(method = c("GBM",
                             "GBM",
                             "GBM"),
                  metric = c("auc",
                             "acc",
                             "brier_score"), 
                  val = c(auc_rf,
                          acc_rf,
                          brier_rf)),
       model = NULL)
}

validate_gbm_stump = function(xtrain, xtest, ytrain, ytest){
  gbmod            = gbm(y~., data = data.frame(xtrain,y=as.numeric(factor(ytrain, levels = c("0", "1")))-1), distribution = "bernoulli",
                         n.trees = 500,
                         interaction.depth = 1,
                         shrinkage = 0.05)
  ypred            = predict(gbmod, newdata = data.frame(xtest), type = "response",
                             n.trees = 500)
  auc_rf           = pROC::auc(ytest, ypred)
  acc_rf           = mean(ifelse(ypred > 0.5, 1, 0) ==  ytest)
  brier_rf         = brier_score(ytest, ypred)
  
  list(data.frame(method = c("Boosted Stumps",
                             "Boosted Stumps",
                             "Boosted Stumps"),
                  metric = c("auc",
                             "acc",
                             "brier_score"), 
                  val = c(auc_rf,
                          acc_rf,
                          brier_rf)),
       model = NULL)
}

validate_rpart = function(xtrain, xtest, ytrain, ytest){
  rfmod            = rpart(y~., data = data.frame(xtrain,y=factor(ytrain, levels = c("0", "1"))))
  ypred            = predict(rfmod, data.frame(xtest), type = "prob")[,2]
  auc_rf           = pROC::auc(ytest, ypred)
  acc_rf           = mean(ifelse(ypred > 0.5, 1, 0) ==  ytest)
  brier_rf         = brier_score(ytest, ypred)
  
  list(data.frame(method = c("CART",
                             "CART",
                             "CART"),
                  metric = c("auc",
                             "acc",
                             "brier_score"), 
                  val = c(auc_rf,
                          acc_rf,
                          brier_rf)),
       model = NULL)
}

validate_lasso = function(xtrain, xtest, ytrain, ytest){
  glmmod           = cv.glmnet(x = xtrain, y = ytrain, familiy = "binomial")
  ypred            = predict(glmmod, xtest, type = "response" ,s = "lambda.min")
  auc_rf           = pROC::auc(ytest, ypred)
  acc_rf           = mean(ifelse(ypred > 0.5, 1, 0) ==  ytest)
  brier_rf         = brier_score(ytest, ypred)
  
  list(data.frame(method = c("Lasso",
                             "Lasso",
                             "Lasso"),
                  metric = c("auc",
                             "acc",
                             "brier_score"), val = c(auc_rf, acc_rf, brier_rf)),
       model = NULL)
}

validate_ridge = function(xtrain, xtest, ytrain, ytest){
  glmmod           = cv.glmnet(x = xtrain, y = ytrain, familiy = "binomial", alpha = 0)
  ypred            = predict(glmmod, xtest, type = "response" ,s = "lambda.min")
  auc_rf           = pROC::auc(ytest, ypred)
  acc_rf           = mean(ifelse(ypred > 0.5, 1, 0) ==  ytest)
  brier_rf         = brier_score(ytest, ypred)
  
  list(data.frame(method = c("Ridge",
                             "Ridge",
                             "Ridge"),
                  metric = c("auc",
                             "acc",
                             "brier_score"), val = c(auc_rf, acc_rf, brier_rf)),
       model = NULL)
}

validate_enet = function(xtrain, xtest, ytrain, ytest){
  glmmod           = cv.glmnet(x = xtrain, y = ytrain, familiy = "binomial", alpha = 0.5)
  ypred            = predict(glmmod, xtest, type = "response" ,s = "lambda.min")
  auc_rf           = pROC::auc(ytest, ypred)
  acc_rf           = mean(ifelse(ypred > 0.5, 1, 0) ==  ytest)
  brier_rf         = brier_score(ytest, ypred)
  
  list(data.frame(method = c("ElasticNet",
                             "ElasticNet",
                             "ElasticNet"),
                  metric = c("auc",
                             "acc",
                             "brier_score"), val = c(auc_rf, acc_rf, brier_rf)),
       model = NULL)
}

validate_glm = function(xtrain, xtest, ytrain, ytest){
  glmmod             = glm(y~., data = data.frame(xtrain,y=factor(ytrain, levels = c("0", "1"))), family = binomial())
  ypred              = predict(glmmod, data.frame(xtest), type = "response")
  auc_rf             = pROC::auc(ytest, ypred)
  acc_rf             = mean(ifelse(ypred > 0.5, 1, 0) ==  ytest)
  brier_rf           = brier_score(ytest, ypred)
  
  list(data.frame(method = c("glm",
                             "glm",
                             "glm"),
                  metric = c("auc",
                             "acc", 
                             "brier_score"),
                  val = c(auc_rf, 
                          acc_rf,
                          brier_rf)),
       model = NULL)
}

brier_score = function(y,ypred){
y = as.numeric(as.factor(y))-1
mean((ypred-y)^2)  
}

methods = list(
  CART              = validate_rpart,
  rf                = validate_rf,
  Lasso             = validate_lasso,
  Ridge             = validate_ridge,
  ElasticNet        = validate_enet,
  glm               = validate_glm,
  GBM               = validate_gbm,
  BoostedStumps     = validate_gbm_stump
               )
mean_impute = function(x){
  for(j in 1:ncol(x)){
    x[is.na(x[,j]),j] = mean(x[,j], na.rm=T)
  }
  x
}
datasets = list.files("datasets_class/")
#length(datasets)
#length(datasets)
for(d in 1:length(datasets)){
  data        = read.csv(paste0("datasets_class/", datasets[d]))
  dataname    = strsplit(datasets[d], "[.]")[[1]][1]
  X           = data[,-ncol(data)]
  X           = apply(X, 2, as.numeric)
  # X           = 0.4*scale(X)
  X           = mean_impute(X)
  colnames(X) = paste0("feature_",1:ncol(X))
  y           = data[,ncol(data)]
  y           = as.numeric(as.factor(y))-1
  set.seed(1987)
  ids         = sample(1:nrow(X))
  fold        = rep(1:10, length.out = nrow(X))
    for(l in 1:length(methods)){
      res         = list()
      rule_models = list()
      method      = methods[[l]]
      methodname  = names(methods)[l]
      ### foreach!
      for(i  in 1:10){
        xtrain           = X[ids[fold != i],]
        xtest            = X[ids[fold == i],]
        ytrain           = y[ids[fold != i]]
        ytest            = y[ids[fold == i]]
        temp             = method(xtrain, xtest, ytrain, ytest)
        res[[i]]         = data.frame(fold = i, temp[[1]])
      }
      aggregate          = do.call('rbind', res)
      aggregate$dataset  = dataname
      write.csv(aggregate, paste0("results/", dataname, "_", methodname, ".csv"), row.names = F)
  }
}
